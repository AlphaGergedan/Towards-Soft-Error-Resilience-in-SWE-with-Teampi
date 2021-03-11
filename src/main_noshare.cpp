#include <mpi.h>
#include <teaMPI.h>
#include <unistd.h>

#include <algorithm>
#include <climits>
#include <csetjmp>
#include <fstream>
#include <iostream>
#include <sstream>

#include "blocks/DimSplitMPIOverdecomp.hpp"
#include "io/Reader.hpp"
#include "io/Writer.hpp"
#include "scenarios/LoadNetCDFScenario.hpp"
#include "scenarios/simple_scenarios.hpp"
#include "tools/Args.hpp"
#include "types/Boundary.hpp"

// some data needs to be global, otherwise the checkpoint callbacks cannot
// access it.
jmp_buf jumpBuffer{};
int ranksPerTeam{1};
std::string restartNameInput{""};
std::vector<std::string> backupMetadataNames{};
std::string outputTeamName{""};
std::string backupTeamName{""};
std::string outputNameInput{""};
std::string backupNameInput{""};

unsigned int numCheckpoints{0};

float t(0.0F);
std::vector<std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>> simulationBlocks{};
SWE_Scenario* scenario{nullptr};

void createCheckpointCallback(std::vector<int> failedTeams) {
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    const int myTeam = TMPI_GetTeamNumber();
    std::printf("Rank %i of Team %i writing checkpoint for restoration.\n", myRank, myTeam);

    // write a checkpoint for every block
    for (int i = 0; i < simulationBlocks.size(); i++) {
        simulationBlocks[i]->createCheckpoint(t, backupMetadataNames[i], 0);
    }

    int send = 1;
    for (const auto& failedTeam : failedTeams) {
        MPI_Send(&send, 1, MPI_INT, TMPI_TeamToWorldRank(myRank, failedTeam), 0, TMPI_GetWorldComm());
    }
}

void loadCheckpointCallback(int reloadTeam) {
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    const int myTeam = TMPI_GetTeamNumber();
    std::printf("Rank %i of Team %i loading checkpoint from Team %i.\n", myRank, myTeam, reloadTeam);
    MPI_Status status;
    int recv_buf;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Recv(&recv_buf, 1, MPI_INT, TMPI_TeamToWorldRank(myRank, reloadTeam), 0, TMPI_GetWorldComm(), &status);
    simulationBlocks.clear();
    delete scenario;
    scenario = nullptr;
    restartNameInput = backupNameInput + "_t" + std::to_string(reloadTeam);
    longjmp(jumpBuffer, 1);
}

std::array<int, 4> getNeighbours(int localBlockPositionX, int localBlockPositionY, int blockCountX, int blockCountY,
                                 int myRank) {
    std::array<int, 4> myNeighbours;
    myNeighbours[BND_LEFT] = (localBlockPositionX > 0) ? myRank - blockCountY : -1;
    myNeighbours[BND_RIGHT] = (localBlockPositionX < blockCountX - 1) ? myRank + blockCountY : -1;
    myNeighbours[BND_BOTTOM] = (localBlockPositionY > 0) ? myRank - 1 : -1;
    myNeighbours[BND_TOP] = (localBlockPositionY < blockCountY - 1) ? myRank + 1 : -1;
    return myNeighbours;
}

int main(int argc, char** argv) {
    // Define command line arguments
    tools::Args args;

    args.addOption("simulation-duration", 't', "Time in seconds to simulate");
    args.addOption("checkpoint-interval", 'c', "Seconds to wait between checkpoints");
    args.addOption("resolution-x", 'x', "Number of simulated cells in x-direction");
    args.addOption("resolution-y", 'y', "Number of simulated cells in y-direction");
    args.addOption("output-basepath", 'o', "Output base file name");
    args.addOption("backup-basepath", 'b', "Output base file name");
    args.addOption("restart-basepath", 'r', "Restart base file name", tools::Args::Required, false);
    args.addOption("decomp-factor", 'd', "Blocks per rank", tools::Args::Required, false);

    // Parse command line arguments
    tools::Args::Result ret = args.parse(argc, argv);
    switch (ret) {
    case tools::Args::Error:
        return 1;
    case tools::Args::Help:
        return 0;
    case tools::Args::Success:
        break;
    }

    // Read in command line arguments
    float simulationDuration{args.getArgument<float>("simulation-duration")};
    clock_t checkpointInterval{args.getArgument<clock_t>("checkpoint-interval")};
    int nxRequested{args.getArgument<int>("resolution-x")};
    int nyRequested{args.getArgument<int>("resolution-y")};
    outputNameInput = args.getArgument<std::string>("output-basepath");
    backupNameInput = args.getArgument<std::string>("backup-basepath");
    if (args.isSet("restart-basepath")) {
        restartNameInput = args.getArgument<std::string>("restart-basepath");
    }

    int blocksPerRank = 1;
    if (args.isSet("decomp-factor")) {
        blocksPerRank = args.getArgument<int>("decomp-factor");
    }

    // init teaMPI
    std::function<void(std::vector<int>)> create(createCheckpointCallback);
    std::function<void(int)> load(loadCheckpointCallback);
    TMPI_SetCreateCheckpointCallback(&create);
    TMPI_SetLoadCheckpointCallback(&load);
    TMPI_SetErrorHandlingStrategy(TMPI_WarmSpareErrorHandler);

    // init MPI
    int myRank;
    int provided;
    int requested = MPI_THREAD_MULTIPLE;
    if (setjmp(jumpBuffer) == 0) {
        MPI_Init_thread(&argc, &argv, requested, &provided);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &ranksPerTeam);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    const int myTeam{TMPI_GetTeamNumber()};
    outputTeamName = outputNameInput + "_t" + std::to_string(myTeam);
    backupTeamName = backupNameInput + "_t" + std::to_string(myTeam);

    // Print status
    char hostname[HOST_NAME_MAX];
    gethostname(hostname, HOST_NAME_MAX);

    std::printf("Rank %i of Team %i spawned at %s\n", myRank, myTeam, hostname);
    int totalBlocks = blocksPerRank * ranksPerTeam;

    // number of SWE-Blocks in x- and y-direction
    int blockCountY = std::sqrt(totalBlocks);
    while (totalBlocks % blockCountY != 0) blockCountY--;
    int blockCountX = totalBlocks / blockCountY;

    int startPoint = myRank * blocksPerRank;

    float simulationStart{0.0f};

    if (restartNameInput == "") {
        scenario = new SWE_RadialBathymetryDamBreakScenario{};
        int widthScenario = scenario->getBoundaryPos(BND_RIGHT) - scenario->getBoundaryPos(BND_LEFT);
        int heightScenario = scenario->getBoundaryPos(BND_TOP) - scenario->getBoundaryPos(BND_BOTTOM);

        float dxSimulation = (float)widthScenario / nxRequested;
        float dySimulation = (float)heightScenario / nyRequested;

        for (int i = startPoint; i < startPoint + blocksPerRank; i++) {
            auto myBlock = i;
            int localBlockPositionX = myBlock / blockCountY;
            int localBlockPositionY = myBlock % blockCountY;

            // compute local number of cells for each SWE_Block w.r.t. the
            // simulation domain (particularly not the original scenario
            // domain, which might be finer in resolution) (blocks at the
            // domain boundary are assigned the "remainder" of cells)

            int nxBlockSimulation = (nxRequested) / blockCountX;
            int nxRemainderSimulation = nxRequested - (blockCountX - 1) * nxBlockSimulation;
            int nyBlockSimulation = nyRequested / blockCountY;
            int nyRemainderSimulation = nyRequested - (blockCountY - 1) * nyBlockSimulation;

            int nxLocal = (localBlockPositionX < blockCountX - 1) ? nxBlockSimulation : nxRemainderSimulation;
            int nyLocal = (localBlockPositionY < blockCountY - 1) ? nyBlockSimulation : nyRemainderSimulation;

            // Compute the origin of the local simulation block w.r.t. the
            // original scenario domain.
            float localOriginX =
                scenario->getBoundaryPos(BND_LEFT) + localBlockPositionX * dxSimulation * nxBlockSimulation;
            float localOriginY =
                scenario->getBoundaryPos(BND_BOTTOM) + localBlockPositionY * dySimulation * nyBlockSimulation;
            std::string outputFileName = generateBaseFileName(outputTeamName, localBlockPositionX, localBlockPositionY);
            std::string backupFileName = generateBaseFileName(backupTeamName, localBlockPositionX, localBlockPositionY);

            backupMetadataNames.push_back(backupFileName + "_metadata");

            simulationBlocks.push_back(
                std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>(new SWE_DimensionalSplittingMPIOverdecomp(
                    nxLocal, nyLocal, dxSimulation, dySimulation, localOriginX, localOriginY, 0, outputFileName,
                    backupFileName, true, false)));
        }

        for (int i = startPoint; i < startPoint + blocksPerRank; i++) {
            auto myRank = i;
            int localBlockPositionX = myRank / blockCountY;
            int localBlockPositionY = myRank % blockCountY;
            std::array<int, 4> myNeighbours =
                getNeighbours(localBlockPositionX, localBlockPositionY, blockCountX, blockCountY, myRank);

            int refinedNeighbours[4];
            int realNeighbours[4];
            std::array<std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>, 4> neighbourBlocks;
            std::array<BoundaryType, 4> boundaries;

            for (int j = 0; j < 4; j++) {
                if (myNeighbours[j] >= startPoint && myNeighbours[j] < (startPoint + blocksPerRank)) {
                    refinedNeighbours[j] = -2;
                    realNeighbours[j] = myNeighbours[j];
                    neighbourBlocks[j] = simulationBlocks[myNeighbours[j] - startPoint];
                    boundaries[j] = CONNECT_WITHIN_RANK;
                } else if (myNeighbours[j] == -1) {
                    boundaries[j] = scenario->getBoundaryType((Boundary)j);
                    refinedNeighbours[j] = -1;
                    realNeighbours[j] = -1;
                } else {
                    realNeighbours[j] = myNeighbours[j];
                    refinedNeighbours[j] = myNeighbours[j] / blocksPerRank;
                    boundaries[j] = CONNECT;
                }
            }

            simulationBlocks[i - startPoint]->initScenario(*scenario, boundaries.data());
            simulationBlocks[i - startPoint]->connectNeighbourLocalities(refinedNeighbours);
            simulationBlocks[i - startPoint]->connectNeighbours(realNeighbours);
            simulationBlocks[i - startPoint]->connectLocalNeighbours(neighbourBlocks);
            simulationBlocks[i - startPoint]->setRank(myRank);
            simulationBlocks[i - startPoint]->setDuration(simulationDuration);
            simulationBlocks[i - startPoint]->writer->initMetadataFile(
                backupMetadataNames[i - startPoint], simulationDuration, ranksPerTeam,
                simulationBlocks[i - startPoint]->nx - 2, simulationBlocks[i - startPoint]->ny - 2, 0,
                std::vector<BoundaryType>(boundaries.begin(), boundaries.end()),
                {simulationBlocks[i - startPoint]->originX,
                 simulationBlocks[i - startPoint]->originX +
                     (simulationBlocks[i - startPoint]->nx - 2) * simulationBlocks[i - startPoint]->dx,
                 simulationBlocks[i - startPoint]->originY,
                 simulationBlocks[i - startPoint]->originY +
                     (simulationBlocks[i - startPoint]->ny - 2) * simulationBlocks[i - startPoint]->dy});
        }
    } else {
        std::vector<SWE_Scenario*> scenarios{};
        for (int i = startPoint; i < startPoint + blocksPerRank; i++) {
            auto myBlock = i;
            int localBlockPositionX = myBlock / blockCountY;
            int localBlockPositionY = myBlock % blockCountY;

            io::Reader reader{restartNameInput, outputTeamName,      myRank,
                              ranksPerTeam,     localBlockPositionX, localBlockPositionY};

            int nxLocal = reader.getGridSizeX();
            int nyLocal = reader.getGridSizeY();
            scenario = reader.getScenario();
            simulationDuration = scenario->endSimulation();
            scenarios.push_back(scenario);
            simulationStart = reader.getCurrentTime();

            int widthScenario = scenario->getBoundaryPos(BND_RIGHT) - scenario->getBoundaryPos(BND_LEFT);
            int heightScenario = scenario->getBoundaryPos(BND_TOP) - scenario->getBoundaryPos(BND_BOTTOM);

            float dxSimulation = (float)widthScenario / nxRequested;
            float dySimulation = (float)heightScenario / nyRequested;

            // Compute the origin of the local simulation block w.r.t. the
            // original scenario domain.
            float localOriginX = scenario->getBoundaryPos(BND_LEFT);
            float localOriginY = scenario->getBoundaryPos(BND_BOTTOM);
            std::string outputFileName = generateBaseFileName(outputTeamName, localBlockPositionX, localBlockPositionY);
            std::string backupFileName = generateBaseFileName(backupTeamName, localBlockPositionX, localBlockPositionY);

            backupMetadataNames.push_back(backupFileName + "_metadata");

            simulationBlocks.push_back(
                std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>(new SWE_DimensionalSplittingMPIOverdecomp(
                    nxLocal, nyLocal, dxSimulation, dySimulation, localOriginX, localOriginY, 0, outputFileName,
                    backupFileName, true, true)));
        }

        for (int i = startPoint; i < startPoint + blocksPerRank; i++) {
            auto myBlock = i;
            int localBlockPositionX = myBlock / blockCountY;
            int localBlockPositionY = myBlock % blockCountY;
            std::array<int, 4> myNeighbours =
                getNeighbours(localBlockPositionX, localBlockPositionY, blockCountX, blockCountY, myBlock);

            scenario = scenarios[myBlock - startPoint];

            int refinedNeighbours[4];
            int realNeighbours[4];
            std::array<std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>, 4> neighbourBlocks;
            std::array<BoundaryType, 4> boundaries;

            for (int j = 0; j < 4; j++) {
                if (myNeighbours[j] >= startPoint && myNeighbours[j] < (startPoint + blocksPerRank)) {
                    refinedNeighbours[j] = -2;
                    realNeighbours[j] = myNeighbours[j];
                    neighbourBlocks[j] = simulationBlocks[myNeighbours[j] - startPoint];
                    boundaries[j] = CONNECT_WITHIN_RANK;
                } else if (myNeighbours[j] == -1) {
                    boundaries[j] = scenario->getBoundaryType((Boundary)j);
                    refinedNeighbours[j] = -1;
                    realNeighbours[j] = -1;
                } else {
                    realNeighbours[j] = myNeighbours[j];
                    refinedNeighbours[j] = myNeighbours[j] / blocksPerRank;
                    boundaries[j] = CONNECT;
                }
            }

            simulationBlocks[i - startPoint]->initScenario(*scenario, boundaries.data());
            simulationBlocks[i - startPoint]->connectNeighbourLocalities(refinedNeighbours);
            simulationBlocks[i - startPoint]->connectNeighbours(realNeighbours);
            simulationBlocks[i - startPoint]->connectLocalNeighbours(neighbourBlocks);
            simulationBlocks[i - startPoint]->setRank(myBlock);
            simulationBlocks[i - startPoint]->setDuration(simulationDuration);
            simulationBlocks[i - startPoint]->writer->initMetadataFile(
                backupMetadataNames[i - startPoint], simulationDuration, ranksPerTeam,
                simulationBlocks[i - startPoint]->nx - 2, simulationBlocks[i - startPoint]->ny - 2, 0,
                std::vector<BoundaryType>(boundaries.begin(), boundaries.end()),
                {simulationBlocks[i - startPoint]->originX,
                 simulationBlocks[i - startPoint]->originX +
                     (simulationBlocks[i - startPoint]->nx - 2) * simulationBlocks[i - startPoint]->dx,
                 simulationBlocks[i - startPoint]->originY,
                 simulationBlocks[i - startPoint]->originY +
                     (simulationBlocks[i - startPoint]->ny - 2) * simulationBlocks[i - startPoint]->dy});
        }
    }
    for (auto& block : simulationBlocks) block->sendBathymetry();
    for (auto& block : simulationBlocks) block->recvBathymetry();

    std::vector<float> timesteps;
    // Simulated time
    t = simulationStart;

    float timestep;

    int teamsCount = TMPI_GetInterTeamCommSize();
    std::vector<size_t> myBlockOrder{};
    for (int i = myTeam; i < blocksPerRank; i += teamsCount) {
        myBlockOrder.push_back(i);
    }
    for (int i = 0; i < blocksPerRank; i++) {
        if ((i - myTeam) % teamsCount != 0) {
            myBlockOrder.push_back(i);
        }
    }

    // Write initial timestep
    for (int i = 0; i < simulationBlocks.size(); i++) {
        // TODO Uncoment once we test with output
        // simulationBlocks[i]->writeTimestep(simulationStart);
    }

    // simulate until end of simulation
    while (t <= simulationDuration) {
        double startTime{MPI_Wtime()};
        double elapsedTime{0.0};
        // simulate until the checkpoint is reached
        while (elapsedTime < checkpointInterval) {
            if (t >= simulationDuration) {
                break;
            }

            // exchange boundaries between blocks
            for (int i = 0; i < simulationBlocks.size(); i++) {
                simulationBlocks[i]->setGhostLayer();
            }
            for (int i = 0; i < simulationBlocks.size(); i++) {
                simulationBlocks[i]->receiveGhostLayer();
            }

            for (int i = 0; i < simulationBlocks.size(); i++) {
                simulationBlocks[i]->computeNumericalFluxes();
            }

            // determine max possible timestep
            timesteps.clear();
            for (auto& block : simulationBlocks) timesteps.push_back(block->maxTimestep);
            float minTimestep = *std::min_element(timesteps.begin(), timesteps.end());
            MPI_Allreduce(&minTimestep, &timestep, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
            for (auto& block : simulationBlocks) block->maxTimestep = timestep;

            for (int i = 0; i < simulationBlocks.size(); i++) {
                simulationBlocks[i]->updateUnknowns(timestep);
            }
            t += timestep;

            elapsedTime = MPI_Wtime() - startTime;
            MPI_Bcast(&elapsedTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            if (elapsedTime >= 1.0f && TMPI_GetTeamNumber() == 1 && restartNameInput == "") {
                std::abort();
            }
        }

        // check for failures in teaMPI world comm
        std::cout << "HEARTBEAT: Rank " << myRank << "; Team " << myTeam << '\n';
        MPI_Allreduce(nullptr, nullptr, 0, MPI_INT, MPI_MIN, MPI_COMM_SELF);
        // printf("Returned from heartbeat rank: %d, team %d\n",
        // l_mpiRank, l_teamNumber); Only write timestep when simluation
        // has finished
        if (t >= simulationDuration) {
            std::printf("Rank %i of Team %i writing final checkpoint\n", myRank, myTeam);
            for (int i = 0; i < simulationBlocks.size(); i++) {
                simulationBlocks[i]->createCheckpoint(t, backupMetadataNames[i], 0);
            }
        }
    }

    for (auto& block : simulationBlocks) {
        block->freeMpiType();
    }

    MPI_Finalize();
}
