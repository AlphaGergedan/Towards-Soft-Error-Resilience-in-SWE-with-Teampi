/**
 * @file src/tolerance/swe_softRes_and_hardRes_wTaskSharing.cpp
 *
 * @brief soft and hard error resilience with task sharing
 *
 * Checks for admissability of the computations (also see validateAdmissability
 * in src/blocks/DimSplitMPIOverdecomp.cpp) and only share the results if they
 * are admissable. If they are not, processes of this team and the team above
 * ((myTeam+1)%numberOfTeams) goes into a recovery mode, and demand checkpoint
 * from the team (myTeam+1)%numberOfTeams. After the checkpoints written, the
 * team with the soft error can load the checkpoint and continue its computations
 * in an admissable state. This should (hopefully) not break hard error resilience
 * using warm spares TODO
 *
 *
 * Here is a short pseudo-code for the computation loop:
 *
 *      int nextTeam = (myTeam + 1) % numberOfTeams;
 *      while ( t < simulationDuration ) {
 *          // start heartbeat
 *          while ( timeSinceLastHeartbeat < heartbeatInterval && t < simulationDuration ) {
 *              [...]
 *              for (blocks) {
 *                  if (isPrimary) {
 *                      // report admissability, load checkpoint if not admissable
 *                      // share results
 *                  }
 *                  else {
 *                      // receive report, create checkpoint and wait if replica reports SDC
 *                      // receive secondary blocks
 *                  }
 *              }
 *          }
 *          // end heartbeat
 *      }
 *
 * TODO adding more admissability criteria
 * TODO warm spares integration for also activating hard error resilience
 */


#include <mpi.h>
#include <string>
#include <teaMPI.h>
#include <unistd.h>

#include <algorithm>
#include <climits>
#include <csetjmp>
#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>
#include <thread>

#include "blocks/DimSplitMPIOverdecomp.hpp"
#include "io/Reader.hpp"
#include "io/Writer.hpp"
#include "scenarios/LoadNetCDFScenario.hpp"
#include "scenarios/simple_scenarios.hpp"
#include "tools/Args.hpp"
#include "tools/ftLogger.hpp"
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

float t(0.0F);
std::vector<std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>> simulationBlocks{};
SWE_Scenario* scenario{nullptr};
bool hasRecovered{false};
bool recoveredFromSDC{false};

void createCheckpointCallback(std::vector<int> failedTeams)
{
    int myRankInTeam;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    const int myTeam = TMPI_GetTeamNumber();
    std::printf("Rank %i of Team %i writing checkpoint for restoration.\n", myRankInTeam, myTeam);

    // write a checkpoint for every block
    for (int i = 0; i < simulationBlocks.size(); i++)
    {
        simulationBlocks[i]->createCheckpoint(t, backupMetadataNames[i], 0);
    }

    int send = 1;
    for (const auto& failedTeam : failedTeams)
    {
        MPI_Send(&send, 1, MPI_INT, TMPI_TeamToWorldRank(myRankInTeam, failedTeam), 0, TMPI_GetWorldComm());
    }
}

void loadCheckpointCallback(int reloadTeam)
{
    int myRankInTeam;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    const int myTeam = TMPI_GetTeamNumber();
    std::printf("Rank %i of Team %i loading checkpoint from Team %i.\n", myRankInTeam, myTeam, reloadTeam);
    MPI_Status status;
    int recv_buf;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    MPI_Recv(&recv_buf, 1, MPI_INT, TMPI_TeamToWorldRank(myRankInTeam, reloadTeam), 0, TMPI_GetWorldComm(), &status);
    simulationBlocks.clear();
    delete scenario;
    scenario = nullptr;
    restartNameInput = backupNameInput + "_t" + std::to_string(reloadTeam);
    hasRecovered = true;
    longjmp(jumpBuffer, 1);
}

/* TODO currently in development ... only difference that hasRecovered flag is not set !!! */
void loadCheckpointCallback_AFTER_SDC(int reloadTeam)
{
    int myRankInTeam;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    const int myTeam = TMPI_GetTeamNumber();
    std::printf("Rank %i of Team %i loading checkpoint from Team %i.\n", myRankInTeam, myTeam, reloadTeam);
    MPI_Status status;
    int recv_buf;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    MPI_Recv(&recv_buf, 1, MPI_INT, TMPI_TeamToWorldRank(myRankInTeam, reloadTeam), 0, TMPI_GetWorldComm(), &status);
    simulationBlocks.clear();
    delete scenario;
    scenario = nullptr;
    restartNameInput = backupNameInput + "_t" + std::to_string(reloadTeam);
    //hasRecovered = true;
    recoveredFromSDC = true;
    longjmp(jumpBuffer, 1);
}


std::array<int, 4> getNeighbours(int localBlockPositionX,
                                 int localBlockPositionY,
                                 int blockCountX,
                                 int blockCountY,
                                 int myRank)
{
    std::array<int, 4> myNeighbours;
    myNeighbours[BND_LEFT] = (localBlockPositionX > 0) ? myRank - blockCountY : -1;
    myNeighbours[BND_RIGHT] = (localBlockPositionX < blockCountX - 1) ? myRank + blockCountY : -1;
    myNeighbours[BND_BOTTOM] = (localBlockPositionY > 0) ? myRank - 1 : -1;
    myNeighbours[BND_TOP] = (localBlockPositionY < blockCountY - 1) ? myRank + 1 : -1;
    return myNeighbours;
}

int main(int argc, char** argv)
{
    // Define command line arguments
    tools::Args args;

    args.addOption("simulation-duration", 't', "Time in seconds to simulate");
    args.addOption("heartbeat-interval", 'i', "Wall-clock time in seconds to wait between heartbeats", args.Required, true);
    args.addOption("resolution-x", 'x', "Number of simulated cells in x-direction");
    args.addOption("resolution-y", 'y', "Number of simulated cells in y-direction");
    args.addOption("decomp-factor",
                   'd',
                   "Split each rank into \"TEAMS\" * \"decomp-factor\" blocks",
                   tools::Args::Required,
                   false);
    args.addOption("output-basepath", 'o', "Output base file name");
    args.addOption("backup-basepath", 'b', "Output base file name");
    args.addOption("restart-basepath", 'r', "Restart base file name", tools::Args::Required, false);
    args.addOption("write-output", 'w', "Write output using netcdf writer to the specified output file", args.No, false);
    // TODO make bitflip at random time ?
    args.addOption("inject-bitflip", 'f', "Injects a bit-flip to the first rank right after the simulation time reaches the given time", args.Required, false);
    args.addOption("kill-rank", 'k', "Kills the rank 0 of team 0 at the specified simulation time", args.Required, false);
    args.addOption("verbose", 'v', "Let the simulation produce more output, default: No", args.No, false);

    // Parse command line arguments
    tools::Args::Result ret = args.parse(argc, argv);
    switch (ret)
    {
    case tools::Args::Error:
        return 1;
    case tools::Args::Help:
        return 0;
    case tools::Args::Success:
        break;
    }

    /* Arguments to read */
    float simulationDuration;
    clock_t heartbeatInterval;

    int nxRequested;
    int nyRequested;

    // TODO this must be 1 for now, fix this in the future
    unsigned int decompFactor = 1;
    bool writeOutput = false;
    double bitflip_at = -1.f;
    double kill_rank = -1.f;
    bool verbose = false;

    /* Read in command line arguments */
    simulationDuration = args.getArgument<float>("simulation-duration");
    heartbeatInterval = args.getArgument<clock_t>("heartbeat-interval");
    nxRequested = args.getArgument<int>("resolution-x");
    nyRequested = args.getArgument<int>("resolution-y");
    if (args.isSet("decomp-factor")) {
        decompFactor = args.getArgument<unsigned int>("decomp-factor");
        if (decompFactor == 0)
        {
            decompFactor = 1;
        }
    }
    outputNameInput = args.getArgument<std::string>("output-basepath");
    backupNameInput = args.getArgument<std::string>("backup-basepath");
    if (args.isSet("restart-basepath")) {
        restartNameInput = args.getArgument<std::string>("restart-basepath");
    }
    writeOutput = args.isSet("write-output");
    if (args.isSet("inject-bitflip")) bitflip_at = args.getArgument<double>("inject-bitflip");
    if (args.isSet("kill-rank")) kill_rank = args.getArgument<double>("kill-rank");
    verbose = args.isSet("verbose");

    // init teaMPI
    std::function<void(std::vector<int>)> create(createCheckpointCallback);
    std::function<void(int)> load(loadCheckpointCallback);
    TMPI_SetCreateCheckpointCallback(&create);
    TMPI_SetLoadCheckpointCallback(&load);
    TMPI_SetErrorHandlingStrategy(TMPI_WarmSpareErrorHandler);

    // init MPI
    int myRankInTeam;
    if (setjmp(jumpBuffer) == 0)
    {
        // MPI_Init_thread(&argc, &argv, requested, &provided);
        MPI_Init(&argc, &argv);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &ranksPerTeam);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    const int myTeam{TMPI_GetTeamNumber()};
    int numTeams = TMPI_GetInterTeamCommSize();
    unsigned int blocksPerRank = numTeams * decompFactor;
    outputTeamName = outputNameInput + "_t" + std::to_string(myTeam);
    backupTeamName = backupNameInput + "_t" + std::to_string(myTeam);
    /* in recovery, the team name is already added to restart name in load function */

    /* Print status */
    tools::FtLogger::ft_logger.ft_print_spawnStatus(myRankInTeam, myTeam);

    int totalBlocks = blocksPerRank * ranksPerTeam;

    /* number of SWE-Blocks in x- and y-direction */
    int blockCountY = std::sqrt(totalBlocks);
    while (totalBlocks % blockCountY != 0) blockCountY--;
    int blockCountX = totalBlocks / blockCountY;

    int startPoint = myRankInTeam * blocksPerRank;

    float simulationStart{0.0f};

    /* if not loading from a checkpoint */
    if (restartNameInput == "")
    {
        scenario = new SWE_RadialBathymetryDamBreakScenario{};
        int widthScenario = scenario->getBoundaryPos(BND_RIGHT) - scenario->getBoundaryPos(BND_LEFT);
        int heightScenario = scenario->getBoundaryPos(BND_TOP) - scenario->getBoundaryPos(BND_BOTTOM);

        float dxSimulation = (float)widthScenario / nxRequested;
        float dySimulation = (float)heightScenario / nyRequested;

        for (int currentBlockNr = startPoint; currentBlockNr < startPoint + blocksPerRank; currentBlockNr++)
        {
            int localBlockPositionX = currentBlockNr / blockCountY;
            int localBlockPositionY = currentBlockNr % blockCountY;

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

            std::string outputTeamPosName = genTeamPosName(outputTeamName, localBlockPositionX, localBlockPositionY);
            std::string backupTeamPosName = genTeamPosName(backupTeamName, localBlockPositionX, localBlockPositionY);

            backupMetadataNames.push_back(backupTeamPosName + "_metadata");

            simulationBlocks.push_back(std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>(
                new SWE_DimensionalSplittingMPIOverdecomp(nxLocal,
                                                          nyLocal,
                                                          dxSimulation,
                                                          dySimulation,
                                                          localOriginX,
                                                          localOriginY,
                                                          0,
                                                          outputTeamPosName,
                                                          backupTeamPosName,
                                                          true,
                                                          false)));
        }

        for (int currentBlockNr = startPoint; currentBlockNr < startPoint + blocksPerRank; currentBlockNr++)
        {
            int localBlockPositionX = currentBlockNr / blockCountY;
            int localBlockPositionY = currentBlockNr % blockCountY;
            std::array<int, 4> myNeighbours =
                getNeighbours(localBlockPositionX, localBlockPositionY, blockCountX, blockCountY, currentBlockNr);

            int refinedNeighbours[4];
            int realNeighbours[4];
            std::array<std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>, 4> neighbourBlocks;
            std::array<BoundaryType, 4> boundaries;

            for (int j = 0; j < 4; j++)
            {
                if (myNeighbours[j] >= startPoint && myNeighbours[j] < (startPoint + blocksPerRank))
                {
                    refinedNeighbours[j] = -2;
                    realNeighbours[j] = myNeighbours[j];
                    neighbourBlocks[j] = simulationBlocks[myNeighbours[j] - startPoint];
                    boundaries[j] = CONNECT_WITHIN_RANK;
                }
                else if (myNeighbours[j] == -1)
                {
                    boundaries[j] = scenario->getBoundaryType((Boundary)j);
                    refinedNeighbours[j] = -1;
                    realNeighbours[j] = -1;
                }
                else
                {
                    realNeighbours[j] = myNeighbours[j];
                    refinedNeighbours[j] = myNeighbours[j] / blocksPerRank;
                    boundaries[j] = CONNECT;
                }
            }

            simulationBlocks[currentBlockNr - startPoint]->initScenario(*scenario, boundaries.data());
            simulationBlocks[currentBlockNr - startPoint]->connectNeighbourLocalities(refinedNeighbours);
            simulationBlocks[currentBlockNr - startPoint]->connectNeighbours(realNeighbours);
            simulationBlocks[currentBlockNr - startPoint]->connectLocalNeighbours(neighbourBlocks);
            simulationBlocks[currentBlockNr - startPoint]->setRank(currentBlockNr);
            simulationBlocks[currentBlockNr - startPoint]->setDuration(simulationDuration);
            simulationBlocks[currentBlockNr - startPoint]->writer->initMetadataFile(
                backupMetadataNames[currentBlockNr - startPoint],
                simulationDuration,
                ranksPerTeam,
                simulationBlocks[currentBlockNr - startPoint]->nx,
                simulationBlocks[currentBlockNr - startPoint]->ny,
                0,
                std::vector<BoundaryType>(boundaries.begin(), boundaries.end()),
                {simulationBlocks[currentBlockNr - startPoint]->originX,
                 simulationBlocks[currentBlockNr - startPoint]->originX +
                     (simulationBlocks[currentBlockNr - startPoint]->nx) *
                         simulationBlocks[currentBlockNr - startPoint]->dx,
                 simulationBlocks[currentBlockNr - startPoint]->originY,
                 simulationBlocks[currentBlockNr - startPoint]->originY +
                     (simulationBlocks[currentBlockNr - startPoint]->ny) *
                         simulationBlocks[currentBlockNr - startPoint]->dy});
        }
    }
    else /* loading from a checkpoint */
    {
        std::vector<SWE_Scenario*> scenarios{};
        for (int currentBlockNr = startPoint; currentBlockNr < startPoint + blocksPerRank; currentBlockNr++)
        {
            int localBlockPositionX = currentBlockNr / blockCountY;
            int localBlockPositionY = currentBlockNr % blockCountY;

            io::Reader reader{restartNameInput,
                              outputTeamName,
                              myRankInTeam,
                              ranksPerTeam,
                              localBlockPositionX,
                              localBlockPositionY};

            int nxLocal = reader.getGridSizeX();
            int nyLocal = reader.getGridSizeY();
            scenario = reader.getScenario();
            simulationDuration = scenario->endSimulation();
            scenarios.push_back(scenario);
            simulationStart = reader.getCurrentTime();

            int widthScenario = scenario->getBoundaryPos(BND_RIGHT) - scenario->getBoundaryPos(BND_LEFT);
            int heightScenario = scenario->getBoundaryPos(BND_TOP) - scenario->getBoundaryPos(BND_BOTTOM);

            float dxSimulation = (float)widthScenario / nxLocal;
            float dySimulation = (float)heightScenario / nyLocal;

            // Compute the origin of the local simulation block w.r.t. the
            // original scenario domain.
            float localOriginX = scenario->getBoundaryPos(BND_LEFT);
            float localOriginY = scenario->getBoundaryPos(BND_BOTTOM);
            std::string outputTeamPosName = genTeamPosName(outputTeamName, localBlockPositionX, localBlockPositionY);
            std::string backupTeamPosName = genTeamPosName(backupTeamName, localBlockPositionX, localBlockPositionY);

            backupMetadataNames.push_back(backupTeamPosName + "_meta");
            simulationBlocks.push_back(std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>(
                new SWE_DimensionalSplittingMPIOverdecomp(nxLocal,
                                                          nyLocal,
                                                          dxSimulation,
                                                          dySimulation,
                                                          localOriginX,
                                                          localOriginY,
                                                          0,
                                                          outputTeamPosName,
                                                          backupTeamPosName,
                                                          true,
                                                          true)));
        }

        for (int currentBlockNr = startPoint; currentBlockNr < startPoint + blocksPerRank; currentBlockNr++)
        {
            int localBlockPositionX = currentBlockNr / blockCountY;
            int localBlockPositionY = currentBlockNr % blockCountY;
            std::array<int, 4> myNeighbours =
                getNeighbours(localBlockPositionX, localBlockPositionY, blockCountX, blockCountY, currentBlockNr);

            scenario = scenarios[currentBlockNr - startPoint];

            int refinedNeighbours[4];
            int realNeighbours[4];
            std::array<std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>, 4> neighbourBlocks;
            std::array<BoundaryType, 4> boundaries;

            for (int j = 0; j < 4; j++)
            {
                if (myNeighbours[j] >= startPoint && myNeighbours[j] < (startPoint + blocksPerRank))
                {
                    refinedNeighbours[j] = -2;
                    realNeighbours[j] = myNeighbours[j];
                    neighbourBlocks[j] = simulationBlocks[myNeighbours[j] - startPoint];
                    boundaries[j] = CONNECT_WITHIN_RANK;
                }
                else if (myNeighbours[j] == -1)
                {
                    boundaries[j] = scenario->getBoundaryType((Boundary)j);
                    refinedNeighbours[j] = -1;
                    realNeighbours[j] = -1;
                }
                else
                {
                    realNeighbours[j] = myNeighbours[j];
                    refinedNeighbours[j] = myNeighbours[j] / blocksPerRank;
                    boundaries[j] = CONNECT;
                }
            }

            simulationBlocks[currentBlockNr - startPoint]->initScenario(*scenario, boundaries.data());
            simulationBlocks[currentBlockNr - startPoint]->connectNeighbourLocalities(refinedNeighbours);
            simulationBlocks[currentBlockNr - startPoint]->connectNeighbours(realNeighbours);
            simulationBlocks[currentBlockNr - startPoint]->connectLocalNeighbours(neighbourBlocks);
            simulationBlocks[currentBlockNr - startPoint]->setRank(currentBlockNr);
            simulationBlocks[currentBlockNr - startPoint]->setDuration(simulationDuration);
            simulationBlocks[currentBlockNr - startPoint]->writer->initMetadataFile(
                backupMetadataNames[currentBlockNr - startPoint],
                simulationDuration,
                ranksPerTeam,
                simulationBlocks[currentBlockNr - startPoint]->nx,
                simulationBlocks[currentBlockNr - startPoint]->ny,
                0,
                std::vector<BoundaryType>(boundaries.begin(), boundaries.end()),
                {simulationBlocks[currentBlockNr - startPoint]->originX,
                 simulationBlocks[currentBlockNr - startPoint]->originX +
                     (simulationBlocks[currentBlockNr - startPoint]->nx) *
                         simulationBlocks[currentBlockNr - startPoint]->dx,
                 simulationBlocks[currentBlockNr - startPoint]->originY,
                 simulationBlocks[currentBlockNr - startPoint]->originY +
                     (simulationBlocks[currentBlockNr - startPoint]->ny) *
                         simulationBlocks[currentBlockNr - startPoint]->dy});
            delete scenario;
        }
    }
    for (auto& block : simulationBlocks) block->sendBathymetry();
    for (auto& block : simulationBlocks) block->recvBathymetry();

    std::vector<float> timesteps;
    // Simulated time
    t = simulationStart;

    float timestep;

    std::vector<int> myBlockOrder{};
    // Add my primary blocks first to block order
    for (int i = myTeam; i < blocksPerRank; i += numTeams) { myBlockOrder.push_back(i); }
    // Add all other secondary blocks to block order
    for (int i = 0; i < blocksPerRank; i++)
    {
        if (i % numTeams != myTeam)
        {
            myBlockOrder.push_back(i);
        }
    }

    /* Write zero timestep if not restarting */
    if (simulationStart == 0.f) {
        for (auto& block : simulationBlocks) { block->writeTimestep(0.f); }
    }

    // simulate until end of simulation
    while (t < simulationDuration) {
        // Start Hearbeat
        MPI_Sendrecv(MPI_IN_PLACE,
                     0,
                     MPI_BYTE,
                     MPI_PROC_NULL,
                     1,
                     MPI_IN_PLACE,
                     0,
                     MPI_BYTE,
                     MPI_PROC_NULL,
                     0,
                     MPI_COMM_SELF,
                     MPI_STATUS_IGNORE);

        double timeOfLastHeartbeat{MPI_Wtime()};
        double timeSinceLastHeartbeat{0.0};
        if (verbose) tools::FtLogger::ft_logger.ft_print_HBstart(myTeam, myRankInTeam, timeOfLastHeartbeat, t);

        // simulate until the checkpoint is reached
        while (timeSinceLastHeartbeat < heartbeatInterval && t < simulationDuration)
        {
            if (verbose) tools::FtLogger::ft_logger.ft_print_loop(myTeam, myRankInTeam, timeSinceLastHeartbeat,
                                                                  heartbeatInterval, simulationDuration, t);

            // exchange boundaries between blocks
            for (auto& currentBlock : simulationBlocks) { currentBlock->setGhostLayer(); }
            for (auto& currentBlock : simulationBlocks) { currentBlock->receiveGhostLayer(); }
            const MPI_Comm interTeamComm{TMPI_GetInterTeamComm()};
            if (!hasRecovered && !recoveredFromSDC)
            {
                // Avoid overwriting an old send buffer before everyone reaches this point
                MPI_Barrier(interTeamComm);
            }
            std::vector<MPI_Request> send_reqs(11 * numTeams, MPI_REQUEST_NULL);

            for (int i{0}; i < blocksPerRank; i++)
            {
                // Get a reference to the current block
                const int& currentBlockNr{myBlockOrder[i]};
                auto& currentBlock = *simulationBlocks[currentBlockNr];

                // Size of the update fields incl. ghost layer
                const int fieldSizeX{(currentBlock.nx + 2) * (currentBlock.ny + 2)};
                const int fieldSizeY{(currentBlock.nx + 1) * (currentBlock.ny + 2)};

                // The first [decompFactor] blocks are the ones we always compute ourselves
                if (i < decompFactor) {

                    currentBlock.computeNumericalFluxes();

                    /* inject bitflip if desired */
                    if (bitflip_at >= 0  && t > bitflip_at && myTeam == 0 && myRankInTeam == 0) {
                        currentBlock.injectRandomBitflip();
                        //currentBlock.injectRandomBitflip_intoData();
                        //currentBlock.injectRandomBitflip_intoUpdates();

                        /* prevent any other bitflip */
                        bitflip_at = -1.f;
                    }

                    /* check for soft errors */
                    int destTeam = (myTeam + 1) % numTeams;
                    MPI_Request reqSend;
                    int sendBuf = 0; // TODO char could also be used
                    int admissable = currentBlock.validateAdmissability(t);

//--------------------TODO maybe move the block to a header, we have the same block below again.. -----------------------

                    /* handle soft errors */
                    // soft error only in updates, recomputing..
                    if (admissable == 1) {
                        currentBlock.computeNumericalFluxes();

                        /* check if the silent error is still present */
                        if (currentBlock.validateAdmissability(t) != 0) {
                            std::cout << "-- TEAM " << myTeam << ", Rank " << myRankInTeam << " Warning : SDC detected and cannot be fixed..\n"
                                      << "             Warning the rest of my team and the reload team " << destTeam
                                      << std::endl;
                            sendBuf = -1;
                        } else {
                            std::cout << "-- TEAM " << myTeam << ", Rank " << myRankInTeam << " Warning : SDC detected but it is fixed.."
                                      << std::endl;
                            sendBuf = 1;
                        }
                    }
                    // soft error in data, turn on the recovery mode
                    else if (admissable == 2) {
                        std::cout << "-- TEAM " << myTeam << ", Rank " << myRankInTeam << " Warning : SDC detected but cannot be fixed..\n"
                                  << "             Warning the rest of my team and the reload team " << destTeam
                                  << std::endl;
                        sendBuf = -1;
                    }
                    else {
                        //std::cout << "-- TEAM " << myTeam << ", Rank " << myRankInTeam << " : No SDC is detected"
                                  //<< std::endl;
                        if (verbose) tools::FtLogger::ft_logger.ft_SDC_notDetected(myTeam, myRankInTeam);
                        sendBuf = 1;
                    }

                    /* report to my replica TODO after allreduce, only rank 0 can send this */
                    MPI_Send(&sendBuf, 1, MPI_INT, destTeam, 30, interTeamComm);

                    int teamCheck = 0;

                    /* report to my team */
                    MPI_Allreduce(&sendBuf, &teamCheck, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

                    //for (int i = 0; i < ranksPerTeam; i++) {
                        //if (i != myRankInTeam) {
                            // TODO try to send with Isend here
                            //MPI_Isend(&sendBuf, 1, MPI_INT, i, 20, MPI_COMM_WORLD, &reqSend);
                            //MPI_Request_free(&reqSend);
                            //MPI_Send(&sendBuf, 1, MPI_INT, i, 20, MPI_COMM_WORLD);
                        //}
                    //}

                    // TODO debugging
                    if(sendBuf == 0) {
                        std::cout << "UNKNOWN ERROR ??" << std::endl;
                        MPI_Abort(MPI_COMM_WORLD, MPI_ERR_UNKNOWN);
                    }

                        /* TODO we dont need this ?
                        wait for the checkpoint to load from destTeam
                        std::cout << "\t***Waiting for replica " << destTeam << " to write checkpoint"
                                  << std::endl;
                        //int temp2;
                        // tag 12 for emergency soft error recv
                        //MPI_Recv(&temp2, 1, MPI_INT, destTeam, 40, interTeamComm, MPI_STATUS_IGNORE);
                        */

                    /* load if soft error is present TODO we don't have to wait here? */
                    if(teamCheck == -1) {
                        /* TODO wait for destTeam to write checkpoint, possibly we don't need to wait (there is a MPI_Recv in loadCheckpoint call) */
                        std::cout << "\t*** team " << myTeam << ", rank " << myRankInTeam << "; received from my team, from rank " << i << " the buffer: "
                                  << teamCheck << " : WARNING SOFT ERROR IS PRESENT IN MY TEAM !! loading checkpoint..." << std::endl;
                        loadCheckpointCallback_AFTER_SDC(destTeam); // TODO change this name
                    }

//--------------------TODO maybe move the block to a header, we have the same block below again.. -----------------------


                    // TODO TESTING REPORT CHECK LOCATION
        //-----------------------------------------------------------------
                    // check for soft errors in your own team TODO remove this after debugging
                    for (int i = 0; i < ranksPerTeam; i++) {
                        if (i != myRankInTeam) {
                            //int teamCheck;
                            //MPI_Recv(&teamCheck, 1, MPI_INT, i, 20, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                            // TODO debugging
                            //std::cout << "\t*** team " << myTeam << ", rank " << myRankInTeam << "; received from my team, from rank " << i << " the buffer: "
                                      //<< teamCheck << std::endl;

                            /* if soft error is present */
                            //if (teamCheck == -1) {
                                ///* TODO wait for destTeam to write checkpoint, possibly we don't need to wait (there is a MPI_Recv in loadCheckpoint call) */
                                //std::cout << "\t*** team " << myTeam << ", rank " << myRankInTeam << "; received from my team, from rank " << i << " the buffer: "
                                          //<< teamCheck << " : WARNING SOFT ERROR IS PRESENT IN MY TEAM !! loading.." << std::endl;
                                //loadCheckpointCallback((myTeam + 1) % numTeams);
                            //}
                        }
                    }

        //---------------------------------------------------------------------------------------

                    //std::cout << "Team " << myTeam << ": Block" << currentBlockNr << " calculated timestep "
                              //<< currentBlock.maxTimestep << '\n';
                    if (verbose) tools::FtLogger::ft_logger.ft_block_calculatingTask(myTeam, myRankInTeam, currentBlockNr, currentBlock.maxTimestep);

                    // Send to all other teams
                    for (int destTeam{0}; destTeam < numTeams; destTeam++)
                    {
                        if (destTeam == myTeam || hasRecovered == true)
                        {
                            // Do not send to myself
                            // Skip send if I recovered
                        }
                        else
                        {
                            //std::cout << "Team " << myTeam << ": Sending t=" << t << " to Team " << destTeam << '\n';
                            if (verbose) tools::FtLogger::ft_logger.ft_block_sending(myTeam, myRankInTeam, t, destTeam);

                            /* current timestep */
                            MPI_Isend(&t, 1, MPI_FLOAT, destTeam, 1, interTeamComm, &send_reqs[11 * destTeam]);
                            MPI_Request_free(&send_reqs[11 * destTeam]);

                            /* current block */
                            MPI_Isend(&myBlockOrder[i],
                                      1,
                                      MPI_INT,
                                      destTeam,
                                      2,
                                      interTeamComm,
                                      &send_reqs[11 * destTeam + 1]);
                            MPI_Request_free(&send_reqs[11 * destTeam + 1]);

                            /* h netupdates Left */
                            MPI_Isend(currentBlock.hNetUpdatesLeft.getRawPointer(),
                                      fieldSizeX,
                                      MPI_FLOAT,
                                      destTeam,
                                      2,
                                      interTeamComm,
                                      &send_reqs[11 * destTeam + 2]);
                            MPI_Request_free(&send_reqs[11 * destTeam + 2]);

                            /* h netupdates Right */
                            MPI_Isend(currentBlock.hNetUpdatesRight.getRawPointer(),
                                      fieldSizeX,
                                      MPI_FLOAT,
                                      destTeam,
                                      3,
                                      interTeamComm,
                                      &send_reqs[11 * destTeam + 3]);
                            MPI_Request_free(&send_reqs[11 * destTeam + 3]);

                            /* hu netupdates Left */
                            MPI_Isend(currentBlock.huNetUpdatesLeft.getRawPointer(),
                                      fieldSizeX,
                                      MPI_FLOAT,
                                      destTeam,
                                      4,
                                      interTeamComm,
                                      &send_reqs[11 * destTeam + 4]);
                            MPI_Request_free(&send_reqs[11 * destTeam + 4]);

                            /* hu netupdates Right*/
                            MPI_Isend(currentBlock.huNetUpdatesRight.getRawPointer(),
                                      fieldSizeX,
                                      MPI_FLOAT,
                                      destTeam,
                                      5,
                                      interTeamComm,
                                      &send_reqs[11 * destTeam + 5]);
                            MPI_Request_free(&send_reqs[11 * destTeam + 5]);

                            /* h netupdates below */
                            MPI_Isend(currentBlock.hNetUpdatesBelow.getRawPointer(),
                                      fieldSizeY,
                                      MPI_FLOAT,
                                      destTeam,
                                      6,
                                      interTeamComm,
                                      &send_reqs[11 * destTeam + 6]);
                            MPI_Request_free(&send_reqs[11 * destTeam + 6]);

                            /* h netupdates above */
                            MPI_Isend(currentBlock.hNetUpdatesAbove.getRawPointer(),
                                      fieldSizeY,
                                      MPI_FLOAT,
                                      destTeam,
                                      7,
                                      interTeamComm,
                                      &send_reqs[11 * destTeam + 7]);
                            MPI_Request_free(&send_reqs[11 * destTeam + 7]);

                            /* hv netupdates below */
                            MPI_Isend(currentBlock.hvNetUpdatesBelow.getRawPointer(),
                                      fieldSizeY,
                                      MPI_FLOAT,
                                      destTeam,
                                      8,
                                      interTeamComm,
                                      &send_reqs[11 * destTeam + 8]);
                            MPI_Request_free(&send_reqs[11 * destTeam + 8]);

                            /* hv netupdates above */
                            MPI_Isend(currentBlock.hvNetUpdatesAbove.getRawPointer(),
                                      fieldSizeY,
                                      MPI_FLOAT,
                                      destTeam,
                                      9,
                                      interTeamComm,
                                      &send_reqs[11 * destTeam + 9]);
                            MPI_Request_free(&send_reqs[11 * destTeam + 9]);

                            /* max time step */
                            MPI_Isend(&(currentBlock.maxTimestep),
                                      1,
                                      MPI_FLOAT,
                                      destTeam,
                                      10,
                                      interTeamComm,
                                      &send_reqs[11 * destTeam + 10]);
                            MPI_Request_free(&send_reqs[11 * destTeam + 10]);
                        }
                    } // end of for (int destTeam{0}; destTeam < numTeams; destTeam++)

                } // end of if (i < decompFactor), primary blocks
                else {


        //-----------------------------------------------------------------
                    // TODO debugging test for Soft error ping from previous team

                    /* attention, % is not the modulo operator, it is the remainder operator */
                    int prevTeam = abs((myTeam - 1) % numTeams);
                    int replicaCheck;
                    int softErr;

                    // TODO we can remove this after debugging + teams>2 problem fix
                    int counter = 0;

                    while (true) {
                        /* check the replica if it has soft error */
                        //if (counter == 1) std::cout << "\t### calling MPI_Recv " << std::endl;
                        MPI_Recv(&replicaCheck, 1, MPI_INT, prevTeam, 30, interTeamComm, MPI_STATUS_IGNORE);
                        //if (counter == 1) std::cout << "\t### MPI_Recv finished! " << std::endl;
                        //std::cout << "\t*** team " << myTeam << ", rank " << myRankInTeam << "; received from replica from team " << prevTeam << " the buffer: "
                                  //<< replicaCheck << std::endl;

                        //std::cout << "\t*** team " << myTeam << ", rank " << myRankInTeam << "; replicaCheck Allreduce is called now !" << std::endl;

                        /* inform the team about the replicas */
                        MPI_Allreduce(&replicaCheck, &softErr, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

                        // TODO debugging
                        //std::cout << "\t*** team " << myTeam << ", rank " << myRankInTeam << "; replicaCheck Allreduce result is : " << softErr << std::endl;

                        /* if no SDC is detected from replica TEAM */
                        if (softErr == 1) break;

                        std::cout << "\t+++ TEAM " << myTeam << ": Replica " << myRankInTeam << ", team " << prevTeam << " got soft error"
                                  << std::endl;
                        std::cout << "\t+++ writing checkpoint for ya!"
                                  << std::endl;

                        std::vector<int> failedTeams;
                        failedTeams.push_back(prevTeam);
                        createCheckpointCallback(failedTeams);
                        counter++;
                    }

                    // TODO TESTING REPORT CHECK LOCATION

                    int fieldSizeX = (currentBlock.nx + 2) * (currentBlock.ny + 2);
                    int fieldSizeY = (currentBlock.nx + 1) * (currentBlock.ny + 2);
                    std::vector<MPI_Request> reqs(11, MPI_REQUEST_NULL);
                    std::vector<MPI_Status> stats(11);
                    int source_rank = currentBlockNr % numTeams;
                    float recv_t;
                    int recv_blockNum;
                    if (!hasRecovered)
                    {
                        MPI_Irecv(&recv_t, 1, MPI_FLOAT, source_rank, 1, interTeamComm, &reqs[0]);
                        MPI_Irecv(&recv_blockNum, 1, MPI_INT, source_rank, 2, interTeamComm, &reqs[1]);
                        MPI_Irecv(currentBlock.hNetUpdatesLeft.getRawPointer(),
                                  fieldSizeX,
                                  MPI_FLOAT,
                                  source_rank,
                                  2,
                                  TMPI_GetInterTeamComm(),
                                  &reqs[2]);
                        MPI_Irecv(currentBlock.hNetUpdatesRight.getRawPointer(),
                                  fieldSizeX,
                                  MPI_FLOAT,
                                  source_rank,
                                  3,
                                  TMPI_GetInterTeamComm(),
                                  &reqs[3]);
                        MPI_Irecv(currentBlock.huNetUpdatesLeft.getRawPointer(),
                                  fieldSizeX,
                                  MPI_FLOAT,
                                  source_rank,
                                  4,
                                  TMPI_GetInterTeamComm(),
                                  &reqs[4]);
                        MPI_Irecv(currentBlock.huNetUpdatesRight.getRawPointer(),
                                  fieldSizeX,
                                  MPI_FLOAT,
                                  source_rank,
                                  5,
                                  TMPI_GetInterTeamComm(),
                                  &reqs[5]);
                        MPI_Irecv(currentBlock.hNetUpdatesBelow.getRawPointer(),
                                  fieldSizeY,
                                  MPI_FLOAT,
                                  source_rank,
                                  6,
                                  TMPI_GetInterTeamComm(),
                                  &reqs[6]);
                        MPI_Irecv(currentBlock.hNetUpdatesAbove.getRawPointer(),
                                  fieldSizeY,
                                  MPI_FLOAT,
                                  source_rank,
                                  7,
                                  TMPI_GetInterTeamComm(),
                                  &reqs[7]);
                        MPI_Irecv(currentBlock.hvNetUpdatesBelow.getRawPointer(),
                                  fieldSizeY,
                                  MPI_FLOAT,
                                  source_rank,
                                  8,
                                  TMPI_GetInterTeamComm(),
                                  &reqs[8]);
                        MPI_Irecv(currentBlock.hvNetUpdatesAbove.getRawPointer(),
                                  fieldSizeY,
                                  MPI_FLOAT,
                                  source_rank,
                                  9,
                                  TMPI_GetInterTeamComm(),
                                  &reqs[9]);
                        MPI_Irecv(&(currentBlock.maxTimestep),
                                  1,
                                  MPI_FLOAT,
                                  source_rank,
                                  10,
                                  TMPI_GetInterTeamComm(),
                                  &reqs[10]);
                    }

                    //std::cout << "T" << myTeam << "R" << myRankInTeam << " MPI_Waitall start" << std::endl;
                    int code = MPI_Waitall(11, reqs.data(), stats.data());
                    //std::cout << "T" << myTeam << "R" << myRankInTeam << " MPI_Waitall end ^^" << std::endl;
                    if (code != MPI_SUCCESS || hasRecovered) {
                        std::cout << "*#*#*#*#*#*# unknown error ?  " << std::endl;
                        assert(false); // TODO debugging right now, you can't enter here sorry :(
                        if (code != MPI_SUCCESS)
                        {
                            std::cout << "Team " << myTeam << ": Error in Waitall for block " << currentBlockNr << ": "
                                      << code << '\n';
                        }
                        currentBlock.computeNumericalFluxes();

                        /* TODO handle this part, we should again check SDC, or
                         *      what should we do ? */

//--------------------TODO maybe move the block to a header, we have the same block above -----------------------

                        /* check for soft errors */

//--------------------TODO maybe move the block to a header, we have the same block above -----------------------

                    }
                    else {
                        // currentBlock.computeNumericalFluxes();
                        //std::cout << "Team " << myTeam << ": Received t=" << recv_t << " from Team " << source_rank
                                  //<< '\n';
                        if (verbose) tools::FtLogger::ft_logger.ft_block_received(myTeam, myRankInTeam, recv_t, source_rank);
                    }
                } // end of else, secondary blocks
            }
            hasRecovered = false;
            recoveredFromSDC = false;

            // determine max possible timestep
            timesteps.clear();
            for (auto& block : simulationBlocks) {
                timesteps.push_back(block->maxTimestep);
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : Single Timestep = " << block->maxTimestep
                          << std::endl;
            }
            float minTimestep = *std::min_element(timesteps.begin(), timesteps.end());

            // TODO debugging ALL REDUCE TIME STEP
            //std::cout << "\t*** team " << myTeam << ", rank " << myRankInTeam << "; Timestep Allreduce is called now !" << std::endl;
            MPI_Allreduce(&minTimestep, &timestep, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);


            std::cout << "T" << myTeam << "R" << myRankInTeam
                      << " : Max Timestep = " << timestep << std::endl;

            for (auto& block : simulationBlocks) { block->maxTimestep = timestep; }
            for (auto& block : simulationBlocks) { block->updateUnknowns(timestep); }
            t += timestep;

            /* write output */
            if (writeOutput) {
                if (verbose) tools::FtLogger::ft_logger.ft_writingTimeStep(myTeam, myRankInTeam, t);
                for (auto& block : simulationBlocks) {block->writeTimestep(t);}
            }

            timeSinceLastHeartbeat = MPI_Wtime() - timeOfLastHeartbeat;
            MPI_Bcast(&timeSinceLastHeartbeat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            //MPI_Barrier(TMPI_GetInterTeamComm());
            if (kill_rank >= 0 && t >= kill_rank &&
                myTeam == 0 && myRankInTeam == 0 && restartNameInput == "") {
                std::cout << "TEAM:0 Rank:0 RETURNS ERROR CODE for hard failure simulation...."
                          << std::endl;
                return 1;
                //std::cout << "ETERNAL SLEEP TIME FOR THE TEAM:0 Rank:0 for hard failure simulation...."
                          //<< std::endl;
                //while(true) sleep(10); TODO also check if this is working
            }
        }

        // End Heartbeat
        //std::cout << "Team " << myTeam << ": HEARTBEAT! Current simulation time is " << t << '\n';
        MPI_Sendrecv(MPI_IN_PLACE,
                     0,
                     MPI_BYTE,
                     MPI_PROC_NULL,
                     -1,
                     MPI_IN_PLACE,
                     0,
                     MPI_BYTE,
                     MPI_PROC_NULL,
                     0,
                     MPI_COMM_SELF,
                     MPI_STATUS_IGNORE);
        if (verbose) tools::FtLogger::ft_logger.ft_print_HBend(myTeam, myRankInTeam, timeSinceLastHeartbeat, t);
    }

    for (auto& block : simulationBlocks) { block->freeMpiType(); }

    std::cout << "\n\n\nT" << myTeam << "R" << myRankInTeam << "calling finalize, yay ^^" << std::endl;
    MPI_Finalize();
}
