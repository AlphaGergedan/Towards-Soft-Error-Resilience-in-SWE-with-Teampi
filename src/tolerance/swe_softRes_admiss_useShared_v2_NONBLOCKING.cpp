/**
 * @file src/tolerance/swe_softRes_admiss_useShared_v2.cpp
 *
 * @brief soft error resilience with admissibility checks and using shared tasks
 *
 * Checks for admissibility of the computations (also see validateAdmissibility
 * in src/blocks/DimSplitMPIOverdecomp.cpp) and only share the results if they
 * are admissible. If they are not, a healthy replica sends its block data to
 * the failed replicas. We use shared tasks immediately, and check
 * them for admissibility in case of SDC during transmission (undetected SDC
 * spreads immediately, but saves computation time for secondary blocks).
 *
 * Here is a short pseudo-code for the computation loop: TODO
 *
 *  while (t < simulationDuration) {
 *
 *    // start Heartbeat
 *
 *    while (timeSinceLastHeartbeat < heartbeatInterval &&
 *           t < simulationDuration) {
 *
 *      for each block in primaryBlocks
 *        // computeNumericalFluxes + validate Admissibility criteria
 *        // report to my replicas with MPI_Send
 *
 *      if (SDC)
 *        while(true)
 *          // receive the reload replica
 *          for each block in corruptedPrimaryBlocks
 *            // receive b,h,hv,hu for the current corrupted block
 *            // computeNumericalFluxes + validate admissibility cirteria
 *            // report to my replicas with MPI_Send
 *          if (allFixed) break
 *
 *      while (true) {
 *        for each block in secondaryBlocks
 *          // receive report from replica
 *
 *        if (NoSDC) break
 *          // the healthy replica with the lowest index is the reload replica
 *          // reload replica sends its index to all the corrupted replicas
 *          for each block corruptedBlocks
 *            // reload replica sends b,h,hv,hu
 *      }
 *
 *      for each block in blocks
 *        if (block is primaryBlock) // send block to replicas
 *        else // receive block from the replicas
 *
 *      for each block in secondaryBlocks
 *        // validate admissibility criteria
 *        // report to my replica which sent this block
 *
 *      if (SDC)
 *        while(true)
 *          for each block in corruptedSecondaryBlocks
 *            // receive the updates again for the current corrupted block from the sender replica
 *            // validate admissibility criteria
 *            // report to my replica which sent this block
 *          if (allFixed) break
 *
 *      while (true) {
 *        for each block in primaryBlocks
 *          // receive report from the replicas that received it
 *
 *        if (NoSDC) break
 *          for each block in corruptedBlocks
 *            // send updates to the corrupted replicas
 *      }
 *    }
 *
 *    // end Heartbeat
 *
 *  } // end of while (t < simulationDuration)
 *
 * @author Atamert Rahma rahma@in.tum.de
 */


#include <limits>
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

int main(int argc, char** argv) {
    double startTime = MPI_Wtime();

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
    args.addOption("restart-basepath", 'r', "Restart base file name", tools::Args::Required, false);
    args.addOption("write-output", 'w', "Write output using netcdf writer to the specified output file", args.No, false);
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
    backupNameInput = "BACKUP_" + outputNameInput;
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

    /* Begin the logger */
    tools::FtLogger ft_logger(myTeam, myRankInTeam);
    //ft_logger.ft_print_spawnStatus();
    // Print status
    char hostname[HOST_NAME_MAX];
    gethostname(hostname, HOST_NAME_MAX);
    std::printf("PID %d, Rank %i of Team %i spawned at %s with start time %f\n", getpid(), myRankInTeam, myTeam, hostname, startTime);

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
        // TODO testing the propagation of SDC
        //scenario = new SWE_SeaAtRestScenario();
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
    for (auto& block : simulationBlocks) block->saveBathymetry();

    std::vector<float> timesteps;
    // Simulated time
    t = simulationStart;

    float timestep;

    /**
     * Contains the block numbers that needs to be calculated on the specific
     * rank. There are 'numTeams * decompFactor * ranksPerTeam' blocks in total
     * when sharing (each block es a complete 'task') so this helps us to order
     * the blocks by their numbers across this rank's replicas to compute and
     * share the results.
     *
     * In the naiv case, we have only one block per rank. So there is no
     * sharing.
     */
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
    if (writeOutput && simulationStart == 0.f) {
        for (auto& block : simulationBlocks) { block->writeTimestep(0.f); }
    }

    const MPI_Comm interTeamComm{TMPI_GetInterTeamComm()};

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
        if (verbose) ft_logger.ft_print_HBstart(timeOfLastHeartbeat, t);

        // simulate until the checkpoint is reached
        while (timeSinceLastHeartbeat < heartbeatInterval && t < simulationDuration) {
            if (verbose) ft_logger.ft_print_loop(timeSinceLastHeartbeat, heartbeatInterval,
                                                 simulationDuration, t);

            // exchange boundaries between blocks
            for (auto& currentBlock : simulationBlocks) { currentBlock->setGhostLayer(); }
            for (auto& currentBlock : simulationBlocks) { currentBlock->receiveGhostLayer(); }

            // Avoid overwriting an old send buffer before everyone reaches this point TODO find a better solution
            MPI_Barrier(interTeamComm);

            /* indicates if the primary block is corrupted */
            unsigned char primaryBlocksCorrupted[decompFactor];
            std::vector<MPI_Request> reportSendReqs(decompFactor * numTeams, MPI_REQUEST_NULL);
            unsigned char SDC = 0;

            /* compute and validate the primary blocks */
            for (unsigned int i = 0; i < decompFactor; i++) {
                // Get a reference to the current block
                const int& currentBlockNr{myBlockOrder[i]};
                auto& currentBlock = *simulationBlocks[currentBlockNr];

                currentBlock.computeNumericalFluxes();

                /* Inject a bitflip at random team and random rank */
                if (bitflip_at >= 0  && t > bitflip_at) {
                    /* Seed the random generator */
                    std::srand (static_cast <unsigned> (time(NULL)));
                    int teamToCorrupt = std::rand() % numTeams;
                    int rankToCorrupt = std::rand() % ranksPerTeam;
                    if (myTeam == teamToCorrupt && myRankInTeam == rankToCorrupt) {
                        std::cout << "T" << myTeam << "R" << myRankInTeam
                                << " : INJECTING A BITFLIP" << std::endl;
                        //currentBlock.injectBigNumber_intoData();
                        currentBlock.injectNaN_intoData();
                        //currentBlock.injectRandomBitflip();
                        //currentBlock.injectRandomBitflip_intoData();
                        //currentBlock.injectRandomBitflip_intoUpdates();
                    }

                    /* prevent any other bitflip */
                    bitflip_at = -1.f;
                }

                /* check for soft errors */
                bool admissible = currentBlock.validateAdmissibility(t);
                if (!admissible) {
                    /* try to fix SDC by recomputing */
                    currentBlock.computeNumericalFluxes();

                    /* check if SDC is still present */
                    admissible = currentBlock.validateAdmissibility(t);
                    if (!admissible) {
                        primaryBlocksCorrupted[i] = 1;
                        SDC = 1;
                    }
                    else primaryBlocksCorrupted[i] = 0;

                }
                else primaryBlocksCorrupted[i] = 0;

                /* report to all replicas */
                for (int destTeam = 0; destTeam < numTeams; destTeam++) {
                    // primaryblock validation starts by tag 100 : TODO Isend integration is easy
                    if (destTeam != myTeam) {
                        MPI_Isend(primaryBlocksCorrupted+i, 1, MPI_BYTE, destTeam,
                                  100, interTeamComm, &reportSendReqs[(decompFactor * destTeam) + i]);
                    }
                }
            }

            /* go into recovery mode to receive */
            while (SDC) {
                MPI_Waitall(decompFactor * numTeams, reportSendReqs.data(), MPI_STATUSES_IGNORE);
                // TODO debugging
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : SDC detected in my primary blocks"
                          << std::endl;
                int reloadReplica;
                /* Receive the reload replica : tag 20 for reload replica */
                MPI_Recv(&reloadReplica, 1, MPI_INT, MPI_ANY_SOURCE, 20, interTeamComm, MPI_STATUS_IGNORE);
                // TODO debugging
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : I will receive from my replica in team " << reloadReplica
                          << std::endl;
                /* for each corrupted block receive b,h,hv,hu */
                for (unsigned int i = 0; i < decompFactor; i++) {
                    if (primaryBlocksCorrupted[i]) {
                        /* Get a reference to the current corrupted block */
                        const int& currentBlockNr = myBlockOrder[i];
                        auto& currentCorruptedBlock = *simulationBlocks[currentBlockNr];
                        /* Size of the arrays b,h,hv,hu */
                        const int dataArraySize = (currentCorruptedBlock.nx + 2) * (currentCorruptedBlock.ny + 2);
                        // TODO debugging
                        std::cout << "T" << myTeam << "R" << myRankInTeam
                                  << " : receiving b,h,hv,hu for my primaryBlock" << currentBlockNr
                                  << std::endl;
                        // primaryBlock recovery by tag 21 : receive order is important! --> b,h,hv,hu
                        MPI_Recv(currentCorruptedBlock.b.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, reloadReplica, 21,
                                 interTeamComm, MPI_STATUS_IGNORE);
                        MPI_Recv(currentCorruptedBlock.h.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, reloadReplica, 21,
                                 interTeamComm, MPI_STATUS_IGNORE);
                        MPI_Recv(currentCorruptedBlock.hv.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, reloadReplica, 21,
                                 interTeamComm, MPI_STATUS_IGNORE);
                        MPI_Recv(currentCorruptedBlock.hu.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, reloadReplica, 21,
                                 interTeamComm, MPI_STATUS_IGNORE);
                        // TODO debugging
                        std::cout << "T" << myTeam << "R" << myRankInTeam
                                  << " : b,h,hv,hu for my primaryBlock" << currentBlockNr
                                  << " are received! Thanks replica " << reloadReplica
                                  << std::endl;
                        /* compute and validate again */
                        currentCorruptedBlock.computeNumericalFluxes();
                        bool admissible = currentCorruptedBlock.validateAdmissibility(t);
                        if (!admissible) assert(false); // TODO we must abort right ?
                        /* problem solved for this corrupted block */
                        else {
                            primaryBlocksCorrupted[i] = 0;
                            // TODO debugging
                            std::cout << "T" << myTeam << "R" << myRankInTeam
                                      << " : problem solved for my primaryBlock" << currentBlockNr
                                      << std::endl;
                            //MPI_Send(primaryBlocksCorrupted+i, 1, MPI_BYTE, reloadReplica, 100, interTeamComm); TODO integrate this later.. we assume we solved the SDC
                        }
                    }
                }
                SDC = 0;
                for (unsigned int i = 0; i < decompFactor; i++)
                    if (primaryBlocksCorrupted[i] == 1) SDC = 1;
            } // TODO up this point we assume we fixed the SDC

            /* indicates if the secondary block is corrupted in replica */
            unsigned char secondaryBlocksCorrupted[blocksPerRank - decompFactor];
            MPI_Request reportRecvReqs[blocksPerRank - decompFactor];
            unsigned char SDC_secondary= 0;
            /* receive report from other replicas */
            for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                /* Get a reference to the current corrupted block */
                const int& currentBlockNr = myBlockOrder[i];
                int source_rank = currentBlockNr % numTeams;
                MPI_Irecv(secondaryBlocksCorrupted+(i-decompFactor), 1, MPI_BYTE,
                          source_rank, 100, interTeamComm, &reportRecvReqs[i-decompFactor]);
            }
            /* we must wait for all the reports */
            MPI_Waitall(blocksPerRank - decompFactor, reportRecvReqs, MPI_STATUSES_IGNORE);
            /* check all the received reports */
            for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                if (secondaryBlocksCorrupted[i-decompFactor] == 1) SDC_secondary = 1;
            }
            /* if an error is present */
            if (SDC_secondary) {
                // TODO debugging
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : SDC reported from other replicas" << std::endl;
                /* recovery mode */
                /* figure out if this replica is the reload replica (lowest rank)
                 * lowest secondary block is the lowest rank */
                bool lowestHealthyReplica = true;
                for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                    /* if this replica has lower index than my team */
                    if ((myBlockOrder[i] % numTeams) < myTeam) {
                        /* if this replica is corrupted */
                        lowestHealthyReplica &= secondaryBlocksCorrupted[i-decompFactor];
                    }
                }
                if (lowestHealthyReplica) {
                    // TODO debugging
                    std::cout << "T" << myTeam << "R" << myRankInTeam
                              << " : I will send the blocks to corrupted replicas " << std::endl;
                    /* send b,h,hv,hu to all the corrupted replicas */
                    for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                        if (secondaryBlocksCorrupted[i-decompFactor]) {
                            /* Get a reference to the current corrupted block */
                            const int& currentBlockNr = myBlockOrder[i];
                            auto& currentSecondaryBlock = *simulationBlocks[currentBlockNr];
                            /* Size of the arrays b,h,hv,hu */
                            const int dataArraySize = (currentSecondaryBlock.nx + 2) * (currentSecondaryBlock.ny + 2);
                            int source_rank = currentBlockNr % numTeams;
                            // TODO debugging
                            std::cout << "T" << myTeam << "R" << myRankInTeam
                                      << " : sending to replica in team " << source_rank
                                      << std::endl;
                            /* send the reload replica by tag 20 */
                            MPI_Send(&myTeam, 1, MPI_INT, source_rank, 20, interTeamComm);
                            /* send b,h,hv,hu*/
                            MPI_Send(currentSecondaryBlock.b.getRawPointer(),
                                     dataArraySize, MPI_FLOAT, source_rank, 21,
                                     interTeamComm);
                            MPI_Send(currentSecondaryBlock.h.getRawPointer(),
                                     dataArraySize, MPI_FLOAT, source_rank, 21,
                                     interTeamComm);
                            MPI_Send(currentSecondaryBlock.hv.getRawPointer(),
                                     dataArraySize, MPI_FLOAT, source_rank, 21,
                                     interTeamComm);
                            MPI_Send(currentSecondaryBlock.hu.getRawPointer(),
                                     dataArraySize, MPI_FLOAT, source_rank, 21,
                                     interTeamComm);
                            // TODO debugging
                            std::cout << "T" << myTeam << "R" << myRankInTeam
                                      << " : sent to replica in team " << source_rank
                                      << std::endl;
                        }
                    }
                }
            } // end of recovery mode
            // TODO we just send the data once and forget about the corrupted replicas
            //      1. other teams should wait in case reload team gets corrupted
            //      2. corrupted teams should send a report that they are ok

            /* Primary Block computation+admissibilityChecks+report are finished */

            std::vector<MPI_Request> send_reqs(11 * numTeams, MPI_REQUEST_NULL);

            /* TODO
             *  In this loop do:
             *    If primary block, it is already computed and
             *                      checked for SDCs, we send it
             *    else it is secondary block, we receive it and
             *                      check it for SDC
             */
            for (unsigned int i = 0; i < blocksPerRank; i++) {
                // Get a reference to the current block
                const int& currentBlockNr{myBlockOrder[i]};
                auto& currentBlock = *simulationBlocks[currentBlockNr];

                // Size of the update fields incl. ghost layer
                const int fieldSizeX{(currentBlock.nx + 2) * (currentBlock.ny + 2)};
                const int fieldSizeY{(currentBlock.nx + 1) * (currentBlock.ny + 2)};

                // The first [decompFactor] blocks are the ones we always compute ourselves
                if (i < decompFactor) {
                    // Send to all other teams
                    for (int destTeam{0}; destTeam < numTeams; destTeam++) {
                        if (destTeam == myTeam || hasRecovered == true)
                        {
                            // Do not send to myself
                            // Skip send if I recovered (only if from hard failure)
                        }
                        else {
                            //std::cout << "Team " << myTeam << ": Sending t=" << t << " to Team " << destTeam << '\n';
                            if (verbose) ft_logger.ft_block_sending(t, destTeam);

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
                    int code = MPI_Waitall(11, reqs.data(), MPI_STATUSES_IGNORE);
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
                    }
                    else {
                        // currentBlock.computeNumericalFluxes();
                        //std::cout << "Team " << myTeam << ": Received t=" << recv_t << " from Team " << source_rank
                                  //<< '\n';
                        //if (verbose) ft_logger.ft_block_received(recv_t, source_rank);
                    }
                } // end of else, secondary blocks
            }

            unsigned char receivedBlocksCorrupted[blocksPerRank - decompFactor];
            MPI_Request reportSendReqs_second[blocksPerRank - decompFactor];
            unsigned char SDC_received = 0;
            /* validate the received secondary blocks */
            for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                // Get a reference to the current block
                const int& currentBlockNr = myBlockOrder[i];
                auto& currentReceivedBlock = *simulationBlocks[currentBlockNr];
                int source_rank = currentBlockNr % numTeams;

                /* TODO we can also inject bitflip here for testing */

                /* check for soft errors */
                bool admissible = currentReceivedBlock.validateAdmissibility(t);
                if (!admissible) {
                    /* TODO the error can either be in the sender or MPI communication */
                    receivedBlocksCorrupted[i-decompFactor] = 1;
                    SDC_received = 1;
                }
                else receivedBlocksCorrupted[i-decompFactor] = 0;

                /* report to the owner of the received block */
                // receivedBlock validation starts by tag 200 : TODO Isend integration is easy
                MPI_Isend(receivedBlocksCorrupted+(i-decompFactor), 1, MPI_BYTE, source_rank,
                          200, interTeamComm, &reportSendReqs_second[i-decompFactor]);
                //MPI_Send(receivedBlocksCorrupted+i, 1, MPI_BYTE, source_rank, 200, interTeamComm);
            }

            /* go into recovery mode for receivedBlocks */
            while (SDC_received) {
                MPI_Waitall(blocksPerRank - decompFactor, reportSendReqs_second, MPI_STATUSES_IGNORE);
                /* we already know the reload replica, owner of the block
                 * TODO you can also handle more SDC cases where the owner also is
                 *      corrupted! This version assumes at this point that the owner
                 *      has validated his primary blocks, and they have no SDC */
                /* for each corrupted block receive b,h,hv,hu */
                for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                    if (receivedBlocksCorrupted[i-decompFactor]) {
                        /* Get a reference to the current corrupted block */
                        const int& currentBlockNr = myBlockOrder[i];
                        auto& currentCorruptedBlock = *simulationBlocks[currentBlockNr];
                        int source_rank = currentBlockNr % numTeams;
                        /* Size of the arrays b,h,hv,hu */
                        const int dataArraySize = (currentCorruptedBlock.nx + 2) * (currentCorruptedBlock.ny + 2);

                        // receivedBlock recovery by tag 22 : receive order is important! --> b,h,hv,hu
                        MPI_Recv(currentCorruptedBlock.b.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, source_rank, 22,
                                 interTeamComm, MPI_STATUS_IGNORE);
                        MPI_Recv(currentCorruptedBlock.h.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, source_rank, 22,
                                 interTeamComm, MPI_STATUS_IGNORE);
                        MPI_Recv(currentCorruptedBlock.hv.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, source_rank, 22,
                                 interTeamComm, MPI_STATUS_IGNORE);
                        MPI_Recv(currentCorruptedBlock.hu.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, source_rank, 22,
                                 interTeamComm, MPI_STATUS_IGNORE);
                        /* compute and validate again */
                        currentCorruptedBlock.computeNumericalFluxes();
                        /* receive the maxtimestep from the owner of this block */
                        MPI_Recv(&currentCorruptedBlock.maxTimestep,
                                 1, MPI_FLOAT, source_rank, 22, interTeamComm, MPI_STATUS_IGNORE);

                        bool admissible = currentCorruptedBlock.validateAdmissibility(t);
                        if (!admissible) assert(false); // TODO we must abort right ?
                        /* problem solved for this corrupted block */
                        else {
                            receivedBlocksCorrupted[i-decompFactor] = 0;
                            //MPI_Send(receivedBlocksCorrupted+i, 1, MPI_BYTE, reloadReplica, 100, interTeamComm); TODO integrate this later.. we assume we solved the SDC
                        }
                    }
                }
                SDC_received = 0;
                for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                    if (receivedBlocksCorrupted[i-decompFactor] == 1) SDC_received = 1;
                }
            } // TODO up this point we assume we fixed the SDC

            /* indicates if the sent block is corrupted. We will receive report from each replica
             *
             * [replicaLowest_0, replica_secondLowest_0, .... , replica_highest_0,
             *  replicaLowest_1, replica_secondLowest_1, .... , replica_highest_1,
             *  ...
             *  replicaLowest_decompFactor, replica_secondLowest_decompFactor, .... , replica_highest_decompFactor]
             */
            unsigned char sentBlocksCorrupted[(numTeams-1)*decompFactor];
            MPI_Request reportRecvReqs_second[(numTeams-1)*decompFactor];
            unsigned char SDC_sent = 0;
            /* receive report for all my primary blocks that I sent */
            int index = 0;
            for (unsigned int i = 0; i < decompFactor; i++) {
                for (int replica = 0; replica < numTeams; replica++) {
                    if (replica != myTeam) {
                        /* receive report */
                        MPI_Irecv(sentBlocksCorrupted+index, 1, MPI_BYTE, replica,
                                  200, interTeamComm, reportRecvReqs_second+index);
                        //MPI_Recv(sentBlocksCorrupted+index, 1, MPI_BYTE, replica, 200,
                                 //interTeamComm, MPI_STATUS_IGNORE);
                        if (sentBlocksCorrupted[index] == 1) SDC_sent = 1;
                        index++;
                    }
                }
            }
            /* we must wait for all the reports */
            MPI_Waitall((numTeams-1) * decompFactor, reportRecvReqs_second, MPI_STATUSES_IGNORE);
            /* check all the received reports */
            for (int i = 0; i < (numTeams-1)*decompFactor; i++) {
                if (sentBlocksCorrupted[i] == 1) SDC_sent = 1;
            }
            /* if error is reported */
            if (SDC_sent) {
                /* recovery mode */
                /* sent the blocks to replicas */
                int index = 0;
                for (unsigned int i = 0; i < decompFactor; i++) {
                    for (int replica = 0; replica < numTeams; replica++) {
                        if (replica != myTeam) {
                            /* send if SDC reported */
                            if (sentBlocksCorrupted[index] == 1) {
                                // Get a reference to the current block
                                const int& currentBlockNr{myBlockOrder[i]};
                                auto& currentSentBlock = *simulationBlocks[currentBlockNr];
                                /* Size of the arrays b,h,hv,hu */
                                const int dataArraySize = (currentSentBlock.nx + 2) * (currentSentBlock.ny + 2);
                                // sentBlock recovery by tag 22 : receive order is important! --> b,h,hv,hu
                                /* send b,h,hv,hu*/
                                MPI_Send(currentSentBlock.b.getRawPointer(),
                                         dataArraySize, MPI_FLOAT, replica, 22,
                                         interTeamComm);
                                MPI_Send(currentSentBlock.h.getRawPointer(),
                                         dataArraySize, MPI_FLOAT, replica, 22,
                                         interTeamComm);
                                MPI_Send(currentSentBlock.hv.getRawPointer(),
                                         dataArraySize, MPI_FLOAT, replica, 22,
                                         interTeamComm);
                                MPI_Send(currentSentBlock.hu.getRawPointer(),
                                         dataArraySize, MPI_FLOAT, replica, 22,
                                         interTeamComm);
                                /* This time we also need to send our maxtimestep */
                                MPI_Send(&currentSentBlock.maxTimestep,
                                         1, MPI_FLOAT, replica, 22,
                                         interTeamComm);
                            }
                            index++;
                        }
                    }
                }
            } // end of recovery mode TODO we assume that the error is fixed after sending 1 time .. maybe validate it ?

            // -------------------- task sharing finished, blocks are validated, we can continue

            hasRecovered = false;
            recoveredFromSDC = false;

            // determine max possible timestep
            timesteps.clear();
            for (auto& block : simulationBlocks) {
                timesteps.push_back(block->maxTimestep);
                //std::cout << "T" << myTeam << "R" << myRankInTeam
                          //<< " : Single Timestep = " << block->maxTimestep
                          //<< std::endl;
            }
            float minTimestep = *std::min_element(timesteps.begin(), timesteps.end());
            float test = *std::min_element(timesteps.begin(), timesteps.end()); // TODO for debugging WE DON'T NEED TO ALLLREDUCE  IF WE DON'T COMPUTE SECOND TASKS OURSELVES
            MPI_Allreduce(&minTimestep, &timestep, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
            //std::cout << "T" << myTeam << "R" << myRankInTeam
                      //<< " : Max Timestep = " << timestep << std::endl;
            assert(test == minTimestep);

            for (auto& block : simulationBlocks) block->maxTimestep = timestep;

            /* redundant saving of the previous results for admissibility checks */
            for (auto& block : simulationBlocks) block->savePreviousData();
            for (auto& block : simulationBlocks) block->updateUnknowns(timestep);

            t += timestep;

            /* write output */
            if (writeOutput) {
                if (verbose) ft_logger.ft_writingTimeStep(t);
                for (auto& block : simulationBlocks) {block->writeTimestep(t);}
            }

            timeSinceLastHeartbeat = MPI_Wtime() - timeOfLastHeartbeat;
            MPI_Bcast(&timeSinceLastHeartbeat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            if (kill_rank >= 0 && t >= kill_rank &&
                myTeam == 0 && myRankInTeam == 0 && restartNameInput == "") {
                std::cout << "Not supported... yet" << std::endl;
                assert(false);
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
        if (verbose) ft_logger.ft_print_HBend(timeSinceLastHeartbeat, t);
    }

    for (auto& block : simulationBlocks) { block->freeMpiType(); }

    double totalTime = MPI_Wtime() - startTime;
    std::cout << "Rank " << myRankInTeam << " from TEAM " << myTeam
              << " total time : " << totalTime << std::endl;

    MPI_Finalize();

    return 0;
}
