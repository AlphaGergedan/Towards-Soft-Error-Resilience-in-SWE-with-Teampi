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


#include <array>
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
        // TODO sea at rest can be used to test the propagation of SDC
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
    for (int i = 0; i < blocksPerRank; i++) {
        if (i % numTeams != myTeam) {
            myBlockOrder.push_back(i);
        }
    }

    /* Write zero timestep if not restarting */
    if (writeOutput && simulationStart == 0.f) {
        for (auto& block : simulationBlocks) { block->writeTimestep(0.f); }
    }

    const MPI_Comm interTeamComm{TMPI_GetInterTeamComm()};

    /* array for flags if blocks are corrupted */
    unsigned char primaryBlocksCorrupted[decompFactor];
    unsigned char receivedBlocksCorrupted[blocksPerRank - decompFactor];
    /* indicates if a block is corrupted in this rank */
    unsigned char SDC = 0;
    /* indicates if a block is corrupted in this rank's replicas */
    unsigned char SDC_inReplica = 0;
    /* indicates which replicas are corrupted with flags */
    unsigned char replicaCorrupted[numTeams];
    /* reports to send to replicas for the received blocks */
    unsigned char receivedBlockReports[numTeams];


    /* number of MPI receives to call in the task sharing */
    size_t totalRecvReqs_taskSharing = 3*(blocksPerRank-decompFactor);
    /* number of MPI sends to call in the task sharing*/
    size_t totalsendReqs_taskSharing = (3*decompFactor) * (numTeams-1);

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

            // Avoid overwriting an old send buffer before everyone reaches this point
            MPI_Barrier(interTeamComm);

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
            }

            /* report to all replicas
             * TODO don't send report for all the blocks everytime, send a single report first! */
            for (int destTeam = 0; destTeam < numTeams; destTeam++) {
                // primaryblock validation starts by tag 100 : TODO Isend integration is easy
                if (destTeam != myTeam) MPI_Send(&SDC, 1, MPI_BYTE, destTeam, 100, interTeamComm);
            }

            /* go into recovery mode to receive if SDC is present */
            while (SDC == 1) {
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
                /* send reports for all my primary blocks */
                for (unsigned int i = 0; i < decompFactor; i++) {
                    MPI_Send(primaryBlocksCorrupted+i, 1, MPI_BYTE, reloadReplica, 100, interTeamComm);
                }
                /* for each corrupted block receive b,h,hv,hu */
                for (unsigned int i = 0; i < decompFactor; i++) {
                    if (primaryBlocksCorrupted[i] == 1) {
                        /* Get a reference to the current corrupted block */
                        const int& currentBlockNr = myBlockOrder[i];
                        auto& currentCorruptedBlock = *simulationBlocks[currentBlockNr];
                        /* Size of the arrays b,h,hv,hu */
                        const int dataArraySize = (currentCorruptedBlock.nx + 2) * (currentCorruptedBlock.ny + 2);
                        // TODO debugging
                        std::cout << "T" << myTeam << "R" << myRankInTeam
                                  << " : receiving b,h,hv,hu for my primaryBlock" << currentBlockNr
                                  << ", block nr. is " << currentBlockNr
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
                                      << ", block nr. is " << currentBlockNr
                                      << std::endl;
                            //MPI_Send(primaryBlocksCorrupted+i, 1, MPI_BYTE, reloadReplica, 100, interTeamComm); TODO integrate this later.. we assume we solved the SDC
                        }
                    }
                }
                SDC = 0;
                for (unsigned int i = 0; i < decompFactor; i++)
                    if (primaryBlocksCorrupted[i] == 1) SDC = 1;
            } // TODO up this point we assume we fixed the SDC

            SDC_inReplica = 0;
            /* receive report from other replicas */
            for (int sourceTeam = 0; sourceTeam < numTeams; sourceTeam++) {
                // TODO Irecv integration
                if (sourceTeam != myTeam) MPI_Recv(replicaCorrupted+sourceTeam, 1, MPI_BYTE, sourceTeam, 100, interTeamComm, MPI_STATUS_IGNORE);
                if (replicaCorrupted[sourceTeam] == 1) SDC_inReplica = 1;
            }
            /* if an error is present */
            if (SDC_inReplica == 1) {
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : SDC reported from other replicas" << std::endl;
                /* recovery mode */
                /* figure out if this replica is the reload replica (lowest rank) */
                bool lowestHealthyReplica = true;
                for (unsigned int i = 0; i < myTeam; i++) {
                    lowestHealthyReplica &= replicaCorrupted[i];
                }

                /* if I am the recovery replica, I will save my replicas! */
                if (lowestHealthyReplica) {
                    /* indicates if the secondary block is corrupted in replica */
                    unsigned char secondaryBlocksCorrupted[blocksPerRank - decompFactor];
                    std::cout << "T" << myTeam << "R" << myRankInTeam
                              << " : I will send the blocks to corrupted replicas " << std::endl;
                    /* iterate the corrupted replicas */
                    for (unsigned int replica = 0; replica < numTeams; replica ++) {
                        /* if this replica has reported SDC */
                        if (replicaCorrupted[replica] == 1) {
                            /* send the reload replica by tag 20 */
                            MPI_Send(&myTeam, 1, MPI_INT, replica, 20, interTeamComm);
                        }
                    }
                    /* iterate through all the blocks */
                    for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                        // Get a reference to the current block
                        const int& currentBlockNr = myBlockOrder[i];
                        int currentReplica = currentBlockNr % numTeams;
                        /* if this replica has reported SDC */
                        if (replicaCorrupted[currentReplica] == 1) {
                            /* receive report for the block */
                            MPI_Recv(secondaryBlocksCorrupted+(i-decompFactor), 1, MPI_BYTE,
                                     currentReplica, 100, interTeamComm, MPI_STATUS_IGNORE);
                            /* send if it is corrupted */
                            if (secondaryBlocksCorrupted[i-decompFactor] == 1) {
                                auto& currentSecondaryBlock = *simulationBlocks[currentBlockNr];
                                /* Size of the arrays b,h,hv,hu */
                                const int dataArraySize = (currentSecondaryBlock.nx + 2) * (currentSecondaryBlock.ny + 2);
                                int source_rank = currentBlockNr % numTeams;
                                std::cout << "T" << myTeam << "R" << myRankInTeam
                                          << " : sending to replica in team " << source_rank
                                          << " the block " << currentBlockNr
                                          << std::endl;
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
                                std::cout << "T" << myTeam << "R" << myRankInTeam
                                          << " : sent to replica in team " << source_rank
                                          << std::endl;
                                secondaryBlocksCorrupted[i-decompFactor] = 0;
                            }
                        }
                    }
                }
                /* reset the flags */
                for (int i = 0; i < numTeams; i++) {
                    replicaCorrupted[i] = 0;
                }
                SDC_inReplica = 0;
            }
/*      ||||||||||||||||||| TODO REMOVE THIS BLOCK IF SURE CURRENT METHOD IS CORRECT ||||||||||||||||

            // receive report from other replicas
            for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                // Get a reference to the current corrupted block
                const int& currentBlockNr = myBlockOrder[i];
                int source_rank = currentBlockNr % numTeams;
                MPI_Recv(secondaryBlocksCorrupted+(i-decompFactor), 1, MPI_BYTE,
                         source_rank, 100, interTeamComm, MPI_STATUS_IGNORE);
                if (secondaryBlocksCorrupted[i-decompFactor] == 1) SDC_secondary = 1;
            }
            // if an error is present
            if (SDC_secondary) {
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : SDC reported from other replicas" << std::endl;
                // recovery mode
                // figure out if this replica is the reload replica (lowest rank)
                // lowest secondary block is the lowest rank
                bool lowestHealthyReplica = true;
                for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                    // if this replica has lower index than my team
                    if ((myBlockOrder[i] % numTeams) < myTeam) {
                        // if this replica is corrupted
                        lowestHealthyReplica &= secondaryBlocksCorrupted[i-decompFactor];
                    }
                }
                if (lowestHealthyReplica) {
                    std::cout << "T" << myTeam << "R" << myRankInTeam
                              << " : I will send the blocks to corrupted replicas " << std::endl;
                    // send b,h,hv,hu to all the corrupted replicas
                    for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                        if (secondaryBlocksCorrupted[i-decompFactor]) {
                            // Get a reference to the current corrupted block
                            const int& currentBlockNr = myBlockOrder[i];
                            auto& currentSecondaryBlock = *simulationBlocks[currentBlockNr];
                            //Size of the arrays b,h,hv,hu
                            const int dataArraySize = (currentSecondaryBlock.nx + 2) * (currentSecondaryBlock.ny + 2);
                            int source_rank = currentBlockNr % numTeams;
                            std::cout << "T" << myTeam << "R" << myRankInTeam
                                      << " : sending to replica in team " << source_rank
                                      << std::endl;
                            // send the reload replica by tag 20
                            MPI_Send(&myTeam, 1, MPI_INT, source_rank, 20, interTeamComm);
                            // send b,h,hv,hu
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
                            std::cout << "T" << myTeam << "R" << myRankInTeam
                                      << " : sent to replica in team " << source_rank
                                      << std::endl;
                        }
                    }
                }
            } // end of recovery mode
*/

            // TODO we just send the data once and forget about the corrupted replicas
            //      1. other teams should wait in case reload team gets corrupted
            //      2. corrupted teams should send a report that they are ok

            /* Primary Block computation+admissibilityChecks+report are finished */

            /* agree on a timestep */
            timesteps.clear();
            /* max timestep in primary blocks */
            for (unsigned int i = 0; i < decompFactor; i++) {
                /* Get a reference to the current block */
                const int& currentBlockNr = myBlockOrder[i];
                auto& currentPrimaryBlock = *simulationBlocks[currentBlockNr];
                timesteps.push_back(currentPrimaryBlock.maxTimestep);
            }
            float minTimestep = *std::min_element(timesteps.begin(), timesteps.end());
            float minTimeStepReplica;
            MPI_Allreduce(&minTimestep, &minTimeStepReplica, 1, MPI_FLOAT, MPI_MIN, interTeamComm);
            MPI_Allreduce(&minTimeStepReplica, &timestep, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
            for (auto& block : simulationBlocks) block->maxTimestep = timestep;

            /* post Irecvs for the secondary blocks */
            std::vector<MPI_Request> recv_reqs(totalRecvReqs_taskSharing, MPI_REQUEST_NULL);
            for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                /* Get a reference to the current secondary block */
                const int& currentBlockNr = myBlockOrder[i];
                auto& currentSecondaryBlock = *simulationBlocks[currentBlockNr];
                // Size of h,hv and hu
                const int dataArraySize = currentSecondaryBlock.dataArraySize;
                int source_rank = currentBlockNr % numTeams;
                // h
                MPI_Irecv(currentSecondaryBlock.h.getRawPointer(),
                          dataArraySize,
                          MPI_FLOAT,
                          source_rank,
                          1,
                          interTeamComm,
                          &recv_reqs[3*(i-decompFactor)]);
                // hv
                MPI_Irecv(currentSecondaryBlock.hv.getRawPointer(),
                          dataArraySize,
                          MPI_FLOAT,
                          source_rank,
                          2,
                          interTeamComm,
                          &recv_reqs[3*(i-decompFactor) + 1]);
                // hu
                MPI_Irecv(currentSecondaryBlock.hu.getRawPointer(),
                          dataArraySize,
                          MPI_FLOAT,
                          source_rank,
                          3,
                          interTeamComm,
                          &recv_reqs[3*(i-decompFactor) + 2]);
            }

            std::vector<MPI_Request> send_reqs(totalsendReqs_taskSharing, MPI_REQUEST_NULL);
            int requestIndex = 0;
            /* updateunknowns for my primary blocks + send them */
            for (unsigned int i = 0; i < decompFactor; i++) {
                const int& currentBlockNr = myBlockOrder[i];
                auto& currentPrimaryBlock = *simulationBlocks[currentBlockNr];
                currentPrimaryBlock.savePreviousData();
                currentPrimaryBlock.updateUnknowns(timestep);
                // Size of h,hv and hu
                const int dataArraySize = currentPrimaryBlock.dataArraySize;

                for (int destTeam = 0; destTeam < numTeams; destTeam++) {
                    if(destTeam != myTeam) {
                        // h
                        MPI_Isend(currentPrimaryBlock.h.getRawPointer(),
                                  dataArraySize,
                                  MPI_FLOAT,
                                  destTeam,
                                  1,
                                  interTeamComm,
                                  &send_reqs[requestIndex]);
                        requestIndex++;
                        // hv
                        MPI_Isend(currentPrimaryBlock.hv.getRawPointer(),
                                  dataArraySize,
                                  MPI_FLOAT,
                                  destTeam,
                                  2,
                                  interTeamComm,
                                  &send_reqs[requestIndex]);
                        requestIndex++;
                        // hu
                        MPI_Isend(currentPrimaryBlock.hu.getRawPointer(),
                                  dataArraySize,
                                  MPI_FLOAT,
                                  destTeam,
                                  3,
                                  interTeamComm,
                                  &send_reqs[requestIndex]);
                        requestIndex++;
                    }
                }
            }

            std::vector<int> checked(blocksPerRank-decompFactor, 0);
            int pendingBlocks = blocksPerRank-decompFactor;
            /* check if a secondary block is received here,
             * and check admissibility if received
             */
            while (pendingBlocks) {
                for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                    /* if not received yet */
                    if (!checked[i-decompFactor]) {
                        int h_received, hv_received, hu_received;
                        MPI_Test(&recv_reqs[3*(i-decompFactor)], &h_received, MPI_STATUS_IGNORE);
                        MPI_Test(&recv_reqs[3*(i-decompFactor) + 1], &hv_received, MPI_STATUS_IGNORE);
                        MPI_Test(&recv_reqs[3*(i-decompFactor) + 2], &hu_received, MPI_STATUS_IGNORE);
                        if (h_received && hv_received && hu_received) {
                            const int& currentBlockNr = myBlockOrder[i];
                            auto& currentSecondaryBlock = *simulationBlocks[currentBlockNr];
                            bool admissible = currentSecondaryBlock.validateAdmissibility_dataArrays(t);
                            if (!admissible) {
                                /* We have a problem. We assume this is a communication error or
                                * an SDC from our side and we will try to receive this block later */
                                receivedBlocksCorrupted[i] = 1;
                                SDC = 1;
                            }
                            else receivedBlocksCorrupted[i] = 0;
                            checked[i-decompFactor] = 1;
                            pendingBlocks--;
                        }
                    }
                }
            }
            /* report to owner replicas of the blocks
             * TODO don't send report for all the blocks everytime, send a single report first! */
            for (int destTeam = 0; destTeam < numTeams; destTeam++) {
                if (destTeam != myTeam) {
                    unsigned char report = 0;
                    for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                        /* this block is from my replica in the destTeam */
                        if (myBlockOrder[i] % numTeams == destTeam) {
                            if (receivedBlocksCorrupted[i] == 1) report = 1;
                        }
                    }
                    receivedBlockReports[destTeam] = report;
                    // receivedBlock validation starts by tag 200 : TODO Isend integration is easy
                    MPI_Send(receivedBlockReports+destTeam, 1, MPI_BYTE, destTeam, 200, interTeamComm);
                }
            }

            /* go into recovery mode to receive if SDC is present */
            while (SDC == 1) {
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : SDC detected in received blocks"
                          << std::endl;
                /* iterate all the replicas that we have received */
                for (int destTeam = 0; destTeam < numTeams; destTeam++) {
                    /* if blocks received from replica destTeam is faulty */
                    if (receivedBlockReports[destTeam] == 1 && destTeam != numTeams) {
                        /* send reports for all received blocks from destTeam */
                        for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                            if (myBlockOrder[i] % numTeams == destTeam) {
                                MPI_Send(receivedBlocksCorrupted+i, 1, MPI_BYTE, destTeam, 200, interTeamComm);
                            }
                        }
                    }
                }
                /* for each corrupted block receive b,h,hv,hu */
                for (unsigned int i = decompFactor; i < blocksPerRank; i++) {
                    if (receivedBlocksCorrupted[i] == 1) {
                        /* Get a reference to the current corrupted block */
                        const int& currentBlockNr = myBlockOrder[i];
                        int reloadReplica = currentBlockNr % numTeams;
                        auto& currentCorruptedBlock = *simulationBlocks[currentBlockNr];
                        /* Size of the arrays b,h,hv,hu */
                        const int dataArraySize = currentCorruptedBlock.dataArraySize;
                        // TODO debugging
                        std::cout << "T" << myTeam << "R" << myRankInTeam
                                  << " : receiving b,h,hv,hu for my secondary block " << currentBlockNr
                                  << ", block nr. is " << currentBlockNr
                                  << std::endl;
                        // receivedBlock recovery by tag 22 : receive order is important! --> b,h,hv,hu
                        MPI_Recv(currentCorruptedBlock.b.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, reloadReplica, 22,
                                 interTeamComm, MPI_STATUS_IGNORE);
                        MPI_Recv(currentCorruptedBlock.h.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, reloadReplica, 22,
                                 interTeamComm, MPI_STATUS_IGNORE);
                        MPI_Recv(currentCorruptedBlock.hv.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, reloadReplica, 22,
                                 interTeamComm, MPI_STATUS_IGNORE);
                        MPI_Recv(currentCorruptedBlock.hu.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, reloadReplica, 22,
                                 interTeamComm, MPI_STATUS_IGNORE);
                        // TODO debugging
                        std::cout << "T" << myTeam << "R" << myRankInTeam
                                  << " : b,h,hv,hu for my secondary Block" << currentBlockNr
                                  << " are received! Thanks replica " << reloadReplica
                                  << std::endl;
                        /* validate again */
                        bool admissible = currentCorruptedBlock.validateAdmissibility_dataArrays(t);
                        if (!admissible) assert(false); // TODO we must abort right ?
                        /* problem solved for this corrupted block */
                        else {
                            receivedBlocksCorrupted[i] = 0;
                            // TODO debugging
                            std::cout << "T" << myTeam << "R" << myRankInTeam
                                      << " : problem solved for my secondary Block" << currentBlockNr
                                      << ", block nr. is " << currentBlockNr
                                      << std::endl;
                            //MPI_Send(primaryBlocksCorrupted+i, 1, MPI_BYTE, reloadReplica, 100, interTeamComm); TODO integrate this later.. we assume we solved the SDC
                        }
                    }
                }
                SDC = 0;
                for (unsigned int i = decompFactor; i < blocksPerRank; i++)
                    if (receivedBlocksCorrupted[i] == 1) SDC = 1;
                /* reset the flags */
                for (int i = 0; i < numTeams; i++) receivedBlockReports[i] = 0;

            } // TODO up this point we assume we fixed the SDC

            SDC_inReplica = 0;
            /* receive report from other replicas */
            for (int sourceTeam = 0; sourceTeam < numTeams; sourceTeam++) {
                // TODO Irecv integration
                if (sourceTeam != myTeam) MPI_Recv(replicaCorrupted+sourceTeam, 1, MPI_BYTE, sourceTeam, 200, interTeamComm, MPI_STATUS_IGNORE);
                if (replicaCorrupted[sourceTeam] == 1) SDC_inReplica = 1;
            }
            /* if an error is present */
            if (SDC_inReplica == 1) {
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : SDC reported from other replicas" << std::endl;
                /* recovery mode, I will send my primary blocks again for the
                 * corrupted blocks in the corrupted replica */
                for (int destTeam = 0; destTeam < numTeams; destTeam++) {
                    if (replicaCorrupted[destTeam] == 1 && destTeam != myTeam) {
                        for (unsigned int i = 0; i < decompFactor; i++) {
                            unsigned char corrupted;
                            MPI_Recv(&corrupted, 1, MPI_BYTE, destTeam, 200, interTeamComm, MPI_STATUS_IGNORE);
                            /* if this block is corrupted in replica destTeam
                             * we send it again */
                            if (corrupted == 1) {
                                const int& currentBlockNr = myBlockOrder[i];
                                auto& currentPrimaryBlock = *simulationBlocks[currentBlockNr];
                                /* Size of the arrays b,h,hv,hu */
                                const int dataArraySize = currentPrimaryBlock.dataArraySize;
                                std::cout << "T" << myTeam << "R" << myRankInTeam
                                          << " : sending to replica in team " << destTeam
                                          << " my block " << currentBlockNr
                                          << std::endl;
                                /* send b,h,hv,hu*/
                                MPI_Send(currentPrimaryBlock.b.getRawPointer(),
                                         dataArraySize, MPI_FLOAT, destTeam, 22,
                                         interTeamComm);
                                MPI_Send(currentPrimaryBlock.h.getRawPointer(),
                                         dataArraySize, MPI_FLOAT, destTeam, 22,
                                         interTeamComm);
                                MPI_Send(currentPrimaryBlock.hv.getRawPointer(),
                                         dataArraySize, MPI_FLOAT, destTeam, 22,
                                         interTeamComm);
                                MPI_Send(currentPrimaryBlock.hu.getRawPointer(),
                                         dataArraySize, MPI_FLOAT, destTeam, 22,
                                         interTeamComm);
                                std::cout << "T" << myTeam << "R" << myRankInTeam
                                          << " : sent to replica in team " << destTeam
                                          << std::endl;
                            }
                        }
                    }
                }
                /* reset the flags */
                for (int i = 0; i < numTeams; i++) {
                    replicaCorrupted[i] = 0;
                }
                SDC_inReplica = 0;
            }
            /* wait for the sending to finish */
            MPI_Waitall(totalsendReqs_taskSharing, send_reqs.data(), MPI_STATUSES_IGNORE);

            // -------------------- task sharing finished, blocks are validated, we can continue

            /* TODO DEBUGGING */
            for (auto& block : simulationBlocks) assert(block->maxTimestep = timestep);
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
