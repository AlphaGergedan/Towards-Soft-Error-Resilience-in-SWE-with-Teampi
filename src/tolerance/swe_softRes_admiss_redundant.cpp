/**
 * @file src/tolerance/swe_softRes_admiss_redundant.cpp
 *
 * @brief Method "Redundant": Soft error resilience with admissibility checks and redundant computation
 *
 * Checks for admissibility of the computations (also see validateAdmissibility
 * in src/blocks/DimSplitMPIOverdecomp.cpp). If they are not admissible, a possibly
 * healthy replica sends its block data to the failed replicas. We don't share any tasks
 * unless a replica reports a possible SDC. This way an undetected SDC cannot be immediately
 * spread across all the replicas. This way we cannot save any computation overhead
 * caused by the replication but we can provide resilience for later detected SDCs,
 * because each replica keeps its results for itself.
 *
 * Here is a short pseudo-code for the computation loop:
 *
 *  while (t < simulationDuration) {
 *
 *    // start Heartbeat
 *
 *    while (timeSinceLastHeartbeat < heartbeatInterval &&
 *           t < simulationDuration) {
 *
 *      for each block in blocks
 *        // computeNumericalFluxes + validate Admissibility criteria
 *
 *      for replica in replicas
 *        // send a single report (if blocks have SDC)
 *
 *      // go into recovery mode to receive if SDC is present
 *
 *      for replica in replicas
 *        // receive a single report (if blocks in replica has SDC)
 *
 *      // go into recovery mode to send the blocks if SDC reported
 *
 *      // agree on a timestep
 *      // agree on timeSinceLastHeartbeat (root bcasts for simplification)
 *    }
 *
 *    // end Heartbeat
 *
 *  } // end of while (t < simulationDuration)
 *
 * @author Atamert Rahma rahma@in.tum.de
 */


#include <teaMPI.h>
#include <string>

#include <climits>
#include <csetjmp>
#include <iostream>

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

/* TODO only for hard error resilience. Currently not checked if hard resilience also working */
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

/* TODO only for hard error resilience. Currently not checked if hard resilience also working */
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
                   "Split each rank into \"decomp-factor\" blocks",
                   tools::Args::Required,
                   false);
    args.addOption("output-basepath", 'o', "Output base file name");
    args.addOption("restart-basepath", 'r', "Restart base file name", tools::Args::Required, false);
    args.addOption("write-output", 'w', "Write output using netcdf writer to the specified output file", args.No, false);
    args.addOption("inject-bitflip", 'f', "Injects a random bit-flip into a random data array in the rank 0 of team 0 right after the simulation time reaches the given time", args.Required, false);
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

    // this must be 1 for this method because decomposition doesn't make any sense
    // if we don't do immediate task sharing. User can still set other values to investigate
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
    /* only primary blocks */
    unsigned int blocksPerRank = decompFactor;
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
    else /* loading from a checkpoint TODO only for hard error resilience. Currently not checked if hard resilience also working */
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

    /* Write zero timestep if not restarting */
    if (writeOutput && simulationStart == 0.f) {
        for (auto& block : simulationBlocks) { block->writeTimestep(0.f); }
    }

    const MPI_Comm interTeamComm{TMPI_GetInterTeamComm()};

    /* array for flags if blocks are corrupted */
    unsigned char primaryBlocksCorrupted[decompFactor];
    /* indicates if a block is corrupted in this rank */
    unsigned char SDC = 0;
    /* indicates if a block is corrupted in this rank's replicas */
    unsigned char SDC_inReplica = 0;
    /* indicates which replicas are corrupted with flags */
    unsigned char replicaCorrupted[numTeams];

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

            /* compute and validate the primary blocks */
            for (unsigned int i = 0; i < decompFactor; i++) {
                // Get a reference to the current block
                auto& currentBlock = *simulationBlocks[i];

                currentBlock.computeNumericalFluxes();

                /* Inject a bitflip at team 0 and rank 0 */
                if (bitflip_at >= 0  && t > bitflip_at) {
                    if (myTeam == 0 && myRankInTeam == 0) {
                        std::cout << "T" << myTeam << "R" << myRankInTeam
                                << " : INJECTING A BITFLIP" << std::endl;
                        currentBlock.injectRandomBitflip();
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

            /* report to all replicas */
            for (int destTeam = 0; destTeam < numTeams; destTeam++) {
                // primaryblock validation
                if (destTeam != myTeam)
                    MPI_Send(&SDC, 1, MPI_BYTE, destTeam,
                             MPI_TAG_REPORT_PRIMARY_BLOCK, interTeamComm);
            }

            /* go into recovery mode to receive */
            while (SDC == 1) {
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : SDC detected in my primary blocks"
                          << std::endl;
                int reloadReplica;
                /* Receive the reload replica : tag 20 for reload replica */
                MPI_Recv(&reloadReplica, 1, MPI_INT, MPI_ANY_SOURCE,
                         MPI_TAG_RECEIVE_RELOAD_REPLICA, interTeamComm, MPI_STATUS_IGNORE);
                std::cout << "T" << myTeam << "R" << myRankInTeam
                          << " : I will receive from my replica in team " << reloadReplica
                          << std::endl;
                /* send reports for all my primary blocks */
                for (unsigned int i = 0; i < decompFactor; i++) {
                    MPI_Send(primaryBlocksCorrupted+i, 1, MPI_BYTE, reloadReplica,
                             MPI_TAG_REPORT_PRIMARY_BLOCK, interTeamComm);
                }
                /* for each corrupted block receive b,h,hv,hu */
                for (unsigned int i = 0; i < decompFactor; i++) {
                    if (primaryBlocksCorrupted[i] == 1) {
                        /* Get a reference to the current corrupted block */
                        auto& currentCorruptedBlock = *simulationBlocks[i];
                        /* Size of the arrays b,h,hv,hu */
                        const int dataArraySize = currentCorruptedBlock.dataArraySize;
                        std::cout << "T" << myTeam << "R" << myRankInTeam
                                  << " : receiving b,h,hv,hu for my block " << i
                                  << std::endl;
                        // primaryBlock recovery by tag 21 : receive order is important! --> b,h,hv,hu
                        MPI_Recv(currentCorruptedBlock.b.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, reloadReplica,
                                 MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                                 interTeamComm, MPI_STATUS_IGNORE);
                        MPI_Recv(currentCorruptedBlock.h.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, reloadReplica,
                                 MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                                 interTeamComm, MPI_STATUS_IGNORE);
                        MPI_Recv(currentCorruptedBlock.hv.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, reloadReplica,
                                 MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                                 interTeamComm, MPI_STATUS_IGNORE);
                        MPI_Recv(currentCorruptedBlock.hu.getRawPointer(),
                                 dataArraySize, MPI_FLOAT, reloadReplica,
                                 MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                                 interTeamComm, MPI_STATUS_IGNORE);
                        std::cout << "T" << myTeam << "R" << myRankInTeam
                                  << " : b,h,hv,hu for my block " << i
                                  << " are received! Thanks replica " << reloadReplica
                                  << std::endl;
                        /* compute and validate again */
                        currentCorruptedBlock.computeNumericalFluxes();
                        bool admissible = currentCorruptedBlock.validateAdmissibility(t);
                        if (!admissible) assert(false); // TODO we must abort right ?
                        /* problem solved for this corrupted block */
                        else {
                            primaryBlocksCorrupted[i] = 0;
                            std::cout << "T" << myTeam << "R" << myRankInTeam
                                      << " : problem solved for my primaryBlock" << i
                                      << std::endl;
                            //MPI_Send(primaryBlocksCorrupted+i, 1, MPI_BYTE, reloadReplica, 100, interTeamComm); TODO we assume we solved the SDC
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
                if (sourceTeam != myTeam) {
                    MPI_Recv(replicaCorrupted+sourceTeam, 1, MPI_BYTE, sourceTeam,
                             MPI_TAG_REPORT_PRIMARY_BLOCK, interTeamComm, MPI_STATUS_IGNORE);
                    if (replicaCorrupted[sourceTeam] == 1) SDC_inReplica = 1;
                }
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
                    /* indicates if a block is corrupted in replica */
                    unsigned char replicaBlocksCorrupted[decompFactor];
                    std::cout << "T" << myTeam << "R" << myRankInTeam
                              << " : I will send the blocks to corrupted replicas " << std::endl;
                    /* iterate the corrupted replicas */
                    for (unsigned int replica = 0; replica < numTeams; replica ++) {
                        /* if this replica has reported SDC */
                        if (replicaCorrupted[replica] == 1) {
                            /* send the reload replica by tag 20 */
                            MPI_Send(&myTeam, 1, MPI_INT, replica,
                                     MPI_TAG_RECEIVE_RELOAD_REPLICA, interTeamComm);
                        }
                    }
                    /* iterate through all the blocks */
                    for (unsigned int i = 0; i < decompFactor; i++) {
                        for (int currentReplica = 0; currentReplica < numTeams; currentReplica++) {
                            if (currentReplica != myTeam) {
                                /* if this replica has reported SDC */
                                if (replicaCorrupted[currentReplica] == 1) {
                                    /* receive report for the block */
                                    MPI_Recv(replicaBlocksCorrupted+i, 1, MPI_BYTE,
                                             currentReplica, MPI_TAG_REPORT_PRIMARY_BLOCK,
                                             interTeamComm, MPI_STATUS_IGNORE);
                                    /* send if it is corrupted */
                                    if (replicaBlocksCorrupted[i] == 1) {
                                        auto& currentBlock = *simulationBlocks[i];
                                        /* Size of the arrays b,h,hv,hu */
                                        const int dataArraySize = currentBlock.dataArraySize;
                                        std::cout << "T" << myTeam << "R" << myRankInTeam
                                                << " : sending to replica in team " << currentReplica
                                                << std::endl;
                                        /* send b,h,hv,hu*/
                                        MPI_Send(currentBlock.b.getRawPointer(),
                                                 dataArraySize, MPI_FLOAT, currentReplica,
                                                 MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                                                interTeamComm);
                                        MPI_Send(currentBlock.h.getRawPointer(),
                                                 dataArraySize, MPI_FLOAT, currentReplica,
                                                 MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                                                 interTeamComm);
                                        MPI_Send(currentBlock.hv.getRawPointer(),
                                                 dataArraySize, MPI_FLOAT, currentReplica,
                                                 MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                                                 interTeamComm);
                                        MPI_Send(currentBlock.hu.getRawPointer(),
                                                 dataArraySize, MPI_FLOAT, currentReplica,
                                                 MPI_TAG_RECOVERY_PRIMARY_BLOCK,
                                                 interTeamComm);
                                        std::cout << "T" << myTeam << "R" << myRankInTeam
                                                << " : sent to replica in team " << currentReplica
                                                << std::endl;
                                        replicaBlocksCorrupted[i] = 0;
                                    }
                                }
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

            // determine max possible timestep
            timesteps.clear();
            for (auto& block : simulationBlocks) {
                timesteps.push_back(block->maxTimestep);
                //std::cout << "T" << myTeam << "R" << myRankInTeam
                          //<< " : Single Timestep = " << block->maxTimestep
                          //<< std::endl;
            }
            float minTimestep = *std::min_element(timesteps.begin(), timesteps.end());
            /* we need to allreduce in the real world comm */
            PMPI_Allreduce(&minTimestep, &timestep, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
            //std::cout << "T" << myTeam << "R" << myRankInTeam
                      //<< " : Max Timestep = " << timestep << std::endl;

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
                //while(true) sleep(10);
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
