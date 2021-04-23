/**
 * @file src/swe_softRes_and_hardRes_woutTaskSharing.cpp
 *
 * @brief hard error resiliency without task sharing and no soft error resilience
 *
 * TODO
 *   - Soft error resilience needs to be added
 *   - You can send hashes using heartbeats
 *
 */



#include <bits/c++config.h>
#include <cstdint>
#include <functional>
#include <memory>
#include <mpi.h>
#include <string>
#include <system_error>
#include <teaMPI.h>
#include <type_traits>
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
#include "types/Boundary.hpp"

#include "tools/sha1.hpp"


/* Size of size_t to decide which MPI_Datatype we need */
#if SIZE_MAX == UCHAR_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_CHAR             /* 1 byte. a little extreme ? */
#elif SIZE_MAX == USHRT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_SHORT            /* 2 bytes */
#elif SIZE_MAX == UINT_MAX
   #define MPI_SIZE_T MPI_UNSIGNED                  /* 4 bytes */
#elif SIZE_MAX == ULONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG             /* 8 bytes */
#elif SIZE_MAX == ULLONG_MAX
   #define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG        /* 8 bytes */
#else
    #error "Size of size_t unknown."
#endif


// some data needs to be global, otherwise the checkpoint callbacks cannot
// access it.
jmp_buf jumpBuffer{};
int ranksPerTeam{1};
std::string restartNameInput{""};
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<std::string> backupMetadataNames{};

/* Since we are doing the first option naively, we will not share any
* blocks. Which means it is best that we have only one block per replica,
* so using only one variable makes more sense TODO*/

std::string backupMetadataName;
//------------------------------------------------------------------------------
std::string outputTeamName{""};
std::string backupTeamName{""};
std::string outputNameInput{""};
std::string backupNameInput{""};

float t(0.0F);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>> simulationBlocks{};

/* Since we are doing the first option naively, we will not share any
* blocks. Which means it is best that we have only one block per replica,
* so using only one variable makes more sense TODO*/

std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp> simulationBlock; /* needs to be global
                                                                         * for checkpointcallbacks
                                                                         * below, which will be
                                                                         * passed to teaMPI
                                                                         **/
//------------------------------------------------------------------------------
SWE_Scenario* scenario{nullptr};
bool hasRecovered{false};

void createCheckpointCallback(std::vector<int> failedTeams)
{
    int myRankInTeam;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    const int myTeam = TMPI_GetTeamNumber();
    std::printf("Rank %i of Team %i writing checkpoint for restoration.\n", myRankInTeam, myTeam);

    // write a checkpoint for every block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* TODO on this naiv soft error checking we have only one simulationBlock ! */
    for (int i = 0; i < simulationBlocks.size(); i++)
    {
        simulationBlocks[i]->createCheckpoint(t, backupMetadataNames[i], 0);
    }
//------------------------------------------------------------------------------

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
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* TODO on this naiv soft error checking we have only one simulationBlock ! */
    simulationBlocks.clear();
//------------------------------------------------------------------------------
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

//******************************************************************************

int main(int argc, char** argv)
{
    // Define command line arguments
    tools::Args args;

    args.addOption("simulation-duration", 't', "Time in seconds to simulate");
    args.addOption("resolution-x", 'x', "Number of simulated cells in x-direction");
    args.addOption("resolution-y", 'y', "Number of simulated cells in y-direction");
    args.addOption("decomp-factor",
                   'd',
                   "Split each rank into \"TEAMS\" * \"decomp-factor\" blocks",
                   tools::Args::Required,
                   false);
    args.addOption("output-basepath", 'o', "Output base file name");
    args.addOption("backup-basepath", 'b', "Output base file name");
    args.addOption("restart-basepath", 'r', "Restart base file name", args.Required, false);
    args.addOption("heartbeat-interval", 'i', "Wall-clock time in seconds to wait between heartbeats");
    args.addOption("hash-method", 'm', "Which hashing to use: (1 | 2), default: 1", args.Required, false);
    args.addOption("verbose", 'v', "Let the simulation produce more output, default: No", args.No, false);

    /* TODO add write option*/

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

    int hashOption;

    bool verbose = false;

    /* Read in command line arguments */
    simulationDuration = args.getArgument<float>("simulation-duration");
    /* changed from checkpoint interval TODO, remove this line if sure is valid */
    heartbeatInterval = args.getArgument<clock_t>("heartbeat-interval");
    nxRequested = args.getArgument<int>("resolution-x");
    nyRequested = args.getArgument<int>("resolution-y");

    if (args.isSet("decomp-factor")) {
        decompFactor = args.getArgument<unsigned int>("decomp-factor");
        if (decompFactor == 0) {
            decompFactor = 1;
        }
    }

    outputNameInput = args.getArgument<std::string>("output-basepath");
    backupNameInput = args.getArgument<std::string>("backup-basepath");
    if (args.isSet("restart-basepath")) {
        restartNameInput = args.getArgument<std::string>("restart-basepath");
    }

    hashOption = args.getArgument<int>("hash-method", 1);
    if (hashOption != 1 && hashOption != 2) {
        std::cout << "Invalid hash method. It has to be either:\n"
                  << "  i^420 or\n"
                  << "  the number of which engineers think it's e"
                  << std::endl;
        return 1;
    }


    /* whether to generate a verbose output */
    verbose = args.isSet("verbose");


    // init teaMPI
    std::function<void(std::vector<int>)> create(createCheckpointCallback);
    std::function<void(int)> load(loadCheckpointCallback);
    TMPI_SetCreateCheckpointCallback(&create);
    TMPI_SetLoadCheckpointCallback(&load);
    TMPI_SetErrorHandlingStrategy(TMPI_WarmSpareErrorHandler);

    // init MPI
    int myRankInTeam;
    int provided;
    int requested = MPI_THREAD_MULTIPLE;
    if (setjmp(jumpBuffer) == 0)
    {
        // MPI_Init_thread(&argc, &argv, requested, &provided);
        MPI_Init(&argc, &argv);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &ranksPerTeam);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    const int myTeam{TMPI_GetTeamNumber()};
    int numTeams = TMPI_GetInterTeamCommSize();
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //unsigned int blocksPerRank = numTeams * decompFactor;

    // TODO check
    /* Since we are doing the first option naively, we will not share any
     * blocks. Which means it is best that we have only one block per replica */
    unsigned int blocksPerRank = 1;
//------------------------------------------------------------------------------
    outputTeamName = outputNameInput + "_t" + std::to_string(myTeam);
    backupTeamName = backupNameInput + "_t" + std::to_string(myTeam);

    // Print status
    char hostname[HOST_NAME_MAX];
    gethostname(hostname, HOST_NAME_MAX);

    std::printf("PID %d, Rank %i of Team %i spawned at %s\n", getpid(), myRankInTeam, myTeam, hostname);
    fflush(stdout);

    int totalBlocks = blocksPerRank * ranksPerTeam;

    // number of SWE-Blocks in x- and y-direction
    int blockCountY = std::sqrt(totalBlocks);
    while (totalBlocks % blockCountY != 0) blockCountY--;
    int blockCountX = totalBlocks / blockCountY;

    int startPoint = myRankInTeam * blocksPerRank;

    float simulationStart{0.0f};

    // if not loading from a checkpoint
    if (restartNameInput == "")
    {
        scenario = new SWE_RadialBathymetryDamBreakScenario{};
        int widthScenario = scenario->getBoundaryPos(BND_RIGHT) - scenario->getBoundaryPos(BND_LEFT);
        int heightScenario = scenario->getBoundaryPos(BND_TOP) - scenario->getBoundaryPos(BND_BOTTOM);

        float dxSimulation = (float)widthScenario / nxRequested;
        float dySimulation = (float)heightScenario / nyRequested;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* TODO on this naiv soft error checking we have only one simulationBlock ! */
        /* Loop over the rank's block number in its team. */
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

            /* Add the block to be calculated by this rank. */
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
//------------------------------------------------------------------------------

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* TODO on this naiv soft error checking we have only one simulationBlock ! */
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

            /* For all my neighbours. */
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

            /* For checkpoints. */
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
//------------------------------------------------------------------------------
    }
    else /* loading from a checkpoint */
    {
        std::vector<SWE_Scenario*> scenarios{};
        /* TODO same sitaution as the if block above.. there is only one block in this rank. */
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
    } /* end of else */


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* TODO on this naiv soft error checking we have only one simulationBlock ! */
    for (auto& block : simulationBlocks) block->sendBathymetry();
    for (auto& block : simulationBlocks) block->recvBathymetry();
//------------------------------------------------------------------------------

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

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
// TODO
//
// THERE IS ONLY ONE BLOCK!

    // Add my primary blocks first to block order
    //for (int i = myTeam; i < blocksPerRank; i += numTeams) { myBlockOrder.push_back(i); }
    // Add all other secondary blocks to block order
    //for (int i = 0; i < blocksPerRank; i++)
    //{
    //    if (i % numTeams != myTeam)
    //    {
    //        myBlockOrder.push_back(i);
    //    }
    //}

    myBlockOrder.push_back(0); /* We have only one block,
                                  but this is conventionally same as main.cpp */

    /* We should have only one block */
    assert(simulationBlocks.size() == 1);

//------------------------------------------------------------------------------

    // Write zero timestep
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* TODO on this naiv soft error checking we have only one simulationBlock ! */
    for (auto& block : simulationBlocks) { block->writeTimestep(0.f); }
//------------------------------------------------------------------------------

    /* For debugging only TODO count the heartbeats */
    int heartbeatCounter = 0;

    // simulate until end of simulation
    while (t < simulationDuration)
    {

        /**
         * Start hearbeat before the task begins.
         * If we post a heartbeat after every time step, a single rank would
         * delay the heartbeats of neighouring ranks because point-to-point
         * messages (like receiveGhostLayer) are required before a new
         * simulation time step can be started.
         */
        MPI_Sendrecv(MPI_IN_PLACE,              /* Send buffer      */
                     0,                         /* Send count       */
                     MPI_BYTE,                  /* Send type        */
                     MPI_PROC_NULL,             /* Destination      */
                     1,                         /* Send tag         */
                     MPI_IN_PLACE,              /* Receive buffer   */
                     0,                         /* Receive count    */
                     MPI_BYTE,                  /* Receive type     */
                     MPI_PROC_NULL,             /* Source           */
                     0,                         /* Receive tag      */
                     MPI_COMM_SELF,             /* Communicator     */
                     MPI_STATUS_IGNORE);        /* Status object    */

        double timeOfLastHeartbeat{MPI_Wtime()};
        double timeSinceLastHeartbeat{0.0};

        /* print heartbeat for debugging */
        std::cout << "\n\t++++++++++++++ " << heartbeatCounter << ". "
                  << "HEARTBEAT START ++++++++++++++\n\t"
                  << "  - team:\t\t" << myTeam << "\n\t"
                  << "  - rank:\t\t" << myRankInTeam << "\n\t"
                  << "at t = " << t << ", at MPI_Wtime() : " << timeOfLastHeartbeat << std::endl;


        /* hash the computation and store the result here for comparison
         * with the other replicas
         *
         * TODO hash with SHA1 implementation or std::hash ?
         *      test this. boost::hash_combine for combining std::hashes?
         */

        /* hash function and hash storage */
        std::hash<std::string> hash_fn;
        std::size_t total_hash = 0; /* initialized as 0 because we want the
                                     * first xor operation with the first hash
                                     * would give us the
                                     * hash itself
                                     */

        /* sha1 hash */
        SHA1 checksum;

        // simulate until the heartbeat interval is reached
        while (timeSinceLastHeartbeat < heartbeatInterval && t < simulationDuration) {

            if (verbose) {
                std::cout << "\n"
                          << "--------------- while (timeSinceLastHeartbeat="
                          << timeSinceLastHeartbeat << " < heartbeatInterval="
                          << heartbeatInterval << " && t="
                          << t << " < simulationDuration="
                          << simulationDuration << ")---------------"
                          << "TEAM " << myTeam << ": \t\tBegin iteration at time " << t
                          << std::endl;
            }

            // exchange boundaries between blocks
            for (auto& currentBlock : simulationBlocks) { currentBlock->setGhostLayer(); }
            for (auto& currentBlock : simulationBlocks) { currentBlock->receiveGhostLayer(); }
            const MPI_Comm interTeamComm{TMPI_GetInterTeamComm()};

            /* TODO: Since there is only on blockPerRank we don't need to wait */
            if (blocksPerRank > 1 && !hasRecovered)
            {
                assert(false);
                std::cout << "\n\n\n\n\n\t\t\t************* ATTENTION **************"
                          << "\t\t\tTHIS BLOCK SHOULD NEVER BE VISITED, NAIV.CPP IS PROBLEMATIC!"
                          << "\n\n\n\n\n\t\t\t************* ATTENTION **************"
                          << std::endl;
                // Avoid overwriting an old send buffer before everyone reaches this point
                MPI_Barrier(interTeamComm);
            }
            /* No sending will be make */
            //std::vector<MPI_Request> send_reqs(11 * numTeams, MPI_REQUEST_NULL);

            /* There is just one block in the naive version! TODO unneccesary loop */
            for (int i{0}; i < blocksPerRank; i++)
            {
                // Get a reference to the current block
                const int& currentBlockNr{myBlockOrder[i]};
                auto& currentBlock = *simulationBlocks[currentBlockNr];

                // Size of the update fields incl. ghost layer
                const int fieldSizeX{(currentBlock.nx + 2) * (currentBlock.ny + 2)};
                const int fieldSizeY{(currentBlock.nx + 1) * (currentBlock.ny + 2)};

                // The first [decompFactor] blocks are the ones we always compute ourselves
                // HINT: There is only one block in the naive version anyway.. TODO
                if (i < decompFactor)
                {
                    currentBlock.computeNumericalFluxes();

                    if (verbose) {
                        std::cout << "Team " << myTeam << " calculated timestep "
                                  << currentBlock.maxTimestep
                                  << std::endl;
                    }

                    /* Pointers to the results to be hashed (with their size on top) */
                    // TODO: try to implement this huge block with a function above
                    //       for better readability

                    // fieldSizeX
                    float* calculated_hNetUpdatesLeft = currentBlock.hNetUpdatesLeft.getRawPointer();
                    float* calculated_hNetUpdatesRight = currentBlock.hNetUpdatesRight.getRawPointer();

                    // fieldSizeX
                    float* calculated_huNetUpdatesLeft = currentBlock.huNetUpdatesLeft.getRawPointer();
                    float* calculated_huNetUpdatesRight = currentBlock.huNetUpdatesRight.getRawPointer();

                    // fieldSizeY
                    float* calculated_hNetUpdatesBelow = currentBlock.hNetUpdatesBelow.getRawPointer();
                    float* calculated_hNetUpdatesAbove = currentBlock.hNetUpdatesAbove.getRawPointer();

                    float* calculated_hvNetUpdatesBelow = currentBlock.hvNetUpdatesBelow.getRawPointer();
                    float* calculated_hvNetUpdatesAbove = currentBlock.hvNetUpdatesAbove.getRawPointer();

                    float* calculated_maxTimeStep = &(currentBlock.maxTimestep);

                    /* Convert them to stringis */
                    std::string str_hLeft((const char*) calculated_hNetUpdatesLeft, sizeof(float) * fieldSizeX);
                    std::string str_hRight((const char*) calculated_hNetUpdatesRight, sizeof(float) * fieldSizeX);

                    std::string str_huLeft((const char*) calculated_huNetUpdatesLeft, sizeof(float) * fieldSizeX);
                    std::string str_huRight((const char*) calculated_huNetUpdatesRight, sizeof(float) * fieldSizeX);

                    std::string str_hBelow((const char*) calculated_hNetUpdatesBelow, sizeof(float) * fieldSizeY);
                    std::string str_hAbove((const char*) calculated_hNetUpdatesAbove, sizeof(float) * fieldSizeY);

                    std::string str_hvBelow((const char*) calculated_hvNetUpdatesBelow, sizeof(float) * fieldSizeY);
                    std::string str_hvAbove((const char*) calculated_hvNetUpdatesAbove, sizeof(float) * fieldSizeY);

                    std::string str_maxTimeStep((const char*) calculated_maxTimeStep, sizeof(float));

                    /* OPTION A: hash them using SHA1 (like redmpi did) TODO use library next time.. */
                    checksum.update(str_hLeft);
                    checksum.update(str_hRight);

                    checksum.update(str_huLeft);
                    checksum.update(str_huRight);

                    checksum.update(str_hBelow);
                    checksum.update(str_hAbove);

                    checksum.update(str_hvBelow);
                    checksum.update(str_hvAbove);

                    checksum.update(str_maxTimeStep);

                    /* OPTION B: hash them using std::hash
                         HINT: TODO xor is used when combining the hash but this
                         might not be ideal, see boost::hash_combine */

                    total_hash ^= hash_fn(str_hLeft);
                    total_hash ^= hash_fn(str_hRight);

                    total_hash ^= hash_fn(str_huLeft);
                    total_hash ^= hash_fn(str_huRight);

                    total_hash ^= hash_fn(str_hBelow);
                    total_hash ^= hash_fn(str_hAbove);

                    total_hash ^= hash_fn(str_hvBelow);
                    total_hash ^= hash_fn(str_hvAbove);

                    total_hash ^= hash_fn(str_maxTimeStep);

                    /* No sending to other teams because we avoid this in naive version*/
                } else { /* This else block should never be visited.. Remove this block if sure everyting working fine TODO */
                    std::cout << "\n\n\n\n\n\t\t\t************* ATTENTION **************"
                              << "\t\t\tTHIS BLOCK SHOULD NEVER BE VISITED, NAIV.CPP IS PROBLEMATIC!"
                              << "\n\n\n\n\n\t\t\t************* ATTENTION **************"
                              << std::endl;
                    assert(false);
                } // end of else
            } // end of for (int i{0}; i < blocksPerRank; i++)

            hasRecovered = false;

            // determine max possible timestep TODO THIS IS NOT NECCESSARY WITH ONE BLOCK
            timesteps.clear();
            for (auto& block : simulationBlocks)
            {
                timesteps.push_back(block->maxTimestep);
                if (verbose) {
                    std::cout << "Team " << myTeam << ": Single Timestep " << block->maxTimestep << std::endl;
                }

            }
            float minTimestep = *std::min_element(timesteps.begin(), timesteps.end());
            MPI_Allreduce(&minTimestep, &timestep, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
            if (verbose) {
                std::cout << "Team " << myTeam << ": Max Timestep " << timestep << std::endl;
            }
            for (auto& block : simulationBlocks) { block->maxTimestep = timestep; }
            for (auto& block : simulationBlocks) { block->updateUnknowns(timestep); }
            t += timestep;

             /* write output */
            if (verbose) {
              std::cout << "\n"
                        << "Rank " << myRankInTeam << " from TEAM " << myTeam
                        << " writing output at t = " << t << std::endl;
            }
            for (auto& block : simulationBlocks) {block->writeTimestep(t);}

            /* Update the elapsed time after sending the heartbeat. */
            timeSinceLastHeartbeat = MPI_Wtime() - timeOfLastHeartbeat;

            /* TODO Why send bcast to the rank's team's 0 here ? */
            MPI_Bcast(&timeSinceLastHeartbeat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            // MPI_Barrier(TMPI_GetInterTeamComm());
            if (myTeam == 1 && myRankInTeam == 0 && restartNameInput == "" && t > 5.f)
            {
                // return 1;
            }
        }

        // End Heartbeat send the computed task to compare hashes
        //
        // std::cout << "Team " << myTeam << ": HEARTBEAT! Current simulation time is " << t << '\n';

        /* end the checksum SHA1*/
        std::string final_checksum = checksum.final();

        /* print heartbeat for debugging */
        std::cout << "\n\t+++++++++ " << heartbeatCounter << ". "
                  << "HEARTBEAT END SENDING HASH ++++++++\n\t"
                  << "  - team:\t\t" << myTeam << "\n\t"
                  << "  - rank:\t\t" << myRankInTeam << "\n\t"
                  << "  - size_t total_hash (" << sizeof(total_hash)
                  << " bytes)" << " : " << total_hash
                  << " ; SHA1 checksum (" << final_checksum.size()
                  << " bytes)" << " : " << final_checksum
                  << " at t = " << t << ", at MPI_Wtime() : "
                  << timeOfLastHeartbeat << std::endl;

        heartbeatCounter++;

        /* Send the hash calculated by std::hash and XORed with size size_t */
        if (hashOption == 1) {

            /* heartbeat end : OPTION B */
            MPI_Sendrecv(&total_hash,              /* Send buffer      */
                        1,                         /* Send count       */
                        MPI_SIZE_T,                /* Send type        */
                        MPI_PROC_NULL,             /* Destination      */
                        -1,                        /* Send tag         */
                        MPI_IN_PLACE,              /* Receive buffer   */
                        0,                         /* Receive count    */
                        MPI_BYTE,                  /* Receive type     */
                        MPI_PROC_NULL,             /* Source           */
                        0,                         /* Receive tag      */
                        MPI_COMM_SELF,             /* Communicator     */
                        MPI_STATUS_IGNORE);        /* Status object    */
        }

        /* Send the has calculated by SHA1, size is 20 byte */
        // TODO: normally size is 20 bytes, but when writing in hex
        //       using strings, we get 40 charactrers. So we are sending
        //       more than we have to, fix this by converting 40 chars
        //       to 20 bytes (read them as hex). Possible option could be
        //       to send 20 bytes as 20 chars again.
        if (hashOption == 2) {

            /* heartbeat end : OPTION A */
            MPI_Sendrecv(final_checksum.c_str(),   /* Send buffer      */
                        final_checksum.size(),     /* Send count: 20 ? */
                        MPI_CHAR,                  /* Send type        */
                        MPI_PROC_NULL,             /* Destination      */
                        -1,                        /* Send tag         */
                        MPI_IN_PLACE,              /* Receive buffer   */
                        0,                         /* Receive count    */
                        MPI_BYTE,                  /* Receive type     */
                        MPI_PROC_NULL,             /* Source           */
                        0,                         /* Receive tag      */
                        MPI_COMM_SELF,             /* Communicator     */
                        MPI_STATUS_IGNORE);        /* Status object    */
        }


        ///* Normal heartbeat */
        //MPI_Sendrecv(MPI_IN_PLACE,              /* Send buffer      */
                     //0,                         /* Send count       */
                     //MPI_BYTE,                  /* Send type        */
                     //MPI_PROC_NULL,             /* Destination      */
                     //-1,                        /* Send tag         */
                     //MPI_IN_PLACE,              /* Receive buffer   */
                     //0,                         /* Receive count    */
                     //MPI_BYTE,                  /* Receive type     */
                     //MPI_PROC_NULL,             /* Source           */
                     //0,                         /* Receive tag      */
                     //MPI_COMM_SELF,             /* Communicator     */
                     //MPI_STATUS_IGNORE);        /* Status object    */


        // printf("Returned from heartbeat rank: %d, team %d\n",
        // l_mpiRank, l_teamNumber); Only write timestep when simluation
        // has finished

        /*
        std::printf("Rank %i of Team %i writing output on timestep %f\n", myRankInTeam, myTeam, t);
        std::cout << "\n-------------------------------------------------"
                  << "Rank " << myRankInTeam << " from TEAM " << myTeam
                  << " writing output at t = " << t << std::endl;
        */

        /* TODO there is only one simulation block.. */
        /*
        for (int i = 0; i < (int) simulationBlocks.size(); i++) {
            simulationBlocks[i]->writeTimestep(t);
            //simulationBlocks[i]->createCheckpoint(t, backupMetadataNames[i], 0);
        }
        */

        // End Heartbeat
        //std::cout << "Team " << myTeam << ": HEARTBEAT! Current simulation time is " << t << '\n';
        //MPI_Sendrecv(MPI_IN_PLACE,
                     //0,
                     //MPI_BYTE,
                     //MPI_PROC_NULL,
                     //-1,
                     //MPI_IN_PLACE,
                     //0,
                     //MPI_BYTE,
                     //MPI_PROC_NULL,
                     //0,
                     //MPI_COMM_SELF,
                     //MPI_STATUS_IGNORE);

        // printf("Returned from heartbeat rank: %d, team %d\n",
        // l_mpiRank, l_teamNumber); Only write timestep when simluation
        // has finished
        /* if (t >= simulationDuration) {
            std::printf("Rank %i of Team %i writing final checkpoint\n", myRankInTeam, myTeam);
            for (int i = 0; i < simulationBlocks.size(); i++) {
                simulationBlocks[i]->createCheckpoint(t, backupMetadataNames[i], 0);
            }
        } */
    }

    for (auto& block : simulationBlocks) { block->freeMpiType(); }
    MPI_Finalize();
}
