/**
 * @file src/tolerance/swe_softRes_hardRes_hashes.cpp
 *
 * @brief hard error resiliency with soft error detection without task sharing
 *
 * TODO description
 *
 *
 *
 *  Here is a short pseudo-code for the computation loop:
 *
 *      Hasher swe_hasher;
 *
 *      // start heartbeat and update timeSinceLastHeartbeat
 *
 *      for (unsigned int i = 0; i < numberOfHashes; i++) {
 *
 *          while (t < sendHashAt[i]) {
 *
 *              if (timeSinceLastHeartbeat >= heartbeatInterval)
 *                  // end heartbeat && start heartbeat and update timeSinceLastHeartbeat
 *
 *              // compute
 *
 *              // update hash
 *
 *              // update timeSinceLastHeartbeat
 *          }
 *
 *          // send hash with heartbeat
 *      }
 *
 *      // end the last heartbeat
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

#include "tools/hasher.hpp"


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
std::string backupMetadataName;
//------------------------------------------------------------------------------
std::string outputTeamName{""};
std::string backupTeamName{""};
std::string outputNameInput{""};
std::string backupNameInput{""};

float t(0.0F);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp> simulationBlock;

/* needs to be global
 * for checkpointcallbacks
 * below, which will be
 * passed to teaMPI
 */

//------------------------------------------------------------------------------
SWE_Scenario* scenario{nullptr};
bool hasRecovered{false};

void createCheckpointCallback(std::vector<int> failedTeams)
{
    int myRankInTeam;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    const int myTeam = TMPI_GetTeamNumber();
    std::printf("Rank %i of Team %i writing checkpoint for restoration.\n", myRankInTeam, myTeam);

    // TODO write a checkpoint for every block

    /* Blocking send to the failed teams to tell them I wrote the checkpoint,
     * you can now reload
     */
    int send = 1;
    for (const auto& failedTeam : failedTeams)
    {
        std::cout << "\n-------------"
                  << "I WROTE THE CHECKPOINT YOU CAN NOW RELOAD\n"
                  << " -> my team : " << myTeam << "\n"
                  << " -> my rank in team : " << myRankInTeam << "\n"
                  << " -> sending message to team : " << failedTeam << "\n"
                  << " -> sending message to world rank : " << TMPI_TeamToWorldRank(myRankInTeam, failedTeam)
                  << "\n-------------" << std::endl;
        MPI_Send(&send, 1, MPI_INT, TMPI_TeamToWorldRank(myRankInTeam, failedTeam), 0, TMPI_GetWorldComm());
    }
}

void loadCheckpointCallback(int reloadTeam)
{
    int myRankInTeam;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    const int myTeam = TMPI_GetTeamNumber();

    MPI_Status status;
    int recv_buf;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    MPI_Recv(&recv_buf, 1, MPI_INT, TMPI_TeamToWorldRank(myRankInTeam, reloadTeam), 0, TMPI_GetWorldComm(), &status);
    std::cout << "\n-------------"
              << "I, FAILED TEAM, HAVE RECEIVED YOUR MESSAGE, I SHALL NOW RELOAD !\n"
              << " -> my team : " << myTeam << "\n"
              << " -> my rank in team : " << myRankInTeam << "\n"
              << " -> received message from team : " << reloadTeam << "\n"
              << " -> received message from world rank : " << TMPI_TeamToWorldRank(myRankInTeam, reloadTeam)
              << "\n-------------" << std::endl;

    std::printf("Rank %i of Team %i loading checkpoint from Team %i.\n", myRankInTeam, myTeam, reloadTeam);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* TODO on this naiv soft error checking we have only one simulationBlock ! */
    //simulationBlocks.clear();
//------------------------------------------------------------------------------
    //delete scenario;
    //scenario = nullptr;
    //restartNameInput = backupNameInput + "_t" + std::to_string(reloadTeam);
    //hasRecovered = true;
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
    args.addOption("output-basepath", 'o', "Output base file name");
//TODO remove this after debugging    args.addOption("backup-basepath", 'b', "Output base file name");
//TODO remove this after debugging    args.addOption("restart-basepath", 'r', "Restart base file name", args.Required, false);
    args.addOption("write-output", 'w', "Write output using netcdf writer to the specified output file", args.No, false);
    args.addOption("heartbeat-interval", 'i', "Wall-clock time in seconds to wait between heartbeats", args.Required, true);
    args.addOption("hash-method", 'm', "Which hashing method to use: ( 0=NONE | 1=stdhash ), default: 1", args.Required, false);
    args.addOption("hash-count", 'c', "Number of total hashes to send to the replica", args.Required, true);
    /* TODO inject-bitflip : make this option random (time) */
    args.addOption("inject-bitflip", 'f', "Injects a bit-flip to the first rank right after the simulation time reaches the given time", args.Required, false);
    args.addOption("kill-rank", 'k', "Kills the rank 0 of team 0 at the specified simulation time", args.Required, false);
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
    double heartbeatInterval;

    int nxRequested;
    int nyRequested;

    bool writeOutput = false;

    int hashOption;
    unsigned int numberOfHashes;

    double bitflip_at = -1.f;

    double kill_rank = -1.f;

    bool verbose = false;

    /* Read in command line arguments */
    simulationDuration = args.getArgument<float>("simulation-duration");
    nxRequested = args.getArgument<int>("resolution-x");
    nyRequested = args.getArgument<int>("resolution-y");
    outputNameInput = args.getArgument<std::string>("output-basepath");
    writeOutput = args.isSet("write-output");
    heartbeatInterval = args.getArgument<double>("heartbeat-interval");
    hashOption = args.getArgument<int>("hash-method", 1);
    if (hashOption != 0 && hashOption != 1 && hashOption != 2) {
        std::cout << "Invalid hash method. It has to be either:\n"
                  << "  0,\n"
                  << "  i^420 or\n"
                  << "  the number of which engineers think it's e"
                  << std::endl;
        return 1;
    }
    numberOfHashes = args.getArgument<unsigned int>("hash-count");
    if (args.isSet("inject-bitflip")) bitflip_at = args.getArgument<double>("inject-bitflip");
    if (args.isSet("kill-rank")) kill_rank = args.getArgument<double>("kill-rank");
    verbose = args.isSet("verbose");

    /* fixed backup name */
    //backupNameInput = args.getArgument<std::string>("backup-basepath");
    backupNameInput = "BACKUP_" + outputNameInput;

    /* Compute when the hash intervals are reached */
    float *sendHashAt = new float[numberOfHashes];
    /* Time delta between sending hashes */
    float sendHashDelta = simulationDuration / numberOfHashes;
    /* The first hash is sent after 0 + delta t */
    sendHashAt[0] = sendHashDelta;
    for (unsigned int i = 1; i < numberOfHashes; i++) {
        sendHashAt[i] = sendHashAt[i - 1] + sendHashDelta;
   }


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
        MPI_Init(&argc, &argv);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &ranksPerTeam);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    const int myTeam{TMPI_GetTeamNumber()};
    int numTeams = TMPI_GetInterTeamCommSize();
    assert(numTeams == 2);
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //unsigned int blocksPerRank = numTeams * decompFactor;

    // TODO check
    /* Since we are doing the first option naively, we will not share any
     * blocks. Which means it is best that we have only one block per replica
     * ==> TODO REMOVE the variables we don't need after debugging */
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

    float simulationStart = 0.f;

    int localBlockPositionX = startPoint / blockCountY;
    int localBlockPositionY = startPoint % blockCountY;

    // if not loading from a checkpoint, TODO debug this
    if (restartNameInput == "")
    {
        scenario = new SWE_RadialBathymetryDamBreakScenario{};
        int widthScenario = scenario->getBoundaryPos(BND_RIGHT) - scenario->getBoundaryPos(BND_LEFT);
        int heightScenario = scenario->getBoundaryPos(BND_TOP) - scenario->getBoundaryPos(BND_BOTTOM);

        float dxSimulation = (float)widthScenario / nxRequested;
        float dySimulation = (float)heightScenario / nyRequested;

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

        backupMetadataName = backupTeamPosName + "_metadata";

        simulationBlock = std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>(
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
                                                          false));

        std::array<int, 4> myNeighbours =
            getNeighbours(localBlockPositionX, localBlockPositionY, blockCountX, blockCountY, startPoint);

        int refinedNeighbours[4];
        int realNeighbours[4];
        std::array<std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>, 4> neighbourBlocks;
        std::array<BoundaryType, 4> boundaries;

        /* For all my neighbours. */
        for (int j = 0; j < 4; j++) {
            if (myNeighbours[j] >= startPoint && myNeighbours[j] < (startPoint + blocksPerRank))
            {
                refinedNeighbours[j] = -2;
                realNeighbours[j] = myNeighbours[j];
                neighbourBlocks[j] = simulationBlock;
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

        simulationBlock->initScenario(*scenario, boundaries.data());
        simulationBlock->connectNeighbourLocalities(refinedNeighbours);
        simulationBlock->connectNeighbours(realNeighbours);
        simulationBlock->connectLocalNeighbours(neighbourBlocks);
        simulationBlock->setRank(startPoint);
        simulationBlock->setDuration(simulationDuration);

        /* For checkpoints. TODO check if we need this for hard failure resilience */
        simulationBlock->writer->initMetadataFile(
            backupMetadataName,
            simulationDuration,
            ranksPerTeam,
            simulationBlock->nx,
            simulationBlock->ny,
            0,
            std::vector<BoundaryType>(boundaries.begin(), boundaries.end()),
                {simulationBlock->originX,
                simulationBlock->originX + (simulationBlock->nx) * simulationBlock->dx,
                simulationBlock->originY,
                simulationBlock->originY + (simulationBlock->ny) * simulationBlock->dy});
    }
    else /* loading from a checkpoint */
    {
        std::vector<SWE_Scenario*> scenarios{};

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

        /* TODO check this! here we use meta suffix instead of metadata suffix
         *      because this is reloaded. But why?
         */
        backupMetadataName = backupTeamPosName + "_meta";
        simulationBlock = std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>(
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
                                                      true));

        std::array<int, 4> myNeighbours =
            getNeighbours(localBlockPositionX, localBlockPositionY, blockCountX, blockCountY, startPoint);

        /* only one scenario, no task sharing */
        scenario = scenarios[0];

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
                neighbourBlocks[j] = simulationBlock;
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

        simulationBlock->initScenario(*scenario, boundaries.data());
        simulationBlock->connectNeighbourLocalities(refinedNeighbours);
        simulationBlock->connectNeighbours(realNeighbours);
        simulationBlock->connectLocalNeighbours(neighbourBlocks);
        simulationBlock->setRank(startPoint);
        simulationBlock->setDuration(simulationDuration);
        simulationBlock->writer->initMetadataFile(
            backupMetadataName,
            simulationDuration,
            ranksPerTeam,
            simulationBlock->nx,
            simulationBlock->ny,
            0,
            std::vector<BoundaryType>(boundaries.begin(), boundaries.end()),
                {simulationBlock->originX,
                simulationBlock->originX + (simulationBlock->nx) * simulationBlock->dx,
                simulationBlock->originY,
                simulationBlock->originY + (simulationBlock->ny) * simulationBlock->dy});
        delete scenario;
    } /* end of else */

    simulationBlock->sendBathymetry();
    simulationBlock->recvBathymetry();

    std::vector<float> timesteps;
    // Simulated time
    t = simulationStart;

    float timestep;

    /* write zero timestep, if not restoring */
    if (writeOutput && simulationStart == 0.f) {
        simulationBlock->writeTimestep(0.f);
    }

//------------------------------------------------------------------------------
//
    /* For debugging, only counts the heartbeats with the
     * hashes to make sure the number is right
     */
    int heartbeatCounter = 0;

    // Size of the update fields incl. ghost layer
    const int fieldSizeX = (simulationBlock->nx + 2) * (simulationBlock->ny + 2);
    const int fieldSizeY = (simulationBlock->nx + 1) * (simulationBlock->ny + 2);

    /* for hashing the calculated updates to detect silent data corruptions */
    tools::Hasher swe_hasher = tools::Hasher(fieldSizeX, fieldSizeY, simulationBlock.get());


    /* start the very first hearbeat */
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

    double timeOfLastHeartbeat = MPI_Wtime();
    double timeSinceLastHeartbeat = 0.f;

    if (verbose) {
        /* print heartbeat for debugging */
        std::cout << "\n\t++++++++++++++ "
                  << "FIRST HEARTBEAT START ++++++++++++++\n\t"
                  << "  - team:\t\t" << myTeam << "\n\t"
                  << "  - rank:\t\t" << myRankInTeam << "\n\t"
                  << "at t = " << t << ", at MPI_Wtime() : "
                  << timeOfLastHeartbeat << std::endl;
    }

    for (unsigned int i = 0; i < numberOfHashes; i++) {

        /* Simulate until the next sending is reached */
        while (t < sendHashAt[i]) {

            if (timeSinceLastHeartbeat >= heartbeatInterval) {

                /* end Heartbeat */
                MPI_Sendrecv(MPI_IN_PLACE,              /* Send buffer      */
                             0,                         /* Send count       */
                             MPI_BYTE,                  /* Send type        */
                             MPI_PROC_NULL,             /* Destination      */
                             -1,                        /* Send tag         */
                             MPI_IN_PLACE,              /* Receive buffer   */
                             0,                         /* Receive count    */
                             MPI_BYTE,                  /* Receive type     */
                             MPI_PROC_NULL,             /* Source           */
                             0,                         /* Receive tag      */
                             MPI_COMM_SELF,             /* Communicator     */
                             MPI_STATUS_IGNORE);        /* Status object    */

                if (verbose) {
                    /* print heartbeat for debugging */
                    std::cout << "\n\t+++++++++ "
                              << "HEARTBEAT END ++++++++\n\t"
                              << "  - team:\t\t" << myTeam << "\n\t"
                              << "  - rank:\t\t" << myRankInTeam << "\n\t"
                              << " at t = " << t << ", at MPI_Wtime() : "
                              << timeOfLastHeartbeat << std::endl;
                }


                /* start Hearbeat */
                /* If we post a heartbeat after every time step, a single rank
                 * would delay the heartbeats of neighouring ranks because
                 * point-to-point messages (like receiveGhostLayer) are required
                 * before a new simulation time step can be started.
                 */
                MPI_Sendrecv(MPI_IN_PLACE,              /* Send buffer      */
                            0,                          /* Send count       */
                            MPI_BYTE,                   /* Send type        */
                            MPI_PROC_NULL,              /* Destination      */
                            1,                          /* Send tag         */
                            MPI_IN_PLACE,               /* Receive buffer   */
                            0,                          /* Receive count    */
                            MPI_BYTE,                   /* Receive type     */
                            MPI_PROC_NULL,              /* Source           */
                            0,                          /* Receive tag      */
                            MPI_COMM_SELF,              /* Communicator     */
                            MPI_STATUS_IGNORE);         /* Status object    */

                timeOfLastHeartbeat = MPI_Wtime();
                timeSinceLastHeartbeat = 0.f;

                if (verbose) {
                    /* print heartbeat for debugging */
                    std::cout << "\n\t++++++++++++++ "
                              << "HEARTBEAT START ++++++++++++++\n\t"
                              << "  - team:\t\t" << myTeam << "\n\t"
                              << "  - rank:\t\t" << myRankInTeam << "\n\t"
                              << "at t = " << t << ", at MPI_Wtime() : "
                              << timeOfLastHeartbeat << std::endl;
                }
            }


            // exchange boundaries between blocks
            simulationBlock->setGhostLayer();
            simulationBlock->receiveGhostLayer();

            if(verbose) {
                std::cout << "calculating.. t = " << t
                          << "\t\t\t\t---- TEAM " << myTeam << " Rank "
                          << myRankInTeam << std::endl;
            }

            /* compute current updates */
            simulationBlock->computeNumericalFluxes();

            /* Inject a bitflip at team 0 at rank 0 */
            if (bitflip_at >= 0  && t > bitflip_at && myTeam == 0 && myRankInTeam == 0) {
                simulationBlock->injectRandomBitflip();
                /* prevent any other bitflip */
                bitflip_at = -1.f;
            }

            /* update the hash */
            if (hashOption == 0) {
                /* don't hash. 0 is for 'no hashing' */
            }
            else if (hashOption == 1) {
                swe_hasher.update_stdHash();
            }
            else {
                std::cout << "Unknown hash method.. something went wrong\n"
                          << std::endl;
                MPI_Abort(TMPI_GetWorldComm(), MPI_ERR_UNKNOWN);
            }

            /* Agree on a timestep */
            timestep = simulationBlock->maxTimestep;
            float agreed_timestep;
            PMPI_Allreduce(&timestep, &agreed_timestep, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

            simulationBlock->maxTimestep = agreed_timestep;
            simulationBlock->updateUnknowns(agreed_timestep);

            t += agreed_timestep;

            if(writeOutput) {
                if(verbose) {
                    std::cout << "--> Writing timestep at: " << t
                              << " by team " << myTeam << ", rank " << myRankInTeam
                              << std::endl;
                }
                simulationBlock->writeTimestep(t);
            }

            /* update the elapsed time after sending the heartbeat */
            timeSinceLastHeartbeat = MPI_Wtime() - timeOfLastHeartbeat;

            /* TODO Why send bcast to the rank's team's 0 here, different heartbeats bad ? */
            MPI_Bcast(&timeSinceLastHeartbeat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            // MPI_Barrier(TMPI_GetInterTeamComm());

            /* kill the rank 0 if specified parameter is set */
            if (kill_rank >= 0 && t >= kill_rank &&
                myTeam == 0 && myRankInTeam == 0 && restartNameInput == "") {
                std::cout << "TEAM:0 Rank:0 RETURNS ERROR CODE for hard failure simulation...."
                          << std::endl;
                return 1;
                //std::cout << "ETERNAL SLEEP TIME FOR THE TEAM:1 Rank:0 for hard failure simulation...."
                          //<< std::endl;
                //while(true) sleep(10); TODO also check if this is working
            }

        } // end of while(t < sendHashAt[i])

        // send hash with heartbeat
        /* Finalize hash computation */
        if (hashOption == 0) {
            /* don't hash. 0 is for 'no hashing' */
        }
        else if (hashOption == 1) {

            size_t total_hash = swe_hasher.finalize_stdHash();

            if(verbose) {
                /* print heartbeat for debugging */
                std::cout << "\n\t+++++++++ " << heartbeatCounter << ". "
                          << "HEARTBEAT HASH ++++++++\n\t"
                          << "  - team:\t\t" << myTeam << "\n\t"
                          << "  - rank:\t\t" << myRankInTeam << "\n\t"
                          << "  - size_t total_hash (" << sizeof(total_hash)
                          << " bytes)" << " : " << total_hash
                          << " at t = " << t << std::endl;
            }

            /* single heartbeat with hash */
            MPI_Sendrecv(&total_hash,               /* Send buffer      */
                         1,                         /* Send count       */
                         MPI_SIZE_T,                /* Send type        */
                         MPI_PROC_NULL,             /* Destination      */
                         0,                         /* Send tag         */
                         MPI_IN_PLACE,              /* Receive buffer   */
                         0,                         /* Receive count    */
                         MPI_BYTE,                  /* Receive type     */
                         MPI_PROC_NULL,             /* Source           */
                         0,                         /* Receive tag      */
                         MPI_COMM_SELF,             /* Communicator     */
                         MPI_STATUS_IGNORE);        /* Status object    */
        }
        else {
            std::cout << "Unknown hash method.. something went wrong\n"
                      << std::endl;
            MPI_Abort(TMPI_GetWorldComm(), MPI_ERR_UNKNOWN);
        }

        heartbeatCounter++;

    } // for (unsigned int i = 0; i < numberOfHashes; i++)

    /* end the very last Heartbeat */
    MPI_Sendrecv(MPI_IN_PLACE,              /* Send buffer      */
                 0,                         /* Send count       */
                 MPI_BYTE,                  /* Send type        */
                 MPI_PROC_NULL,             /* Destination      */
                 -1,                        /* Send tag         */
                 MPI_IN_PLACE,              /* Receive buffer   */
                 0,                         /* Receive count    */
                 MPI_BYTE,                  /* Receive type     */
                 MPI_PROC_NULL,             /* Source           */
                 0,                         /* Receive tag      */
                 MPI_COMM_SELF,             /* Communicator     */
                 MPI_STATUS_IGNORE);        /* Status object    */

    if (verbose) {
        /* print heartbeat for debugging */
        std::cout << "\n\t+++++++++ "
        << "LAST HEARTBEAT END ++++++++\n\t"
        << "  - team:\t\t" << myTeam << "\n\t"
        << "  - rank:\t\t" << myRankInTeam << "\n\t"
        << " at t = " << t << ", at MPI_Wtime() : "
        << timeOfLastHeartbeat << std::endl;
    }


    delete[] sendHashAt;
    delete scenario;

    simulationBlock->freeMpiType();

    MPI_Finalize();
}
