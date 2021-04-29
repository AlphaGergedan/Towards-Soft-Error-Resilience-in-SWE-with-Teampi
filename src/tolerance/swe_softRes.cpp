/**
 * @file src/tolerance/swe_softRes_and_hardRes_woutTaskSharing.cpp
 *
 * @brief hard error resiliency with soft error detection without task sharing
 *
 * TODO
 *   - Implement compare buffer with replicas
 *
 *   - Introduce random bitflips, tests
 *
 */



#include <bits/c++config.h>
#include <bitset>
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
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//std::vector<std::string> backupMetadataNames{};


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
    args.addOption("write-output", 'w', "Write output using netcdf writer to the specified output file", args.No, false);
    args.addOption("hash-method", 'm', "Which hashing method to use: (1 | 2), default: 1", args.Required, false);
    args.addOption("hash-count", 'c', "Number of total hashes to send to the replica", args.Required, true);
    args.addOption("inject-bitflip", 'f', "Injects a bit-flip to the first rank right after the simulation time reaches the given time", args.Required, false);
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

    int nxRequested;
    int nyRequested;

    std::string outputNameInput = "";

    bool writeOutput = false;

    int hashOption;
    unsigned int numberOfHashes;

    double bitflip_at = -1.f;

    bool verbose = false;

    /* Read in command line arguments */
    simulationDuration = args.getArgument<float>("simulation-duration");
    nxRequested = args.getArgument<int>("resolution-x");
    nyRequested = args.getArgument<int>("resolution-y");
    outputNameInput = args.getArgument<std::string>("output-basepath");
    writeOutput = args.isSet("write-output");
    hashOption = args.getArgument<int>("hash-method", 1);
    if (hashOption != 1 && hashOption != 2) {
        std::cout << "Invalid hash method. It has to be either:\n"
                  << "  i^420 or\n"
                  << "  the number of which engineers think it's e"
                  << std::endl;
        return 1;
    }
    numberOfHashes = args.getArgument<unsigned int>("hash-count");
    bitflip_at = args.getArgument<double>("inject-bitflip");
    verbose = args.isSet("verbose");

    /* Simulation time */
    float t(0.0F);

    /* one block per rank TODO remove simulationBlocks variable */
    std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp> simulationBlock;

    clock_t heartbeatInterval;

    std::string outputTeamName = "";

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
    //TMPI_SetErrorHandlingStrategy(TMPI_WarmSpareErrorHandler);
    // no error strategy !

    /*
     * we don't set any error handling strategies for now,
     * TMPI_NoErrorHandler is set by default
     */


    // init MPI
    int myRankInTeam;
    int ranksPerTeam;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ranksPerTeam);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    const int myTeam{TMPI_GetTeamNumber()};
    int numTeams = TMPI_GetInterTeamCommSize();

    // TODO check
    /* Since we are doing the first option naively, we will not share any
     * blocks. Which means it is best that we have only one block per replica */
    unsigned int blocksPerRank = 1;

    outputTeamName = outputNameInput + "_t" + std::to_string(myTeam);

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

    // not loading from a checkpoint
    SWE_Scenario *scenario = new SWE_RadialBathymetryDamBreakScenario{};
    int widthScenario = scenario->getBoundaryPos(BND_RIGHT) - scenario->getBoundaryPos(BND_LEFT);
    int heightScenario = scenario->getBoundaryPos(BND_TOP) - scenario->getBoundaryPos(BND_BOTTOM);

    float dxSimulation = (float)widthScenario / nxRequested;
    float dySimulation = (float)heightScenario / nyRequested;

    int localBlockPositionX = startPoint / blockCountY;
    int localBlockPositionY = startPoint % blockCountY;

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
    std::string backupTeamPosName = genTeamPosName("BACKUP_" + outputTeamName, localBlockPositionX, localBlockPositionY);

    /* Add the block to be calculated by this rank. */
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

    /* TODO on this naiv soft error checking we have only one simulationBlock ! */
    /* Loop over the rank's block number in its team. */
    for (unsigned int currentBlockNr = startPoint; currentBlockNr < startPoint + blocksPerRank; currentBlockNr++) {
        assert(currentBlockNr == startPoint); // TODO only one block

        //backupMetadataNames.push_back(backupTeamPosName + "_metadata"); TODO

    }


    /* TODO on this naiv soft error checking we have only one simulationBlock ! */
    for (unsigned int currentBlockNr = startPoint; currentBlockNr < startPoint + blocksPerRank; currentBlockNr++) {
        int localBlockPositionX = currentBlockNr / blockCountY;
        int localBlockPositionY = currentBlockNr % blockCountY;
        std::array<int, 4> myNeighbours =
            getNeighbours(localBlockPositionX, localBlockPositionY, blockCountX, blockCountY, currentBlockNr);

        int refinedNeighbours[4];
        int realNeighbours[4];
        std::array<std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>, 4> neighbourBlocks;
        std::array<BoundaryType, 4> boundaries;

        /* For all my neighbours. */
        for (int j = 0; j < 4; j++) {
            if (myNeighbours[j] >= startPoint && myNeighbours[j] < (startPoint + blocksPerRank)) {
                refinedNeighbours[j] = -2;
                realNeighbours[j] = myNeighbours[j];
                neighbourBlocks[j] = simulationBlock;
                boundaries[j] = CONNECT_WITHIN_RANK;
            }
            else if (myNeighbours[j] == -1) {
                boundaries[j] = scenario->getBoundaryType((Boundary)j);
                refinedNeighbours[j] = -1;
                realNeighbours[j] = -1;
            }
            else {
                realNeighbours[j] = myNeighbours[j];
                refinedNeighbours[j] = myNeighbours[j] / blocksPerRank;
                boundaries[j] = CONNECT;
            }
        }

        simulationBlock->initScenario(*scenario, boundaries.data());
        simulationBlock->connectNeighbourLocalities(refinedNeighbours);
        simulationBlock->connectNeighbours(realNeighbours);
        simulationBlock->connectLocalNeighbours(neighbourBlocks);
        simulationBlock->setRank(currentBlockNr);
        simulationBlock->setDuration(simulationDuration);
    }


    /* on this soft error detection we have only one simulationBlock ! */
    simulationBlock->sendBathymetry();
    simulationBlock->recvBathymetry();

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

//------------------------------------------------------------------------------

    // Write zero timestep
    if (writeOutput) {
        simulationBlock->writeTimestep(0.f);
    }
//------------------------------------------------------------------------------

    unsigned int heartbeatCounter = 0;

    for (unsigned int i = 0; i < numberOfHashes; i++) {

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



        /* Simulate until the next sending is reached */
        while (t < sendHashAt[i]) {

            // exchange boundaries between blocks
            simulationBlock->setGhostLayer();
            simulationBlock->receiveGhostLayer();

            auto& currentBlock = *simulationBlock;

            // Size of the update fields incl. ghost layer
            const int fieldSizeX{(currentBlock.nx + 2) * (currentBlock.ny + 2)};
            const int fieldSizeY{(currentBlock.nx + 1) * (currentBlock.ny + 2)};

            if(verbose) {
                std::cout << "calculating.. t = " << t
                          << "\t\t\t\t---- TEAM " << myTeam << " Rank "
                          << myRankInTeam << std::endl;
            }

            currentBlock.computeNumericalFluxes();

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// MOVE THIS SECTION TO A HEADER FILE
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

            /* Inject a bitflip at team 0 at rank 0 */
            if (bitflip_at >= 0  && t > bitflip_at && myTeam == 0 && myRankInTeam == 0) {

                /* index of the float we want to corrupt */
                size_t flipAt_float = fieldSizeX / 2;

                std::cout << "\n............Injecting_a_bit_flip.................\n";
                std::cout << "old value : " << calculated_huNetUpdatesLeft[flipAt_float]
                          << "\n";

                /* flip only the first bit with the XOR operation */
                ((unsigned int *)calculated_huNetUpdatesLeft)[flipAt_float] ^= 0x80000000;

                std::cout << "new value : " << calculated_huNetUpdatesLeft[flipAt_float]
                          << "\n"
                          << std::endl;

                /* prevent any other bitflip */
                bitflip_at = -1.f;
            }


            /* Convert them to strings */
            std::string str_hLeft((const char*) calculated_hNetUpdatesLeft, sizeof(float) * fieldSizeX);
            std::string str_hRight((const char*) calculated_hNetUpdatesRight, sizeof(float) * fieldSizeX);

            std::string str_huLeft((const char*) calculated_huNetUpdatesLeft, sizeof(float) * fieldSizeX);
            std::string str_huRight((const char*) calculated_huNetUpdatesRight, sizeof(float) * fieldSizeX);

            std::string str_hBelow((const char*) calculated_hNetUpdatesBelow, sizeof(float) * fieldSizeY);
            std::string str_hAbove((const char*) calculated_hNetUpdatesAbove, sizeof(float) * fieldSizeY);

            std::string str_hvBelow((const char*) calculated_hvNetUpdatesBelow, sizeof(float) * fieldSizeY);
            std::string str_hvAbove((const char*) calculated_hvNetUpdatesAbove, sizeof(float) * fieldSizeY);

            std::string str_maxTimeStep((const char*) calculated_maxTimeStep, sizeof(float));

            if (hashOption == 2) {
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
            }
            else if (hashOption == 1) {
                /*
                    OPTION B: hash them using std::hash
                    HINT: TODO xor is used when combining the hash but this
                    might not be ideal, see boost::hash_combine
                */

                total_hash ^= hash_fn(str_hLeft);
                total_hash ^= hash_fn(str_hRight);

                total_hash ^= hash_fn(str_huLeft);
                total_hash ^= hash_fn(str_huRight);

                total_hash ^= hash_fn(str_hBelow);
                total_hash ^= hash_fn(str_hAbove);

                total_hash ^= hash_fn(str_hvBelow);
                total_hash ^= hash_fn(str_hvAbove);

                total_hash ^= hash_fn(str_maxTimeStep);
            }
            else {
                std::cout << "Unknown hash method.. something wrong\n" << std::endl;
                MPI_Abort(TMPI_GetWorldComm(), MPI_ERR_UNKNOWN);
            }



// MOVE THIS SECTION TO A HEADER FILE
//------------------------------------------------------------------------------

            /* Agree on a timestep */
            timestep = simulationBlock->maxTimestep;
            float agreed_timestep;
            MPI_Allreduce(&timestep, &agreed_timestep, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

            // TODO only one block
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
        }

        // TODO send hashes
        // End Heartbeat send the computed task to compare hashes
        //
        // std::cout << "Team " << myTeam << ": HEARTBEAT! Current simulation time is " << t << '\n';
        std::string final_checksum = "";
        if (hashOption == 2) {
            /* end the checksum SHA1*/
            final_checksum = checksum.final();
        }

        if(verbose) {
            /* print heartbeat for debugging */
            std::cout << "\n\t+++++++++ " << heartbeatCounter << ". "
                      << "HEARTBEAT HASH ++++++++\n\t"
                      << "  - team:\t\t" << myTeam << "\n\t"
                      << "  - rank:\t\t" << myRankInTeam << "\n\t"
                      << "  - size_t total_hash (" << sizeof(total_hash)
                      << " bytes)" << " : " << total_hash
                      << " ; SHA1 checksum (" << final_checksum.size()
                      << " bytes)" << " : " << final_checksum
                      << " at t = " << t << std::endl;
        }

        heartbeatCounter++;

        /* Send the hash calculated by std::hash and XORed with size size_t */
        if (hashOption == 1) {

            /* heartbeat end : OPTION B */
            MPI_Sendrecv(&total_hash,              /* Send buffer      */
                        1,                         /* Send count       */
                        MPI_SIZE_T,                /* Send type        */
                        MPI_PROC_NULL,             /* Destination      */
                        0,                         /* Send tag         */               // TODO CHECK IF 0 TAG IS WORKING
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
        else if (hashOption == 2) {

            /* heartbeat end : OPTION A */
            MPI_Sendrecv(final_checksum.c_str(),   /* Send buffer      */
                        final_checksum.size(),     /* Send count: 20 ? */
                        MPI_CHAR,                  /* Send type        */
                        MPI_PROC_NULL,             /* Destination      */
                        0,                         /* Send tag         */               // TODO CHECK IF 0 TAG IS WORKING
                        MPI_IN_PLACE,              /* Receive buffer   */
                        0,                         /* Receive count    */
                        MPI_BYTE,                  /* Receive type     */
                        MPI_PROC_NULL,             /* Source           */
                        0,                         /* Receive tag      */
                        MPI_COMM_SELF,             /* Communicator     */
                        MPI_STATUS_IGNORE);        /* Status object    */
        }
        else {
            std::cout << "Unknown hash method.. something wrong\n" << std::endl;
            MPI_Abort(TMPI_GetWorldComm(), MPI_ERR_UNKNOWN);
        }

    }


    simulationBlock->freeMpiType();
    MPI_Finalize();
    return 0;
}
