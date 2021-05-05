/**
 * @file src/tolerance/swe_softRes.cpp
 *
 * @brief Soft error resilience only with detection using hashes
 *
 * This file tries to integrate hashes to the heartbeats used in
 * the tmpi library, which helps us to compare the results of the
 * replicas.
 *
 * TODO introduce tests
 * TODO fix hashing method with sha1
 *
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
    args.addOption("hash-method", 'm', "Which hashing method to use: ( 0=NONE | 1=stdhash | 2=SHA1 ), default: 1", args.Required, false);
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
    verbose = args.isSet("verbose");

    /* Simulation time */
    float t(0.0F);

    /* one block per rank TODO remove simulationBlocks variable */
    std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp> simulationBlock;

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

    // init MPI
    int myRankInTeam;
    int ranksPerTeam;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ranksPerTeam);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    const int myTeam{TMPI_GetTeamNumber()};
    int numTeams = TMPI_GetInterTeamCommSize();
    assert(numTeams == 2);

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

    /* TODO backups can be removed in this version */
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

    // Size of the update fields incl. ghost layer
    const int fieldSizeX{(simulationBlock->nx + 2) * (simulationBlock->ny + 2)};
    const int fieldSizeY{(simulationBlock->nx + 1) * (simulationBlock->ny + 2)};

    /* for hashing the calculated updates to detect silent data corruptions */
    Hasher swe_hasher = Hasher(
      fieldSizeX, fieldSizeY,
      simulationBlock->hNetUpdatesLeft.getRawPointer(),
      simulationBlock->hNetUpdatesRight.getRawPointer(),
      simulationBlock->huNetUpdatesLeft.getRawPointer(),
      simulationBlock->huNetUpdatesRight.getRawPointer(),
      simulationBlock->hNetUpdatesBelow.getRawPointer(),
      simulationBlock->hNetUpdatesAbove.getRawPointer(),
      simulationBlock->hvNetUpdatesBelow.getRawPointer(),
      simulationBlock->hvNetUpdatesAbove.getRawPointer(),
      &(simulationBlock->maxTimestep));


    for (unsigned int i = 0; i < numberOfHashes; i++) {

        /* Simulate until the next sending is reached */
        while (t < sendHashAt[i]) {

            // exchange boundaries between blocks
            simulationBlock->setGhostLayer();
            simulationBlock->receiveGhostLayer();

            auto& currentBlock = *simulationBlock;

            if(verbose) {
                std::cout << "calculating.. t = " << t
                          << "\t\t\t\t---- TEAM " << myTeam << " Rank "
                          << myRankInTeam << std::endl;
            }

            /* compute current updates */
            currentBlock.computeNumericalFluxes();

            /* Inject a bitflip at team 0 at rank 0 */
            if (bitflip_at >= 0  && t > bitflip_at && myTeam == 0 && myRankInTeam == 0) {

                /* index of the float we want to corrupt */
                size_t flipAt_float = fieldSizeX / 2;

                float *calculated_huNetUpdatesLeft = simulationBlock->huNetUpdatesLeft.getRawPointer();

                std::cout << "\n............Injecting..a..bit..flip.................\n"
                          << "old value : " <<     calculated_huNetUpdatesLeft[flipAt_float]
                          << "\n...............DATA..CORRUPTED......................\n"
                          << "\n";

                /* flip only the first bit with the XOR operation */
                ((unsigned int *)calculated_huNetUpdatesLeft)[flipAt_float] ^= 0x80000000;

                std::cout << "new value : " << calculated_huNetUpdatesLeft[flipAt_float]
                          << "\n"
                          << std::endl;

                /* prevent any other bitflip */
                bitflip_at = -1.f;
            }

            /* compute the hash */
            if (hashOption == 0) {
                /* don't hash. 0 is for 'no hashing' */
            }
            else if (hashOption == 1) {
                swe_hasher.update_stdHash();
            }
            else if (hashOption == 2) {
                swe_hasher.update_SHA1();
            }
            else {
                std::cout << "Unknown hash method.. something went wrong\n"
                          << std::endl;
                MPI_Abort(TMPI_GetWorldComm(), MPI_ERR_UNKNOWN);
            }

            /* Agree on a timestep */
            timestep = simulationBlock->maxTimestep;
            float agreed_timestep;
            MPI_Allreduce(&timestep, &agreed_timestep, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

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
        } // end of t < sendHashAt[i]

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
        else if (hashOption == 2) {

            /* get 20 bytes from SHA1 hasher
             * SHA1 has a fixed output size of 120 bits == 20 bytes */
            unsigned char *total_hash = swe_hasher.finalize_SHA1();
            std::cout << "\nsorry.. but SHA1 hashing is still in development :(\n"
                      << std::endl;

            assert(false);

            if(verbose) {
                /* print heartbeat for debugging */
                std::cout << "\n\t+++++++++ " << heartbeatCounter << ". "
                          << "HEARTBEAT HASH ++++++++\n\t"
                          << "  - team:\t\t" << myTeam << "\n\t"
                          << "  - rank:\t\t" << myRankInTeam << "\n\t"
                          << "  - unsigned char total_hash (" << 20
                          << " bytes)" << " : " << total_hash // TODO print only 20 bytes here
                          << " at t = " << t << std::endl;
            }

            /* single heartbeat with hash */
            MPI_Sendrecv(total_hash,                /* Send buffer      */
                         20,                        /* Send count       */
                         MPI_CHAR,                  /* Send type        */
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

    simulationBlock->freeMpiType();
    MPI_Finalize();
    return 0;
}
