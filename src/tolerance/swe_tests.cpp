/**
 * @file src/tolerance/swe_tests.cpp
 *
 * @brief Tests for fault tolerant swe applications
 *
 * TODO add tests
 * TODO add more tests
 * TODO description
 *
 * @author Atamert Rahma
 */


#include <functional>
#include <memory>
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

/* -- GLOBAL DECLARATIONS -- */
bool SDC_injected = false;
size_t testsPassed = 0;

std::string lightShade = "\u2591\u2591\u2591\u2591\u2591\u2591\u2591\u2591\u2591\u2591\u2591\u2591\u2591\u2591\u2591";
std::string doubleVertical = "\u2551";
std::string doubleHorizontal = "\u2550";
std::string markPassed = "\u2714";
std::string markFailed = "\u2716";
std::string markRightArrow = "\u279c";
std::string color_blue = "\033[1;34m";
std::string color_green = "\033[1;32m";
std::string color_red = "\033[1;31m";
std::string color_norm = "\033[0m";

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

/* TODO currently in development ... only difference is that hasRecovered flag is not set !!! */
void loadCheckpointCallback2_softErrorRecovery(int reloadTeam) {
    int myRankInTeam;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    const int myTeam = TMPI_GetTeamNumber();
    std::printf("Rank %i of Team %i loading checkpoint from Team %i.\n", myRankInTeam, myTeam, reloadTeam);
    MPI_Status status;
    int recv_buf;
    MPI_Recv(&recv_buf, 1, MPI_INT, TMPI_TeamToWorldRank(myRankInTeam, reloadTeam), 0, TMPI_GetWorldComm(), &status);
    simulationBlocks.clear();
    delete scenario;
    scenario = nullptr;
    restartNameInput = "TEST_" + std::to_string(testsPassed+1) + "_" +
                     backupNameInput + "_t" + std::to_string(reloadTeam);
    //hasRecovered = true;
    recoveredFromSDC = true;
    longjmp(jumpBuffer, 1);
}
void createCheckpointCallback(std::vector<int> failedTeams) {
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


std::array<int, 4> getNeighbours(int localBlockPositionX,
                                 int localBlockPositionY,
                                 int blockCountX,
                                 int blockCountY,
                                 int myRank) {
    std::array<int, 4> myNeighbours;
    myNeighbours[BND_LEFT] = (localBlockPositionX > 0) ? myRank - blockCountY : -1;
    myNeighbours[BND_RIGHT] = (localBlockPositionX < blockCountX - 1) ? myRank + blockCountY : -1;
    myNeighbours[BND_BOTTOM] = (localBlockPositionY > 0) ? myRank - 1 : -1;
    myNeighbours[BND_TOP] = (localBlockPositionY < blockCountY - 1) ? myRank + 1 : -1;
    return myNeighbours;
}

/* Data sturcture for easier syntax on arguments passed by tests */
struct TestArguments {
    float simulationDuration;
    clock_t heartbeatInterval;
    int nxRequested;
    int nyRequested;
    unsigned int decompFactor;
    bool writeOutput;
    double bitflip_at;          /* simulation time to inject the bitflip */

    /* Constructor */
    TestArguments(float simulationDuration_, clock_t heartbeatInterval_,
                  int nxRequested_, int nyRequested_, unsigned int decompFactor_,
                  bool writeOutput_, double bitflip_at_) {
        simulationDuration = simulationDuration_;
        heartbeatInterval = heartbeatInterval_;
        nxRequested = nxRequested_;
        nyRequested = nyRequested_;
        decompFactor = decompFactor_;
        writeOutput = writeOutput_;
        bitflip_at = bitflip_at_;
    }
};


//------------------------------------------------------------------------------


/**
 * Runs the simulation hard coded, and injects
 * a desired bitflip at a random location.
 * Location can be either updates, data or
 * both.
 *
 * Information usually obtained over the
 * command line is hardcoded here, using
 * \f$nxRequested = nyRequested = 500\f$.
 *
 * @param simulationDuration
 * @param heartbeatInterval
 * @param nxRequested
 * @param nyRequested
 * @param decompFactor
 * @param writeOutput
 * @param bitflip_at
 * @param bitflipLocation 0 for updates, 1 for data, 2 for random
 * @param bitflipType 0 for NaN, 1 for negative water height
 * @param rankToCorrupt
 *
 * @return value of the corrupted number
 */
int run(TestArguments *args, int bitflipLocation, int bitflipType, int rankToCorrupt) {

    /* read the hard coded user arguments */
    float simulationDuration = args->simulationDuration;
    clock_t heartbeatInterval = args->heartbeatInterval;
    int nxRequested = args->nxRequested;
    int nyRequested = args->nyRequested;
    unsigned int decompFactor = args->decompFactor;
    bool writeOutput = args->writeOutput;
    double bitflip_at = args->bitflip_at;

    int myRankInTeam;
    int worldRank;
    MPI_Comm_size(MPI_COMM_WORLD, &ranksPerTeam);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    PMPI_Comm_rank(TMPI_GetWorldComm(), &worldRank);

    const int myTeam{TMPI_GetTeamNumber()};
    int numTeams = TMPI_GetInterTeamCommSize();
    unsigned int blocksPerRank = numTeams * decompFactor;

    outputTeamName.clear();
    outputTeamName = "TEST_" + std::to_string(testsPassed+1) + "_" +
                     outputNameInput + "_t" + std::to_string(myTeam);
    backupTeamName.clear();
    backupTeamName = "TEST_" + std::to_string(testsPassed+1) + "_" +
                     backupNameInput + "_t" + std::to_string(myTeam);
    /* in recovery, the team name is already added to restart name in load function */

    /* Begin the logger */
    tools::FtLogger ft_logger(myTeam, myRankInTeam);
    //ft_logger.ft_print_spawnStatus();

    int totalBlocks = blocksPerRank * ranksPerTeam;

    /* number of SWE-Blocks in x- and y-direction */
    int blockCountY = std::sqrt(totalBlocks);
    while (totalBlocks % blockCountY != 0) blockCountY--;
    int blockCountX = totalBlocks / blockCountY;

    int startPoint = myRankInTeam * blocksPerRank;

    float simulationStart{0.0f};

    /* if not loading from a checkpoint */
    if (restartNameInput == "") {
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
    else { /* loading from a checkpoint */
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
        //ft_logger.ft_print_HBstart(timeOfLastHeartbeat, t);

        // simulate until the checkpoint is reached
        while (timeSinceLastHeartbeat < heartbeatInterval && t < simulationDuration) {
            //std::cout << "T" << myTeam << "R" << myRankInTeam << " t = " << t << std::endl; TODO remove after debugging
            // exchange boundaries between blocks
            for (auto& currentBlock : simulationBlocks) { currentBlock->setGhostLayer(); }
            for (auto& currentBlock : simulationBlocks) { currentBlock->receiveGhostLayer(); }
            const MPI_Comm interTeamComm{TMPI_GetInterTeamComm()};
            if (!hasRecovered && !recoveredFromSDC)
            {
                // TODO debugging printing again
                //std::cout << "T" << myTeam << "R" << myRankInTeam << " Barrier in the comp loop.. " << std::endl;
                // Avoid overwriting an old send buffer before everyone reaches this point
                MPI_Barrier(interTeamComm);
                //std::cout << "T" << myTeam << "R" << myRankInTeam << " Barrier in the comp loop.. FINISHED" << std::endl;
            }

            /* 0 for no SDC, 1 for SDC */
            unsigned char reportFlag = 0;

            /* compute and check the primary blocks for SDCs
             * and set the reportFlag accordingly
             */
            for (unsigned int i = 0; i < decompFactor; i++) {
                // Get a reference to the current block
                const int& currentBlockNr{myBlockOrder[i]};
                auto& currentBlock = *simulationBlocks[currentBlockNr];

                currentBlock.computeNumericalFluxes();

                /* TODO inject bitflip if desired */
                if (!SDC_injected && bitflip_at >= 0  && t > bitflip_at && rankToCorrupt == worldRank) {
                    /* bitflip location at updates */
                    if (bitflipLocation == 0) {
                        /* inject NaN */
                        if (bitflipType == 0) {
                            currentBlock.injectNaN_intoUpdates();
                        }
                        /* negative water height */
                        else if (bitflipType == 1) {
                            currentBlock.injectNegativeWaterHeight_intoUpdates();
                        }
                        /* random injection */
                        else {
                            currentBlock.injectRandomBitflip_intoUpdates();
                        }
                    }
                    /* bitflip location at data */
                    else if (bitflipLocation == 1) {
                        /* inject NaN */
                        if (bitflipType == 0) {
                            currentBlock.injectNaN_intoData();
                        }
                        /* negative water height */
                        else if (bitflipType == 1) {
                            currentBlock.injectNegativeWaterHeight_intoData();
                        }
                        /* bathymetry change */
                        else if (bitflipType == 2) {
                            currentBlock.injectBathymetryChange_intoData();
                        }
                        /* random injection */
                        else {
                            currentBlock.injectRandomBitflip_intoData();
                        }
                    }
                    /* random location and random bitflip */
                    else {
                        currentBlock.injectRandomBitflip();
                    }

                    /* prevent any other bitflip */
                    bitflip_at = -1.f;
                }

                /* check for soft errors */
                int admissable = currentBlock.validateAdmissability(t);

                /* handle soft errors
                 * if SDC detected in updates only, then recompute */
                if (admissable == 0) {
                    //ft_logger.ft_SDC_notDetected();
                }
                else if (admissable == 1) {
                    currentBlock.computeNumericalFluxes();

                    /* check if SDC is still present */
                    if (currentBlock.validateAdmissability(t) != 0) {
                        std::cout << "-- TEAM " << myTeam << ", Rank " << myRankInTeam << " Warning : SDC detected but cannot be fixed..\n"
                                  << "             Warning the rest of my team and all my replicas"
                                  << std::endl;
                        reportFlag = 1;
                    } else {
                        std::cout << "-- TEAM " << myTeam << ", Rank " << myRankInTeam << " Warning : SDC detected but it is fixed.."
                                  << std::endl;
                        reportFlag = 1;
                    }
                }
                // soft error in data, turn on the recovery mode
                else if (admissable == 2) {
                    std::cout << "-- TEAM " << myTeam << ", Rank " << myRankInTeam << " Warning : SDC detected but cannot be fixed..\n"
                              << "             Warning the rest of my team and all my replicas"
                              << std::endl;
                    reportFlag = 1;
                }
                else {
                    std::cout << "Unknown error (admissability)" << std::endl;
                    assert(false);
                }

                //if (verbose) ft_logger.ft_block_calculatingTask(currentBlockNr, currentBlock.maxTimestep);
            } // primary block computation + admissability checks for SDCs are finished

            unsigned char teamCheck = 0;
            /* Report to my team */
            MPI_Allreduce(&reportFlag, &teamCheck, 1, MPI_BYTE, MPI_MAX, MPI_COMM_WORLD);

            /* Rank 0 of my team reports to its replicas */
            if (myRankInTeam == 0) {
                for (int destTeam = 0; destTeam < numTeams; destTeam++) {
                    if (destTeam != myTeam)
                      MPI_Send(&teamCheck, 1, MPI_BYTE, destTeam, 30, interTeamComm);
                }
            }

            /* load if SDC reported in my team */
            if(teamCheck == 1) {
                // TODO for debugging
                //std::cout << "\t*** team " << myTeam << ", rank " << myRankInTeam
                          //<< "; received SDC report from my Team, teamCheck = "
                          //<< static_cast<int>(teamCheck) << " : WARNING! SOFT ERROR IS PRESENT IN MY TEAM !!"
                          //<< std::endl;

                int reloadTeam;
                /* Receive the reload team */
                MPI_Recv(&reloadTeam, 1, MPI_INT, MPI_ANY_SOURCE, 40,
                         interTeamComm, MPI_STATUS_IGNORE);
                // TODO for debugging
                //std::cout << "\t*** team " << myTeam << ", rank " << myRankInTeam
                          //<< "; loading checkpoint from team " << reloadTeam
                          //<< "..." << std::endl;
                loadCheckpointCallback2_softErrorRecovery(reloadTeam);
            }

            /* Report flags to be received from the other teams */
            unsigned char receivedReportFlag[numTeams];
            /* Flag to set to indicate that we have already received from that team */
            unsigned char replicaCheck[numTeams];
            for (int i = 0; i < numTeams; i++) {
                receivedReportFlag[i] = 0;
                replicaCheck[i] = 0;
            }
            while (true) {
                bool noSDC = true;
                /* Receive report from all my replicas */
                for (int sourceTeam = 0; sourceTeam < numTeams; sourceTeam++) {
                    /* don't receive if already received */
                    if (!replicaCheck[sourceTeam] && sourceTeam != myTeam) {
                        if (myRankInTeam == 0) {
                            MPI_Recv(receivedReportFlag+sourceTeam, 1, MPI_BYTE,
                                    sourceTeam, 30, interTeamComm,
                                    MPI_STATUS_IGNORE);
                        }
                        /* Inform my team about the other teams */
                        MPI_Bcast(receivedReportFlag+sourceTeam, 1, MPI_BYTE, 0,
                                  MPI_COMM_WORLD);

                        /* if no error is reported */
                        if (!receivedReportFlag[sourceTeam])
                            replicaCheck[sourceTeam]++;
                        else
                            noSDC = false;
                    }
                }

                /* Break if no SDC is reported */
                if (noSDC) break;

                /* If myTeam is the team with the smallest index, that has not
                 * been corrupted, then write checkpoint for all the corrupted
                 * teams
                 */
                bool writeCheckpoint = true;
                for (int i = 0; i < myTeam; i++) writeCheckpoint &= receivedReportFlag[i] == 1;

                /* Extract the teams with SDC */
                if (writeCheckpoint) {
                    std::vector<int> failedTeams;
                    for (int i = 0; i < numTeams; i++) {
                        if (receivedReportFlag[i] == 1) {
                            // TODO for debugging
                            //std::cout << "\t+++ TEAM " << myTeam << ": Replica "
                                      //<< myRankInTeam << ", team " << i
                                      //<< " got soft error\n"
                                      //<< "\t+++ writing checkpoint for ya!"
                                      //<< std::endl;
                            failedTeams.push_back(i);
                            /* Send myTeam to let the failed team know which
                             * team is writing checkpoint for it */
                            MPI_Send(&myTeam, 1, MPI_INT, i, 40, interTeamComm);
                        }
                    }
                    createCheckpointCallback(failedTeams);
                }
            }

            /* Primary Block computation + Admissability checks + SDC Reports
             * are finished */

            std::vector<MPI_Request> send_reqs(11 * numTeams, MPI_REQUEST_NULL);

            /*  In this loop do:
             *    If primary block, it is already computed and
             *                      checked for SDCs, we send it
             *    else it is secondary block, we receive it
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
                    int code = MPI_Waitall(11, reqs.data(), stats.data());
                    //std::cout << "T" << myTeam << "R" << myRankInTeam << " MPI_Waitall end ^^" << std::endl;
                    if (code != MPI_SUCCESS || hasRecovered) {
                        std::cout << "*#*#*#*#*#*# unknown error ?  " << std::endl;
                        assert(false);
                        if (code != MPI_SUCCESS) {
                            std::cout << "Team " << myTeam << ": Error in Waitall for block " << currentBlockNr << ": "
                                      << code << '\n';
                        }
                        currentBlock.computeNumericalFluxes();
                    }
                } // end of else, secondary blocks
            }
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
            MPI_Allreduce(&minTimestep, &timestep, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
            //std::cout << "T" << myTeam << "R" << myRankInTeam
                      //<< " : Max Timestep = " << timestep << std::endl;

            for (auto& block : simulationBlocks) { block->maxTimestep = timestep; }
            for (auto& block : simulationBlocks) { block->updateUnknowns(timestep); }
            t += timestep;

            timeSinceLastHeartbeat = MPI_Wtime() - timeOfLastHeartbeat;
            MPI_Bcast(&timeSinceLastHeartbeat, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            // TODO for debugging
            //std::cout << "T" << myTeam << "R" << myRankInTeam << " ALIVE " << std::endl;
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
        //ft_logger.ft_print_HBend(timeSinceLastHeartbeat, t);
    }
    // TODO too much debugging printing, remove them later
    //std::cout << "T" << myTeam << "R" << myRankInTeam << " ACTUALLY FINISHED" << std::endl;
    //std::cout << "T" << myTeam << "R" << myRankInTeam << " clearing simulationBlocks " << std::endl;
    for (auto& block : simulationBlocks) { block->freeMpiType(); block.reset(); }
    /* delete the pointers as well */
    simulationBlocks.clear();
    //std::cout << "T" << myTeam << "R" << myRankInTeam << " clearing scenario" << std::endl;
    if (!SDC_injected) delete scenario; scenario = nullptr;
    //std::cout << "T" << myTeam << "R" << myRankInTeam << " clearing t" << std::endl;
    t = 0.f;
    //std::cout << "T" << myTeam << "R" << myRankInTeam << " clearing backupMetadataNames" << std::endl;
    for (auto& data : backupMetadataNames) { data.clear(); }
    backupMetadataNames.clear(); outputTeamName.clear(); backupTeamName.clear();
    //std::cout << "T" << myTeam << "R" << myRankInTeam << " clearing restartnameinput" << std::endl;
    restartNameInput.clear();

    // TODO return meaningful things like 3.14.. JOKING! HAH! but you can return TestResults struct for example..
    return 0;
}

/* -- TESTS -- */

        /* TODO inject bitflip into udpates
         *      inject a very big value / small value
         */

/* NaN injection */
int TEST_bitflipIntoUpdates_1(TestArguments *args, int rankToCorrupt) {
    return run(args, 0, 0, rankToCorrupt);
}

/* Negative water height injection */
int TEST_bitflipIntoUpdates_2(TestArguments *args, int rankToCorrupt) {
    return run(args, 0, 1, rankToCorrupt);
}

/* NaN injection */
int TEST_bitflipIntoData_1(TestArguments *args, int rankToCorrupt) {
    return run(args, 1, 0, rankToCorrupt);
}

/* Negative water height injection */
int TEST_bitflipIntoData_2(TestArguments *args, int rankToCorrupt) {
    return run(args, 1, 1, rankToCorrupt);
}

/* Bathymetry change injection */
int TEST_bitflipIntoData_3(TestArguments *args, int rankToCorrupt) {
    return run(args, 1, 2, rankToCorrupt);
}


        /* TODO inject bitflip into data
         *      inject a very big / small value
         */

//------------------------------------------------------------------------------

/**
 * TODO Description
 *
 * TODO introduce injection information and different ranks / teams injecting SDCs, not just rank 0
 * TODO introduce TestResults struct to return the corrupted result and check if it is
 *      really corrected, and return FAILED if correction failed!
 */
void run_all_tests(TestArguments *args,
        std::vector<std::function<int(TestArguments*,int)>*> tests_data,
        std::vector<std::function<int(TestArguments*,int)>*> tests_updates) {
    int worldRank;
    int worldSize;
    int myTeam;
    int myRankInTeam;
    int jumped = 0; // TODO can be removed after debugging ?
    int rankToCorrupt;
    PMPI_Comm_rank(TMPI_GetWorldComm(), &worldRank);
    PMPI_Comm_size(TMPI_GetWorldComm(), &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    myTeam = TMPI_GetTeamNumber();

    /* Seed the random generator */
    std::srand (static_cast <unsigned> (time(NULL)));
    rankToCorrupt = 0; // TODO change this after debugging

    /* iterate through all the tests and run them */

    /* tests for SDCs in data arrays */
    for (std::function<int(TestArguments*,int)> *test : tests_data) {
        if (worldRank == 0) std::cout << color_blue << "\n" << markRightArrow << " TEST " << (testsPassed+1) << " IS RUNNING..\n" << color_norm << std::endl;
        //rankToCorrupt = std::rand() % worldSize; TODO do this after debugging
        SDC_injected = false; jumped = 0;
        PMPI_Barrier(TMPI_GetWorldComm()); /* wait here for the other ranks */

        jumped = setjmp(jumpBuffer);
        if (jumped == 1) SDC_injected = true;

        (*test)(args,rankToCorrupt);

        if(worldRank == 0) std::cout << color_green << "\n-- TEST " << (testsPassed+1) << " HAVE PASSED!\t\t" << markPassed << "\n" << color_norm << std::endl;
        testsPassed++;
    }

    std::cout << "\n\n DATA ARRAY TESTS PASSED ! NOW TESTING WITH UPDATE ARRAYS.." << std::endl;

    /* tests for SDCs in update arrays */
    for (std::function<int(TestArguments*,int)> *test : tests_updates) {
        if (worldRank == 0) std::cout << color_blue << "\n" << markRightArrow << " TEST " << (testsPassed+1) << " IS RUNNING..\n" << color_norm << std::endl;
        //rankToCorrupt = std::rand() % worldSize; TODO do this after debugging
        SDC_injected = false; jumped = 0;
        PMPI_Barrier(TMPI_GetWorldComm()); /* wait here for the other ranks */

        jumped = setjmp(jumpBuffer);
        if (jumped == 1) SDC_injected = true;

        (*test)(args,rankToCorrupt);

        if(worldRank == 0) std::cout << color_green << "\n-- TEST " << (testsPassed+1) << " HAVE PASSED!\t\t" << markPassed << "\n" << color_norm << std::endl;
        testsPassed++;
    }
}


/* Initializes MPI and runs the tests
 *
 * @return result Returns 0 if all tests have passed
 */
int main(int argc, char** argv) {

    /* Some variables for tests */
    int result = 0;
    size_t testsPassed = 0;

    /* set string names */
    outputNameInput = "swe_tests";
    backupNameInput = "BACKUP_swe_tests";

    /* Test parameters can be changed here */
    float simulationDuration = 5;
    clock_t heartbeatInterval = 1;
    int nxRequested = 500; int nyRequested = 500;
    unsigned int decompFactor = 1;
    bool writeOutput = false;
    double bitflip_at = 2.5f;

    TestArguments *args = new TestArguments(simulationDuration,
                                            heartbeatInterval,
                                            nxRequested, nyRequested,
                                            decompFactor,
                                            writeOutput,
                                            bitflip_at);

    /* all the test cases */
    std::vector<std::function<int(TestArguments*,int)>*> tests_data;
    std::vector<std::function<int(TestArguments*,int)>*> tests_updates;

    /* add more tests cases here */
    // test functions for SDC in the data arrays
    std::function<int(TestArguments*,int)> t_data_1(TEST_bitflipIntoData_1);
    std::function<int(TestArguments*,int)> t_data_2(TEST_bitflipIntoData_2);
    std::function<int(TestArguments*,int)> t_data_3(TEST_bitflipIntoData_3);
    // test functions for SDC in the update arrays
    std::function<int(TestArguments*,int)> t_updates_1(TEST_bitflipIntoUpdates_1);
    std::function<int(TestArguments*,int)> t_updates_2(TEST_bitflipIntoUpdates_2);
    // deploy test functions for testing!
    tests_data.push_back(&t_data_1);        // NaN injection
    tests_data.push_back(&t_data_2);        // Negative water height injection
    tests_data.push_back(&t_data_3);        // Bathymetry change injection
    tests_updates.push_back(&t_updates_1);  // NaN injection
    // TODO add more tests : discrete maximum criteria (very big or very small values)
    //                       negative water height on updates...

    /*************** START INITIALIZING ***************/
    /* Init teaMPI
     * No hard error resilience tests in this test file */
    TMPI_SetErrorHandlingStrategy(TMPI_NoErrorHandler);
    /* Init MPI */
    MPI_Init(&argc, &argv);
    /* TODO also print something fancy */
    // TODO change this so that we can have all the tests in one for loop

    // TODO move these variables to global or something.. they are used everywhere !
    int worldRank;
    PMPI_Comm_rank(TMPI_GetWorldComm(), &worldRank);

    /*************** RUN ALL THE TESTS  ***************/
    if (worldRank == 0) {
        std::cout << "\n" << color_green << lightShade
                  << " STARTING TESTS "
                  << lightShade << color_norm << std::endl;
    }
    run_all_tests(args, tests_data, tests_updates);
    MPI_Finalize();
    /*************** TEST RUNS FINISHED ***************/
    if (worldRank == 0) {
        std::cout << "\n" << color_green << lightShade
                  << " TESTS ARE FINISHED ! "
                  << lightShade << color_norm << std::endl;
    }

    delete args;
    return result;
}

