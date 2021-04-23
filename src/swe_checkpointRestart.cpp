/**
 * @file src/swe_checkpointRestart.cpp
 *
 * @brief baseline model for resilience, checkpoint/restart.
 *
 * TODO
 *   - remove dependencies with teaMPI
 *
 *   - use checkpointsLeft field. This field was not used in other
 *     failure mitigation methods because our checkpoints were real
 *     time dependent.
 *
 *  - remove unneccesary blocks
 *
 *  - don't allow -r and -c at the same time (also x,y and t)
 *
 *  - tested with 1, 2 and 4 processes but with 6 processes:
 *    ftagree:agreement (ERA) Error: Agreement ID (0.0).1 was started by some processor, but I never completed to it
 *
 */



#include <memory>
#include <mpi.h>
#include <string>
//#include <teaMPI.h>
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

// some data needs to be global, otherwise the checkpoint callbacks cannot
// access it.
jmp_buf jumpBuffer{};
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::vector<std::string> backupMetadataNames{};

/* Since we are doing the first option naively, we will not share any
* blocks. Which means it is best that we have only one block per replica,
* so using only one variable makes more sense TODO*/

std::string backupMetadataName;
//------------------------------------------------------------------------------

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
    args.addOption("checkpoint-count", 'c', "Number of total checkpoints");
    args.addOption("resolution-x", 'x', "Number of simulated cells in x-direction");
    args.addOption("resolution-y", 'y', "Number of simulated cells in y-direction");
    args.addOption("output-basepath", 'o', "Output base file name");

    /* TODO: Remove this backup and restart options as you improve the checkpointing */
    args.addOption("backup-basepath", 'b', "Output base file name");
    args.addOption("restart-basepath", 'r', "Restart base file name", tools::Args::Required, false);

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

    /* Parameters to read */
    float simulationDuration;
    int numberOfCheckpoints;

    int nxRequested, nyRequested;

    // TODO remove backup name. Start the simulation with a unique name so that
    //      you know the simulation is not finished and restart from that output
    //      Example : use the suffix _CONT_
    //      of course with this we can remove restart name as well
    //      Example : indicate that you are checkpointing with a flag
    std::string outputNameInput = "";
    std::string backupNameInput = "";
    std::string restartNameInput = "";

    // TODO since we are naively writing checkpoints, we have only one team.
    //      no replication is necessary so no team name is necessary.
    std::string outputTeamName{""};
    std::string backupTeamName{""};

    // Read in command line arguments
    simulationDuration = args.getArgument<float>("simulation-duration");
    numberOfCheckpoints = args.getArgument<int>("checkpoint-count");
    assert(numberOfCheckpoints > 0);

    nxRequested = args.getArgument<int>("resolution-x");
    nyRequested = args.getArgument<int>("resolution-y");
    outputNameInput = args.getArgument<std::string>("output-basepath");
    backupNameInput = args.getArgument<std::string>("backup-basepath");
    if (args.isSet("restart-basepath")) {
        restartNameInput = args.getArgument<std::string>("restart-basepath");
    }

    /* init MPI */
    int ranksPerTeam = 1;
    int myRankInTeam;
    if (setjmp(jumpBuffer) == 0) {
        MPI_Init(&argc, &argv);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &ranksPerTeam);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    const int myTeam = 0;
    int numTeams = 1;
    unsigned int blocksPerRank = 1;

    // TODO : remove this later, we have only one team
    outputTeamName = outputNameInput + "_t" + std::to_string(myTeam);
    backupTeamName = backupNameInput + "_t" + std::to_string(myTeam);

    /* also add team to the restart name TODO fix when file not found */
    if (args.isSet("restart-basepath")) {
        restartNameInput = restartNameInput + "_t" + std::to_string(myTeam);
    }

    // Print status
    char hostname[HOST_NAME_MAX];
    gethostname(hostname, HOST_NAME_MAX);

    std::printf("Rank %i of Team %i spawned at %s\n", myRankInTeam, myTeam, hostname);


    /* int totalBlocks = blocksPerRank * ranksPerTeam; */
    int totalBlocks = ranksPerTeam; // actually total number of global ranks
                                    // (there is only one team)

    // number of SWE-Blocks in x- and y-direction
    int blockCountY = std::sqrt(totalBlocks);
    while (totalBlocks % blockCountY != 0) blockCountY--;
    int blockCountX = totalBlocks / blockCountY;

    /* We have only one team, meaning blocksPerRank equals to 1*/
    /* int startPoint = myRankInTeam * blocksPerRank; */
    int startPoint = myRankInTeam;

    float simulationStart = 0.f;

    // if not loading from a checkpoint
    if (restartNameInput == "") {
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

            /* only in classical checkpointing */
            numberOfCheckpoints = reader.getRemainingCheckpoints();

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

    /* Compute when the checkpoints are reached */
    float *checkpointInstantOfTime = new float[numberOfCheckpoints];
    /* Time delta between checkpoints */
    float checkpointTimeDelta = simulationDuration / numberOfCheckpoints;
    /* The first checkpoint is reached after 0 + delta t */
    checkpointInstantOfTime[0] = checkpointTimeDelta;
    for (int i = 1; i < numberOfCheckpoints; i++) {
        checkpointInstantOfTime[i] = checkpointInstantOfTime[i - 1] + checkpointTimeDelta;
    }

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

//------------------------------------------------------------------------------

    // Write zero timestep
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/* TODO on this naiv soft error checking we have only one simulationBlock ! */
    for (auto& block : simulationBlocks) { block->writeTimestep(0.f); }
//------------------------------------------------------------------------------


    /* Simulate start TODO */

    std::cout << "+---------------------+\n"
              << "| Starting Simulation |\n"
              << "+---------------------+\n"
              << " # We are at t = " << t << ", total: " << simulationDuration
              << "\t\t\t" << "---- Rank " << myRankInTeam
              << std::endl;

    /* We should only have one block per rank */
    assert(simulationBlocks.size() == 1);

    for (int i = 0; i < numberOfCheckpoints; i++) {

        /* Simulate until the checkpoint is reached */
        while (t < checkpointInstantOfTime[i]) {

            // exchange boundaries between blocks TODO only one block
            for (auto& currentBlock : simulationBlocks) { currentBlock->setGhostLayer(); }
            for (auto& currentBlock : simulationBlocks) { currentBlock->receiveGhostLayer(); }

            // Get a reference to the current block TODO there is only one block
            const int& currentBlockNr = 0;
            auto& currentBlock = *simulationBlocks[currentBlockNr];

            std::cout << "calculating.. t = " << t << "\t\t/ "
                      << simulationDuration;

            currentBlock.computeNumericalFluxes();

            hasRecovered = false;

            /* Agree on a timestep */
            timestep = simulationBlocks[0]->maxTimestep;
            float agreed_timestep;
            MPI_Allreduce(&timestep, &agreed_timestep, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

            // TODO only one block
            for (auto& block : simulationBlocks) { block->maxTimestep = agreed_timestep; }
            for (auto& block : simulationBlocks) { block->updateUnknowns(agreed_timestep); }

            t += agreed_timestep;

            std::cout << "  --> Writing timestep at: " << t
                      << "\t\t\t\t---- Rank " << myRankInTeam /* equals global rank with one team */
                      << std::endl;

            /* TODO if write */
            for (auto& block : simulationBlocks) {block->writeTimestep(t);}

        }

        /* checkpoints left */
        int numberOfCheckpointsLeft = (numberOfCheckpoints - 1) - i;

        // TODO only one block
        for (int i = 0; i < (int) simulationBlocks.size(); i++) {
            assert(i == 0);
            std::cout << "      --- CHECKPOINTING AT : " << t << " ---" << std::endl;
            simulationBlocks[i]->createCheckpoint(t, backupMetadataNames[i], numberOfCheckpointsLeft);
        }


    }

    for (auto& block : simulationBlocks) { block->freeMpiType(); }
    MPI_Finalize();
}
