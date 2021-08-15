/**
 * @file src/tolerance/swe_noRes.cpp
 *
 * @brief Method "NoRes": No resilience for benchmarking
 *
 * @author Atamert Rahma rahma@in.tum.de
 *
 * No error resilience for benchmarking. However this can provide
 * a naive soft error detection if we run the application twice,
 * and even soft error resilience (detection + correction) if we
 * run the application 3 times (assuming that we would have at least
 * 2 equal solutions). This can be used as a baseline model.
 *
 * Here is a short pseudo-code of the computation loop:
 *
 *      while (t < simulationDuration) {
 *        compute
 *        agreeOnTimestep
 *        updateUnknowns
 *      }
 */


#include <memory>
#include <mpi.h>
#include <climits>
#include <iostream>
#include <unistd.h>

#include "blocks/DimSplitMPIOverdecomp.hpp"
#include "io/Reader.hpp"
#include "io/Writer.hpp"
#include "scenarios/simple_scenarios.hpp"
#include "tools/Args.hpp"
#include "types/Boundary.hpp"


//------------------------------------------------------------------------------

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

//******************************************************************************

int main(int argc, char** argv) {
    double startTime = MPI_Wtime();

    // Define command line arguments
    tools::Args args;

    args.addOption("simulation-duration", 't', "Time in seconds to simulate");
    args.addOption("resolution-x", 'x', "Number of simulated cells in x-direction");
    args.addOption("resolution-y", 'y', "Number of simulated cells in y-direction");
    args.addOption("output-basepath", 'o', "Output base file name");
    args.addOption("write-output", 'w', "Write output using netcdf writer to the specified output file", args.No, false);
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

    /* Parameters to read */
    float simulationDuration;
    int nxRequested, nyRequested;
    std::string outputNameInput = "";
    std::string outputTeamName = "";

    /* whether to write output */
    bool writeOutput = false;
    bool verbose = false;

    // Read in command line arguments
    simulationDuration = args.getArgument<float>("simulation-duration");
    nxRequested = args.getArgument<int>("resolution-x");
    nyRequested = args.getArgument<int>("resolution-y");
    outputNameInput = args.getArgument<std::string>("output-basepath");
    writeOutput = args.isSet("write-output");
    verbose = args.isSet("verbose");

    /* one block per rank */
    std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp> simulationBlock;

    /* Simulation time */
    float t = 0.f;

    /* init MPI */
    int ranksPerTeam = 1;
    int myRankInTeam;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ranksPerTeam);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInTeam);
    const int myTeam = 0;
    int numTeams = 1;
    unsigned int blocksPerRank = 1;

    outputTeamName = outputNameInput + "_t" + std::to_string(myTeam);

    // Print status
    char hostname[HOST_NAME_MAX];
    gethostname(hostname, HOST_NAME_MAX);

    std::printf("PID %d, Rank %i of Team %i spawned at %s with start time %f\n", getpid(), myRankInTeam, myTeam, hostname, startTime);

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

    SWE_Scenario *scenario = new SWE_RadialBathymetryDamBreakScenario{};
    int widthScenario = scenario->getBoundaryPos(BND_RIGHT) - scenario->getBoundaryPos(BND_LEFT);
    int heightScenario = scenario->getBoundaryPos(BND_TOP) - scenario->getBoundaryPos(BND_BOTTOM);

    float dxSimulation = (float)widthScenario / nxRequested;
    float dySimulation = (float)heightScenario / nyRequested;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    /* Loop over the rank's block number in its team. */
    int localBlockPositionX = startPoint / blockCountY;
    int localBlockPositionY = startPoint % blockCountY;

    //compute local number of cells for each SWE_Block w.r.t. the
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

    simulationBlock = std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>(
        new SWE_DimensionalSplittingMPIOverdecomp(nxLocal,
                                                  nyLocal,
                                                  dxSimulation,
                                                  dySimulation,
                                                  localOriginX,
                                                  localOriginY,
                                                  0,
                                                  outputTeamPosName,
                                                  "",
                                                  true, // always write for checkpointing
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
            // neighbourBlocks[j] = simulationBlocks[myNeighbours[j] - startPoint];
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

//------------------------------------------------------------------------------
    simulationBlock->sendBathymetry();
    simulationBlock->recvBathymetry();

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

//------------------------------------------------------------------------------

    // Write zero timestep
    if (writeOutput) simulationBlock->writeTimestep(0.f);

    if (verbose) {
    std::cout << "+---------------------+\n"
              << "| Starting Simulation |\n"
              << "+---------------------+\n"
              << " # We are at t = " << t << ", total: " << simulationDuration
              << "\t\t\t" << "---- Rank " << myRankInTeam
              << std::endl;
    }

    auto& block = *simulationBlock;

    while (t < simulationDuration) {
        // exchange boundaries between blocks
        block.setGhostLayer();
        block.receiveGhostLayer();

        if (verbose) {
            std::cout << "calculating.. t = " << t << "\t\t/ "
                      << simulationDuration;
        }

        if(!writeOutput) {
            std::cout << "\t\t\t\t---- Rank " << myRankInTeam << std::endl;
        }

        block.computeNumericalFluxes();

        /* Agree on a timestep */
        timestep = block.maxTimestep;
        float agreed_timestep;
        MPI_Allreduce(&timestep, &agreed_timestep, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

        block.maxTimestep = agreed_timestep;
        block.updateUnknowns(agreed_timestep);

        t += agreed_timestep;

        if (writeOutput) {
            if (verbose) {
                std::cout << "  --> Writing timestep at: " << t
                        << "\t\t\t\t---- Rank " << myRankInTeam /* equals global rank with one team */
                        << std::endl;
            }
            block.writeTimestep(t);
        }
    }

    simulationBlock->freeMpiType();
    double totalTime = MPI_Wtime() - startTime;
    std::cout << "Rank " << myRankInTeam << " from TEAM " << myTeam
              << " total time : " << totalTime << std::endl;

    MPI_Finalize();
    return 0;
}
