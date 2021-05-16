/**
 * @file
 * This file is part of SWE.
 *
 * @author Jurek Olden
 * (jurek.olden
 * AT in.tum.de)
 *
 * @section LICENSE
 *
 * SWE is free
 * software: you can
 * redistribute it and/or modify
 * it under the terms of
 * the GNU General Public
 * License as published by
 * the Free Software
 * Foundation, either version 3 of
 * the License, or
 * (at your option) any
 * later version.
 *
 * SWE is
 * distributed in the hope that it will be
 * useful,
 * but WITHOUT ANY WARRANTY;
 * without even the implied warranty of

 * * MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the
 * GNU
 * General Public License for more details.

 * *
 * You should have received a
 * copy of the GNU General Public License
 *
 * along with SWE.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 *
 *
 * @section DESCRIPTION
 *
 * Implementation of the SWE_Block abstract class that
 * uses dimensional
 * splitting.
 * It extends the computational domain to two
 * dimensions by
 * decomposing 2D updates
 * to updates on the x- and y-axis.
 *

 */

#ifndef SWE_DIM_SPLIT_MPI_OVERDECOMP_HPP
#define SWE_DIM_SPLIT_MPI_OVERDECOMP_HPP

#include <mpi.h>
#include <time.h>

#include <ctime>

#include "blocks/Block.hpp"
#include "io/Writer.hpp"
#include "scenarios/Scenario.hpp"
#include "types/Float2DNative.hpp"
#if SOLVER == 0
    #include "solvers/HLLEFun.hpp"
#endif

class SWE_DimensionalSplittingMPIOverdecomp : public SWE_Block<Float2DNative>
{
public:
    // Constructor/Destructor
    SWE_DimensionalSplittingMPIOverdecomp() = default;
    SWE_DimensionalSplittingMPIOverdecomp(int cellCountHorizontal,
                                          int cellCountVertical,
                                          float cellSizeHorizontal,
                                          float cellSizeVertical,
                                          float originX,
                                          float originY,
                                          bool localTimestepping,
                                          std::string name,
                                          std::string backupName,
                                          bool write,
                                          bool existingFile);
    ~SWE_DimensionalSplittingMPIOverdecomp(){};

    // Interface methods
    void setGhostLayer();
    void receiveGhostLayer();
    void computeNumericalFluxes();

    void updateUnknowns(float dt);

    SWE_DimensionalSplittingMPIOverdecomp* left;
    SWE_DimensionalSplittingMPIOverdecomp* right;
    SWE_DimensionalSplittingMPIOverdecomp* bottom;
    SWE_DimensionalSplittingMPIOverdecomp* top;
    void setLeft(SWE_DimensionalSplittingMPIOverdecomp* argLeft);
    void setRight(SWE_DimensionalSplittingMPIOverdecomp* argRight);
    void freeMpiType();
    void freeWriter();

    void connectNeighbourLocalities(int neighbourRankId[]);
    void connectLocalNeighbours(std::array<std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>, 4> neighbourBlocks);

    int neighbourLocality[4];

    void writeTimestep(float timestep);
    void createCheckpoint(float t, std::string backupMetadataName, int checkpointsLeft);

    /* validates physical and numerical admissability criteria */
    bool validateAdmissability(float timestep);

    struct blockData_s
    {
        float t;
        int blockNr;

        float* hNetUpdatesLeft;
        float* hNetUpdatesRight;

        float* huNetUpdatesLeft;
        float* huNetUpdatesRight;

        float* hNetUpdatesBelow;
        float* hNetUpdatesAbove;

        float* hvNetUpdatesBelow;
        float* hvNetUpdatesAbove;
    };
    MPI_Datatype blockData_t;

    std::shared_ptr<blockData_s> getBlockDataBuffer(float t, int blockNr);
    MPI_Datatype getBlockMPIType();

    MPI_Datatype HORIZONTAL_BOUNDARY;
    std::shared_ptr<io::Writer> writer;
    bool write;
#if WAVE_PROPAGATION_SOLVER == 0
    //! Hybrid solver (f-wave + augmented)
    // solver::Hybrid<float> solver;
    solver::HLLEFun<float> solver;
#elif WAVE_PROPAGATION_SOLVER == 1
    //! F-wave Riemann solver
    solver::FWave<float> solver;
#elif WAVE_PROPAGATION_SOLVER == 2
    //! Approximate Augmented Riemann solver
    solver::AugRie<float> solver;
#endif

    // Temporary values after x-sweep and before y-sweep
    Float2DNative hStar;
    Float2DNative huStar;

    // net updates per cell
    Float2DNative hNetUpdatesLeft;
    Float2DNative hNetUpdatesRight;

    Float2DNative huNetUpdatesLeft;
    Float2DNative huNetUpdatesRight;

    Float2DNative hNetUpdatesBelow;
    Float2DNative hNetUpdatesAbove;

    Float2DNative hvNetUpdatesBelow;
    Float2DNative hvNetUpdatesAbove;

    void sendBathymetry();
    void recvBathymetry();
};

#endif // SWE_DIM_SPLIT_MPI_OVERDECOMP_HPP
