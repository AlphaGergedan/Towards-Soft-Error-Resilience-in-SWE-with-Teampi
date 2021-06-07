/**
 * @file
 * This file is part of an SWE fork created for the
 *
 * Tsunami-Simulation Bachelor Lab Course.
 *
 * @author Philipp Samfass,
 *
 * Alexander Pöppl  (samfass@in.tum.de poeppl@in.tum.de(
 *
 * @section
 * LICENSE

 * *
 * SWE is free software: you can redistribute it and/or modify

 * * it under
 * the terms of the GNU General Public License as published by
 *
 * the Free
 * Software Foundation, either version 3 of the License, or
 * (at
 * your option)
 * any later version.
 *
 * SWE is distributed in the hope that
 * it will be
 * useful,
 * but WITHOUT ANY WARRANTY; without even the implied
 * warranty of
 *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
 * the
 * GNU General
 * Public License for more details.
 *
 * You should have
 * received a copy of the
 * GNU General Public License
 * along with SWE.  If
 * not, see
 * <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *

 * *
 * Implementation of SWE_DimensionalSplittingMPIOverdecomp.hh
 *
 */
#include "blocks/DimSplitMPIOverdecomp.hpp"

#include <bitset>
#include <cmath>
#include <limits>
#include <math.h>
#include <mpi.h>
#include <ostream>
#include <unistd.h>

#include <algorithm>
#include <cassert>

/*
 * Constructor of a SWE_DimensionalSplittingMPIOverdecomp Block.
 * Computational domain is [1,...,nx]*[1,...,ny]
 * Ghost layer consists of two additional rows and columns
 *
 * State variables h, hu, hv and b are defined on the whole grid (including ghost layer)
 * Net updates coming from above/below/left/right are defined for each cell.
 *
 * Net updates are computed on all rows first, then on all columns, the total net updates are then composed
 * from the two 1D solutions.
 *
 * This strategy only works, if the timestep chosen w.r.t. to the maximum horizontal wave speeds
 * also satisfies the CFL-condition in y-direction.
 *
 * @param l_nx Size of the computational domain in x-direction
 * @param l_ny Size of the computational domain in y-direction
 * @param l_dx Cell width
 * @param l_dy Cell height
 */
SWE_DimensionalSplittingMPIOverdecomp::SWE_DimensionalSplittingMPIOverdecomp(
    int nx, int ny, float dx, float dy, float originX, float originY, bool localTimestepping, std::string outputName,
    std::string backupName, bool write, bool existingFile)
    : /*
       * Important note concerning grid allocations:
       * Since index shifts all over the place are bug-prone and maintenance unfriendly,
       * an index of [x][y] is at the actual position x,y on the actual grid.
       * This implies that the allocation size in any direction might be larger than the number of values needed.
       * So if, for instance, array[x][y] needs to hold values in the domain [1,a][1,b],
       * it will be allocated with size (a+1, b+1) instead of (a, b).
       * array[0][0] is then unused.
       */

      // Initialize grid metadata using the base class constructor
      SWE_Block(nx, ny, dx, dy, originX, originY, localTimestepping),
      write(write),
      // intermediate state Q after x-sweep
      hStar(nx + 1, ny + 2),
      huStar(nx + 1, ny + 2),

      /*
       * Temporary storage for the net updates per grid cell during a sweep.
       * There are four update values per cell:
       * Left-going wave from the right edge, analogue for the left edge.
       * Down-going wave from the top edge, analogue for the bottom edge
       */

      // For the x-sweep
      hNetUpdatesLeft(nx + 2, ny + 2),
      hNetUpdatesRight(nx + 2, ny + 2),

      huNetUpdatesLeft(nx + 2, ny + 2),
      huNetUpdatesRight(nx + 2, ny + 2),

      // For the y-sweep
      hNetUpdatesBelow(nx + 1, ny + 2),
      hNetUpdatesAbove(nx + 1, ny + 2),

      hvNetUpdatesBelow(nx + 1, ny + 2),
      hvNetUpdatesAbove(nx + 1, ny + 2),

      // for SDC detection
      b_replica(nx + 2, ny + 2),
      h_prev(nx + 2, ny + 2),
      hv_prev(nx + 2, ny + 2),
      hu_prev(nx + 2, ny + 2) {
    MPI_Type_vector(nx, 1, ny + 2, MPI_FLOAT, &HORIZONTAL_BOUNDARY);
    MPI_Type_commit(&HORIZONTAL_BOUNDARY);
    if (write) {
        writer = io::Writer::createWriterInstance(outputName, backupName, b, {{1, 1, 1, 1}}, nx, ny, dx, dy, 0.0F, 0.0F,
                                                  originX, originY, 1, existingFile);
    }
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* set the array sizes */
    dataArraySize = ((nx + 2) * (ny + 2));
    fieldSizeX = ((nx + 2) * (ny + 2));
    fieldSizeY = ((nx + 1) * (ny + 2));

    /*// Modeled after
    https://gitlab.lrz.de/exahype/ExaHyPE-Engine/-/blob/master/ExaHyPE/exahype/records/ADERDGCellDescription.cpp#L384
    *accessed on 12/12/2020* int Attributes = 10; int blocklen[Attributes] = {1, 1, nx, nx, nx, nx, ny, ny, ny, ny};
    MPI_Aint disp[Attributes];
    MPI_Type_create_struct(Attributes, blocklen, )*/
}

void SWE_DimensionalSplittingMPIOverdecomp::writeTimestep(float time) {
    if (write) {
        writer->writeTimeStep(h, hu, hv, time);
    }
}

void SWE_DimensionalSplittingMPIOverdecomp::createCheckpoint(float time, std::string backupMetadataName,
                                                             int checkpointsLeft) {
    // TODO do we really need to write timestep here? If writing outputs
    //      at the same time, this might conflict. Alternative would be
    //      to write output before writing checkpoint (even while benchmarking)
//    double t1 = MPI_Wtime();
    writeTimestep(time);
//    double t2 = MPI_Wtime();
    writer->updateMetadataFile(backupMetadataName, time, checkpointsLeft);
//    double t3 = MPI_Wtime();
    writer->commitBackup();
//    double t4 = MPI_Wtime();
//    std::cout << "\t\t-- writeTimestep(time) in " << t2-t1 << " seconds" << std::endl;
//    std::cout << "\t\t-- writer->updateMetadataFile(backupMetadataName, time, checkpointsLeft) in " << t3-t2 << " seconds" << std::endl;
//    std::cout << "\t\t-- writer->commitBackup() in " << t4-t3 << " seconds" << std::endl;
    std::cout << "Wrote timestep with time " << time << '\n';
}

void SWE_DimensionalSplittingMPIOverdecomp::connectLocalNeighbours(
    std::array<std::shared_ptr<SWE_DimensionalSplittingMPIOverdecomp>, 4> neighbourBlocks) {
    for (int i = 0; i < 4; i++) {
        if (boundaryType[i] == CONNECT_WITHIN_RANK) {
            switch (i) {
            case BND_LEFT:
                left = neighbourBlocks[i].get();
                break;
            case BND_RIGHT:
                right = neighbourBlocks[i].get();
                break;
            case BND_BOTTOM:
                bottom = neighbourBlocks[i].get();
                break;
            case BND_TOP:
                top = neighbourBlocks[i].get();
                break;
            }
        }
    }
}

void SWE_DimensionalSplittingMPIOverdecomp::connectNeighbourLocalities(int p_neighbourRankId[]) {
    for (int i = 0; i < 4; i++) {
        neighbourLocality[i] = p_neighbourRankId[i];
    }
}

void SWE_DimensionalSplittingMPIOverdecomp::freeMpiType() { MPI_Type_free(&HORIZONTAL_BOUNDARY); }

int getTag(int rank, int tag) {
    // return (tag*100000) + rank;
    // max tag is 32767
    // return (tag*1000)
    return (tag * 6553) + (rank % 511);
}

void SWE_DimensionalSplittingMPIOverdecomp::recvBathymetry() {
    /***********
     * RECEIVE *
     **********/

    MPI_Request recvReqs[4];
    MPI_Status stati[4];

    if (boundaryType[BND_LEFT] == CONNECT) {
        int startIndex = 1;
        MPI_Irecv(b.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourLocality[BND_LEFT],
                  getTag(myRank, MPI_TAG_OUT_B_RIGHT), MPI_COMM_WORLD, &recvReqs[BND_LEFT]);
    } else {
        recvReqs[BND_LEFT] = MPI_REQUEST_NULL;
    }

    if (boundaryType[BND_RIGHT] == CONNECT) {
        int startIndex = (nx + 1) * (ny + 2) + 1;
        MPI_Irecv(b.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourLocality[BND_RIGHT],
                  getTag(myRank, MPI_TAG_OUT_B_LEFT), MPI_COMM_WORLD, &recvReqs[BND_RIGHT]);
    } else {
        recvReqs[BND_RIGHT] = MPI_REQUEST_NULL;
    }

    if (boundaryType[BND_BOTTOM] == CONNECT) {
        for (int i = 1; i < nx + 1; i++) {
            MPI_Irecv(&b[i][0], 1, MPI_FLOAT, neighbourLocality[BND_BOTTOM], getTag(myRank, MPI_TAG_OUT_B_TOP),
                      MPI_COMM_WORLD, &recvReqs[BND_BOTTOM]);
        }
    } else {
        recvReqs[BND_BOTTOM] = MPI_REQUEST_NULL;
    }

    if (boundaryType[BND_TOP] == CONNECT) {
        for (int i = 1; i < nx + 1; i++) {
            MPI_Irecv(&b[i][ny + 1], 1, MPI_FLOAT, neighbourLocality[BND_TOP], getTag(myRank, MPI_TAG_OUT_B_BOTTOM),
                      MPI_COMM_WORLD, &recvReqs[BND_TOP]);
        }
    } else {
        recvReqs[BND_TOP] = MPI_REQUEST_NULL;
    }

    MPI_Waitall(4, recvReqs, stati);
}

void SWE_DimensionalSplittingMPIOverdecomp::sendBathymetry() {
    if (boundaryType[BND_RIGHT] == CONNECT_WITHIN_RANK) {
        for (int i = 1; i < ny + 1; i++) {
            b[nx + 1][i] = right->getBathymetry()[1][i];
        }
    }
    if (boundaryType[BND_LEFT] == CONNECT_WITHIN_RANK) {
        for (int i = 1; i < ny + 1; i++) {
            b[0][i] = left->getBathymetry()[nx][i];
        }
    }
    if (boundaryType[BND_TOP] == CONNECT_WITHIN_RANK) {
        for (int i = 1; i < nx + 1; i++) {
            b[i][ny + 1] = top->getBathymetry()[i][1];
        }
    }
    if (boundaryType[BND_BOTTOM] == CONNECT_WITHIN_RANK) {
        for (int i = 1; i < nx + 1; i++) {
            b[i][0] = bottom->getBathymetry()[i][ny];
        }
    }
    // The requests generated by the Isends are immediately freed, since we will
    // wait on the requests generated by the corresponding receives
    MPI_Request req;
    /*********
     * SEND *
     ********/
    if (boundaryType[BND_LEFT] == CONNECT) {
        int startIndex = ny + 2 + 1;
        MPI_Isend(b.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourLocality[BND_LEFT],
                  getTag(neighbourRankId[BND_LEFT], MPI_TAG_OUT_B_LEFT), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);
    }
    if (boundaryType[BND_RIGHT] == CONNECT) {
        int startIndex = nx * (ny + 2) + 1;
        MPI_Isend(b.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourLocality[BND_RIGHT],
                  getTag(neighbourRankId[BND_RIGHT], MPI_TAG_OUT_B_RIGHT), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);
    }
    if (boundaryType[BND_BOTTOM] == CONNECT) {
        for (int i = 1; i < nx + 1; i++) {
            MPI_Isend(&b[i][1], 1, MPI_FLOAT, neighbourLocality[BND_BOTTOM],
                      getTag(neighbourRankId[BND_BOTTOM], MPI_TAG_OUT_B_BOTTOM), MPI_COMM_WORLD, &req);
            MPI_Request_free(&req);
        }
    }
    if (boundaryType[BND_TOP] == CONNECT) {
        for (int i = 1; i < nx + 1; i++) {
            MPI_Isend(&b[i][ny], 1, MPI_FLOAT, neighbourLocality[BND_TOP],
                      getTag(neighbourRankId[BND_TOP], MPI_TAG_OUT_B_TOP), MPI_COMM_WORLD, &req);
            MPI_Request_free(&req);
        }
    }
}

void SWE_DimensionalSplittingMPIOverdecomp::setGhostLayer() {
    // Apply appropriate conditions for OUTFLOW/WALL boundaries
    SWE_Block::applyBoundaryConditions();
    if (boundaryType[BND_RIGHT] == CONNECT_WITHIN_RANK && isReceivable(BND_RIGHT)) {
        borderTimestep[BND_RIGHT] = right->getTotalLocalTimestep();
        for (int i = 1; i < ny + 1; i++) {
            bufferH[nx + 1][i] = right->getWaterHeight()[1][i];
            bufferHu[nx + 1][i] = right->getMomentumHorizontal()[1][i];
            bufferHv[nx + 1][i] = right->getMomentumVertical()[1][i];
        }
    }
    if (boundaryType[BND_LEFT] == CONNECT_WITHIN_RANK && isReceivable(BND_LEFT)) {
        borderTimestep[BND_LEFT] = left->getTotalLocalTimestep();
        for (int i = 1; i < ny + 1; i++) {
            bufferH[0][i] = left->getWaterHeight()[left->nx][i];
            bufferHu[0][i] = left->getMomentumHorizontal()[left->nx][i];
            bufferHv[0][i] = left->getMomentumVertical()[left->nx][i];
        }
    }
    if (boundaryType[BND_TOP] == CONNECT_WITHIN_RANK && isReceivable(BND_TOP)) {
        borderTimestep[BND_TOP] = top->getTotalLocalTimestep();
        for (int i = 1; i < nx + 1; i++) {
            bufferH[i][ny + 1] = top->getWaterHeight()[i][1];
            bufferHu[i][ny + 1] = top->getMomentumHorizontal()[i][1];
            bufferHv[i][ny + 1] = top->getMomentumVertical()[i][1];
        }
    }
    if (boundaryType[BND_BOTTOM] == CONNECT_WITHIN_RANK && isReceivable(BND_BOTTOM)) {
        borderTimestep[BND_BOTTOM] = bottom->getTotalLocalTimestep();
        for (int i = 1; i < nx + 1; i++) {
            bufferH[i][0] = bottom->getWaterHeight()[i][bottom->ny];
            bufferHu[i][0] = bottom->getMomentumHorizontal()[i][bottom->ny];
            bufferHv[i][0] = bottom->getMomentumVertical()[i][bottom->ny];
        }
    }
    MPI_Request req;
    float totalLocalTimestep = getTotalLocalTimestep();
    if (boundaryType[BND_LEFT] == CONNECT && isSendable(BND_LEFT)) {
        int startIndex = ny + 2 + 1;

        MPI_Isend(h.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourLocality[BND_LEFT],
                  getTag(neighbourRankId[BND_LEFT], MPI_TAG_OUT_H_LEFT), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);

        MPI_Isend(hu.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourLocality[BND_LEFT],
                  getTag(neighbourRankId[BND_LEFT], MPI_TAG_OUT_HU_LEFT), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);

        MPI_Isend(hv.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourLocality[BND_LEFT],
                  getTag(neighbourRankId[BND_LEFT], MPI_TAG_OUT_HV_LEFT), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);

        MPI_Isend(&totalLocalTimestep, 1, MPI_FLOAT, neighbourLocality[BND_LEFT],
                  getTag(neighbourRankId[BND_LEFT], MPI_TAG_TIMESTEP_LEFT), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);
    }
    if (boundaryType[BND_RIGHT] == CONNECT && isSendable(BND_RIGHT)) {
        int startIndex = nx * (ny + 2) + 1;

        MPI_Isend(h.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourLocality[BND_RIGHT],
                  getTag(neighbourRankId[BND_RIGHT], MPI_TAG_OUT_H_RIGHT), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);

        MPI_Isend(hu.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourLocality[BND_RIGHT],
                  getTag(neighbourRankId[BND_RIGHT], MPI_TAG_OUT_HU_RIGHT), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);

        MPI_Isend(hv.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourLocality[BND_RIGHT],
                  getTag(neighbourRankId[BND_RIGHT], MPI_TAG_OUT_HV_RIGHT), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);

        MPI_Isend(&totalLocalTimestep, 1, MPI_FLOAT, neighbourLocality[BND_RIGHT],
                  getTag(neighbourRankId[BND_RIGHT], MPI_TAG_TIMESTEP_RIGHT), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);
    }
    if (boundaryType[BND_BOTTOM] == CONNECT && isSendable(BND_BOTTOM)) {
        // int code =
        MPI_Isend(&h[1][1], 1, HORIZONTAL_BOUNDARY, neighbourLocality[BND_BOTTOM],
                  getTag(neighbourRankId[BND_BOTTOM], MPI_TAG_OUT_H_BOTTOM), MPI_COMM_WORLD, &req);
        // if(code != MPI_SUCCESS)
        //	printf("%d: No success %d\n", myRank, code);
        MPI_Request_free(&req);

        MPI_Isend(&hu[1][1], 1, HORIZONTAL_BOUNDARY, neighbourLocality[BND_BOTTOM],
                  getTag(neighbourRankId[BND_BOTTOM], MPI_TAG_OUT_HU_BOTTOM), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);

        MPI_Isend(&hv[1][1], 1, HORIZONTAL_BOUNDARY, neighbourLocality[BND_BOTTOM],
                  getTag(neighbourRankId[BND_BOTTOM], MPI_TAG_OUT_HV_BOTTOM), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);
        // printf("%d: Sent to bottom %d, %f at %f\n", myRank,
        // neighbourRankId[BND_BOTTOM], h[1][1], originX);

        MPI_Isend(&totalLocalTimestep, 1, MPI_FLOAT, neighbourLocality[BND_BOTTOM],
                  getTag(neighbourRankId[BND_BOTTOM], MPI_TAG_TIMESTEP_BOTTOM), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);
    }
    if (boundaryType[BND_TOP] == CONNECT && isSendable(BND_TOP)) {
        MPI_Isend(&h[1][ny], 1, HORIZONTAL_BOUNDARY, neighbourLocality[BND_TOP],
                  getTag(neighbourRankId[BND_TOP], MPI_TAG_OUT_H_TOP), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);

        MPI_Isend(&hu[1][ny], 1, HORIZONTAL_BOUNDARY, neighbourLocality[BND_TOP],
                  getTag(neighbourRankId[BND_TOP], MPI_TAG_OUT_HU_TOP), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);

        MPI_Isend(&hv[1][ny], 1, HORIZONTAL_BOUNDARY, neighbourLocality[BND_TOP],
                  getTag(neighbourRankId[BND_TOP], MPI_TAG_OUT_HV_TOP), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);
        // printf("%d: Sent to top %d, %f at %f\n", myRank,
        // neighbourRankId[BND_TOP], h[1][ny], originX);

        MPI_Isend(&totalLocalTimestep, 1, MPI_FLOAT, neighbourLocality[BND_TOP],
                  getTag(neighbourRankId[BND_TOP], MPI_TAG_TIMESTEP_TOP), MPI_COMM_WORLD, &req);
        MPI_Request_free(&req);
    }
    /*printf("%d: send   %d:%d %d:%d %d:%d %d:%d \n", myRank,
            neighbourRankId[BND_LEFT],getTag(neighbourRankId[BND_LEFT], MPI_TAG_TIMESTEP_LEFT)
        ,   neighbourRankId[BND_RIGHT],getTag(neighbourRankId[BND_RIGHT], MPI_TAG_TIMESTEP_RIGHT)
            ,neighbourRankId[BND_BOTTOM], getTag(neighbourRankId[BND_BOTTOM], MPI_TAG_TIMESTEP_BOTTOM)
            ,neighbourRankId[BND_TOP],getTag(neighbourRankId[BND_TOP], MPI_TAG_TIMESTEP_TOP));
*/
    assert(h.getRows() == ny + 2);
    assert(hu.getRows() == ny + 2);
    assert(hv.getRows() == ny + 2);
    assert(h.getCols() == nx + 2);
    assert(hu.getCols() == nx + 2);
    assert(hv.getCols() == nx + 2);

    /*********
     * SEND *
     ********/
}

void SWE_DimensionalSplittingMPIOverdecomp::receiveGhostLayer() {
    /***********
     * RECEIVE *
     **********/

    // The requests generated by the Isends are immediately freed, since we will
    // wait on the requests generated by the corresponding receives

    // std::cout << myRank << "| start recv\n";

    MPI_Status status;
    // 4 Boundaries times 3 arrays (h, hu, hv) means 12 requests
    MPI_Request recvReqs[16];
    MPI_Status stati[16];

    if (boundaryType[BND_LEFT] == CONNECT && isReceivable(BND_LEFT)) {
        int startIndex = 1;

        MPI_Irecv(bufferH.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourLocality[BND_LEFT],
                  getTag(myRank, MPI_TAG_OUT_H_RIGHT), MPI_COMM_WORLD, &recvReqs[0]);
        MPI_Irecv(bufferHu.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourLocality[BND_LEFT],
                  getTag(myRank, MPI_TAG_OUT_HU_RIGHT), MPI_COMM_WORLD, &recvReqs[1]);
        MPI_Irecv(bufferHv.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourLocality[BND_LEFT],
                  getTag(myRank, MPI_TAG_OUT_HV_RIGHT), MPI_COMM_WORLD, &recvReqs[2]);
        MPI_Irecv(&borderTimestep[BND_LEFT], 1, MPI_FLOAT, neighbourLocality[BND_LEFT],
                  getTag(myRank, MPI_TAG_TIMESTEP_RIGHT), MPI_COMM_WORLD, &recvReqs[3]);
    } else {
        recvReqs[0] = MPI_REQUEST_NULL;
        recvReqs[1] = MPI_REQUEST_NULL;
        recvReqs[2] = MPI_REQUEST_NULL;
        recvReqs[3] = MPI_REQUEST_NULL;
    }

    if (boundaryType[BND_RIGHT] == CONNECT && isReceivable(BND_RIGHT)) {
        int startIndex = (nx + 1) * (ny + 2) + 1;
        MPI_Irecv(bufferH.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourLocality[BND_RIGHT],
                  getTag(myRank, MPI_TAG_OUT_H_LEFT), MPI_COMM_WORLD, &recvReqs[4]);
        MPI_Irecv(bufferHu.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourLocality[BND_RIGHT],
                  getTag(myRank, MPI_TAG_OUT_HU_LEFT), MPI_COMM_WORLD, &recvReqs[5]);
        MPI_Irecv(bufferHv.getRawPointer() + startIndex, ny, MPI_FLOAT, neighbourLocality[BND_RIGHT],
                  getTag(myRank, MPI_TAG_OUT_HV_LEFT), MPI_COMM_WORLD, &recvReqs[6]);
        MPI_Irecv(&borderTimestep[BND_RIGHT], 1, MPI_FLOAT, neighbourLocality[BND_RIGHT],
                  getTag(myRank, MPI_TAG_TIMESTEP_LEFT), MPI_COMM_WORLD, &recvReqs[7]);
    } else {
        recvReqs[4] = MPI_REQUEST_NULL;
        recvReqs[5] = MPI_REQUEST_NULL;
        recvReqs[6] = MPI_REQUEST_NULL;
        recvReqs[7] = MPI_REQUEST_NULL;
    }

    if (boundaryType[BND_BOTTOM] == CONNECT && isReceivable(BND_BOTTOM)) {
        MPI_Irecv(&bufferH[1][0], 1, HORIZONTAL_BOUNDARY, neighbourLocality[BND_BOTTOM],
                  getTag(myRank, MPI_TAG_OUT_H_TOP), MPI_COMM_WORLD, &recvReqs[8]);
        MPI_Irecv(&bufferHu[1][0], 1, HORIZONTAL_BOUNDARY, neighbourLocality[BND_BOTTOM],
                  getTag(myRank, MPI_TAG_OUT_HU_TOP), MPI_COMM_WORLD, &recvReqs[9]);
        MPI_Irecv(&bufferHv[1][0], 1, HORIZONTAL_BOUNDARY, neighbourLocality[BND_BOTTOM],
                  getTag(myRank, MPI_TAG_OUT_HV_TOP), MPI_COMM_WORLD, &recvReqs[10]);
        MPI_Irecv(&borderTimestep[BND_BOTTOM], 1, MPI_FLOAT, neighbourLocality[BND_BOTTOM],
                  getTag(myRank, MPI_TAG_TIMESTEP_TOP), MPI_COMM_WORLD, &recvReqs[11]);
    } else {
        recvReqs[8] = MPI_REQUEST_NULL;
        recvReqs[9] = MPI_REQUEST_NULL;
        recvReqs[10] = MPI_REQUEST_NULL;
        recvReqs[11] = MPI_REQUEST_NULL;
    }

    if (boundaryType[BND_TOP] == CONNECT && isReceivable(BND_TOP)) {
        MPI_Irecv(&bufferH[1][ny + 1], 1, HORIZONTAL_BOUNDARY, neighbourLocality[BND_TOP],
                  getTag(myRank, MPI_TAG_OUT_H_BOTTOM), MPI_COMM_WORLD, &recvReqs[12]);
        MPI_Irecv(&bufferHu[1][ny + 1], 1, HORIZONTAL_BOUNDARY, neighbourLocality[BND_TOP],
                  getTag(myRank, MPI_TAG_OUT_HU_BOTTOM), MPI_COMM_WORLD, &recvReqs[13]);
        MPI_Irecv(&bufferHv[1][ny + 1], 1, HORIZONTAL_BOUNDARY, neighbourLocality[BND_TOP],
                  getTag(myRank, MPI_TAG_OUT_HV_BOTTOM), MPI_COMM_WORLD, &recvReqs[14]);
        MPI_Irecv(&borderTimestep[BND_TOP], 1, MPI_FLOAT, neighbourLocality[BND_TOP],
                  getTag(myRank, MPI_TAG_TIMESTEP_BOTTOM), MPI_COMM_WORLD, &recvReqs[15]);
    } else {
        recvReqs[12] = MPI_REQUEST_NULL;
        recvReqs[13] = MPI_REQUEST_NULL;
        recvReqs[14] = MPI_REQUEST_NULL;
        recvReqs[15] = MPI_REQUEST_NULL;
    }

    int code = MPI_Waitall(16, recvReqs, stati);
    if (code != MPI_SUCCESS) {
        printf("%d: No success %d  %d:%d:%d:%d:\n", myRank, code, getTag(myRank, MPI_TAG_TIMESTEP_RIGHT),
               getTag(myRank, MPI_TAG_TIMESTEP_LEFT), getTag(myRank, MPI_TAG_TIMESTEP_TOP),
               getTag(myRank, MPI_TAG_TIMESTEP_BOTTOM));
    }

    checkAllGhostlayers();
}

/**
 * Compute net updates for the block.
 * The member variable #maxTimestep will be updated with the
 * maximum allowed time step size
 */
void SWE_DimensionalSplittingMPIOverdecomp::computeNumericalFluxes() {
    if (!allGhostlayersInSync()) return;
    // maximum (linearized) wave speed within one iteration
    float maxWaveSpeed = (float)0.;
    /***************************************************************************************
     * compute the net-updates for the vertical edges
     **************************************************************************************/

    for (int i = 1; i < nx + 2; i++) {
        const int ny_end = ny + 1;

        for (int j = 1; j < ny_end; ++j) {
            float maxEdgeSpeed;

            solver.computeNetUpdates(h[i - 1][j], h[i][j], hu[i - 1][j], hu[i][j], b[i - 1][j], b[i][j],
                                     hNetUpdatesLeft[i - 1][j - 1], hNetUpdatesRight[i - 1][j - 1],
                                     huNetUpdatesLeft[i - 1][j - 1], huNetUpdatesRight[i - 1][j - 1], maxEdgeSpeed);
            maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
        }
    }

    /***************************************************************************************
     * compute the net-updates for the horizontal edges
     **************************************************************************************/

    for (int i = 1; i < nx + 1; i++) {
        const int ny_end = ny + 2;

        for (int j = 1; j < ny_end; j++) {
            float maxEdgeSpeed;

            solver.computeNetUpdates(h[i][j - 1], h[i][j], hv[i][j - 1], hv[i][j], b[i][j - 1], b[i][j],
                                     hNetUpdatesBelow[i - 1][j - 1], hNetUpdatesAbove[i - 1][j - 1],
                                     hvNetUpdatesBelow[i - 1][j - 1], hvNetUpdatesAbove[i - 1][j - 1], maxEdgeSpeed);

            // update the maximum wave speed
            maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
            // maxTestSpeed = std::max (maxTestSpeed, maxEdgeSpeed);
        }
    }

    if (maxWaveSpeed > 0.00001) {
        maxTimestep = std::min(dx / maxWaveSpeed, dy / maxWaveSpeed);

        maxTimestep *= (float).4;  // CFL-number = .5
    } else {
        // might happen in dry cells
        maxTimestep = std::numeric_limits<float>::max();
    }
}

/**
 * Updates the unknowns with the already computed net-updates.
 *
 * @param dt time step width used in the update. The timestep has to be equal to maxTimestep calculated by
 * computeNumericalFluxes(), since this is the step width used for the intermediary updates after the x-sweep.
 */
void SWE_DimensionalSplittingMPIOverdecomp::updateUnknowns(float dt) {
    if (!allGhostlayersInSync()) return;
    // update cell averages with the net-updates
    dt = maxTimestep;
    for (int i = 1; i < nx + 1; i++) {
        const int ny_end = ny + 1;

        for (int j = 1; j < ny_end; j++) {
            h[i][j] -= dt / dx * (hNetUpdatesRight[i - 1][j - 1] + hNetUpdatesLeft[i][j - 1]) +
                       dt / dy * (hNetUpdatesAbove[i - 1][j - 1] + hNetUpdatesBelow[i - 1][j]);
            hu[i][j] -= dt / dx * (huNetUpdatesRight[i - 1][j - 1] + huNetUpdatesLeft[i][j - 1]);
            hv[i][j] -= dt / dy * (hvNetUpdatesAbove[i - 1][j - 1] + hvNetUpdatesBelow[i - 1][j]);

            if (h[i][j] < 0) {
                // TODO: dryTol
#ifndef NDEBUG
                // Only print this warning when debug is enabled
                // Otherwise we cannot vectorize this loop
                if (h[i][j] < -0.1) {
                    std::cerr << "Warning, negative height: (i,j)=(" << i << "," << j << ")=" << h[i][j] << std::endl;
                    std::cerr << "         b: " << b[i][j] << std::endl;
                }
#endif  // NDEBUG \
       // zero (small) negative depths
                h[i][j] = hu[i][j] = hv[i][j] = 0.;
            } else if (h[i][j] < 0.1)
                hu[i][j] = hv[i][j] = 0.;  // no water, no speed!
        }
    }
    iterationNumber++;
}

/**
 * Redundantly saves bathymetry data to check for SDC during computation
 */
void SWE_DimensionalSplittingMPIOverdecomp::saveBathymetry() {
    /* Copy the bathymetry data into redundant storage for SDC detection */
    std::memcpy(b_replica.getRawPointer(), b.getRawPointer(), dataArraySize * sizeof(float));
}

/**
 * Save previous data arrays h, hv and hu to later check for DMP admissibility.
 * This function should be called before the data arrays has been updated
 */
void SWE_DimensionalSplittingMPIOverdecomp::savePreviousData() {
    /* Copy h */
    std::memcpy(h_prev.getRawPointer(), h.getRawPointer(), dataArraySize * sizeof(float));
    /* Copy hv */
    std::memcpy(hv_prev.getRawPointer(), hv.getRawPointer(), dataArraySize * sizeof(float));
    /* Copy hu */
    std::memcpy(hu_prev.getRawPointer(), hu.getRawPointer(), dataArraySize * sizeof(float));
}

/**
 * Validates the following physical and numerical admissibility
 * criteria, should be called after the computeNumericalFluxes
 * but before the updateUnknowns method:
 *   1. Physical Admissibility
 *      - Bathymetry is unchanged
 *      - No negative water height
 *   2. Numerical Admissibility
 *      - No NaN values
 *      - Discrete Maximum Principle (DMP)
 *
 * Warning: there is still a chance for silent data corruptions, even
 *          if the results are admissible. Adding more admissibility
 *          checks can reduce this probability even more
 *          This check only validates the update and data arrays stored in this block
 *
 * @param timeStep current timestep of the simulation
 * @return admissible returns true if admissible
 */
bool SWE_DimensionalSplittingMPIOverdecomp::validateAdmissibility(float timestep) {
    /* return value */
    bool admissible = true;

    /* admissibility */
    bool updatesAdmissible = true;
    bool dataAdmissible = true;

    float dt = maxTimestep;

    admissible &= !std::isnan(maxTimestep);
    /* loop over the computation domain only */
    for (int i = 1; i < nx + 1; i++) {
        for (int j = 1; j < ny + 2; j++) {
            /**********************************************
              1. NUMERICAL ADMISSIBILITY : no NaN values
             **********************************************/
            updatesAdmissible &= !std::isnan(hNetUpdatesAbove[i][j]);
            updatesAdmissible &= !std::isnan(hNetUpdatesBelow[i][j]);
            updatesAdmissible &= !std::isnan(hvNetUpdatesAbove[i][j]);
            updatesAdmissible &= !std::isnan(hvNetUpdatesBelow[i][j]);
            updatesAdmissible &= !std::isnan(hNetUpdatesLeft[i][j]);
            updatesAdmissible &= !std::isnan(hNetUpdatesRight[i][j]);
            updatesAdmissible &= !std::isnan(huNetUpdatesLeft[i][j]);
            updatesAdmissible &= !std::isnan(huNetUpdatesRight[i][j]);
            dataAdmissible &= !std::isnan(b[i][j]);
            dataAdmissible &= !std::isnan(h[i][j]);
            dataAdmissible &= !std::isnan(hv[i][j]);
            dataAdmissible &= !std::isnan(hu[i][j]);

            /**********************************************
              2. NUMERICAL ADMISSIBILITY : DMP
             **********************************************/
            /* get the neigbours max and min */
            float h_max = std::max({h_prev[i][j-1], h_prev[i][j+1], h_prev[i-1][j], h_prev[i+1][j]});
            float h_min = std::min({h_prev[i][j-1], h_prev[i][j+1], h_prev[i-1][j], h_prev[i+1][j]});
            float hv_max = std::max({hv_prev[i][j-1], hv_prev[i][j+1], hv_prev[i-1][j], hv_prev[i+1][j]});
            float hv_min = std::min({hv_prev[i][j-1], hv_prev[i][j+1], hv_prev[i-1][j], hv_prev[i+1][j]});
            float hu_max = std::max({hu_prev[i][j-1], hu_prev[i][j+1], hu_prev[i-1][j], hu_prev[i+1][j]});
            float hu_min = std::min({hu_prev[i][j-1], hu_prev[i][j+1], hu_prev[i-1][j], hu_prev[i+1][j]});

            /* Relaxation factor d */
            float d = 100.f;

            /*
             * check Discrete Maximum Principle holding only if we already have
             * a previous solution. Previous solution must be saved by calling
             * the function savePreviousData by the application before writing
             * the results with updateUnknowns
             *
             * TODO add DMP checks to update arrays as well
             *      alternative: save the previous updates and try to 'guess' the next update interval
             */
            /* check the current data arrays */
            if (iterationNumber > 0) {
                //bool cond1 = h_min - d <= h[i][j] && h[i][j] <= h_max + d;
                //bool cond2 = hv_min - d <= hv[i][j] && hv[i][j] <= hv_max + d;
                //bool cond3 = hu_min - d <= hu[i][j] && hu[i][j] <= hu_max + d;
                //// TODO for debugging
                //if (!cond1) {
                    //std::cout << "\n\nh neighbours : [" << h_prev[i][j-1] << ", "
                              //<< h_prev[i][j+1] << ", " << h_prev[i-1][j] << ", "
                              //<< h_prev[i+1][j] << "]\n"
                              //<< "h_min-d  < h < h_max+d\n" << h_min-d << " < " << h[i][j] << " < " << h_max+d
                              //<< std::endl;
                //}
                //if (!cond2) {
                    //std::cout << "\n\nhv neighbours : [" << hv_prev[i][j-1] << ", "
                              //<< hv_prev[i][j+1] << ", " << hv_prev[i-1][j] << ", "
                              //<< hv_prev[i+1][j] << "]\n"
                              //<< "hv_min-d < hv < hv_max+d\n" << hv_min-d << " < " << hv[i][j] << " < " << hv_max+d
                              //<< std::endl;
                //}
                //if (!cond3) {
                    //std::cout << "\n\nhu neighbours : [" << hu_prev[i][j-1] << ", "
                              //<< hu_prev[i][j+1] << ", " << hu_prev[i-1][j] << ", "
                              //<< hu_prev[i+1][j] << "]\n"
                              //<< "hu_min-d < hu < hu_max+d\n" << hu_min-d << " < " << hu[i][j] << " < " << hu_max+d
                              //<< std::endl;
                //}

                dataAdmissible &= h_min - d <= h[i][j] && h[i][j] <= h_max + d;
                dataAdmissible &= hv_min - d <= hv[i][j] && hv[i][j] <= hv_max + d;
                dataAdmissible &= hu_min - d <= hu[i][j] && hu[i][j] <= hu_max + d;
            }

            /******************************************************
              1. PHYSICAL ADMISSIBILITY : bathymetry is unchanged
             ******************************************************/
            dataAdmissible &= b[i][j] == b_replica[i][j];


            /******************************************************
              2. PHYSICAL ADMISSIBILITY : no negative water height
             ******************************************************/
            dataAdmissible &= h[i][j] >= 0;
        }
    }
    /* Check NaNs in the ghost layers */

    /* loop over the ghost layers with indices
     * i = 0 and i = nx + 1 */
    for (int j = 0; j < ny + 2; j++) {
        // access i = 0
        updatesAdmissible &= !std::isnan(hNetUpdatesAbove[0][j]);
        updatesAdmissible &= !std::isnan(hNetUpdatesBelow[0][j]);
        updatesAdmissible &= !std::isnan(hvNetUpdatesAbove[0][j]);
        updatesAdmissible &= !std::isnan(hvNetUpdatesBelow[0][j]);
        updatesAdmissible &= !std::isnan(hNetUpdatesLeft[0][j]);
        updatesAdmissible &= !std::isnan(hNetUpdatesRight[0][j]);
        updatesAdmissible &= !std::isnan(huNetUpdatesLeft[0][j]);
        updatesAdmissible &= !std::isnan(huNetUpdatesRight[0][j]);
        dataAdmissible &= !std::isnan(b[0][j]);
        dataAdmissible &= !std::isnan(h[0][j]);
        dataAdmissible &= !std::isnan(hv[0][j]);
        dataAdmissible &= !std::isnan(hu[0][j]);
        // access i = nx + 1
        updatesAdmissible &= !std::isnan(hNetUpdatesLeft[nx + 1][j]);
        updatesAdmissible &= !std::isnan(hNetUpdatesRight[nx + 1][j]);
        updatesAdmissible &= !std::isnan(huNetUpdatesLeft[nx + 1][j]);
        updatesAdmissible &= !std::isnan(huNetUpdatesRight[nx + 1][j]);
        dataAdmissible &= !std::isnan(b[nx + 1][j]);
        dataAdmissible &= !std::isnan(h[nx + 1][j]);
        dataAdmissible &= !std::isnan(hv[nx + 1][j]);
        dataAdmissible &= !std::isnan(hu[nx + 1][j]);
    }
    /* loop over the ghost layers with indices
     * j = 0 and j = nx + 1 */
    for (int i = 1; i < nx + 1; i++) {
        // access j = 0
        updatesAdmissible &= !std::isnan(hNetUpdatesAbove[i][0]);
        updatesAdmissible &= !std::isnan(hNetUpdatesBelow[i][0]);
        updatesAdmissible &= !std::isnan(hvNetUpdatesAbove[i][0]);
        updatesAdmissible &= !std::isnan(hvNetUpdatesBelow[i][0]);
        updatesAdmissible &= !std::isnan(hNetUpdatesLeft[i][0]);
        updatesAdmissible &= !std::isnan(hNetUpdatesRight[i][0]);
        updatesAdmissible &= !std::isnan(huNetUpdatesLeft[i][0]);
        updatesAdmissible &= !std::isnan(huNetUpdatesRight[i][0]);
        dataAdmissible &= !std::isnan(b[i][0]);
        dataAdmissible &= !std::isnan(h[i][0]);
        dataAdmissible &= !std::isnan(hv[i][0]);
        dataAdmissible &= !std::isnan(hu[i][0]);
        // access j = ny + 1
        updatesAdmissible &= !std::isnan(hNetUpdatesAbove[i][ny + 1]);
        updatesAdmissible &= !std::isnan(hNetUpdatesBelow[i][ny + 1]);
        updatesAdmissible &= !std::isnan(hvNetUpdatesAbove[i][ny + 1]);
        updatesAdmissible &= !std::isnan(hvNetUpdatesBelow[i][ny + 1]);
        updatesAdmissible &= !std::isnan(hNetUpdatesLeft[i][ny + 1]);
        updatesAdmissible &= !std::isnan(hNetUpdatesRight[i][ny + 1]);
        updatesAdmissible &= !std::isnan(huNetUpdatesLeft[i][ny + 1]);
        updatesAdmissible &= !std::isnan(huNetUpdatesRight[i][ny + 1]);
        dataAdmissible &= !std::isnan(b[i][ny + 1]);
        dataAdmissible &= !std::isnan(h[i][ny + 1]);
        dataAdmissible &= !std::isnan(hv[i][ny + 1]);
        dataAdmissible &= !std::isnan(hu[i][ny + 1]);
    }
    /* is data admissible */
    admissible = updatesAdmissible && dataAdmissible;
    return admissible ;
}

//------------------------------------------------------------------------------


/**
 * Injects a random bit flip into one of the following arrays:

 *  --> b, h, hv, hu, hNetUpdatesLeft, hNetUpdatesRight,
 *      huNetUpdatesLeft, huNetUpdatesRight, hNetUpdatesAbove,
 *      hNetUpdatesBelow, hvNetUpdatesAbove, hvNetUpdatesBelow
 *
 * The array, element and bit to corrupt is selected randomly.
 */
void SWE_DimensionalSplittingMPIOverdecomp::injectRandomBitflip() {
    float* data_arrays[12] = {
        /* arrays with size (nx+2)*(ny+2) */
        b.getRawPointer(), h.getRawPointer(),
        hv.getRawPointer(), hu.getRawPointer(),
        hNetUpdatesLeft.getRawPointer(), hNetUpdatesRight.getRawPointer(),
        huNetUpdatesLeft.getRawPointer(), huNetUpdatesRight.getRawPointer(),
        /* arrays with size (nx+1)*(ny+2) */
        hNetUpdatesAbove.getRawPointer(), hNetUpdatesBelow.getRawPointer(),
        hvNetUpdatesAbove.getRawPointer(), hvNetUpdatesBelow.getRawPointer()};
    unsigned int arraySize;

    /* randomly select the data array */
    srand (static_cast <unsigned> (time(NULL)));
    int rand_index = std::rand() % 12;
    arraySize = (rand_index < 8) ? fieldSizeX : fieldSizeY;

    /* randomly select the float index and bit to flip */
    int rand_float = std::rand() % arraySize;
    int rand_bit = std::rand() % 32;

    /* float to flip */
    float target = (data_arrays[rand_index])[rand_float];
    std::bitset<32> *targetFloat = reinterpret_cast<std::bitset<32> *>(&target);

    float oldValue = (data_arrays[rand_index])[rand_float];
    /* flip the bit and write it */
    targetFloat->flip(rand_bit); (data_arrays[rand_index])[rand_float] = target;
    float newValue = (data_arrays[rand_index])[rand_float];
    assert(target == newValue);

    print_injectionRandom(rand_index, rand_float, rand_bit, oldValue, newValue);
}


/**
 * Injects a random bit flip into one of the following arrays:
 *
 *  --> hNetUpdatesLeft, hNetUpdatesRight, huNetUpdatesLeft,
 *      huNetUpdatesRight, hNetUpdatesAbove, hNetUpdatesBelow,
 *      hvNetUpdatesAbove, hvNetUpdatesBelow
 */
void SWE_DimensionalSplittingMPIOverdecomp::injectRandomBitflip_intoUpdates() {
    float* data_arrays[8] = {
        /* arrays with size (nx+2)*(ny+2) */
        hNetUpdatesLeft.getRawPointer(), hNetUpdatesRight.getRawPointer(),
        huNetUpdatesLeft.getRawPointer(), huNetUpdatesRight.getRawPointer(),
        /* arrays with size (nx+1)*(ny+2) */
        hNetUpdatesAbove.getRawPointer(), hNetUpdatesBelow.getRawPointer(),
        hvNetUpdatesAbove.getRawPointer(), hvNetUpdatesBelow.getRawPointer()};
    unsigned int arraySize;

    /* randomly select the data array */
    srand (static_cast <unsigned> (time(NULL)));
    int rand_index = std::rand() % 8;
    arraySize = (rand_index < 4) ? fieldSizeX : fieldSizeY;

    /* randomly select the float index and bit to flip */
    int rand_float = std::rand() % arraySize;
    int rand_bit = std::rand() % 32;

    /* float to flip */
    float target = (data_arrays[rand_index])[rand_float];
    std::bitset<32> *targetFloat = reinterpret_cast<std::bitset<32> *>(&target);

    float oldValue = (data_arrays[rand_index])[rand_float];
    /* flip the bit and write it */
    targetFloat->flip(rand_bit); (data_arrays[rand_index])[rand_float] = target;
    float newValue = (data_arrays[rand_index])[rand_float];
    assert(target == newValue);

    print_injectionIntoUpdates(rand_index, rand_float, rand_bit, oldValue, newValue);
}


/**
 * Injects a random bit flip into one of the following arrays:
 *
 *  --> b, h, hv, hu
 *
 * The array, element and bit to corrupt is selected randomly.
 */
void SWE_DimensionalSplittingMPIOverdecomp::injectRandomBitflip_intoData() {
    float* data_arrays[4] = { b.getRawPointer(), h.getRawPointer(),
                              hv.getRawPointer(), hu.getRawPointer()};
    unsigned int arraySize = dataArraySize;

    /* randomly select the data array */
    srand (static_cast <unsigned> (time(NULL)));
    int rand_index = std::rand() % 4;

    /* randomly select the float index and bit to flip */
    int rand_float = std::rand() % arraySize;
    int rand_bit = std::rand() % 32;

    /* this was for the propagation test */
    ///* float to flip */
    //float target = (data_arrays[3])[200];
    //std::bitset<32> *targetFloat = reinterpret_cast<std::bitset<32> *>(&target);
//
    //float oldValue = (data_arrays[3])[200];
    ///* flip the bit and write it */
    //targetFloat->flip(15); (data_arrays[3])[200] = target;
    //float newValue = (data_arrays[3])[200];
    //assert(target == newValue);

    /* float to flip */
    float target = (data_arrays[rand_index])[rand_float];
    std::bitset<32> *targetFloat = reinterpret_cast<std::bitset<32> *>(&target);

    float oldValue = (data_arrays[rand_index])[rand_float];
    /* flip the bit and write it */
    targetFloat->flip(rand_bit); (data_arrays[rand_index])[rand_float] = target;
    float newValue = (data_arrays[rand_index])[rand_float];
    assert(target == newValue);

    print_injectionIntoData(rand_index, rand_float, rand_bit, oldValue, newValue);
}

// ------------ FOR TESTS  --------------- //

/**
 * Injects a NaN into one of the following arrays:
 *
 *  --> b, h, hv, hu
 *
 * The array and its element is selected randomly.
 */
void SWE_DimensionalSplittingMPIOverdecomp::injectNaN_intoData() {
    float* data_arrays[4] = { b.getRawPointer(), h.getRawPointer(),
                              hv.getRawPointer(), hu.getRawPointer()};
    unsigned int arraySize = dataArraySize;

    /* randomly select the data array */
    srand (static_cast <unsigned> (time(NULL)));
    int rand_index = std::rand() % 4;

    /* randomly select the float index */
    int rand_float = std::rand() % arraySize;

    float oldValue = data_arrays[rand_index][rand_float];
    /* inject NaN */
    data_arrays[rand_index][rand_float] = NAN;
    float newValue = data_arrays[rand_index][rand_float];
    assert(std::isnan(newValue));

    print_injectionIntoData(rand_index, rand_float, -1, oldValue, newValue);
}

/**
 * Injects negative water height into the array h
 *
 * The element in the array is selected randomly.
 */
void SWE_DimensionalSplittingMPIOverdecomp::injectNegativeWaterHeight_intoData() {
    float* data_array = h.getRawPointer();
    unsigned int arraySize = dataArraySize;

    /* randomly select the float index */
    srand (static_cast <unsigned> (time(NULL)));
    int rand_float = std::rand() % arraySize;

    /* float to flip */
    float target = data_array[rand_float];
    std::bitset<32> *targetFloat = reinterpret_cast<std::bitset<32> *>(&target);

    float oldValue = data_array[rand_float];
    /* flip the sign bit (leftmost) and write it */
    targetFloat->flip(31); data_array[rand_float] = target;
    float newValue = data_array[rand_float];
    assert(oldValue == -newValue);
    assert(target == newValue);
    print_injectionIntoData(1, rand_float, 31, oldValue, newValue);
}

/**
 * Injects random bitflip into the array b
 *
 * The element in the array and the bit is selected randomly.
 */
void SWE_DimensionalSplittingMPIOverdecomp::injectBathymetryChange_intoData() {
    float* data_array = b.getRawPointer();
    unsigned int arraySize = dataArraySize;

    /* randomly select the float index and bit to flip */
    int rand_float = std::rand() % arraySize;
    int rand_bit = std::rand() % 32;

    /* float to flip */
    float target = data_array[rand_float];
    std::bitset<32> *targetFloat = reinterpret_cast<std::bitset<32> *>(&target);

    float oldValue = data_array[rand_float];
    /* flip the bit and write it */
    targetFloat->flip(rand_bit); data_array[rand_float] = target;
    float newValue = data_array[rand_float];
    assert(target == newValue);

    print_injectionIntoData(0, rand_float, rand_bit, oldValue, newValue);
}

/**
 * Injects infinity into one of the following arrays:
 *
 *  --> b, h, hv, hu
 *
 * The array and its element is selected randomly.
 */
void SWE_DimensionalSplittingMPIOverdecomp::injectInf_intoData() {
    float* data_arrays[4] = { b.getRawPointer(), h.getRawPointer(),
                              hv.getRawPointer(), hu.getRawPointer()};
    unsigned int arraySize = dataArraySize;

    /* randomly select the data array */
    srand (static_cast <unsigned> (time(NULL)));
    int rand_index = std::rand() % 4;

    /* randomly select the float index */
    int rand_float = std::rand() % arraySize;

    float oldValue = data_arrays[rand_index][rand_float];
    /* inject Inf */
    data_arrays[rand_index][rand_float] = INFINITY;
    float newValue = data_arrays[rand_index][rand_float];
    assert(newValue == INFINITY);

    print_injectionIntoData(rand_index, rand_float, -1, oldValue, newValue);
}

/**
 * Injects -infinity into one of the following arrays:
 *
 *  --> b, h, hv, hu
 *
 * The array and its element is selected randomly.
 */
void SWE_DimensionalSplittingMPIOverdecomp::injectnInf_intoData() {
    float* data_arrays[4] = { b.getRawPointer(), h.getRawPointer(),
                              hv.getRawPointer(), hu.getRawPointer()};
    unsigned int arraySize = dataArraySize;

    /* randomly select the data array */
    srand (static_cast <unsigned> (time(NULL)));
    int rand_index = std::rand() % 4;

    /* randomly select the float index */
    int rand_float = std::rand() % arraySize;

    float oldValue = data_arrays[rand_index][rand_float];
    /* inject -Inf */
    data_arrays[rand_index][rand_float] = -INFINITY;
    float newValue = data_arrays[rand_index][rand_float];
    assert(newValue == -INFINITY);

    print_injectionIntoData(rand_index, rand_float, -1, oldValue, newValue);
}

/**
 * Injects a big number into one of the following arrays:
 *
 *  --> b, h, hv, hu
 *
 * The array and its element is selected randomly.
 */
void SWE_DimensionalSplittingMPIOverdecomp::injectBigNumber_intoData() {
    float* data_arrays[4] = { b.getRawPointer(), h.getRawPointer(),
                              hv.getRawPointer(), hu.getRawPointer()};
    unsigned int arraySize = dataArraySize;

    /* randomly select the data array */
    srand (static_cast <unsigned> (time(NULL)));
    int rand_index = std::rand() % 4;

    ///* randomly select the float index */
    //int rand_float = std::rand() % arraySize;
//
    //float oldValue = data_arrays[3][200];
    ///* inject big number */
    //data_arrays[3][200] = 4294967296.f;
    //float newValue = data_arrays[3][200];
    //assert(newValue == 4294967296.f);

    /* randomly select the float index */
    int rand_float = std::rand() % arraySize;

    float oldValue = data_arrays[rand_index][rand_float];
    /* inject big number */
    data_arrays[rand_index][rand_float] = 4294967296.f;
    float newValue = data_arrays[rand_index][rand_float];
    assert(newValue == 4294967296.f);


    print_injectionIntoData(rand_index, rand_float, -1, oldValue, newValue);
}

/**
 * Injects a small number into one of the following arrays:
 *
 *  --> b, h, hv, hu
 *
 * The array and its element is selected randomly.
 */
void SWE_DimensionalSplittingMPIOverdecomp::injectSmallNumber_intoData() {
    float* data_arrays[4] = { b.getRawPointer(), h.getRawPointer(),
                              hv.getRawPointer(), hu.getRawPointer()};
    unsigned int arraySize = dataArraySize;

    /* randomly select the data array */
    srand (static_cast <unsigned> (time(NULL)));
    int rand_index = std::rand() % 4;

    /* randomly select the float index */
    int rand_float = std::rand() % arraySize;

    float oldValue = data_arrays[rand_index][rand_float];
    /* inject small number */
    data_arrays[rand_index][rand_float] = -4294967296.f;
    float newValue = data_arrays[rand_index][rand_float];
    assert(newValue == -4294967296.f);

    print_injectionIntoData(rand_index, rand_float, -1, oldValue, newValue);
}

//-------------------------------------------------------


/**
 * Injects a NaN into one of the following arrays:
 *
 *  --> hNetUpdatesLeft, hNetUpdatesRight, huNetUpdatesLeft,
 *      huNetUpdatesRight, hNetUpdatesAbove, hNetUpdatesBelow,
 *      hvNetUpdatesAbove, hvNetUpdatesBelow
 *
 * The array and its element is selected randomly.
 */
void SWE_DimensionalSplittingMPIOverdecomp::injectNaN_intoUpdates() {
    float* data_arrays[8] = {
        /* arrays with size (nx+2)*(ny+2) */
        hNetUpdatesLeft.getRawPointer(), hNetUpdatesRight.getRawPointer(),
        huNetUpdatesLeft.getRawPointer(), huNetUpdatesRight.getRawPointer(),
        /* arrays with size (nx+1)*(ny+2) */
        hNetUpdatesAbove.getRawPointer(), hNetUpdatesBelow.getRawPointer(),
        hvNetUpdatesAbove.getRawPointer(), hvNetUpdatesBelow.getRawPointer()};
    unsigned int arraySize;

    /* randomly select the data array */
    srand (static_cast <unsigned> (time(NULL)));
    int rand_index = std::rand() % 8;
    arraySize = (rand_index < 4) ? fieldSizeX : fieldSizeY;

    /* randomly select the float index */
    int rand_float = std::rand() % arraySize;

    float oldValue = data_arrays[rand_index][rand_float];
    /* inject NaN */
    data_arrays[rand_index][rand_float] = NAN;
    float newValue = data_arrays[rand_index][rand_float];
    assert(std::isnan(newValue));

    print_injectionIntoUpdates(rand_index, rand_float, -1, oldValue, newValue);
}

/**
 * Injects Inf into one of the following arrays:
 *
 *  --> hNetUpdatesLeft, hNetUpdatesRight, huNetUpdatesLeft,
 *      huNetUpdatesRight, hNetUpdatesAbove, hNetUpdatesBelow,
 *      hvNetUpdatesAbove, hvNetUpdatesBelow
 *
 * The array and its element is selected randomly.
 */
void SWE_DimensionalSplittingMPIOverdecomp::injectInf_intoUpdates() {
    float* data_arrays[8] = {
        /* arrays with size (nx+2)*(ny+2) */
        hNetUpdatesLeft.getRawPointer(), hNetUpdatesRight.getRawPointer(),
        huNetUpdatesLeft.getRawPointer(), huNetUpdatesRight.getRawPointer(),
        /* arrays with size (nx+1)*(ny+2) */
        hNetUpdatesAbove.getRawPointer(), hNetUpdatesBelow.getRawPointer(),
        hvNetUpdatesAbove.getRawPointer(), hvNetUpdatesBelow.getRawPointer()};
    unsigned int arraySize;

    /* randomly select the data array */
    srand (static_cast <unsigned> (time(NULL)));
    int rand_index = std::rand() % 8;
    arraySize = (rand_index < 4) ? fieldSizeX : fieldSizeY;

    /* randomly select the float index */
    int rand_float = std::rand() % arraySize;

    float oldValue = data_arrays[rand_index][rand_float];
    /* inject Inf */
    data_arrays[rand_index][rand_float] = INFINITY;
    float newValue = data_arrays[rand_index][rand_float];
    assert(newValue == INFINITY);

    print_injectionIntoUpdates(rand_index, rand_float, -1, oldValue, newValue);
}

/**
 * Injects -Inf into one of the following arrays:
 *
 *  --> hNetUpdatesLeft, hNetUpdatesRight, huNetUpdatesLeft,
 *      huNetUpdatesRight, hNetUpdatesAbove, hNetUpdatesBelow,
 *      hvNetUpdatesAbove, hvNetUpdatesBelow
 *
 * The array and its element is selected randomly.
 */
void SWE_DimensionalSplittingMPIOverdecomp::injectnInf_intoUpdates() {
    float* data_arrays[8] = {
        /* arrays with size (nx+2)*(ny+2) */
        hNetUpdatesLeft.getRawPointer(), hNetUpdatesRight.getRawPointer(),
        huNetUpdatesLeft.getRawPointer(), huNetUpdatesRight.getRawPointer(),
        /* arrays with size (nx+1)*(ny+2) */
        hNetUpdatesAbove.getRawPointer(), hNetUpdatesBelow.getRawPointer(),
        hvNetUpdatesAbove.getRawPointer(), hvNetUpdatesBelow.getRawPointer()};
    unsigned int arraySize;

    /* randomly select the data array */
    srand (static_cast <unsigned> (time(NULL)));
    int rand_index = std::rand() % 8;
    arraySize = (rand_index < 4) ? fieldSizeX : fieldSizeY;

    /* randomly select the float index */
    int rand_float = std::rand() % arraySize;

    float oldValue = data_arrays[rand_index][rand_float];
    /* inject -Inf */
    data_arrays[rand_index][rand_float] = -INFINITY;
    float newValue = data_arrays[rand_index][rand_float];
    assert(newValue == -INFINITY);

    print_injectionIntoUpdates(rand_index, rand_float, -1, oldValue, newValue);
}

/**
 * Injects a big number into one of the following arrays:
 *
 *  --> hNetUpdatesLeft, hNetUpdatesRight, huNetUpdatesLeft,
 *      huNetUpdatesRight, hNetUpdatesAbove, hNetUpdatesBelow,
 *      hvNetUpdatesAbove, hvNetUpdatesBelow
 *
 * The array and its element is selected randomly.
 */
void SWE_DimensionalSplittingMPIOverdecomp::injectBigNumber_intoUpdates() {
    float* data_arrays[8] = {
        /* arrays with size (nx+2)*(ny+2) */
        hNetUpdatesLeft.getRawPointer(), hNetUpdatesRight.getRawPointer(),
        huNetUpdatesLeft.getRawPointer(), huNetUpdatesRight.getRawPointer(),
        /* arrays with size (nx+1)*(ny+2) */
        hNetUpdatesAbove.getRawPointer(), hNetUpdatesBelow.getRawPointer(),
        hvNetUpdatesAbove.getRawPointer(), hvNetUpdatesBelow.getRawPointer()};
    unsigned int arraySize;

    /* randomly select the data array */
    srand (static_cast <unsigned> (time(NULL)));
    int rand_index = std::rand() % 8;
    arraySize = (rand_index < 4) ? fieldSizeX : fieldSizeY;

    /* randomly select the float index */
    int rand_float = std::rand() % arraySize;

    float oldValue = data_arrays[rand_index][rand_float];
    /* inject big number */
    data_arrays[rand_index][rand_float] = 4294967296.f;
    float newValue = data_arrays[rand_index][rand_float];
    assert(newValue == 4294967296.f);

    print_injectionIntoUpdates(rand_index, rand_float, -1, oldValue, newValue);
}

/**
 * Injects a small number into one of the following arrays:
 *
 *  --> hNetUpdatesLeft, hNetUpdatesRight, huNetUpdatesLeft,
 *      huNetUpdatesRight, hNetUpdatesAbove, hNetUpdatesBelow,
 *      hvNetUpdatesAbove, hvNetUpdatesBelow
 *
 * The array and its element is selected randomly.
 */
void SWE_DimensionalSplittingMPIOverdecomp::injectSmallNumber_intoUpdates() {
    float* data_arrays[8] = {
        /* arrays with size (nx+2)*(ny+2) */
        hNetUpdatesLeft.getRawPointer(), hNetUpdatesRight.getRawPointer(),
        huNetUpdatesLeft.getRawPointer(), huNetUpdatesRight.getRawPointer(),
        /* arrays with size (nx+1)*(ny+2) */
        hNetUpdatesAbove.getRawPointer(), hNetUpdatesBelow.getRawPointer(),
        hvNetUpdatesAbove.getRawPointer(), hvNetUpdatesBelow.getRawPointer()};
    unsigned int arraySize;

    /* randomly select the data array */
    srand (static_cast <unsigned> (time(NULL)));
    int rand_index = std::rand() % 8;
    arraySize = (rand_index < 4) ? fieldSizeX : fieldSizeY;

    /* randomly select the float index */
    int rand_float = std::rand() % arraySize;

    float oldValue = data_arrays[rand_index][rand_float];
    /* inject small number */
    data_arrays[rand_index][rand_float] = -4294967296.f;
    float newValue = data_arrays[rand_index][rand_float];
    assert(newValue == -4294967296.f);

    print_injectionIntoUpdates(rand_index, rand_float, -1, oldValue, newValue);
}


//------------------------------------------------------------------------------
/* Printer Implementations */

void SWE_DimensionalSplittingMPIOverdecomp::print_injectionIntoData(int rand_index, int rand_float, int rand_bit, float oldValue, float newValue) {
    std::cout.precision(std::numeric_limits<float>::max_digits10);
    std::cout << "\n............Injecting..a..bit..flip.................\n"
              << "  Corruption at array index " << rand_index << " of [b, h, hv, hu]\n"
              << "             at element index " << rand_float << "\n"
              << "             at bit no. " << rand_bit << " (indexed from right)\n"
              << "  --> old value : " << oldValue << "\n"
              << "  --> new value : " << newValue << "\n"
              << std::endl;
}

void SWE_DimensionalSplittingMPIOverdecomp::print_injectionIntoUpdates(int rand_index, int rand_float, int rand_bit, float oldValue, float newValue) {
    std::cout.precision(std::numeric_limits<float>::max_digits10);
    std::cout << "\n............Injecting..a..bit..flip.................\n"
              << "  Corruption at array index " << rand_index << " of "
              << "[hNetUpdatesLeft, hNetUpdatesRight, huNetUpdatesLeft, huNetUpdatesRight, "
              << "hNetUpdatesAbove, hNetUpdatesBelow, hvNetUpdatesAbove, hvNetUpdatesBelow]\n"
              << "             at element index " << rand_float << "\n"
              << "             at bit no. " << rand_bit << " (indexed from right)\n"
              << "  --> old value : " << oldValue << "\n"
              << "  --> new value : " << newValue << "\n"
              << std::endl;
}

void SWE_DimensionalSplittingMPIOverdecomp::print_injectionRandom(int rand_index, int rand_float, int rand_bit, float oldValue, float newValue) {
    std::cout.precision(std::numeric_limits<float>::max_digits10);
    std::cout << "\n............Injecting..a..bit..flip.................\n"
              << "  Corruption at array index " << rand_index << " of [b, h, hv, hu, "
              << "hNetUpdatesLeft, hNetUpdatesRight, huNetUpdatesLeft, huNetUpdatesRight, "
              << "hNetUpdatesAbove, hNetUpdatesBelow, hvNetUpdatesAbove, hvNetUpdatesBelow]\n"
              << "             at element index " << rand_float << "\n"
              << "             at bit no. " << rand_bit << " (indexed from right)\n"
              << "  --> old value : " << oldValue << "\n"
              << "  --> new value : " << newValue << "\n"
              << std::endl;
}
//------------------------------------------------------------------------------


MPI_Datatype SWE_DimensionalSplittingMPIOverdecomp::getBlockMPIType() { return blockData_t; }
