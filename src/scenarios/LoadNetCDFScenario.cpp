#include "LoadNetCDFScenario.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <sstream>

// Reads the netCDF files and creates scenario from it
/*TODO the values in the outer ghost-layer are not stored in the files
 *but are needed for creating the scenario. I just use the adjacent values although that is not 100%
 *correct. See last function and comment.
 */
SWE_LoadNetCdfScenario::SWE_LoadNetCdfScenario(std::string& i_file,
                                               float i_endTime,
                                               std::vector<BoundaryType>& i_boundaryTypes,
                                               std::vector<float>& i_boundaryPositions)
    : endTime(i_endTime),
      boundaryTypes(i_boundaryTypes),
      boundaryPositions(i_boundaryPositions)
{
    std::string file = i_file + ".nc";
    int error = nc_open(file.c_str(), NC_NOWRITE, &dataFile);

    if (error != NC_NOERR)
    {
        fprintf(stderr, "Reading error: %s\n", nc_strerror(error));
        assert(false);
        return;
    }
    // TODO Error Handling
    nc_inq_varid(dataFile, "time", &timeVar);
    nc_inq_varid(dataFile, "x", &xVar);
    nc_inq_varid(dataFile, "y", &yVar);
    nc_inq_varid(dataFile, "h", &hVar);
    nc_inq_varid(dataFile, "hu", &huVar);
    nc_inq_varid(dataFile, "hv", &hvVar);
    nc_inq_varid(dataFile, "b", &bVar);

    nc_inq_dimid(dataFile, "time", &timeDim);
    nc_inq_dimid(dataFile, "x", &xDim);
    nc_inq_dimid(dataFile, "y", &yDim);

    nc_inq_dimlen(dataFile, timeDim, &numTimesteps);
    nc_inq_dimlen(dataFile, xDim, &xLen);
    nc_inq_dimlen(dataFile, yDim, &yLen);

    size_t pos = numTimesteps - 1;
    size_t index[] = {0};
    float tempX1, tempX2, tempY1, tempY2;
    nc_get_var1_float(dataFile, xVar, index, &tempX1);
    nc_get_var1_float(dataFile, yVar, index, &tempY1);
    index[0] = 1;
    nc_get_var1_float(dataFile, xVar, index, &tempX2);
    nc_get_var1_float(dataFile, yVar, index, &tempY2);

    offsetX = tempX1; // center of first cell
    offsetY = tempY1; // center of first cell
    dX = tempX2 - tempX1;
    dY = tempY2 - tempY1;

    float* x = new float[xLen];
    float* y = new float[yLen];

    h = new std::vector<float>(xLen * yLen);
    hu = new std::vector<float>(xLen * yLen);
    hv = new std::vector<float>(xLen * yLen);
    b = new std::vector<float>(xLen * yLen);

    size_t start[] = {pos, 0, 0}; // time, y, x
    size_t count[] = {1, yLen, xLen};

    nc_get_vara_float(dataFile, yVar, &start[1], &count[1], y);
    nc_get_vara_float(dataFile, xVar, &start[2], &count[2], x);
    nc_get_vara_float(dataFile, bVar, &start[1], &count[1], b->data());
    nc_get_vara_float(dataFile, hVar, start, count, h->data());
    nc_get_vara_float(dataFile, huVar, start, count, hu->data());
    nc_get_vara_float(dataFile, hvVar, start, count, hv->data());

    xVec.resize(xLen);
    yVec.resize(yLen);
    std::copy(x, x + xLen, xVec.begin());
    std::copy(y, y + yLen, yVec.begin());

    delete (x);
    delete (y);
}

SWE_LoadNetCdfScenario::~SWE_LoadNetCdfScenario()
{
    // delete(x);
    // delete(y);

    delete (h);
    delete (hu);
    delete (hv);
    delete (b);
}

float SWE_LoadNetCdfScenario::getWaterHeight(float x, float y) { return h->at(toValidCoords(x, y)); }

float SWE_LoadNetCdfScenario::getVeloc_u(float x, float y)
{
    if ((h->at(toValidCoords(x, y)) != 0.0f))
    {
        return hu->at(toValidCoords(x, y)) / (h->at(toValidCoords(x, y)));
    }
    return 0.0f;
}

float SWE_LoadNetCdfScenario::getVeloc_v(float x, float y)
{
    if ((h->at(toValidCoords(x, y)) != 0.0f))
    {
        return hv->at(toValidCoords(x, y)) / (h->at(toValidCoords(x, y)));
    }
    return 0.0f;
}

float SWE_LoadNetCdfScenario::getMomentum_u(float x, float y) { return hu->at(toValidCoords(x, y)); }

float SWE_LoadNetCdfScenario::getMomentum_v(float x, float y) { return hv->at(toValidCoords(x, y)); }

float SWE_LoadNetCdfScenario::getBathymetry(float x, float y) { return b->at(toValidCoords(x, y)); }

float SWE_LoadNetCdfScenario::waterHeightAtRest() { return 0.0; }

float SWE_LoadNetCdfScenario::endSimulation() { return endTime; }

BoundaryType SWE_LoadNetCdfScenario::getBoundaryType(Boundary edge) { return boundaryTypes.at(static_cast<int>(edge)); }

float SWE_LoadNetCdfScenario::getBoundaryPos(Boundary edge) { return boundaryPositions.at(static_cast<int>(edge)); }

// For border regions this code just returns clostest known value.
// Will only happen when loadig bathemetry but better solution should be found
int SWE_LoadNetCdfScenario::toValidCoords(float x, float y)
{
    x -= offsetX; // Abstand zur 1. Zelle
    y -= offsetY;
    int indexX = std::round((x - .5f * dX) / (dX));
    int indexY = std::round((y - .5f * dY) / (dY));
    indexX = std::max(0, std::min(indexX, (int)xLen - 1));
    indexY = std::max(0, std::min(indexY, (int)yLen - 1));

    return (indexY * xLen + indexX) % (xLen * yLen);
}
