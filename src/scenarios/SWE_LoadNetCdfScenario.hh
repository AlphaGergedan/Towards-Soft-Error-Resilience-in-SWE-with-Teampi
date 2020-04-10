#ifndef SWELOADNETCDFSCENARIO_HPP
#define  SWELOADNETCDFSCENARIO_HPP

#include <vector>
#include <string>
#include <netcdf.h>

#include "SWE_Scenario.hh"
#include "tools/help.hh"


class SWE_LoadNetCdfScenario : public SWE_Scenario{ 
    public:

    SWE_LoadNetCdfScenario(std::string &file, float endTime, 
                           std::vector<BoundaryType> &boundaryTypes,
                           std::vector<float> &boundaryPositions);

    ~SWE_LoadNetCdfScenario();

    float getWaterHeight(float x, float y);
    float getVeloc_u(float x, float y);
    float getVeloc_v(float x, float y);
    float getBathymetry(float x, float y);
    
    float waterHeightAtRest();

    float endSimulation();
    
    BoundaryType getBoundaryType(BoundaryEdge edge);
    float getBoundaryPos(BoundaryEdge edge);

    private:
        int validCoords(float x, float y, int *indexX, int *indexY);

    private:
        float endTime;
        int dataFile;
        size_t numTimesteps, xLen, yLen;
        int timeDim, xDim, yDim;
        int timeVar, xVar, yVar;
        int hVar, huVar, hvVar;
        int bVar;
        std::vector<float> xVec, yVec;
        Float2D *h, *hu, *hv, *b;
        std::string fileName;
        std::vector<BoundaryType> boundaryTypes;
        std::vector<float> boundaryPositions;
    
    
};

#endif