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

    virtual ~SWE_LoadNetCdfScenario();

    virtual float getWaterHeight(float x, float y);
    virtual float getVeloc_u(float x, float y);
    virtual float getVeloc_v(float x, float y);
    virtual float getBathymetry(float x, float y);
    virtual float getMomentum_u(float x, float y);
    virtual float getMomentum_v(float x, float y);

    virtual float waterHeightAtRest();

    virtual float endSimulation();
    
    virtual BoundaryType getBoundaryType(BoundaryEdge edge);
    virtual float getBoundaryPos(BoundaryEdge edge);

    private:
        int toValidCoords(float x, float y);

    private:
        float endTime, dX, dY;
        int dataFile;
        size_t numTimesteps, xLen, yLen;
        int timeDim, xDim, yDim;
        int timeVar, xVar, yVar;
        int hVar, huVar, hvVar;
        int bVar;
        std::vector<float> xVec, yVec;
        std::vector<float> *h, *hu, *hv, *b;
        std::string fileName;
        std::vector<BoundaryType> boundaryTypes;
        std::vector<float> boundaryPositions;
    
    
};

#endif