#include <algorithm>
#include <cassert>
#include <cmath>

#include <sstream>

#include "SWE_LoadNetCdfScenario.hh"
#include "tools/Logger.hh"



SWE_LoadNetCdfScenario::SWE_LoadNetCdfScenario(std::string &i_file, float i_endTime, 
                           std::vector<BoundaryType> &i_boundaryTypes,
                           std::vector<float> &i_boundaryPositions):endTime(i_endTime),
                           boundaryTypes(i_boundaryTypes),
                           boundaryPositions(i_boundaryPositions){


    tools::Logger::logger.printString("Loading NetCdf");
    std::string file = i_file + ".nc";
    int error = nc_open(file.c_str(), NC_NOWRITE, &dataFile);
    if(error != NC_NOERR){
        assert(false);
        return;
    }
    //TODO Error Handling
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
    size_t index[] = {pos};
    nc_get_var1_float(dataFile, timeVar, index, &endTime);
    index[0] = 0;
    nc_get_var1_float(dataFile, xVar, index, &dX);
    nc_get_var1_float(dataFile, yVar, index, &dY);

    tools::Logger::logger.printString("dX: " + std::to_string(dX));
    tools::Logger::logger.printString("dY: " + std::to_string(dY));
    float *x = new float[xLen];
    float *y = new float[yLen];

    h = new Float2D(xLen, yLen);
    hu = new Float2D(xLen, yLen);
    hv = new Float2D(xLen, yLen);
    b = new Float2D(xLen, yLen);

    size_t index3D[] = {pos, 0, 0};
    size_t count3D[] = {1, yLen, yLen};

    nc_get_vara_float(dataFile, yVar, &index3D[2], &count3D[2], y);
    count3D[2] = xLen;
    nc_get_vara_float(dataFile, xVar, &index3D[2], &count3D[2], x);
    nc_get_vara_float(dataFile, bVar, &index3D[1], &count3D[1], b->elemVector());
    nc_get_vara_float(dataFile, hVar, index3D, count3D, h->elemVector());
    nc_get_vara_float(dataFile, huVar, index3D, count3D, hu->elemVector());
    nc_get_vara_float(dataFile, hvVar, index3D, count3D, hv->elemVector());

    xVec.resize(xLen);
    yVec.resize(yLen);
    std::copy(x, x+xLen, xVec.begin());
    std::copy(y, y+yLen, yVec.begin());

    for(size_t i = 0; i < xVec.size(); i++){
        std::stringstream xStream;
        xStream << std::setprecision(15) << xVec.at(i);
        //tools::Logger::logger.printString("x: " + xStream.str());
    }

    for(size_t i = 0; i < yVec.size(); i++){
        std::stringstream yStream;
        yStream << std::setprecision(15) << yVec.at(i);
        //tools::Logger::logger.printString("y: " + yStream.str());
    }

    delete(x);
    delete(y);
    
}

SWE_LoadNetCdfScenario::~SWE_LoadNetCdfScenario(){
    //delete(x);
    //delete(y);
    
    
    delete(h);
    delete(hu);
    delete(hv);
    delete(b);
}



float SWE_LoadNetCdfScenario::getWaterHeight(float x, float y){
    int indexX, indexY;
    if(toValidCoords(x,y,&indexX, &indexY)){
        return (*h)[indexX][indexY];
    }
    return 0.0;
}

float SWE_LoadNetCdfScenario::getVeloc_u(float x, float y){
    int indexX, indexY;
    if(toValidCoords(x,y,&indexX, &indexY)){
        return (*hu)[indexX][indexY];
    }
    return 0.0;
}
    
float SWE_LoadNetCdfScenario::getVeloc_v(float x, float y){
    int indexX, indexY;
    if(toValidCoords(x,y,&indexX, &indexY)){
        return (*hv)[indexX][indexY];
    }
    return 0.0;
}
    
float SWE_LoadNetCdfScenario::getBathymetry(float x, float y){
    int indexX, indexY;
    if(toValidCoords(x,y,&indexX, &indexY)){
        return (*b)[indexX][indexY];
    }
    return 0.0;
}
    
float SWE_LoadNetCdfScenario::waterHeightAtRest(){
    return 0.0;
}

float SWE_LoadNetCdfScenario::endSimulation(){
    return endTime;
}
    
BoundaryType SWE_LoadNetCdfScenario::getBoundaryType(BoundaryEdge edge){
    return boundaryTypes.at(static_cast<int>(edge));
}

float SWE_LoadNetCdfScenario::getBoundaryPos(BoundaryEdge edge){
    return boundaryPositions.at(static_cast<int>(edge));
}

//Change vector and find to something better
int SWE_LoadNetCdfScenario::toValidCoords(float x, float y, int *indexX, int *indexY){
    std::stringstream xStream, yStream, ixStream, iyStream;
    xStream << std::setprecision(15) << x;
    yStream << std::setprecision(15) << y;
    //tools::Logger::logger.printString("accessing x:" + xStream.str());
    //tools::Logger::logger.printString("accessing y:" + yStream.str());
    *indexX = std::round((x-dX)/(2*dX));
    *indexY = std::round((y-dY)/(2*dY));
    if(*indexX >= (int) xLen || *indexY >= (int) yLen){
        tools::Logger::logger.printString("this shouldnt happen");
    }

    if(*indexX < 0 || *indexY < 0){ tools::Logger::logger.printString("this shouldnt happen");
         tools::Logger::logger.printString("this shouldnt happen2");
    }
    ixStream << std::setprecision(15) << *indexX;
    iyStream << std::setprecision(15) << *indexY;

    //tools::Logger::logger.printString("returning x:" + ixStream.str());
    //tools::Logger::logger.printString("returning y:" + iyStream.str());
    return 1;

}