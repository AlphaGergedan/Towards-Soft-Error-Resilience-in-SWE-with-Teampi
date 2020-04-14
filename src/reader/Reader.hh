#ifndef READER_HH_
#define READER_HH_

#include <string>
#include <vector>

#include "scenarios/SWE_LoadNetCdfScenario.hh"

namespace io{
    class Reader;
}

class io::Reader{
    public:
        Reader(std::string backupFilename, int rank, int mpiSize, int blockPosX, int blockPosY);
        ~Reader();
        float getRemainingTime();
        int getGridSizeX();
        int getGridSizeY();
        int getRemainingCheckpoints();
        float getCurrentTime();
        SWE_Scenario *getScenario();

    protected:
        std::string backupFilename;
        SWE_LoadNetCdfScenario *scenario;
        int rank, mpiCurrentSize, mpiExpectedSize;
        int blockPosX, blockPosY;
        int remainingCheckpoints;
        int sizeX;
        int sizeY;
        float totalTime;
        float currentTime;
        std::vector<BoundaryType> boundaryTypes;
        std::vector<float> boudaryPositions;

    protected:
        void readMetadataFile(std::string filename);
        
};

#endif