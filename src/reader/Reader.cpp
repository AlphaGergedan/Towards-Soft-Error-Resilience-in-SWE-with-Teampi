#include "Reader.hh"
#include "tools/help.hh"
#include "tools/Logger.hh"
#include <fstream>

//Netcdf file extension hardcoded
//Reads corresponding files and creates new scenatrio from it

io::Reader::Reader(std::string i_backupFilename, std::string i_outputFilename, int i_rank, int i_mpiSizeCurrent, int i_blockPosX, int i_blockPosY):
                    backupFilename(i_backupFilename), outputFilename(i_outputFilename),
                    rank(i_rank), mpiCurrentSize(i_mpiSizeCurrent),
                    blockPosX(i_blockPosX), blockPosY(i_blockPosY){
    
    readMetadataFile(i_backupFilename + "_metadata");

    std::string backup = generateBaseFileName(backupFilename, blockPosX, blockPosY);
    std::string output = generateBaseFileName(outputFilename, blockPosX, blockPosY);
    tools::Logger::logger.printString(backup +" " + output);
    std::ifstream src(backup + ".nc", std::ios::binary);
	std::ofstream out(output + ".nc", std::ios::binary);
	out << src.rdbuf();
    out.close();
    scenario = new SWE_LoadNetCdfScenario(backup, totalTime, boundaryTypes, boudaryPositions);

    
}

io::Reader::~Reader(){
    
}

//Todo Errorhandling when file not found, sorry :(
//Reads the Metadata file and sets corresponding values
void io::Reader::readMetadataFile(std::string filename){
    mpiExpectedSize = std::stoi(readConfigureFileValue(filename, "ranks"));
    sizeX = std::stoi(readConfigureFileValue(filename, "grid_size_x"));
    sizeY = std::stoi(readConfigureFileValue(filename, "grid_size_y"));
    totalTime = std::stof(readConfigureFileValue(filename, "total_time"));
    currentTime = std::stof(readConfigureFileValue(filename, "current_time"));
    remainingCheckpoints = std::stoi(readConfigureFileValue(filename, "checkpoints_unfinished"));
    for(int i = 0; i <= 3; i++){
        std::string boundaryType = "boundary_type_" + std::to_string(i);
        std::string boudaryPosition = "boundary_position_" + std::to_string(i);

        BoundaryType type = static_cast<BoundaryType> (std::stoi(readConfigureFileValue(filename, boundaryType)));
        float pos = std::stof(readConfigureFileValue(filename, boudaryPosition));

        boundaryTypes.push_back(type);
        boudaryPositions.push_back(pos);
    }
    
}

 float io::Reader::getRemainingTime(){
        return totalTime - currentTime;
}
        
float io::Reader::getCurrentTime(){
    return currentTime;
}
int io::Reader::getRemainingCheckpoints(){
    return remainingCheckpoints;
}

SWE_Scenario *io::Reader::getScenario(){
    return scenario;
}

int io::Reader::getGridSizeX(){
    return sizeX;
}

int io::Reader::getGridSizeY(){
    return sizeY;
}