#include "io/Writer.hpp"

#include <memory>
#include <vector>

#include "cassert"
#include "io/NetCDFWriter.hpp"
#include "scenarios/Scenario.hpp"

std::shared_ptr<io::Writer> io::Writer::createWriterInstance(std::string& outputName, std::string& backupName,
                                                             const Float2D& bathymetry,
                                                             const BoundarySize& boundarySize, int nX, int nY, float dX,
                                                             float dY, float offsetX, float offsetY, float originX,
                                                             float originY, int flush, bool existingFile) {
    // construct a NetCdfWriter
    auto writer = std::make_shared<io::NetCdfWriter>(outputName, backupName, bathymetry, boundarySize, nX, nY, dX, dY,
                                                     existingFile, originX, originY, flush);
    return writer;
}

void io::Writer::initMetadataFile(std::string meatadataName, float totalTime, int ranks, int gridSizeX, int gridSizeY,
                                  int numCheckpoints, const std::vector<BoundaryType>& boundaryTypes,
                                  const std::vector<float>& boundaryPositions) {
    assert(boundaryTypes.size() == 4);
    assert(boundaryPositions.size() == 4);
    std::ofstream metadataFile(meatadataName);
    metadataFile << "total_time=" << std::to_string(totalTime) << "\n";
    metadataFile << "current_time="
                 << "0.0"
                 << "\n";
    metadataFile << "ranks=" << std::to_string(ranks) << " \n";
    metadataFile << "grid_size_x=" << std::to_string(gridSizeX) << "\n";
    metadataFile << "grid_size_y=" << std::to_string(gridSizeY) << "\n";
    for (int i = 0; i <= 3; i++) {
        metadataFile << "boundary_type_" << std::to_string(i) << "=" << std::to_string(boundaryTypes[i]) << "\n";
        metadataFile << "boundary_position_" << std::to_string(i) << "=" << std::to_string(boundaryPositions[i])
                     << "\n";
    }

    metadataFile << "checkpoints_unfinished=" << std::to_string(numCheckpoints) << std::endl;
    metadataFile.close();
}

void io::Writer::updateMetadataFile(std::string metadataName, float currentTime, int checkpointsLeft) {
    replaceConfigureFileValue(metadataName, "current_time", std::to_string(currentTime));
    replaceConfigureFileValue(metadataName, "checkpoints_unfinished", std::to_string(checkpointsLeft));
}
