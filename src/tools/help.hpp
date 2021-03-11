/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader, Kaveh Rahnema
 * @author Sebastian Rettenberger
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * TODO
 */

#ifndef SWE_HELP_HPP
#define SWE_HELP_HPP

#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

//-------- Methods for Visualistion of Results --------

/**
 * generate output filenames for the single-SWE_Block version
 * (for serial and OpenMP-parallelised versions that use only a
 *  single SWE_Block - one output file is generated per checkpoint)
 *
 *  @deprecated
 */
inline std::string generateFileName(std::string outputTeamName, int timeStep)
{
    std::ostringstream FileName;
    FileName << outputTeamName << timeStep << ".vtk";
    return FileName.str();
};

/**
 * Generates an output file name for a multiple SWE_Block version based on the ordering of the blocks.
 *
 * @param i_baseName base name of the output.
 * @param i_blockPositionX position of the SWE_Block in x-direction.
 * @param i_blockPositionY position of the SWE_Block in y-direction.
 * @param i_fileExtension file extension of the output file.
 * @return
 *
 * @deprecated
 */
inline std::string generateFileName(std::string i_baseName,
                                    int i_blockPositionX,
                                    int i_blockPositionY,
                                    std::string i_fileExtension = ".nc")
{
    std::ostringstream l_fileName;

    l_fileName << i_baseName << "_" << i_blockPositionX << i_blockPositionY << i_fileExtension;
    return l_fileName.str();
};

/**
 * generate output filename for the multiple-SWE_Block version
 * (for serial and parallel (OpenMP and MPI) versions that use
 *  multiple SWE_Blocks - for each block, one output file is
 *  generated per checkpoint)
 *
 *  @deprecated
 */
inline std::string generateFileName(std::string outputTeamName,
                                    int timeStep,
                                    int block_X,
                                    int block_Y,
                                    std::string i_fileExtension = ".vts")
{
    std::ostringstream FileName;
    FileName << outputTeamName << "_" << block_X << "_" << block_Y << "_" << timeStep << i_fileExtension;
    return FileName.str();
};

/**
 * Generates an output file name for a multiple SWE_Block version based on the ordering of the blocks.
 *
 * @param i_baseName base name of the output.
 * @param i_blockPositionX position of the SWE_Block in x-direction.
 * @param i_blockPositionY position of the SWE_Block in y-direction.
 *
 * @return the output filename <b>without</b> timestep information and file extension
 */
inline std::string genTeamPosName(std::string const& i_baseName, int i_blockPositionX, int i_blockPositionY)
{
    std::ostringstream l_fileName;

    l_fileName << i_baseName << '_' << i_blockPositionX << '_' << i_blockPositionY;
    return l_fileName.str();
}

/**
 * generate output filename for the ParaView-Container-File
 * (to visualize multiple SWE_Blocks per checkpoint)
 */
inline std::string generateContainerFileName(std::string outputTeamName, int timeStep)
{
    std::ostringstream FileName;
    FileName << outputTeamName << "_" << timeStep << ".pvts";
    return FileName.str();
};

inline int replaceConfigureFileValue(std::string filename, std::string key, std::string newValue)
{
    std::string tempName = filename + "_temp";
    std::ifstream fileIn;
    std::ofstream fileOut;
    int found = 0;

    fileIn.open(filename);
    fileOut.open(tempName);
    while (!fileIn.eof())
    {
        std::string line;
        std::string currentKey;
        size_t pos;

        std::getline(fileIn, line);
        pos = line.find("=");
        if (pos == std::string::npos)
        {
            fileOut << line;
            continue;
        }
        currentKey = line.substr(0, pos);
        if (key == currentKey)
        {
            found++;
            fileOut << currentKey << "=" << newValue << "\n";
        }
        else
        {
            fileOut << line << "\n";
        }
    }

    std::remove(filename.c_str());
    std::rename(tempName.c_str(), filename.c_str());
    return found;
}

inline std::string readConfigureFileValue(const std::string& filename, const std::string& key)
{
    std::ifstream fileIn;
    fileIn.open(filename);
    while (!fileIn.eof())
    {
        std::string line;
        std::string currentKey;
        size_t pos;

        std::getline(fileIn, line);
        pos = line.find("=");

        if (pos == std::string::npos)
            continue;

        currentKey = line.substr(0, pos);
        if (key == currentKey)
        {
            return line.substr(pos + 1);
        }
    }
    return "";
}

#endif // SWE_HELP_HPP
