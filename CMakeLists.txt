cmake_minimum_required(VERSION 3.13)

project(swe_with_teaMPI_and_fault_tolerance
    DESCRIPTION "SWE with teaMPI with different fault tolerance techniques."
    LANGUAGES CXX
)

# teampi uses 11 originally
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#*******************************************************************************
# solver options
option(SOLVER "Wave Propagation Solver" "hlle")
set(SOLVERS "fwave, hlle, augrie")

# This displays a combo box in CMake GUIs
set_property(CACHE SOLVER PROPERTY STRINGS ${SOLVERS})

option(USE_DEBUG "Debug compile mode" OFF)
option(USE_PROFILING "Compiles for profiling" OFF)

#*******************************************************************************
# sources

# Baseline for hard error resilience
add_executable(swe_checkpointRestart)                   # --> baseline version for hard resilience

# Hard error resilience with task sharing to reduce redundancy
add_executable(swe_hardRes_wTaskSharing)                # --> old version

# Combination of hard and soft error resilience using warmspares + hashes
add_executable(swe_softRes_hardRes_hashes) # --> soft + hard error resilience

########################### SOFT RESILIENCE METHODS #############################

# METHOD 1   : BASELINE, NO ERROR RESILIENCE (NoRes)
add_executable(swe_noRes)                               # --> baseline model for soft error resilience

# METHOD 2   : SOFT ERROR DETECTION USING HASHES (Hashes)
add_executable(swe_softRes_hashes)                             # --> soft error detection

# METHOD 3 : SOFT ERROR RESILIENCE WITH ADMISSIBILITY CHECKS AND TASK SHARING + USING SHARED TASKS IMMEDIATELY (Sharing)
add_executable(swe_softRes_admiss_useShared)

# METHOD 4   : SOFT ERROR RESILIENCE WITH ADMISSIBILITY CHECKS AND TASK SHARING + ONLY USING SHARED TASKS IF SDC IS DETECTED (Redundant)
add_executable(swe_softRes_admiss_redundant)


target_sources(swe_noRes PRIVATE src/tolerance/swe_noRes.cpp)
target_sources(swe_hardRes_wTaskSharing PRIVATE src/constants.hpp)

target_sources(swe_hardRes_wTaskSharing PRIVATE src/tolerance/swe_hardRes_wTaskSharing.cpp)
target_sources(swe_hardRes_wTaskSharing PRIVATE src/constants.hpp)

target_sources(swe_softRes_hashes PRIVATE src/tolerance/swe_softRes_hashes.cpp)
target_sources(swe_softRes_hashes PRIVATE src/constants.hpp)

target_sources(swe_softRes_hardRes_hashes PRIVATE src/tolerance/swe_softRes_hardRes_hashes.cpp)
target_sources(swe_softRes_hardRes_hashes PRIVATE src/constants.hpp)

target_sources(swe_softRes_admiss_useShared PRIVATE src/tolerance/swe_softRes_admiss_useShared.cpp)
target_sources(swe_softRes_admiss_useShared PRIVATE src/constants.hpp)

target_sources(swe_softRes_admiss_redundant PRIVATE src/tolerance/swe_softRes_admiss_redundant.cpp)
target_sources(swe_softRes_admiss_redundant PRIVATE src/constants.hpp)

target_sources(swe_checkpointRestart PRIVATE src/tolerance/swe_checkpointRestart.cpp)
target_sources(swe_checkpointRestart PRIVATE src/constants.hpp)


# Each subdirectory has it's own CMakeLists.txt to target their sources
add_subdirectory(src/types)
add_subdirectory(src/tools)
add_subdirectory(src/io)
add_subdirectory(src/scenarios)
add_subdirectory(src/solvers)
add_subdirectory(src/blocks)


# Include the src directory
target_include_directories(swe_noRes PRIVATE src)
target_include_directories(swe_hardRes_wTaskSharing PRIVATE src)
target_include_directories(swe_softRes_hashes PRIVATE src)
target_include_directories(swe_softRes_hardRes_hashes PRIVATE src)
target_include_directories(swe_softRes_admiss_useShared PRIVATE src)
target_include_directories(swe_softRes_admiss_redundant PRIVATE src)
target_include_directories(swe_checkpointRestart PRIVATE src)

#*******************************************************************************
# libraries

# NetCDF
include(cmake/FindNetCDF.cmake)
if(NOT NETCDF_FOUND)
  message(FATAL_ERROR "NetCDF could not be found!")
endif()
target_include_directories(swe_noRes SYSTEM PRIVATE NETCDF_INCLUDE_DIRS)
target_link_libraries(swe_noRes PRIVATE ${NETCDF_LIBRARIES})

target_include_directories(swe_hardRes_wTaskSharing SYSTEM PRIVATE NETCDF_INCLUDE_DIRS)
target_link_libraries(swe_hardRes_wTaskSharing PRIVATE ${NETCDF_LIBRARIES})

target_include_directories(swe_softRes_hashes SYSTEM PRIVATE NETCDF_INCLUDE_DIRS)
target_link_libraries(swe_softRes_hashes PRIVATE ${NETCDF_LIBRARIES})

target_include_directories(swe_softRes_hardRes_hashes SYSTEM PRIVATE NETCDF_INCLUDE_DIRS)
target_link_libraries(swe_softRes_hardRes_hashes PRIVATE ${NETCDF_LIBRARIES})

target_include_directories(swe_softRes_admiss_useShared SYSTEM PRIVATE NETCDF_INCLUDE_DIRS)
target_link_libraries(swe_softRes_admiss_useShared PRIVATE ${NETCDF_LIBRARIES})

target_include_directories(swe_softRes_admiss_redundant SYSTEM PRIVATE NETCDF_INCLUDE_DIRS)
target_link_libraries(swe_softRes_admiss_redundant PRIVATE ${NETCDF_LIBRARIES})

target_include_directories(swe_checkpointRestart SYSTEM PRIVATE NETCDF_INCLUDE_DIRS)
target_link_libraries(swe_checkpointRestart PRIVATE ${NETCDF_LIBRARIES})

# ULFM
# Set this to the directory, that contains the "bin" folder which contains ULFM's MPI executables
#set(MPI_HOME "/home/${USERNAME}/ulfm2/build/bin" CACHE PATH "Base directory for MPI" FORCE)
# or specify -DMPI_HOME="build-directory" while running cmake
find_package(MPI REQUIRED)
if (MPI_FOUND)
    target_link_libraries(swe_noRes PUBLIC MPI::MPI_CXX)
    target_link_libraries(swe_hardRes_wTaskSharing PUBLIC MPI::MPI_CXX)
    target_link_libraries(swe_softRes_hashes PUBLIC MPI::MPI_CXX)
    target_link_libraries(swe_softRes_hardRes_hashes PUBLIC MPI::MPI_CXX)
    target_link_libraries(swe_softRes_admiss_useShared PUBLIC MPI::MPI_CXX)
    target_link_libraries(swe_softRes_admiss_redundant PUBLIC MPI::MPI_CXX)
    target_link_libraries(swe_checkpointRestart PUBLIC MPI::MPI_CXX)
else()
    message(FATAL_ERROR "MPI not found!")
endif()


#*******************************************************************************
# teaMPI library

add_library(tmpi SHARED
    teampi-soft-error-resilience/lib/CommStats.cpp
    teampi-soft-error-resilience/lib/RankControl.cpp
    teampi-soft-error-resilience/lib/Wrapper.cpp
    teampi-soft-error-resilience/lib/Rank.cpp
    teampi-soft-error-resilience/lib/teaMPI.cpp
    teampi-soft-error-resilience/lib/Timing.cpp
    teampi-soft-error-resilience/lib/ErrorHandling/WarmSpareStrategy.cpp
    teampi-soft-error-resilience/lib/ErrorHandling/RespawnProcStrategy.cpp
    teampi-soft-error-resilience/lib/ErrorHandling/KillTeamStrategy.cpp
)

# Link tmpi, NoRes and swe_checkpointRestart don't need it
target_link_libraries(swe_hardRes_wTaskSharing PRIVATE tmpi)
target_link_libraries(swe_softRes_hashes PRIVATE tmpi)
target_link_libraries(swe_softRes_hardRes_hashes PRIVATE tmpi)
target_link_libraries(swe_softRes_admiss_useShared PRIVATE tmpi)
target_link_libraries(swe_softRes_admiss_redundant PRIVATE tmpi)

# include path
target_include_directories(swe_hardRes_wTaskSharing PRIVATE teampi-soft-error-resilience/include)
target_include_directories(swe_softRes_hashes PRIVATE teampi-soft-error-resilience/include)
target_include_directories(swe_softRes_hardRes_hashes PRIVATE teampi-soft-error-resilience/include)
target_include_directories(swe_softRes_admiss_useShared PRIVATE teampi-soft-error-resilience/include)
target_include_directories(swe_softRes_admiss_redundant PRIVATE teampi-soft-error-resilience/include)

#*******************************************************************************
# Warnings & Flags

if (USE_DEBUG)
    target_compile_options(tmpi PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>
    )
    target_compile_options(swe_noRes PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>
    )
    target_compile_options(swe_hardRes_wTaskSharing PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>
    )
    target_compile_options(swe_softRes_hashes PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>
    )
    target_compile_options(swe_softRes_hardRes_hashes PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>
    )
    target_compile_options(swe_softRes_admiss_useShared  PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>
    )
    target_compile_options(swe_softRes_admiss_redundant PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>
    )
    target_compile_options(swe_checkpointRestart PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>
    )
elseif (USE_PROFILING)
    target_compile_options(tmpi PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-g -O2>
    )
    target_compile_options(swe_noRes PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-g -O2>
    )
    target_compile_options(swe_hardRes_wTaskSharing PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-g -O2>
    )
    target_compile_options(swe_softRes_hashes PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-g -O2>
    )
    target_compile_options(swe_softRes_hardRes_hashes PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-g -O2>
    )
    target_compile_options(swe_softRes_admiss_useShared PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-g -O2>
    )
    target_compile_options(swe_softRes_admiss_redundant PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-g -O2>
    )
    target_compile_options(swe_checkpointRestart PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-g -O2>
    )
else () # release
    target_compile_options(tmpi PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>: -O3 >
    )
    target_compile_options(swe_noRes PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>: -O3 >
    )
    target_compile_options(swe_hardRes_wTaskSharing PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>: -O3 >
    )
    target_compile_options(swe_softRes_hashes PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>: -O3 >
    )
    target_compile_options(swe_softRes_hardRes_hashes PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>: -O3 >
    )
    target_compile_options(swe_softRes_admiss_useShared PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>: -O3 >
    )
    target_compile_options(swe_softRes_admiss_redundant PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>: -O3 >
    )
    target_compile_options(swe_checkpointRestart PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>: -O3 >
    )
endif ()

#*******************************************************************************
# Export compiler commands for use in YouCompleteMe (optional)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
