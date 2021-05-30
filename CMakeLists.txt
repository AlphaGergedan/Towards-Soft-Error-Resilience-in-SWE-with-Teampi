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
# There is only one solver, but having the option is nice
set(SOLVER "hlle" CACHE STRING "Riemann solver")

# This displays a combo box in CMake GUIs
set_property(CACHE SOLVER PROPERTY STRINGS ${SOLVERS})

option(USE_DEBUG "Debug compile mode" OFF)

#*******************************************************************************
# sources

# Hard error resilience with task sharing to reduce redundancy
add_executable(swe_hardRes_wTaskSharing)                # --> old version

# Soft error resilience by comparing hashes between redundant replicas
add_executable(swe_softRes)                             # --> soft error detection

# Combination of hard and soft error resilience using warmspares + hashes
add_executable(swe_softRes_and_hardRes_woutTaskSharing) # --> soft + hard error resilience

# Baseline for hard error resilience
add_executable(swe_checkpointRestart)                   # --> baseline version for hard resilience

# Combination of hard and soft error resilience using warmspares + prechecks + task sharing
add_executable(swe_softRes_and_hardRes_wTaskSharing)    # --> soft + hard error resilience

# tests TODO improve this with unit tests (gtest might be good)
add_executable(swe_tests)


target_sources(swe_hardRes_wTaskSharing PRIVATE src/main.cpp)
target_sources(swe_hardRes_wTaskSharing PRIVATE src/constants.hpp)

target_sources(swe_softRes PRIVATE src/tolerance/swe_softRes.cpp)
target_sources(swe_softRes PRIVATE src/constants.hpp)

target_sources(swe_softRes_and_hardRes_woutTaskSharing PRIVATE src/tolerance/swe_softRes_and_hardRes_woutTaskSharing.cpp)
target_sources(swe_softRes_and_hardRes_woutTaskSharing PRIVATE src/constants.hpp)

target_sources(swe_softRes_and_hardRes_wTaskSharing PRIVATE src/tolerance/swe_softRes_and_hardRes_wTaskSharing.cpp)
target_sources(swe_softRes_and_hardRes_wTaskSharing PRIVATE src/constants.hpp)

target_sources(swe_checkpointRestart PRIVATE src/tolerance/swe_checkpointRestart.cpp)
target_sources(swe_checkpointRestart PRIVATE src/constants.hpp)

target_sources(swe_tests PRIVATE src/tolerance/swe_tests.cpp)
target_sources(swe_tests PRIVATE src/constants.hpp)


# Each subdirectory has it's own CMakeLists.txt to target their sources
add_subdirectory(src/types)
add_subdirectory(src/tools)
add_subdirectory(src/io)
add_subdirectory(src/scenarios)
add_subdirectory(src/solvers)
add_subdirectory(src/blocks)


# Include the src directory
target_include_directories(swe_hardRes_wTaskSharing PRIVATE src)

target_include_directories(swe_softRes PRIVATE src)

target_include_directories(swe_softRes_and_hardRes_woutTaskSharing PRIVATE src)

target_include_directories(swe_softRes_and_hardRes_wTaskSharing PRIVATE src)

target_include_directories(swe_checkpointRestart PRIVATE src)

target_include_directories(swe_tests PRIVATE src)


#*******************************************************************************
# libraries

# NetCDF
include(cmake/FindNetCDF.cmake)
if(NOT NETCDF_FOUND)
  message(FATAL_ERROR "NetCDF could not be found!")
endif()
target_include_directories(swe_hardRes_wTaskSharing SYSTEM PRIVATE NETCDF_INCLUDE_DIRS)
target_link_libraries(swe_hardRes_wTaskSharing PRIVATE ${NETCDF_LIBRARIES})

target_include_directories(swe_softRes SYSTEM PRIVATE NETCDF_INCLUDE_DIRS)
target_link_libraries(swe_softRes PRIVATE ${NETCDF_LIBRARIES})

target_include_directories(swe_softRes_and_hardRes_woutTaskSharing SYSTEM PRIVATE NETCDF_INCLUDE_DIRS)
target_link_libraries(swe_softRes_and_hardRes_woutTaskSharing PRIVATE ${NETCDF_LIBRARIES})

target_include_directories(swe_softRes_and_hardRes_wTaskSharing SYSTEM PRIVATE NETCDF_INCLUDE_DIRS)
target_link_libraries(swe_softRes_and_hardRes_wTaskSharing PRIVATE ${NETCDF_LIBRARIES})

target_include_directories(swe_checkpointRestart SYSTEM PRIVATE NETCDF_INCLUDE_DIRS)
target_link_libraries(swe_checkpointRestart PRIVATE ${NETCDF_LIBRARIES})

target_include_directories(swe_tests SYSTEM PRIVATE NETCDF_INCLUDE_DIRS)
target_link_libraries(swe_tests PRIVATE ${NETCDF_LIBRARIES})

# ULFM
# Set this to the directory, that contains the "bin" folder which contains ULFM's MPI executables
#set(MPI_HOME "/home/${USERNAME}/ssd/lib/ULFM/ulfm2/build/bin" CACHE PATH "Base directory for MPI" FORCE)
find_package(MPI REQUIRED)
if (MPI_FOUND)
    target_link_libraries(swe_hardRes_wTaskSharing PUBLIC MPI::MPI_CXX)
    target_link_libraries(swe_softRes PUBLIC MPI::MPI_CXX)
    target_link_libraries(swe_softRes_and_hardRes_woutTaskSharing PUBLIC MPI::MPI_CXX)
    target_link_libraries(swe_softRes_and_hardRes_wTaskSharing PUBLIC MPI::MPI_CXX)
    target_link_libraries(swe_checkpointRestart PUBLIC MPI::MPI_CXX) # checkpointing doesn't need ulfm
    target_link_libraries(swe_tests PUBLIC MPI::MPI_CXX) # checkpointing doesn't need ulfm
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

# Link tmpi, swe_checkpointRestart doesn't need it
target_link_libraries(swe_hardRes_wTaskSharing PRIVATE tmpi)
target_link_libraries(swe_softRes PRIVATE tmpi)
target_link_libraries(swe_softRes_and_hardRes_woutTaskSharing PRIVATE tmpi)
target_link_libraries(swe_softRes_and_hardRes_wTaskSharing PRIVATE tmpi)
target_link_libraries(swe_tests PRIVATE tmpi)

# include path
target_include_directories(swe_hardRes_wTaskSharing PRIVATE teampi-soft-error-resilience/include)
target_include_directories(swe_softRes PRIVATE teampi-soft-error-resilience/include)
target_include_directories(swe_softRes_and_hardRes_woutTaskSharing PRIVATE teampi-soft-error-resilience/include)
target_include_directories(swe_softRes_and_hardRes_wTaskSharing PRIVATE teampi-soft-error-resilience/include)
target_include_directories(swe_tests PRIVATE teampi-soft-error-resilience/include)


#*******************************************************************************
# add unit tests for fault tolerant swe
#add_subdirectory(tests)

#*******************************************************************************
# Warnings & Flags

if (USE_DEBUG)
    target_compile_options(tmpi PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>
    )
    target_compile_options(swe_hardRes_wTaskSharing PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>
    )
    target_compile_options(swe_softRes PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>
    )
    target_compile_options(swe_softRes_and_hardRes_woutTaskSharing PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>
    )
    target_compile_options(swe_softRes_and_hardRes_wTaskSharing PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>
    )
    target_compile_options(swe_checkpointRestart PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>
    )
else ()
    target_compile_options(tmpi PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>: -O3 >
    )
    target_compile_options(swe_hardRes_wTaskSharing PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>: -O3 >
    )
    target_compile_options(swe_softRes PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>: -O3 >
    )
    target_compile_options(swe_softRes_and_hardRes_woutTaskSharing PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>: -O3 >
    )
    target_compile_options(swe_softRes_and_hardRes_wTaskSharing PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>: -O3 >
    )
    target_compile_options(swe_checkpointRestart PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>: -O3 >
    )
endif ()

target_compile_options(swe_tests PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>: -O3 >
)


#*******************************************************************************
# Export compiler commands for use in YouCompleteMe
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
