cmake_minimum_required(VERSION 3.13)

project(swe-fault-resilience
    DESCRIPTION "SWE with different fault tolerance techniques."
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

#*******************************************************************************
# sources

# More executables can be added for other failure mitigation methods
add_executable(swe-fault-resilience)

target_sources(swe-fault-resilience PRIVATE src/main.cpp)
target_sources(swe-fault-resilience PRIVATE src/constants.hpp)

# Each subdirectory has it's own CMakeLists.txt to target their sources
add_subdirectory(src/types)
add_subdirectory(src/tools)
add_subdirectory(src/io)
add_subdirectory(src/scenarios)
add_subdirectory(src/solvers)
add_subdirectory(src/blocks)

# Include the src directory
target_include_directories(swe-fault-resilience PRIVATE src)

#*******************************************************************************
# libraries

# NetCDF
include(cmake/FindNetCDF.cmake)
if(NOT NETCDF_FOUND)
  message(FATAL_ERROR "NetCDF could not be found!")
endif()
target_include_directories(swe-fault-resilience SYSTEM PRIVATE NETCDF_INCLUDE_DIRS)
target_link_libraries(swe-fault-resilience PRIVATE ${NETCDF_LIBRARIES})

# ULFM
# Set this to the directory, that contains the "bin" folder which contains ULFM's MPI executables
#set(MPI_HOME "/home/${USERNAME}/ssd/lib/ULFM/ulfm2/build/bin" CACHE PATH "Base directory for MPI" FORCE)
find_package(MPI REQUIRED)
if (MPI_FOUND)
    target_link_libraries(swe-fault-resilience PUBLIC MPI::MPI_CXX)
else()
    message(FATAL_ERROR "MPI not found!")
endif()

# teaMPI: TODO delete this section after you configured teaMPI submodule
# Set this to the directory containing libtmpi
#target_include_directories(swe-fault-resilience SYSTEM PRIVATE teampi-soft-error-resilience/lib)
#find_library(teaMPI tmpi HINTS teampi-soft-error-resilience/lib)
#if(NOT teaMPI)
#  message(FATAL_ERROR "teaMPI could not be found!")
#endif()
#target_link_libraries(swe-fault-resilience PRIVATE ${teaMPI})

#*******************************************************************************
# teaMPI library

#project(tmpi
#    DESCRIPTION "see teaMPI TODO")
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
#install(TARGETS tmpi
#    LIBRARY DESTINATION teampi-soft-error-resilience/lib
#    PUBLIC_HEADER DESTINATION teampi-soft-error-resilience/include
#    )

# Link tmpi
target_link_libraries(swe-fault-resilience PRIVATE tmpi)

# include path
target_include_directories(swe-fault-resilience PRIVATE teampi-soft-error-resilience/include)


#*******************************************************************************
# Warnings

#target_compile_options(tmpi PUBLIC
#$<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>)

target_compile_options(swe-fault-resilience PUBLIC
$<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter -g -pedantic>
)

#*******************************************************************************

# Export compiler commands for use in YouCompleteMe
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

