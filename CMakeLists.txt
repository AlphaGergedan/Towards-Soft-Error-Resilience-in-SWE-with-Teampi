cmake_minimum_required(VERSION 3.0)

project(swe CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# There is only one solver, but having the option is nice
set(SOLVER "hlle" CACHE STRING "Riemann solver")
# This displays a combo box in CMake GUIs
set_property(CACHE SOLVER PROPERTY STRINGS ${SOLVERS})

add_executable(swe src/main.cpp)

target_sources(swe PRIVATE src/constants.hpp)

target_include_directories(swe PRIVATE src)

add_subdirectory(src/types)
add_subdirectory(src/tools)
add_subdirectory(src/io)
add_subdirectory(src/scenarios)
add_subdirectory(src/solvers)
add_subdirectory(src/blocks)

# NetCDF
include(FindNetCDF.cmake)
if(NOT NETCDF_FOUND)
  message(FATAL_ERROR "NetCDF could not be found!")
endif()
target_include_directories(swe SYSTEM PRIVATE NETCDF_INCLUDE_DIRS)
target_link_libraries(swe PRIVATE ${NETCDF_LIBRARIES})

# ULFM
# Set this to the directory, that contains the "bin" folder which contains ULFM's MPI executables
# set(MPI_HOME "/home/username/opt" CACHE PATH "Base directory for MPI" FORCE)
find_package(MPI REQUIRED)
target_link_libraries(swe PUBLIC MPI::MPI_CXX)

# teaMPI
# Set this to the directory containing libtmpi
target_include_directories(swe SYSTEM PRIVATE ../teaMPI/lib)
find_library(teaMPI tmpi HINTS ../teaMPI/lib)
if(NOT teaMPI)
  message(FATAL_ERROR "teaMPI could not be found!")
endif()
target_link_libraries(swe PRIVATE ${teaMPI})
