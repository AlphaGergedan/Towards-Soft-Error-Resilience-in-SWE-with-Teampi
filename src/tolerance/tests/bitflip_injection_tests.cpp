/**
 * @file tolerance/tests/bitflip_injection_tests.cpp
 * @brief Implementation of bitflip_injection_tests.hpp
 *
 * @author Atamert Rahma rahma@in.tum.de
 */

#include "tolerance/tests/bitflip_injection_tests.hpp"
#include <iostream>


/**
 * test function to run the bitflips
 *
 * Runs the simulation hard coded, and injects a desired bitflip at a random location.
 * Location can be either updates, data or both. It is given as parameter bitflipLocation.
 * Bitflip type can also be selected if desired, using the parameter bitflipType
 *
 * Information usually obtained over the command line can be provided with TestArguments args
 * parameter.
 *
 * @param bitflipLocation 0 for updates, 1 for data, 2 for random
 * @param bitflipType 0 for NaN, 1 Infinity, 2 -Infinity, 3 big number, 4 small number, 5 negative h, 6 change in b
 * @param args Simulation parameters which are normally obtained from the command line as user input
 * @param rankToCorrupt world rank of the processor to corrupt
 */
static std::function<void(FT_tests::TestArguments *args, int bitflipLocation, int bitflipType, int rankToCorrupt)> *run = nullptr;

void FT_tests::setRun(std::function<void(FT_tests::TestArguments *args,
                                    int bitflipLocation,
                                    int bitflipType,
                                    int rankToCorrupt)> *currentRun) {
    run = currentRun;
}

/* NaN injection */
void FT_tests::TEST_bitflipIntoUpdates_1(FT_tests::TestArguments *args, int rankToCorrupt) {
    (*run)(args, 0, 0, rankToCorrupt);
}

/* Inf injection */
void FT_tests::TEST_bitflipIntoUpdates_2(FT_tests::TestArguments *args, int rankToCorrupt) {
    (*run)(args, 0, 1, rankToCorrupt);
}

/* -Inf injection */
void FT_tests::TEST_bitflipIntoUpdates_3(FT_tests::TestArguments *args, int rankToCorrupt) {
    (*run)(args, 0, 2, rankToCorrupt);
}

/* Big number injection for DMP */
void FT_tests::TEST_bitflipIntoUpdates_4(FT_tests::TestArguments *args, int rankToCorrupt) {
    (*run)(args, 0, 3, rankToCorrupt);
}

/* Small number injection for DMP */
void FT_tests::TEST_bitflipIntoUpdates_5(FT_tests::TestArguments *args, int rankToCorrupt) {
    (*run)(args, 0, 4, rankToCorrupt);
}

//------------------------------------------------------------------------

/* NaN injection */
void FT_tests::TEST_bitflipIntoData_1(FT_tests::TestArguments *args, int rankToCorrupt) {
    (*run)(args, 1, 0, rankToCorrupt);
}

/* Inf injection */
void FT_tests::TEST_bitflipIntoData_2(FT_tests::TestArguments *args, int rankToCorrupt) {
    (*run)(args, 1, 1, rankToCorrupt);
}

/* -Inf injection */
void FT_tests::TEST_bitflipIntoData_3(FT_tests::TestArguments *args, int rankToCorrupt) {
    (*run)(args, 1, 2, rankToCorrupt);
}

/* Big number injection for DMP */
void FT_tests::TEST_bitflipIntoData_4(FT_tests::TestArguments *args, int rankToCorrupt) {
    (*run)(args, 1, 3, rankToCorrupt);
}

/* Small number injection for DMP */
void FT_tests::TEST_bitflipIntoData_5(FT_tests::TestArguments *args, int rankToCorrupt) {
    (*run)(args, 1, 4, rankToCorrupt);
}

/* Negative water height injection */
void FT_tests::TEST_bitflipIntoData_6(FT_tests::TestArguments *args, int rankToCorrupt) {
    (*run)(args, 1, 5, rankToCorrupt);
}

/* Bathymetry change injection */
void FT_tests::TEST_bitflipIntoData_7(FT_tests::TestArguments *args, int rankToCorrupt) {
    (*run)(args, 1, 6, rankToCorrupt);
}
