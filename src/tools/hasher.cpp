/**
 * @file src/tools/hasher.cpp
 *
 * @brief Implementation of src/tools/hasher.hpp
 *
 * @author Atamert Rahma rahma@in.tum.de
 */


#include "tools/hasher.hpp"


/**
 * Constructor: save the pointers to the updates to hash
 *              them later.
 *
 * @param fieldSizeX
 * @param fieldSizeY
 * @param currentBlock Current dimensional splitting block that we want to hash
 */
tools::Hasher::Hasher(int fieldSizeX, int fieldSizeY, SWE_DimensionalSplittingMPIOverdecomp *currentBlock) {

    Hasher::fieldSizeX = fieldSizeX;
    //Hasher::fieldSizeY = fieldSizeY; // wedo not use this parameter

    /* Data arrays */
    Hasher::h = currentBlock->getWaterHeight().getRawPointer();
    Hasher::hv = currentBlock->getMomentumVertical().getRawPointer();
    Hasher::hu = currentBlock->getMomentumHorizontal().getRawPointer();
    Hasher::b = currentBlock->getBathymetry().getRawPointer();

    Hasher::maxTimeStep = &(currentBlock->maxTimestep);

    /* initialized as 0 because we want the first xor operation with the
     * first hash to give the hash itself
     */
    Hasher::total_hash = 0;
}


/**
 * Returns the final hash value calculated by the standard library
 *
 * @return r Final hash value
 */
size_t tools::Hasher::finalize_stdHash() {
    size_t r = total_hash;
    total_hash = 0;
    return r;
}

/* hash with std::hash */
void tools::Hasher::update_stdHash() {
    /* update the strings */
    updateStrings();

    /* update the hash */
    total_hash ^= hash_fn(str_h);
    total_hash ^= hash_fn(str_hv);
    total_hash ^= hash_fn(str_hu);
    total_hash ^= hash_fn(str_b);

    total_hash ^= hash_fn(str_maxTimeStep);
}

/* Tries to improve hashing by not converting them to strings
 * However this method seems to be slower than string hashing..
 * Compiler seems to optimize the hasher::update_stdHash() better
 *
 * Warning: this method seems to be slower than update_stdHash().
 */
void tools::Hasher::update_stdHash_float() {
    /* update the hash using float hasher  */
    int i = 0;
    while (i < fieldSizeX) {
        total_hash ^= hash_fn_float(h[i]);
        total_hash ^= hash_fn_float(hv[i]);
        total_hash ^= hash_fn_float(hu[i]);
        total_hash ^= hash_fn_float(b[i]);
        i++;
    }
    /* lastly hash the max time step */
    total_hash ^= hash_fn_float(*maxTimeStep);
}

/**
 * Converts the float arrays to strings to be hashed later
 * creating a single string seems to be slower than xor operations
 */
void tools::Hasher::updateStrings() {
    str_h = std::string((const char*) h, sizeof(float) * fieldSizeX);
    str_hv = std::string((const char*) hv, sizeof(float) * fieldSizeX);
    str_hu = std::string((const char*) hu, sizeof(float) * fieldSizeX);
    str_b = std::string((const char*) b, sizeof(float) * fieldSizeX);
    str_maxTimeStep = std::string((const char*) maxTimeStep, sizeof(float));
}
