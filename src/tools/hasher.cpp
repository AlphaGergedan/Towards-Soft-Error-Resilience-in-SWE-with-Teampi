/**
 * @file src/tools/hasher.cpp
 *
 * @brief Implementation of src/tools/hasher.hpp
 *
 * TODO XOR operation is used when combining hashes but see boost::hash_combine
 *      as it might be better/faster
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
    Hasher::fieldSizeY = fieldSizeY;

    /* Update arrays */
    Hasher::hNetUpdatesLeft = currentBlock->hNetUpdatesLeft.getRawPointer();
    Hasher::hNetUpdatesRight = currentBlock->hNetUpdatesRight.getRawPointer();

    Hasher::huNetUpdatesLeft = currentBlock->huNetUpdatesLeft.getRawPointer();
    Hasher::huNetUpdatesRight = currentBlock->huNetUpdatesRight.getRawPointer();

    Hasher::hNetUpdatesBelow = currentBlock->hNetUpdatesBelow.getRawPointer();
    Hasher::hNetUpdatesAbove = currentBlock->hNetUpdatesAbove.getRawPointer();

    Hasher::hvNetUpdatesBelow = currentBlock->hvNetUpdatesBelow.getRawPointer();
    Hasher::hvNetUpdatesAbove = currentBlock->hvNetUpdatesAbove.getRawPointer();

    Hasher::maxTimeStep = &(currentBlock->maxTimestep);

    /* Data arrays */
    Hasher::h = currentBlock->getWaterHeight().getRawPointer();
    Hasher::hv = currentBlock->getMomentumVertical().getRawPointer();
    Hasher::hu = currentBlock->getMomentumHorizontal().getRawPointer();
    Hasher::b = currentBlock->getBathymetry().getRawPointer();

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
    total_hash ^= hash_fn(str_hLeft);
    total_hash ^= hash_fn(str_hRight);
    total_hash ^= hash_fn(str_huLeft);
    total_hash ^= hash_fn(str_huRight);
    total_hash ^= hash_fn(str_hBelow);
    total_hash ^= hash_fn(str_hAbove);
    total_hash ^= hash_fn(str_hvBelow);
    total_hash ^= hash_fn(str_hvAbove);
    total_hash ^= hash_fn(str_maxTimeStep);

    total_hash ^= hash_fn(str_h);
    total_hash ^= hash_fn(str_hv);
    total_hash ^= hash_fn(str_hu);
    total_hash ^= hash_fn(str_b);
}

/* Tries to improve hashing by not converting them to strings
 * However this method seems to be slower than string hashing..
 * Compiler is able to optimize the hasher::update_stdHash() better
 *
 * Warning: this method seems to be slower than update_stdHash(). Do not use this one */
void tools::Hasher::update_stdHash_float() {
    /* update the hash using float hasher  */
    int i = 0;
    while (i < fieldSizeY) {
        /* fieldSizeX length arrays */
        total_hash ^= hash_fn_float(hNetUpdatesLeft[i]);
        total_hash ^= hash_fn_float(hNetUpdatesRight[i]);
        total_hash ^= hash_fn_float(huNetUpdatesLeft[i]);
        total_hash ^= hash_fn_float(huNetUpdatesRight[i]);
        total_hash ^= hash_fn_float(h[i]);
        total_hash ^= hash_fn_float(hv[i]);
        total_hash ^= hash_fn_float(hu[i]);
        total_hash ^= hash_fn_float(b[i]);
        /* fieldSizeY length arrays */
        total_hash ^= hash_fn_float(hNetUpdatesBelow[i]);
        total_hash ^= hash_fn_float(hNetUpdatesAbove[i]);
        total_hash ^= hash_fn_float(hvNetUpdatesBelow[i]);
        total_hash ^= hash_fn_float(hvNetUpdatesAbove[i]);

        i++;
    }
    while (i < fieldSizeX) {
        /* fieldSizeX length arrays */
        total_hash ^= hash_fn_float(hNetUpdatesLeft[i]);
        total_hash ^= hash_fn_float(hNetUpdatesRight[i]);
        total_hash ^= hash_fn_float(huNetUpdatesLeft[i]);
        total_hash ^= hash_fn_float(huNetUpdatesRight[i]);
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
    str_hLeft = std::string((const char*) hNetUpdatesLeft, sizeof(float) * fieldSizeX);
    str_hRight = std::string((const char*) hNetUpdatesRight, sizeof(float) * fieldSizeX);
    str_huLeft = std::string((const char*) huNetUpdatesLeft, sizeof(float) * fieldSizeX);
    str_huRight = std::string((const char*) huNetUpdatesRight, sizeof(float) * fieldSizeX);

    str_hBelow = std::string((const char*) hNetUpdatesBelow, sizeof(float) * fieldSizeY);
    str_hAbove = std::string((const char*) hNetUpdatesAbove, sizeof(float) * fieldSizeY);
    str_hvBelow = std::string((const char*) hvNetUpdatesBelow, sizeof(float) * fieldSizeY);
    str_hvAbove = std::string((const char*) hvNetUpdatesAbove, sizeof(float) * fieldSizeY);

    str_maxTimeStep = std::string((const char*) maxTimeStep, sizeof(float));

    str_h = std::string((const char*) h, sizeof(float) * fieldSizeX);
    str_hv = std::string((const char*) hv, sizeof(float) * fieldSizeX);
    str_hu = std::string((const char*) hu, sizeof(float) * fieldSizeX);
    str_b = std::string((const char*) b, sizeof(float) * fieldSizeX);
}
