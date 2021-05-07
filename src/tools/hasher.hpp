/**
 * @file src/tools/hasher.cpp
 *
 * @brief Small header for hashing calculated variables to be used in tolerance
 *
 * TODO improve hashing methods, converting to string might not be useful
 *
 */



#ifndef HASHER_HPP
#define HASHER_HPP

#include "tools/sha1.hpp"
#include <string.h>
#include <assert.h>

/* Hash the net updates calculated from SWE applications */
class Hasher {

public:

    /* Constructor */
    Hasher(int fieldSizeX, int fieldSizeY,
           float* hNetUpdatesLeft, float* hNetUpdatesRight,
           float* huNetUpdatesLeft, float* huNetUpdatesRight,
           float* hNetUpdatesBelow, float* hNetUpdatesAbove,
           float* hvNetUpdatesBelow, float* hvNetUpdatesAbove,
           float* maxTimeStep);

    /* hashes with std::hash */
    void update_stdHash();
    size_t finalize_stdHash();

    /* hashes with SHA-1 */
    void update_SHA1();
    unsigned char* finalize_SHA1();


private:

    /* Size of the update fields incl. ghost layer */
    int fieldSizeX, fieldSizeY;

    /* h udpates to hash, with size fieldSizeX */
    float *hNetUpdatesLeft, *hNetUpdatesRight;

    /* hu updates to hash, with size fieldSizeX */
    float *huNetUpdatesLeft, *huNetUpdatesRight;

    /* h updates to hash, with size fieldSizeY */
    float *hNetUpdatesBelow, *hNetUpdatesAbove;

    /* hv updates to hash, with size fieldSizeY */
    float *hvNetUpdatesBelow, *hvNetUpdatesAbove;

    /* max time step to hash, with size 1 */
    float* maxTimeStep;

    /* Strings values to hash
     * TODO Conversion is expensive.. strings are also expensive.
     *      change to std::hash<float> or something if possbible
     */
    std::string str_hLeft, str_hRight, str_huLeft, str_huRight, str_hBelow,
        str_hAbove, str_hvBelow, str_hvAbove, str_maxTimeStep;

    /* hash function and hash storage
     *
     * TODO calculate this hash using a float hash function,
     *      this way we must convert them to strings first !
     */
    std::hash<std::string> hash_fn;
    std::size_t total_hash;

    /* sha1 hash */
    SHA1 checksum;

    /* Converts the datas to strings for easy hashing.. should be avoided for
     * performance reasons TODO */
    void updateStrings();
};


/**
 * Constructor: save the pointers to the updates to hash
 *              them later.
 *
 * @param fieldSizeX TODO
 * @param fieldSizeY
 * @param hNetUpdatesLeft
 * @param hNetUpdatesRight
 * @param huNetUpdatesLeft
 * @param huNetUpdatesRight
 * @param hNetUpdatesBelow
 * @param hNetUpdatesAbove
 * @param hvNetUpdatesBelow
 * @param hvNetUpdatesAbove
 * @param maxTimeStep
 */
inline Hasher::Hasher(int fieldSizeX, int fieldSizeY,
                      float* hNetUpdatesLeft, float* hNetUpdatesRight,
                      float* huNetUpdatesLeft, float* huNetUpdatesRight,
                      float* hNetUpdatesBelow, float* hNetUpdatesAbove,
                      float* hvNetUpdatesBelow, float* hvNetUpdatesAbove,
                      float* maxTimeStep) {

    Hasher::fieldSizeX = fieldSizeX;
    Hasher::fieldSizeY = fieldSizeY;

    Hasher::hNetUpdatesLeft = hNetUpdatesLeft;
    Hasher::hNetUpdatesRight = hNetUpdatesRight;

    Hasher::huNetUpdatesLeft = huNetUpdatesLeft;
    Hasher::huNetUpdatesRight = huNetUpdatesRight;

    Hasher::hNetUpdatesBelow = hNetUpdatesBelow;
    Hasher::hNetUpdatesAbove = hNetUpdatesAbove;

    Hasher::hvNetUpdatesBelow = hvNetUpdatesBelow;
    Hasher::hvNetUpdatesAbove = hvNetUpdatesAbove;

    Hasher::maxTimeStep = maxTimeStep;

    /* initialized as 0 because we want the first xor operation with the
     * first hash to give the hash itself
     */
    Hasher::total_hash = 0;
}

/**
 * TODO
 *
 * @return r
 */
inline size_t Hasher::finalize_stdHash() {
    size_t r = total_hash;
    total_hash = 0;
    return r;
}

/**
 * TODO
 * https://stackoverflow.com/questions/3221170/how-to-turn-a-hex-string-into-an-unsigned-char-array
 *
 * @return
 */
inline unsigned char* Hasher::finalize_SHA1() {

    /* returns the final SHA1 hash as string */
    std::string final_sha1 = checksum.final();
    assert(final_sha1.size() == 40);
    const char *s = (const char*) final_sha1.data();

    /* we will return this */
    unsigned char resultOfsha1[20];

    for (int i = 0; i < 4; i++) {
        // TODO
    }
    return resultOfsha1;
}

/* hash with std::hash */
inline void Hasher::update_stdHash() {

    /* update the strings */
    updateStrings();

    /* TODO is this better than concat the strings then xor?*/
    total_hash ^= hash_fn(Hasher::str_hLeft);
    total_hash ^= hash_fn(Hasher::str_hRight);
    total_hash ^= hash_fn(Hasher::str_huLeft);
    total_hash ^= hash_fn(Hasher::str_huRight);
    total_hash ^= hash_fn(str_hBelow);
    total_hash ^= hash_fn(str_hAbove);
    total_hash ^= hash_fn(str_hvBelow);
    total_hash ^= hash_fn(str_hvAbove);
    total_hash ^= hash_fn(str_maxTimeStep);
}


/**
 * hash with SHA-1
 * TODO
 */
inline void Hasher::update_SHA1() {

    /* update the strings */
    updateStrings();

    /* concat the strings */

    /* hash the strings TODO don't hash them one by one ! */
    checksum.update(str_hLeft);
    checksum.update(str_hRight);
    checksum.update(str_huLeft);
    checksum.update(str_huRight);
    checksum.update(str_hBelow);
    checksum.update(str_hAbove);
    checksum.update(str_hvBelow);
    checksum.update(str_hvAbove);
    checksum.update(str_maxTimeStep);
}


/**
 *
 * TODO
 *
 *
 */
inline void Hasher::updateStrings() {
    str_hLeft = std::string((const char*) Hasher::hNetUpdatesLeft, sizeof(float) * fieldSizeX);
    Hasher::str_hRight = std::string((const char*) Hasher::hNetUpdatesRight, sizeof(float) * fieldSizeX);
    Hasher::str_huLeft = std::string((const char*) Hasher::huNetUpdatesLeft, sizeof(float) * fieldSizeX);
    Hasher::str_huRight = std::string((const char*) Hasher::huNetUpdatesRight, sizeof(float) * fieldSizeX);

    Hasher::str_hBelow = std::string((const char*) Hasher::hNetUpdatesBelow, sizeof(float) * fieldSizeY);
    Hasher::str_hAbove = std::string((const char*) Hasher::hNetUpdatesAbove, sizeof(float) * fieldSizeY);
    Hasher::str_hvBelow = std::string((const char*) Hasher::hvNetUpdatesBelow, sizeof(float) * fieldSizeY);
    Hasher::str_hvAbove = std::string((const char*) Hasher::hvNetUpdatesAbove, sizeof(float) * fieldSizeY);

    Hasher::str_maxTimeStep = std::string((const char*) Hasher::maxTimeStep, sizeof(float));
}


#endif // HASHER_HPP
