/**
 * @file src/tools/hasher.hpp
 *
 * @brief Small header for hashing calculated variables to be used in tolerance
 *
 * Generates and finalizes hashes for swe blocks. The pointers should not be
 * empty and the hasher should be calles after the compute net updates.
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
    void update_stdHash_float();
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
    std::hash<float> hash_fn_float;
    std::size_t total_hash;

    /* sha1 hash */
    SHA1 checksum;

    /* Converts the datas to strings for easy hashing.. should be avoided for
     * performance reasons TODO */
    void updateStrings();
};


#endif // HASHER_HPP
