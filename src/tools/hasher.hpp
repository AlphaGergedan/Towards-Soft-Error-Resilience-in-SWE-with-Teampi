/**
 * @file src/tools/hasher.hpp
 *
 * @brief Small header for hashing calculated variables to be used in tolerance
 *
 * Generates and finalizes hashes for swe blocks. The pointers should not be
 * empty and the hasher should be called after the compute net-updates.
 *
 * @author Atamert Rahma rahma@in.tum.de
 */


#ifndef HASHER_HPP
#define HASHER_HPP

#include "blocks/DimSplitMPIOverdecomp.hpp"
#include <string.h>
#include <assert.h>

/* Hash the net updates calculated from SWE applications */
namespace tools {
  class Hasher {

    public:

        /* Constructor */
        Hasher(int fieldSizeX, int fieldSizeY, SWE_DimensionalSplittingMPIOverdecomp *currentBlock);

        /* hashes with std::hash */
        void update_stdHash();
        void update_stdHash_float();
        size_t finalize_stdHash();

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
        float *maxTimeStep;

        /* data arrays, we also need to hash them */
        float *h, *hv, *hu, *b;

        /* Strings values to hash, converting to strings may sound expensive but
        * float hasher seems not be optimized very well by compiler */
        std::string str_hLeft, str_hRight, str_huLeft, str_huRight, str_hBelow,
            str_hAbove, str_hvBelow, str_hvAbove, str_maxTimeStep,
            str_h, str_hu, str_hv, str_b;

        /* hash function and hash storage */
        std::hash<std::string> hash_fn;
        std::hash<float> hash_fn_float; // Seems to be slower !
        std::size_t total_hash;

        /* Converts the datas to strings for easy hashing */
        void updateStrings();

  }; // end of class FtLogger
} // end of namepace tools


#endif // HASHER_HPP
