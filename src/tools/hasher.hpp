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

        /* Constructor
         * int fieldSizeX: should be given as (nx+2)*(ny+2) including the ghost layers
         * int fieldSizeY: should be given as (nx+1)*(ny+2)
         * currentBlock: corresponding block
         */
        Hasher(int fieldSizeX, int fieldSizeY, SWE_DimensionalSplittingMPIOverdecomp *currentBlock);

        /* hashes with std::hash */
        void update_stdHash();
        void update_stdHash_float();
        size_t finalize_stdHash();

    private:

        /* Size of the fields of the given block */
        int fieldSizeX, fieldSizeY;

        /* max time step to hash, with size 1 */
        float *maxTimeStep;

        /* data arrays, we also need to hash them */
        float *h, *hv, *hu, *b;

        /* Strings values to hash, converting to strings may sound expensive but
         * float hasher seems not be optimized very well by compiler */
        std::string str_maxTimeStep, str_h, str_hu, str_hv, str_b;

        /* hash function and hash storage */
        std::hash<std::string> hash_fn;
        std::hash<float> hash_fn_float; // Seems to be slower !
        std::size_t total_hash;

        /* Converts the datas to strings for easy hashing */
        void updateStrings();

  }; // end of class FtLogger
} // end of namepace tools


#endif // HASHER_HPP
