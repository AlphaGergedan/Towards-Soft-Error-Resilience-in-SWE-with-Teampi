/**
 * @file tolerance/tests/bitflip_injection_tests.hpp
 * @brief Bitflip injection cases for tests
 *
 * User should set a run first. It can then be used
 * to test various error injection cases including
 * injection into data or update arrays and injection
 * of different types of errors like NaNs or negative
 * water height.
 *
 * See tolerance/swe_tests.cpp to see which cases
 * are tested for which methods.
 *
 * @author Atamert Rahma rahma@in.tum.de
 */

#ifndef FT_TESTS_HPP
#define FT_TESTS_HPP

#include <cstdlib>
#include <functional>


namespace FT_tests {

  /* Data sturcture for easier syntax on arguments passed by tests */
  struct TestArguments {
      float simulationDuration;
      clock_t heartbeatInterval;
      int nxRequested;
      int nyRequested;
      unsigned int decompFactor;
      bool writeOutput;
      unsigned int numberOfHashes;
      double bitflip_at;          /* simulation time to inject the bitflip */

      /* Constructor */
      TestArguments(float simulationDuration_, clock_t heartbeatInterval_,
                    int nxRequested_, int nyRequested_, unsigned int decompFactor_,
                    bool writeOutput_, double bitflip_at_) {
          simulationDuration = simulationDuration_;
          heartbeatInterval = heartbeatInterval_;
          nxRequested = nxRequested_;
          nyRequested = nyRequested_;
          decompFactor = decompFactor_;
          writeOutput = writeOutput_;
          bitflip_at = bitflip_at_;
      }
  };

  /* sets the main run function */
  void setRun(std::function<void(FT_tests::TestArguments *args,
                                 int bitflipLocation,
                                 int bitflipType,
                                 int rankToCorrupt)> *run);

  /* ------------ TESTS ------------ */

  /* NaN injection */
  void TEST_bitflipIntoUpdates_1(TestArguments *args, int rankToCorrupt);

  /* Inf injection */
  void TEST_bitflipIntoUpdates_2(TestArguments *args, int rankToCorrupt);

  /* -Inf injection */
  void TEST_bitflipIntoUpdates_3(TestArguments *args, int rankToCorrupt);

  /* Big number injection for DMP */
  void TEST_bitflipIntoUpdates_4(TestArguments *args, int rankToCorrupt);

  /* Small number injection for DMP */
  void TEST_bitflipIntoUpdates_5(TestArguments *args, int rankToCorrupt);

  //------------------------------------------------------------------------

  /* NaN injection */
  void TEST_bitflipIntoData_1(TestArguments *args, int rankToCorrupt);

  /* Inf injection */
  void TEST_bitflipIntoData_2(TestArguments *args, int rankToCorrupt);

  /* -Inf injection */
  void TEST_bitflipIntoData_3(TestArguments *args, int rankToCorrupt);

  /* Big number injection for DMP */
  void TEST_bitflipIntoData_4(TestArguments *args, int rankToCorrupt);

  /* Small number injection for DMP */
  void TEST_bitflipIntoData_5(TestArguments *args, int rankToCorrupt);

  /* Negative water height injection */
  void TEST_bitflipIntoData_6(TestArguments *args, int rankToCorrupt);

  /* Bathymetry change injection */
  void TEST_bitflipIntoData_7(TestArguments *args, int rankToCorrupt);

}

#endif
