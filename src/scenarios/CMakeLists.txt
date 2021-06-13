target_sources(swe_noRes PRIVATE
  Scenario.hpp
  simple_scenarios.hpp
  LoadNetCDFScenario.cpp
  LoadNetCDFScenario.hpp)

target_sources(swe_hardRes_wTaskSharing PRIVATE
  Scenario.hpp
  simple_scenarios.hpp
  LoadNetCDFScenario.cpp
  LoadNetCDFScenario.hpp)

target_sources(swe_softRes_hashes PRIVATE
  Scenario.hpp
  simple_scenarios.hpp
  LoadNetCDFScenario.cpp
  LoadNetCDFScenario.hpp)

target_sources(swe_softRes_hardRes_hashes PRIVATE
  Scenario.hpp
  simple_scenarios.hpp
  LoadNetCDFScenario.cpp
  LoadNetCDFScenario.hpp)

target_sources(swe_softRes_admiss_useShared_v1 PRIVATE
  Scenario.hpp
  simple_scenarios.hpp
  LoadNetCDFScenario.cpp
  LoadNetCDFScenario.hpp)

target_sources(swe_softRes_admiss_useShared_v2 PRIVATE
  Scenario.hpp
  simple_scenarios.hpp
  LoadNetCDFScenario.cpp
  LoadNetCDFScenario.hpp)

# TODO testing out new algorithm NONBLOCKING
target_sources(swe_softRes_admiss_useShared_v2_NONBLOCKING PRIVATE
  Scenario.hpp
  simple_scenarios.hpp
  LoadNetCDFScenario.cpp
  LoadNetCDFScenario.hpp)

target_sources(swe_softRes_admiss_redundant PRIVATE
  Scenario.hpp
  simple_scenarios.hpp
  LoadNetCDFScenario.cpp
  LoadNetCDFScenario.hpp)

target_sources(swe_checkpointRestart PRIVATE
  Scenario.hpp
  simple_scenarios.hpp
  LoadNetCDFScenario.cpp
  LoadNetCDFScenario.hpp)

target_sources(swe_tests PRIVATE
  Scenario.hpp
  simple_scenarios.hpp
  LoadNetCDFScenario.cpp
  LoadNetCDFScenario.hpp)
