target_sources(swe_hardRes_wTaskSharing PRIVATE
  Scenario.hpp
  simple_scenarios.hpp
  LoadNetCDFScenario.cpp
  LoadNetCDFScenario.hpp)

target_sources(swe_softRes PRIVATE
  Scenario.hpp
  simple_scenarios.hpp
  LoadNetCDFScenario.cpp
  LoadNetCDFScenario.hpp)

target_sources(swe_softRes_and_hardRes_woutTaskSharing PRIVATE
  Scenario.hpp
  simple_scenarios.hpp
  LoadNetCDFScenario.cpp
  LoadNetCDFScenario.hpp)

target_sources(swe_checkpointRestart PRIVATE
  Scenario.hpp
  simple_scenarios.hpp
  LoadNetCDFScenario.cpp
  LoadNetCDFScenario.hpp)
