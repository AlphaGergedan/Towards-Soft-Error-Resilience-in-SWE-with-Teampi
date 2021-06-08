target_sources(swe_noRes PRIVATE
  Block.hpp
  DimSplitMPIOverdecomp.cpp
  DimSplitMPIOverdecomp.hpp)

target_sources(swe_hardRes_wTaskSharing PRIVATE
  Block.hpp
  DimSplitMPIOverdecomp.cpp
  DimSplitMPIOverdecomp.hpp)

target_sources(swe_softRes_hashes PRIVATE
  Block.hpp
  DimSplitMPIOverdecomp.cpp
  DimSplitMPIOverdecomp.hpp)

target_sources(swe_softRes_hardRes_hashes PRIVATE
  Block.hpp
  DimSplitMPIOverdecomp.cpp
  DimSplitMPIOverdecomp.hpp)

target_sources(swe_softRes_admiss_useShared_v1 PRIVATE
  Block.hpp
  DimSplitMPIOverdecomp.cpp
  DimSplitMPIOverdecomp.hpp)

# TODO testing out new algorithm
target_sources(swe_softRes_admiss_useShared_v2 PRIVATE
  Block.hpp
  DimSplitMPIOverdecomp.cpp
  DimSplitMPIOverdecomp.hpp)

target_sources(swe_checkpointRestart PRIVATE
  Block.hpp
  DimSplitMPIOverdecomp.cpp
  DimSplitMPIOverdecomp.hpp)

target_sources(swe_tests PRIVATE
  Block.hpp
  DimSplitMPIOverdecomp.cpp
  DimSplitMPIOverdecomp.hpp)
