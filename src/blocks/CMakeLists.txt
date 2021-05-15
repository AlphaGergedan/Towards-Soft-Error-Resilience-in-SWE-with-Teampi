target_sources(swe_hardRes_wTaskSharing PRIVATE
  Block.hpp
  DimSplitMPIOverdecomp.cpp
  DimSplitMPIOverdecomp.hpp)

target_sources(swe_softRes PRIVATE
  Block.hpp
  DimSplitMPIOverdecomp.cpp
  DimSplitMPIOverdecomp.hpp)

target_sources(swe_softRes_and_hardRes_woutTaskSharing PRIVATE
  Block.hpp
  DimSplitMPIOverdecomp.cpp
  DimSplitMPIOverdecomp.hpp)

target_sources(swe_softRes_and_hardRes_wTaskSharing PRIVATE
  Block.hpp
  DimSplitMPIOverdecomp.cpp
  DimSplitMPIOverdecomp.hpp)

target_sources(swe_checkpointRestart PRIVATE
  Block.hpp
  DimSplitMPIOverdecomp.cpp
  DimSplitMPIOverdecomp.hpp)
