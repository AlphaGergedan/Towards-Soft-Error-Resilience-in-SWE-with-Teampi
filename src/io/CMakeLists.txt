target_sources(swe_hardRes_wTaskSharing PRIVATE
  Writer.cpp
  Writer.hpp
  NetCDFWriter.cpp
  NetCDFWriter.hpp
  Reader.cpp
  Reader.hpp)

target_sources(swe_softRes_and_hardRes_woutTaskSharing PRIVATE
  Writer.cpp
  Writer.hpp
  NetCDFWriter.cpp
  NetCDFWriter.hpp
  Reader.cpp
  Reader.hpp)

target_sources(swe_checkpointRestart PRIVATE
  Writer.cpp
  Writer.hpp
  NetCDFWriter.cpp
  NetCDFWriter.hpp
  Reader.cpp
  Reader.hpp)
