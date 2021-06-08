target_sources(swe_noRes PRIVATE
  Writer.cpp
  Writer.hpp
  NetCDFWriter.cpp
  NetCDFWriter.hpp
  Reader.cpp
  Reader.hpp)


target_sources(swe_hardRes_wTaskSharing PRIVATE
  Writer.cpp
  Writer.hpp
  NetCDFWriter.cpp
  NetCDFWriter.hpp
  Reader.cpp
  Reader.hpp)

target_sources(swe_softRes_hashes PRIVATE
  Writer.cpp
  Writer.hpp
  NetCDFWriter.cpp
  NetCDFWriter.hpp
  Reader.cpp
  Reader.hpp)

target_sources(swe_softRes_hardRes_hashes PRIVATE
  Writer.cpp
  Writer.hpp
  NetCDFWriter.cpp
  NetCDFWriter.hpp
  Reader.cpp
  Reader.hpp)

target_sources(swe_softRes_admiss_useShared_v1 PRIVATE
  Writer.cpp
  Writer.hpp
  NetCDFWriter.cpp
  NetCDFWriter.hpp
  Reader.cpp
  Reader.hpp)

# testing out new algorithm
target_sources(swe_softRes_admiss_useShared_v2 PRIVATE
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

target_sources(swe_tests PRIVATE
  Writer.cpp
  Writer.hpp
  NetCDFWriter.cpp
  NetCDFWriter.hpp
  Reader.cpp
  Reader.hpp)
