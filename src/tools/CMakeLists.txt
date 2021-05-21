target_sources(swe_hardRes_wTaskSharing PRIVATE
  help.hpp
  Args.hpp
  ftLogger.hpp
  ftLogger.cpp)

target_sources(swe_softRes PRIVATE
  help.hpp
  Args.hpp
  ftLogger.hpp
  ftLogger.cpp
  hasher.cpp
  hasher.hpp
  sha1.hpp)

target_sources(swe_softRes_and_hardRes_woutTaskSharing PRIVATE
  help.hpp
  Args.hpp
  ftLogger.hpp
  ftLogger.cpp
  hasher.cpp
  hasher.hpp
  sha1.hpp)

target_sources(swe_softRes_and_hardRes_wTaskSharing PRIVATE
  help.hpp
  Args.hpp
  ftLogger.hpp
  ftLogger.cpp)

target_sources(swe_checkpointRestart PRIVATE
  help.hpp
  Args.hpp
  ftLogger.hpp
  ftLogger.cpp)
