
target_sources(swe_noRes PRIVATE
  help.hpp
  Args.hpp)

target_sources(swe_hardRes_wTaskSharing PRIVATE
  help.hpp
  Args.hpp
  ftLogger.hpp
  ftLogger.cpp)

target_sources(swe_softRes_hashes PRIVATE
  help.hpp
  Args.hpp
  ftLogger.hpp
  ftLogger.cpp
  hasher.cpp
  hasher.hpp)

target_sources(swe_softRes_hardRes_hashes PRIVATE
  help.hpp
  Args.hpp
  ftLogger.hpp
  ftLogger.cpp
  hasher.cpp
  hasher.hpp)

#target_sources(swe_softRes_admiss_useShared_v1 PRIVATE
#help.hpp
#Args.hpp
#ftLogger.hpp
#ftLogger.cpp)

#target_sources(swe_softRes_admiss_useShared_v2 PRIVATE
#help.hpp
#Args.hpp
#ftLogger.hpp
#ftLogger.cpp)

#target_sources(swe_softRes_admiss_useShared_v2_NONBLOCKING PRIVATE
#help.hpp
#Args.hpp
#ftLogger.hpp
#ftLogger.cpp)

target_sources(swe_softRes_admiss_useShared PRIVATE
  help.hpp
  Args.hpp
  ftLogger.hpp
  ftLogger.cpp
  Reports.hpp
  Reports.cpp)

target_sources(swe_softRes_admiss_redundant PRIVATE
  help.hpp
  Args.hpp
  ftLogger.hpp
  ftLogger.cpp
  Reports.hpp
  Reports.cpp)

target_sources(swe_checkpointRestart PRIVATE
  help.hpp
  Args.hpp
  ftLogger.hpp
  ftLogger.cpp)

#target_sources(swe_tests PRIVATE
#help.hpp
#Args.hpp
#ftLogger.hpp
#ftLogger.cpp
#hasher.cpp
#hasher.hpp)
