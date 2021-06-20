target_sources(swe_noRes PRIVATE
  Boundary.hpp
  Float2D.hpp
  Float2DNative.hpp
  Float2DBuffer.hpp)

target_sources(swe_hardRes_wTaskSharing PRIVATE
  Boundary.hpp
  Float2D.hpp
  Float2DNative.hpp
  Float2DBuffer.hpp)

target_sources(swe_softRes_hashes PRIVATE
  Boundary.hpp
  Float2D.hpp
  Float2DNative.hpp
  Float2DBuffer.hpp)

target_sources(swe_softRes_hardRes_hashes PRIVATE
  Boundary.hpp
  Float2D.hpp
  Float2DNative.hpp
  Float2DBuffer.hpp)

#target_sources(swe_softRes_admiss_useShared_v1 PRIVATE
#Boundary.hpp
#Float2D.hpp
#Float2DNative.hpp
#Float2DBuffer.hpp)

#target_sources(swe_softRes_admiss_useShared_v2 PRIVATE
#Boundary.hpp
#Float2D.hpp
#Float2DNative.hpp
#Float2DBuffer.hpp)

#target_sources(swe_softRes_admiss_useShared_v2_NONBLOCKING PRIVATE
#Boundary.hpp
#Float2D.hpp
#Float2DNative.hpp
#Float2DBuffer.hpp)

target_sources(swe_softRes_admiss_useShared PRIVATE
  Boundary.hpp
  Float2D.hpp
  Float2DNative.hpp
  Float2DBuffer.hpp)

target_sources(swe_softRes_admiss_redundant PRIVATE
  Boundary.hpp
  Float2D.hpp
  Float2DNative.hpp
  Float2DBuffer.hpp)

target_sources(swe_checkpointRestart PRIVATE
  Boundary.hpp
  Float2D.hpp
  Float2DNative.hpp
  Float2DBuffer.hpp)

#target_sources(swe_tests PRIVATE
#Boundary.hpp
#Float2D.hpp
#Float2DNative.hpp
#Float2DBuffer.hpp)
