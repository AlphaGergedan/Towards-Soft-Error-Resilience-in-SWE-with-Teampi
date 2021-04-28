target_sources(swe_hardRes_wTaskSharing PRIVATE
  Boundary.hpp
  Float2D.hpp
  Float2DNative.hpp
  Float2DBuffer.hpp)

target_sources(swe_softRes PRIVATE
  Boundary.hpp
  Float2D.hpp
  Float2DNative.hpp
  Float2DBuffer.hpp)

target_sources(swe_softRes_and_hardRes_woutTaskSharing PRIVATE
  Boundary.hpp
  Float2D.hpp
  Float2DNative.hpp
  Float2DBuffer.hpp)

target_sources(swe_checkpointRestart PRIVATE
  Boundary.hpp
  Float2D.hpp
  Float2DNative.hpp
  Float2DBuffer.hpp)
