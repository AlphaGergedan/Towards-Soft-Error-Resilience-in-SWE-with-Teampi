if(SOLVER STREQUAL "hlle")
  set(SOLVER_HEADER "HLLEFun.hpp")
  target_compile_definitions(swe_hardRes_wTaskSharing PRIVATE SOLVER=0)
  target_compile_definitions(swe_softRes_and_hardRes_woutTaskSharing PRIVATE SOLVER=0)
  target_compile_definitions(swe_checkpointRestart PRIVATE SOLVER=0)
else()
message(FATAL_ERROR "SOLVER must be one of: ${SOLVERS}")
endif()

message(STATUS "Using Riemann solver: ${SOLVER}")

target_sources(swe_hardRes_wTaskSharing PRIVATE ${SOLVER_HEADER})
target_sources(swe_softRes_and_hardRes_woutTaskSharing PRIVATE ${SOLVER_HEADER})
target_sources(swe_checkpointRestart PRIVATE ${SOLVER_HEADER})
