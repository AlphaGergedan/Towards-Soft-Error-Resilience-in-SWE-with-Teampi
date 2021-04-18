if(SOLVER STREQUAL "hlle")
  set(SOLVER_HEADER "HLLEFun.hpp")
  target_compile_definitions(swe-fault-resilience PRIVATE SOLVER=0)
else()
message(FATAL_ERROR "SOLVER must be one of: ${SOLVERS}")
endif()

message(STATUS "Using Riemann solver: ${SOLVER}")

target_sources(swe-fault-resilience PRIVATE ${SOLVER_HEADER})
