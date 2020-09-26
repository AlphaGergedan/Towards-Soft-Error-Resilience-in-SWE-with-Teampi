if (NetCDF_FOUND AND ENABLE_NETCDF)
   set(SELECTED_WRITER_COMPILATION_UNIT SWE_LoadNetCdfScenario.cpp)
   set(SELECTED_WRITER_HEADER SWE_LoadNetCdfScenario.hh)
endif()

target_sources(swe 
    PUBLIC
    ${SELECTED_WRITER_HEADER}
    SWE_Scenario.hh
    SWE_simple_scenarios.hh
    
    PRIVATE
    ${SELECTED_WRITER_COMPILATION_UNIT}
)
