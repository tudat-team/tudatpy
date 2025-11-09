add_executable(ComaAnalysis_debug_solar_longitude main.cpp)

target_link_libraries(ComaAnalysis_debug_solar_longitude
    tudat_aerodynamics
    tudat_basic_astrodynamics
    tudat_basic_mathematics
    tudat_interpolators
    tudat_simulation_environment_setup
    ${Boost_LIBRARIES}
)

set_target_properties(ComaAnalysis_debug_solar_longitude PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
