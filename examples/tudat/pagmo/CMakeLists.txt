

### Add Tudat pagmo examples
add_subdirectory(problems)

## Disable the CEC2013 example application on Windows as CEC2013 is not supported by MinGW
#if( NOT MINGW )
#    add_executable(application_PagmoCEC2013Comparison "${SRCROOT}/cec2013OptimizerComparison.cpp")
#    setup_executable_target(application_PagmoCEC2013Comparison "${SRCROOT}")
#    target_link_libraries(application_PagmoCEC2013Comparison pagmo pagmo2_library_example_problems tudat_input_output pthread ${Boost_LIBRARIES} )
#endif()


TUDAT_ADD_EXECUTABLE(application_PagmoHimmelblauOptimization
        "himmelblauOptimization.cpp"
        tudat_pagmo_example_problems tudat_input_output)

    TUDAT_ADD_EXECUTABLE(application_PagmoMgaTransferExample
            "mgaTransferExample.cpp"
            tudat_pagmo_example_problems tudat_mission_segments ${Tudat_PROPAGATION_LIBRARIES} )

#add_executable(application_PagmoHimmelblauOptimization "${SRCROOT}/himmelblauOptimization.cpp")
#setup_executable_target(application_PagmoHimmelblauOptimization "${SRCROOT}")
#target_link_libraries(application_PagmoHimmelblauOptimization pagmo pagmo2_library_example_problems tudat_input_output pthread ${Boost_LIBRARIES} )

#add_executable(application_PagmoEarthMarsTransferExample "${SRCROOT}/earthMarsTransferExample.cpp")
#setup_executable_target(application_PagmoEarthMarsTransferExample "${SRCROOT}")
#target_link_libraries(application_PagmoEarthMarsTransferExample pagmo pagmo2_library_example_problems ${Tudat_PROPAGATION_LIBRARIES}
#     pthread ${Boost_LIBRARIES} )

# add_executable(application_PagmoMgaTransferExample "${SRCROOT}/mgaTransferExample.cpp")
# setup_executable_target(application_PagmoMgaTransferExample "${SRCROOT}")
# target_link_libraries(application_PagmoMgaTransferExample pagmo pagmo2_library_example_problems ${Tudat_PROPAGATION_LIBRARIES}
#      pthread ${Boost_LIBRARIES} )


# add_executable(application_PagmoMultiObjectiveEarthMarsTransferExample "${SRCROOT}/multiObjectiveEarthMarsTransferExample.cpp")
# setup_executable_target(application_PagmoMultiObjectiveEarthMarsTransferExample "${SRCROOT}")
# target_link_libraries(application_PagmoMultiObjectiveEarthMarsTransferExample pagmo pagmo2_library_example_problems ${Tudat_PROPAGATION_LIBRARIES}
#      pthread ${Boost_LIBRARIES} )

#add_executable(application_PagmoPropagationTargetingExample "${SRCROOT}/propagationTargetingExample.cpp")
#setup_executable_target(application_PagmoPropagationTargetingExample "${SRCROOT}")
#target_link_libraries(application_PagmoPropagationTargetingExample pagmo pagmo2_library_example_problems ${Tudat_PROPAGATION_LIBRARIES}
#      pthread ${Boost_LIBRARIES} )

#if( BUILD_WITH_SOCIS_2019 )

# add_executable(application_PagmoHodographicShapingTrajectoryExample "${SRCROOT}/hodographicShapingTrajectoryExample.cpp")
# setup_executable_target(application_PagmoHodographicShapingTrajectoryExample "${SRCROOT}")
# target_link_libraries(application_PagmoHodographicShapingTrajectoryExample tudat_shape_based_methods tudat_low_thrust_trajectories
#     tudat_numerical_quadrature pagmo ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

# add_executable(application_PagmoHodographicShapingFullOptimisation "${SRCROOT}/hodographicShapingFullOptimisationExample.cpp")
# setup_executable_target(application_PagmoHodographicShapingFullOptimisation "${SRCROOT}")
# target_link_libraries(application_PagmoHodographicShapingFullOptimisation tudat_shape_based_methods tudat_low_thrust_trajectories
#     tudat_numerical_quadrature pagmo ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES} )

## # Add low-thrust trajectory design example
## add_executable(application_PagmoSimsFlanaganTrajectoryExample "${SRCROOT}/simsFlanaganTrajectoryExample.cpp")
## setup_executable_target(application_PagmoSimsFlanaganTrajectoryExample "${SRCROOT}")
## target_link_libraries(application_PagmoSimsFlanaganTrajectoryExample pagmo2_library_example_problems tudat_low_thrust_trajectories tudat_shape_based_methods
##     tudat_numerical_quadrature ${Tudat_PROPAGATION_LIBRARIES} pthread pagmo ${Boost_LIBRARIES} )

#endif( )
