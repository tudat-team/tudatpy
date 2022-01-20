#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#
#    References
#


# Create lists of static libraries for ease of use
list(APPEND TUDAT_EXTERNAL_LIBRARIES ${CSpice_LIBRARIES})
list(APPEND TUDAT_EXTERNAL_INTERFACE_LIBRARIES Tudat::tudat_spice_interface)

list(APPEND TUDAT_ITRS_LIBRARIES "")

if (TUDAT_BUILD_WITH_SOFA_INTERFACE)
    list(APPEND TUDAT_EXTERNAL_LIBRARIES ${Sofa_LIBRARIES})
    list(APPEND TUDAT_EXTERNAL_INTERFACE_LIBRARIES Tudat::tudat_sofa_interface)
    list(APPEND TUDAT_ITRS_LIBRARIES Tudat::tudat_earth_orientation)
endif ()

if (TUDAT_BUILD_WITH_NRLMSISE00)
    list(APPEND TUDAT_EXTERNAL_INTERFACE_LIBRARIES ${NRLMSISE00_LIBRARIES})
endif ()


if (TUDAT_BUILD_WITH_JSON_INTERFACE)
    list(APPEND TUDAT_EXTERNAL_LIBRARIES ${nlohmann_json_LIBRARIES})
endif ()

# if (TUDAT_BUILD_WITH_NRLMSISE00)
#     list(APPEND TUDAT_EXTERNAL_LIBRARIES nrlmsise00)
# endif ()

# if (USE_GSL)
#     list(APPEND TUDAT_EXTERNAL_LIBRARIES gsl)
# endif ()

# # Find PaGMO library on local system.
# if( USE_PAGMO )
#   list(APPEND TUDAT_EXTERNAL_LIBRARIES pthread)
# endif( )




list(APPEND Tudat_PROPAGATION_LIBRARIES
        Tudat::tudat_propagation_setup
        Tudat::tudat_shape_based_methods
        Tudat::tudat_low_thrust_trajectories
        Tudat::tudat_environment_setup
        Tudat::tudat_ground_stations
        Tudat::tudat_aerodynamics
        Tudat::tudat_system_models
        Tudat::tudat_geometric_shapes
        Tudat::tudat_relativity
        Tudat::tudat_gravitation
        Tudat::tudat_mission_segments
        Tudat::tudat_electromagnetism
        Tudat::tudat_propulsion
        Tudat::tudat_ephemerides
        ${TUDAT_ITRS_LIBRARIES}
        Tudat::tudat_numerical_integrators
        Tudat::tudat_reference_frames
        Tudat::tudat_statistics
        Tudat::tudat_propagators
        ${TUDAT_EXTERNAL_INTERFACE_LIBRARIES}
        Tudat::tudat_basic_astrodynamics        
        Tudat::tudat_numerical_quadrature
        Tudat::tudat_interpolators
        Tudat::tudat_root_finders
        Tudat::tudat_basic_mathematics
        Tudat::tudat_input_output
        Tudat::tudat_basics
#        ${TUDAT_EXTERNAL_LIBRARIES}
        )

if (TUDAT_BUILD_WITH_ESTIMATION_TOOLS)

    list(APPEND Tudat_ESTIMATION_LIBRARIES
            Tudat::tudat_estimation_setup
            Tudat::tudat_propagation_setup
            Tudat::tudat_environment_setup
            Tudat::tudat_observation_models
            Tudat::tudat_ground_stations
            Tudat::tudat_acceleration_partials
            Tudat::tudat_torque_partials
            Tudat::tudat_observation_partials
            Tudat::tudat_orbit_determination
            Tudat::tudat_estimatable_parameters
            Tudat::tudat_aerodynamics
            Tudat::tudat_system_models
            Tudat::tudat_geometric_shapes
            Tudat::tudat_relativity
            Tudat::tudat_gravitation
            Tudat::tudat_mission_segments
            Tudat::tudat_electromagnetism
            Tudat::tudat_propulsion
            Tudat::tudat_ephemerides
            ${TUDAT_ITRS_LIBRARIES}
            Tudat::tudat_numerical_integrators
            Tudat::tudat_reference_frames
            Tudat::tudat_statistics
            Tudat::tudat_propagators
            ${TUDAT_EXTERNAL_INTERFACE_LIBRARIES}
            Tudat::tudat_basic_astrodynamics
            Tudat::tudat_interpolators
            Tudat::tudat_root_finders
            Tudat::tudat_basic_mathematics
            Tudat::tudat_input_output
            Tudat::tudat_basics
            )

else ()
    list(APPEND Tudat_ESTIMATION_LIBRARIES ${Tudat_PROPAGATION_LIBRARIES})
endif ()

add_library(Tudat_LIBRARIES INTERFACE)
target_link_libraries(Tudat_LIBRARIES INTERFACE ${Tudat_ESTIMATION_LIBRARIES})

#set(Tudat_LIBRARIES ${Tudat_ESTIMATION_LIBRARIES})
