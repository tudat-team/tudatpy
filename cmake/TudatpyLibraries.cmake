#    Copyright (c) 2010-2024, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

# External libraries
list(
    APPEND TUDATPY_EXTERNAL_LIBRARIES
    ${Boost_LIBRARIES}
    ${CSpice_LIBRARIES}
    ${Sofa_LIBRARIES}
    ${NRLMSISE00_LIBRARIES}
    ${nlohmann_json_LIBRARIES}
)


# Internal libraries
list(
    APPEND TUDATPY_INTERNAL_LIBRARIES
    tudat_propagation_setup
    tudat_shape_based_methods
    tudat_low_thrust_trajectories
    tudat_environment_setup
    tudat_ground_stations
    tudat_aerodynamics
    tudat_system_models
    tudat_geometric_shapes
    tudat_relativity
    tudat_gravitation
    tudat_mission_segments
    tudat_electromagnetism
    tudat_propulsion
    tudat_ephemerides
    tudat_numerical_integrators
    tudat_reference_frames
    tudat_statistics
    tudat_propagators
    tudat_basic_astrodynamics
    tudat_numerical_quadrature
    tudat_interpolators
    tudat_root_finders
    tudat_basic_mathematics
    tudat_input_output
    tudat_basics
    tudat_data
    tudat_estimation_setup
    tudat_observation_models
    tudat_acceleration_partials
    tudat_torque_partials
    tudat_observation_partials
    tudat_orbit_determination
    tudat_estimatable_parameters
    tudat_earth_orientation
    # Interfaces
    tudat_spice_interface
)
