 #    Copyright (c) 2010-2018, Delft University of Technology
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
 list(APPEND TUDAT_EXTERNAL_LIBRARIES "")
 list(APPEND TUDAT_EXTERNAL_INTERFACE_LIBRARIES "")
 list(APPEND TUDAT_ITRS_LIBRARIES "")

 if(USE_SOFA)
  list(APPEND TUDAT_EXTERNAL_LIBRARIES sofa)
  list(APPEND TUDAT_EXTERNAL_INTERFACE_LIBRARIES tudat_sofa_interface )
  list(APPEND TUDAT_ITRS_LIBRARIES tudat_earth_orientation )
 endif()

 if(USE_CSPICE)
  list(APPEND TUDAT_EXTERNAL_LIBRARIES cspice)
  list(APPEND TUDAT_EXTERNAL_INTERFACE_LIBRARIES tudat_spice_interface )
 endif()

 if(USE_JSON)
   list(APPEND TUDAT_EXTERNAL_LIBRARIES nlohmann_json)
 endif()

 if(USE_NRLMSISE00)
  list(APPEND TUDAT_EXTERNAL_LIBRARIES nrlmsise00)
 endif()

 if(USE_GSL)
  list(APPEND TUDAT_EXTERNAL_LIBRARIES gsl)
 endif()

 # Find PaGMO library on local system.
 if( USE_PAGMO )
   list(APPEND TUDAT_EXTERNAL_LIBRARIES pthread)
 endif( )

 list(APPEND TUDAT_PROPAGATION_LIBRARIES tudat_trajectory_design tudat_propagation_setup tudat_environment_setup tudat_ground_stations tudat_propagators
     tudat_aerodynamics tudat_system_models tudat_geometric_shapes tudat_relativity tudat_gravitation tudat_mission_segments
     tudat_electro_magnetism tudat_propulsion tudat_ephemerides ${TUDAT_ITRS_LIBRARIES} tudat_numerical_integrators tudat_reference_frames
     tudat_statistics tudat_propagators ${TUDAT_EXTERNAL_INTERFACE_LIBRARIES} tudat_basic_astrodynamics tudat_interpolators tudat_root_finders tudat_filters
     tudat_basic_mathematics tudat_input_output tudat_basics ${TUDAT_EXTERNAL_LIBRARIES})

if( BUILD_WITH_ESTIMATION_TOOLS )

 list(APPEND TUDAT_ESTIMATION_LIBRARIES tudat_trajectory_design tudat_estimation_setup tudat_propagation_setup tudat_environment_setup  tudat_observation_models tudat_ground_stations tudat_acceleration_partials
    tudat_torque_partials  tudat_observation_partials tudat_orbit_determination tudat_estimatable_parameters tudat_propagators
     tudat_aerodynamics tudat_system_models tudat_geometric_shapes tudat_relativity tudat_gravitation tudat_mission_segments
     tudat_electro_magnetism tudat_propulsion tudat_ephemerides ${TUDAT_ITRS_LIBRARIES} tudat_numerical_integrators tudat_reference_frames
     tudat_statistics tudat_propagators ${TUDAT_EXTERNAL_INTERFACE_LIBRARIES} tudat_basic_astrodynamics tudat_interpolators tudat_root_finders tudat_filters
     tudat_basic_mathematics tudat_input_output tudat_basics ${TUDAT_EXTERNAL_LIBRARIES})

else( )
    list(APPEND TUDAT_ESTIMATION_LIBRARIES ${TUDAT_PROPAGATION_LIBRARIES})
endif( )
