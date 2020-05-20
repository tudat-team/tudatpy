#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#
#    Notes
#

# Add subdirectories.
add_subdirectory("aerodynamics")
add_subdirectory("basic_astro")
if (TUDAT_BUILD_WITH_SOFA_INTERFACE)
    add_subdirectory("earth_orientation")
endif ()
add_subdirectory("electromagnetism")
if (TUDAT_BUILD_WITH_SPICE_INTERFACE)
    add_subdirectory("ephemerides")
endif ()
add_subdirectory("gravitation")
add_subdirectory("ground_stations")
if (BUILD_WITH_SOCIS_2019)
    add_subdirectory("low_thrust")
endif ()
add_subdirectory("mission_segments")
if (TUDAT_BUILD_WITH_ESTIMATION_TOOLS)
    add_subdirectory("observation_models")
    add_subdirectory("orbit_determination")
endif ()
add_subdirectory("propagators")
add_subdirectory("propulsion")
add_subdirectory("reference_frames")
add_subdirectory("relativity")
add_subdirectory("system_models")
add_subdirectory("trajectory_design")

#add_subdirectory("${SRCROOT}${BASICASTRODYNAMICSDIR}")
#add_subdirectory("${SRCROOT}${ELECTROMAGNETISMDIR}")
#add_subdirectory("${SRCROOT}${EPHEMERIDESDIR}")
#add_subdirectory("${SRCROOT}${GRAVITATIONDIR}")
#add_subdirectory("${SRCROOT}${MISSIONSEGMENTSDIR}")
#add_subdirectory("${SRCROOT}${REFERENCEFRAMESDIR}")
#if(TUDAT_BUILD_WITH_ESTIMATION_TOOLS )
#    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/observation_models")
#    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/orbit_determination")
#endif( )
#add_subdirectory("${SRCROOT}${PROPAGATORSDIR}")
#add_subdirectory("${SRCROOT}${RELATIVITYDIR}")
#add_subdirectory("${SRCROOT}${SYSTEMMODELSDIR}")
#add_subdirectory("${SRCROOT}${TRAJECTORYDIR}")
#add_subdirectory("${SRCROOT}${PROPULSIONDIR}")
#if( BUILD_WITH_SOCIS_2019 )
#  add_subdirectory("${SRCROOT}${LOWTHRUSTTRAJECTORIESDIR}")
#endif()
#add_subdirectory("${SRCROOT}${GROUNDSTATIONSDIR}")
#if(USE_SOFA)
#  add_subdirectory("${SRCROOT}${EARTHORIENTATIONDIR}")
#endif()

# Get target properties for static libraries.
#get_target_property(AERODYNAMICSSOURCES tudat_aerodynamics SOURCES)
#get_target_property(BASICASTRODYNAMICSSOURCES tudat_basic_astrodynamics SOURCES)
#get_target_property(ELECTROMAGNETISMSOURCES tudat_electromagnetism SOURCES)
#get_target_property(EPHEMERIDESSOURCES tudat_ephemerides SOURCES)
#get_target_property(GRAVITATIONSOURCES tudat_gravitation SOURCES)
#get_target_property(MISSIONSEGMENTSSOURCES tudat_mission_segments SOURCES)
#get_target_property(REFERENCEFRAMESSOURCES tudat_reference_frames SOURCES)
#get_target_property(RELATIVITYSOURCES tudat_relativity SOURCES)
#get_target_property(PROPAGATORSDIR tudat_propagators SOURCES)
#get_target_property(SYSTEMMODELSSOURCES tudat_system_models SOURCES)
#get_target_property(TRAJECTORYSOURCES tudat_trajectory_design SOURCES)
#get_target_property(PROPULSIONSOURCES tudat_propulsion SOURCES)
#get_target_property(GROUNDSTATIONSSOURCES tudat_ground_stations SOURCES)
#get_target_property(LOWTHRUSTTRAJECTORIESDIRSOURCES tudat_low_thrust_trajectories SOURCES)
#
#if(TUDAT_BUILD_WITH_ESTIMATION_TOOLS )
#get_target_property(OBSERVATIONMODELSDIR tudat_observation_models SOURCES)
#get_target_property(ORBITDETERMINATIONSOURCES tudat_orbit_determination SOURCES)
#endif( )
#
#if(USE_SOFA)
#  get_target_property(EARTHORIENTATIONSOURCES tudat_earth_orientation SOURCES)
#endif()
