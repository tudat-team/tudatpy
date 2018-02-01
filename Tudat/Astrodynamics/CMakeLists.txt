 #    Copyright (c) 2010-2018, Delft University of Technology
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
# Define the main sub-directories.
set(AERODYNAMICSDIR "${ASTRODYNAMICSDIR}/Aerodynamics")
set(BASICASTRODYNAMICSDIR "${ASTRODYNAMICSDIR}/BasicAstrodynamics")
set(ELECTROMAGNETISMDIR "${ASTRODYNAMICSDIR}/ElectroMagnetism")
set(EPHEMERIDESDIR "${ASTRODYNAMICSDIR}/Ephemerides")
set(GRAVITATIONDIR "${ASTRODYNAMICSDIR}/Gravitation")
set(GROUNDSTATIONSDIR "${ASTRODYNAMICSDIR}/GroundStations")
set(MISSIONSEGMENTSDIR "${ASTRODYNAMICSDIR}/MissionSegments")
set(REFERENCEFRAMESDIR "${ASTRODYNAMICSDIR}/ReferenceFrames")
set(OBSERVATIONMODELSDIR "${ASTRODYNAMICSDIR}/ObservationModels")
set(ORBITDETERMINATIONDIR "${ASTRODYNAMICSDIR}/OrbitDetermination")
set(PROPAGATORSDIR "${ASTRODYNAMICSDIR}/Propagators")
set(RELATIVITYDIR "${ASTRODYNAMICSDIR}/Relativity")
set(SYSTEMMODELSDIR "${ASTRODYNAMICSDIR}/SystemModels")
set(PROPULSIONDIR "${ASTRODYNAMICSDIR}/Propulsion")
if(USE_SOFA)
  set(EARTHORIENTATIONDIR "${ASTRODYNAMICSDIR}/EarthOrientation")
endif()

# Add subdirectories.
add_subdirectory("${SRCROOT}${AERODYNAMICSDIR}")
add_subdirectory("${SRCROOT}${BASICASTRODYNAMICSDIR}")
add_subdirectory("${SRCROOT}${ELECTROMAGNETISMDIR}")
add_subdirectory("${SRCROOT}${EPHEMERIDESDIR}")
add_subdirectory("${SRCROOT}${GRAVITATIONDIR}")
add_subdirectory("${SRCROOT}${MISSIONSEGMENTSDIR}")
add_subdirectory("${SRCROOT}${REFERENCEFRAMESDIR}")
add_subdirectory("${SRCROOT}${OBSERVATIONMODELSDIR}")
add_subdirectory("${SRCROOT}${ORBITDETERMINATIONDIR}")
add_subdirectory("${SRCROOT}${PROPAGATORSDIR}")
add_subdirectory("${SRCROOT}${RELATIVITYDIR}")
add_subdirectory("${SRCROOT}${SYSTEMMODELSDIR}")
add_subdirectory("${SRCROOT}${PROPULSIONDIR}")
add_subdirectory("${SRCROOT}${GROUNDSTATIONSDIR}")
if(USE_SOFA)
  add_subdirectory("${SRCROOT}${EARTHORIENTATIONDIR}")
endif()

# Get target properties for static libraries.
get_target_property(AERODYNAMICSSOURCES tudat_aerodynamics SOURCES)
get_target_property(BASICASTRODYNAMICSSOURCES tudat_basic_astrodynamics SOURCES)
get_target_property(ELECTROMAGNETISMSOURCES tudat_electro_magnetism SOURCES)
get_target_property(EPHEMERIDESSOURCES tudat_ephemerides SOURCES)
get_target_property(GRAVITATIONSOURCES tudat_gravitation SOURCES)
get_target_property(MISSIONSEGMENTSSOURCES tudat_mission_segments SOURCES)
get_target_property(REFERENCEFRAMESSOURCES tudat_reference_frames SOURCES)
get_target_property(OBSERVATIONMODELSDIR tudat_observation_models SOURCES)
get_target_property(ORBITDETERMINATIONSOURCES tudat_orbit_determination SOURCES)
get_target_property(RELATIVITYSOURCES tudat_relativity SOURCES)
get_target_property(PROPAGATORSDIR tudat_propagators SOURCES)
get_target_property(SYSTEMMODELSSOURCES tudat_system_models SOURCES)
get_target_property(PROPULSIONSOURCES tudat_propulsion SOURCES)
get_target_property(GROUNDSTATIONSSOURCES tudat_ground_stations SOURCES)
if(USE_SOFA)
  get_target_property(EARTHORIENTATIONSOURCES tudat_earth_orientation SOURCES)
endif()
