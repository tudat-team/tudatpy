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
add_subdirectory("ephemerides")
add_subdirectory("gravitation")
add_subdirectory("ground_stations")
add_subdirectory("low_thrust")
add_subdirectory("mission_segments")
if (TUDAT_BUILD_WITH_ESTIMATION_TOOLS)
    add_subdirectory("observation_models")
    add_subdirectory("orbit_determination")
endif ()
add_subdirectory("propagators")
add_subdirectory("propulsion")
add_subdirectory("reference_frames")
add_subdirectory("relativity")
#add_subdirectory("system_models")
