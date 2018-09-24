#    Copyright (c) 2010-2018, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Add source files.
file(GLOB_RECURSE SIMULATION_ESTIMATION_SETUP_HEADERS ${ESTIMATIONSETUPSDIR} ABSOLUTE ${CODEROOT} *.h)
file(GLOB_RECURSE SIMULATION_ESTIMATION_SETUP_SOURCES ${ESTIMATIONSETUPSDIR} ABSOLUTE ${CODEROOT} *.cpp)


# Add static libraries.
add_library(tudat_estimation_setup STATIC ${SIMULATION_ESTIMATION_SETUP_SOURCES} ${SIMULATION_ESTIMATION_SETUP_HEADERS} )
setup_tudat_library_target(tudat_estimation_setup "${SRCROOT}${ESTIMATIONSETUPSDIR}")

