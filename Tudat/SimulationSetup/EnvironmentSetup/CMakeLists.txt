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
file(GLOB_RECURSE SIMULATION_ENVIRONMENT_SETUP_HEADERS ${ENVIRONMENTSETUPDIR} ABSOLUTE ${CODEROOT} *.h)
file(GLOB_RECURSE SIMULATION_ENVIRONMENT_SETUP_SOURCES ${ENVIRONMENTSETUPDIR} ABSOLUTE ${CODEROOT} *.cpp)


# Add static libraries.
add_library(tudat_environment_setup STATIC ${SIMULATION_ENVIRONMENT_SETUP_SOURCES} ${SIMULATION_ENVIRONMENT_SETUP_HEADERS} )
setup_tudat_library_target(tudat_environment_setup "${SRCROOT}${ENVIRONMENTSETUPDIR}")


