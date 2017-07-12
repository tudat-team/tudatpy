 #    Copyright (c) 2010-2017, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

#Set the source files.
set(JSONINTERFACE_SOURCES
    "${SRCROOT}${JSONINTERFACEDIR}/void.cpp"
)

# Set the header files.
set(JSONINTERFACE_HEADERS
    "${SRCROOT}${JSONINTERFACEDIR}/integratorSettings.h"
    "${SRCROOT}${JSONINTERFACEDIR}/simulation.h"
    "${SRCROOT}${JSONINTERFACEDIR}/utilities.h"
)

# Add static libraries.
add_library(tudat_json_interface STATIC ${JSONINTERFACE_SOURCES} ${JSONINTERFACE_HEADERS})
setup_tudat_library_target(tudat_json_interface "${SRCROOT}${JSONINTERFACEDIR}")

# Add unit tests.
# add_executable(test_SpiceInterface "${SRCROOT}${EXTERNALDIR}/SpiceInterface/UnitTests/unitTestSpiceInterface.cpp")
# setup_custom_test_program(test_SpiceInterface "${SRCROOT}${EXTERNALDIR}/SpiceInterface")
# target_link_libraries(test_SpiceInterface tudat_ephemerides tudat_basic_mathematics tudat_spice_interface tudat_basic_astrodynamics ${SPICE_LIBRARIES} ${Boost_LIBRARIES})
