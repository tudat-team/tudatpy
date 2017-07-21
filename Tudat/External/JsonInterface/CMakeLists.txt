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
    "${SRCROOT}${JSONINTERFACEDIR}/jsonInterface.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/keys.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/path.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Body/body.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Body/atmosphere.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Body/ephemeris.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Body/aerodynamics.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Body/gravityField.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Body/rotationModel.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Body/shapeModel.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Body/radiationPressure.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/integrator.cpp"
)

# Set the header files.
set(JSONINTERFACE_HEADERS
    "${SRCROOT}${JSONINTERFACEDIR}/jsonInterface.h"
    "${SRCROOT}${JSONINTERFACEDIR}/keys.h"
    "${SRCROOT}${JSONINTERFACEDIR}/path.h"
    "${SRCROOT}${JSONINTERFACEDIR}/units.h"
    "${SRCROOT}${JSONINTERFACEDIR}/utilities.h"
    "${SRCROOT}${JSONINTERFACEDIR}/simulation.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Body/body.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Body/atmosphere.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Body/ephemeris.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Body/aerodynamics.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Body/gravityField.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Body/rotationModel.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Body/shapeModel.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Body/radiationPressure.h"
    "${SRCROOT}${JSONINTERFACEDIR}/integrator.h"
)

# Add static libraries.
add_library(tudat_json_interface STATIC ${JSONINTERFACE_SOURCES} ${JSONINTERFACE_HEADERS})
setup_tudat_library_target(tudat_json_interface "${SRCROOT}${JSONINTERFACEDIR}")

# Add unit tests.
# add_executable(test_SpiceInterface "${SRCROOT}${EXTERNALDIR}/SpiceInterface/UnitTests/unitTestSpiceInterface.cpp")
# setup_custom_test_program(test_SpiceInterface "${SRCROOT}${EXTERNALDIR}/SpiceInterface")
# target_link_libraries(test_SpiceInterface tudat_ephemerides tudat_basic_mathematics tudat_spice_interface tudat_basic_astrodynamics ${SPICE_LIBRARIES} ${Boost_LIBRARIES})
