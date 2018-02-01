 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

#Set the source files.
set(SPICEINTERFACE_SOURCES
  "${SRCROOT}${EXTERNALDIR}/SpiceInterface/spiceEphemeris.cpp"
  "${SRCROOT}${EXTERNALDIR}/SpiceInterface/spiceRotationalEphemeris.cpp"
  "${SRCROOT}${EXTERNALDIR}/SpiceInterface/spiceInterface.cpp"
)

# Set the header files.
set(SPICEINTERFACE_HEADERS
  "${SRCROOT}${EXTERNALDIR}/SpiceInterface/spiceEphemeris.h"
  "${SRCROOT}${EXTERNALDIR}/SpiceInterface/spiceRotationalEphemeris.h"
  "${SRCROOT}${EXTERNALDIR}/SpiceInterface/spiceInterface.h"
)

# Add static libraries.
add_library(tudat_spice_interface STATIC ${SPICEINTERFACE_SOURCES} ${SPICEINTERFACE_HEADERS})
setup_tudat_library_target(tudat_spice_interface "${SRCROOT}${SPICEINTERFACEDIR}")

# Add unit tests.
add_executable(test_SpiceInterface "${SRCROOT}${EXTERNALDIR}/SpiceInterface/UnitTests/unitTestSpiceInterface.cpp")
setup_custom_test_program(test_SpiceInterface "${SRCROOT}${EXTERNALDIR}/SpiceInterface")
target_link_libraries(test_SpiceInterface tudat_ephemerides tudat_basic_mathematics tudat_spice_interface tudat_basic_astrodynamics ${SPICE_LIBRARIES} ${Boost_LIBRARIES})
