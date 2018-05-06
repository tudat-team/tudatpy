 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

# Set the source files.
set(TORQUE_PARTIALS_SOURCES
  "${SRCROOT}${TORQUEPARTIALSDIR}/inertialTorquePartial.cpp"
  "${SRCROOT}${TORQUEPARTIALSDIR}/constantTorquePartial.cpp"
  "${SRCROOT}${TORQUEPARTIALSDIR}/secondDegreeGravitationalTorquePartial.cpp"
  "${SRCROOT}${TORQUEPARTIALSDIR}/sphericalHarmonicGravitationalTorquePartial.cpp"
)

# Set the header files.
set(TORQUE_PARTIALS_HEADERS
  "${SRCROOT}${TORQUEPARTIALSDIR}/torquePartial.h"
  "${SRCROOT}${TORQUEPARTIALSDIR}/inertiaTensorPartial.h"
  "${SRCROOT}${TORQUEPARTIALSDIR}/inertialTorquePartial.h"
  "${SRCROOT}${TORQUEPARTIALSDIR}/constantTorquePartial.h"
  "${SRCROOT}${TORQUEPARTIALSDIR}/secondDegreeGravitationalTorquePartial.h"
  "${SRCROOT}${TORQUEPARTIALSDIR}/sphericalHarmonicGravitationalTorquePartial.h"
)

# Add static libraries.
add_library(tudat_torque_partials STATIC ${TORQUE_PARTIALS_SOURCES} ${TORQUE_PARTIALS_HEADERS})
setup_tudat_library_target(tudat_torque_partials "${SRCROOT}{TORQUEPARTIALSDIR}")

if(USE_CSPICE)

    # Add unit tests
    add_executable(test_TorquePartials "${SRCROOT}${TORQUEPARTIALSDIR}/UnitTests/unitTestTorquePartials.cpp")
    setup_custom_test_program(test_TorquePartials "${SRCROOT}${TORQUEPARTIALSDIR}")
    target_link_libraries(test_TorquePartials ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

endif( )


