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
    "${SRCROOT}${JSONINTERFACEDIR}/Support/keys.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/modular.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/path.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/valueAccess.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/spice.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/body.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/atmosphere.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/ephemeris.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/aerodynamics.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/gravityField.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/rotationModel.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/shapeModel.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/radiationPressure.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/gravityFieldVariation.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/variable.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/termination.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/acceleration.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/thrust.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/torque.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Mathematics/interpolator.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/export.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/validation.cpp"
)

# Set the header files.
set(JSONINTERFACE_HEADERS
    "${SRCROOT}${JSONINTERFACEDIR}/Support/errorHandling.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/keys.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/modular.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/path.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/units.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/utilities.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/valueAccess.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/valueConversions.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/spice.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/body.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/atmosphere.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/ephemeris.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/aerodynamics.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/gravityField.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/rotationModel.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/shapeModel.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/radiationPressure.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/gravityFieldVariation.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/variable.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/propagator.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/termination.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/referenceFrames.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/acceleration.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/thrust.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/torque.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Mathematics/integrator.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Mathematics/interpolator.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/export.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/validation.h"
    "${SRCROOT}${JSONINTERFACEDIR}/simulation.h"
)

# Add static libraries.
add_library(tudat_json_interface STATIC ${JSONINTERFACE_SOURCES} ${JSONINTERFACE_HEADERS})
setup_tudat_library_target(tudat_json_interface "${SRCROOT}${JSONINTERFACEDIR}")

# Add unit tests.
# add_executable(test_SpiceInterface "${SRCROOT}${EXTERNALDIR}/SpiceInterface/UnitTests/unitTestSpiceInterface.cpp")
# setup_custom_test_program(test_SpiceInterface "${SRCROOT}${EXTERNALDIR}/SpiceInterface")
# target_link_libraries(test_SpiceInterface tudat_ephemerides tudat_basic_mathematics tudat_spice_interface tudat_basic_astrodynamics ${SPICE_LIBRARIES} ${Boost_LIBRARIES})
