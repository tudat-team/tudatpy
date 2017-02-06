 #    Copyright (c) 2010-2017, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

set(ACCELERATIONPARTIALSDIR "${ORBITDETERMINATIONDIR}/AccelerationPartials")
set(ESTIMATABLEPARAMETERSDIR "${ORBITDETERMINATIONDIR}/EstimatableParameters")
set(OBSERVATIONPARTIALSDIR "${ORBITDETERMINATIONDIR}/ObservationPartials")
set(LIGHTTIMECORRECTIONPARTIALSDIR "${ORBITDETERMINATIONDIR}/LightTimeCorrectionPartials")

# Add subdirectories.
add_subdirectory("${SRCROOT}${ACCELERATIONPARTIALSDIR}")
add_subdirectory("${SRCROOT}${ESTIMATABLEPARAMETERSDIR}")
add_subdirectory("${SRCROOT}${OBSERVATIONPARTIALSDIR}")

# Get target properties for static libraries.
get_target_property(ACCELERATIONPARTIALSSOURCES tudat_acceleration_partials SOURCES)
get_target_property(ESTIMATABLEPARAMETERSOURCES tudat_estimatable_parameters SOURCES)
get_target_property(OBSERVATIONPARTIALSSOURCES tudat_observation_partials SOURCES)

# Set the source files.
set(ORBIT_DETERMINATION_SOURCES
  "${SRCROOT}${ORBITDETERMINATIONDIR}/stateDerivativePartial.cpp"
)

# Set the header files.
set(ORBIT_DETERMINATION_HEADERS
  "${SRCROOT}${ORBITDETERMINATIONDIR}/stateDerivativePartial.h"
  "${SRCROOT}${ORBITDETERMINATIONDIR}/orbitDeterminationManager.h"
  "${SRCROOT}${ORBITDETERMINATIONDIR}/podInputOutputTypes.h"
)


# Add static libraries.
add_library(tudat_orbit_determination STATIC ${ORBIT_DETERMINATION_SOURCES} ${ORBIT_DETERMINATION_HEADERS})
setup_tudat_library_target(tudat_orbit_determination "${SRCROOT}{ORBITDETERMINATIONDIR}")

add_executable(test_EstimationFromPosition "${SRCROOT}${ORBITDETERMINATIONDIR}/UnitTests/unitTestEstimationFromIdealData.cpp")
setup_custom_test_program(test_EstimationFromPosition "${SRCROOT}${ORBITDETERMINATIONDIR}")
target_link_libraries(test_EstimationFromPosition ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})


