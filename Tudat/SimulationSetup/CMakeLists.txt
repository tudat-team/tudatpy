#    Copyright (c) 2010-2016, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

set(ESTIMATIONSETUPSDIR "${SIMULATIONSETUPDIR}/EstimationSetup")

# Add source files.
set(SIMULATION_SETUP_SOURCES
  "${SRCROOT}${SIMULATIONSETUPDIR}/createAccelerationModels.cpp"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createAccelerationPartials.cpp"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createStateDerivativePartials.cpp"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createEstimatableParameters.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createAtmosphereModel.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createEphemeris.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createFlightConditions.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createGravityFieldVariations.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createGravityField.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createRotationModel.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createBodies.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createBodyShapeModel.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createGroundStations.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createBodyDeformationModel.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createRadiationPressureInterface.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createThrustModelGuidance.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createMassRateModels.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/defaultBodies.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/body.cpp"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createOneWayRangePartials.cpp"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createAngularPositionPartials.cpp"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createLightTimeCorrectionPartials.cpp"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createPositionPartials.cpp"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createLightTimeCalculator.cpp"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createLightTimeCorrection.cpp"
)

# Add header files.
set(SIMULATION_SETUP_HEADERS 
  "${SRCROOT}${SIMULATIONSETUPDIR}/accelerationSettings.h"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createAccelerationPartials.h"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createStateDerivativePartials.h"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createEstimatableParameters.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createAccelerationModels.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createAtmosphereModel.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createEphemeris.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createFlightConditions.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createGravityFieldVariations.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createGravityField.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createRotationModel.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createBodies.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createBodyShapeModel.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createGroundStations.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createBodyDeformationModel.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createRadiationPressureInterface.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createThrustModelGuidance.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createMassRateModels.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/defaultBodies.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/body.h"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createOneWayRangePartials.h"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createAngularPositionPartials.h"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createPositionPartials.h"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createLightTimeCorrectionPartials.h"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createObservationPartials.h"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createObservationPartials.h"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createObservationModel.h"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createLightTimeCalculator.h"
  "${SRCROOT}${ESTIMATIONSETUPSDIR}/createLightTimeCorrection.h"
)

# Add static libraries.
add_library(tudat_simulation_setup STATIC ${SIMULATION_SETUP_SOURCES} ${SIMULATION_SETUP_HEADERS} )
setup_tudat_library_target(tudat_simulation_setup "${SRCROOT}${SIMULATIONSETUPDIR}")

# Add unit tests.
if(USE_CSPICE)
add_executable(test_EnvironmentSetup "${SRCROOT}${SIMULATIONSETUPDIR}/UnitTests/unitTestEnvironmentModelSetup.cpp")
setup_custom_test_program(test_EnvironmentSetup "${SRCROOT}${SIMULATIONSETUPDIR}/")
target_link_libraries(test_EnvironmentSetup ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_AccelerationModelSetup "${SRCROOT}${SIMULATIONSETUPDIR}/UnitTests/unitTestAccelerationModelSetup.cpp")
setup_custom_test_program(test_AccelerationModelSetup "${SRCROOT}${SIMULATIONSETUPDIR}/")
target_link_libraries(test_AccelerationModelSetup ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
endif()
