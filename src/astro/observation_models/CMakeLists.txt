#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

set(OBSERVABLECORRECTIONSDIR "${OBSERVATIONMODELSDIR}/ObservableCorrections")


# Set the source files.
set(observation_models_SOURCES
  "lightTimeSolution.cpp"
  "observableTypes.cpp"
  "linkTypeDefs.cpp"
  "observationViabilityCalculator.cpp"
  "observationManager.cpp"
  "observationSimulator.cpp"
  "observationModel.cpp"
  "simulateObservations.cpp"
  "corrections/firstOrderRelativisticCorrection.cpp"
)

# Set the header files.
set(observation_models_HEADERS
  "angularPositionObservationModel.h"
  "lightTimeSolution.h"
  "linkTypeDefs.h"
  "observableTypes.h"
  "observationModel.h"
  "observationManager.h"
  "observationSimulator.h"
  "observationBias.h"
  "nWayRangeObservationModel.h"
  "oneWayRangeObservationModel.h"
  "oneWayDopplerObservationModel.h"
  "oneWayDifferencedRangeRateObservationModel.h"
  "simulateObservations.h"
  "positionObservationModel.h"
  "velocityObservationModel.h"
  "twoWayDopplerObservationModel.h"
  "eulerAngleObservationModel.h"
  "tests/testLightTimeCorrections.h"
  "observationViabilityCalculator.h"
  "corrections/lightTimeCorrection.h"
  "corrections/firstOrderRelativisticCorrection.h"
)

TUDAT_ADD_LIBRARY("observation_models"
     "${observation_models_SOURCES}"
     "${observation_models_HEADERS}"
     )



 # Add static libraries.
#add_library(tudat_observation_models STATIC ${OBSERVATION_MODELS_SOURCES} ${OBSERVATION_MODELS_HEADERS})
#setup_tudat_library_target(tudat_observation_models "${SRCROOT}${OBSERVATIONMODELSDIR}")

#if(TUDAT_BUILD_WITH_SPICE_INTERFACE)
#
#    add_executable(test_LightTime "${SRCROOT}${OBSERVATIONMODELSDIR}/tests/unitTestLightTimeSolution.cpp")
#    setup_custom_test_program(test_LightTime "${SRCROOT}${OBSERVATIONMODELSDIR}")
#    target_link_libraries(test_LightTime tudat_observation_models tudat_ephemerides tudat_basic_astrodynamics tudat_spice_interface ${SPICE_LIBRARIES} ${Boost_LIBRARIES})
#
#    add_executable(test_AngularPositionModel "${SRCROOT}${OBSERVATIONMODELSDIR}/tests/unitTestAngularPositionModel.cpp")
#    setup_custom_test_program(test_AngularPositionModel "${SRCROOT}${OBSERVATIONMODELSDIR}")
#    target_link_libraries(test_AngularPositionModel ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#    add_executable(test_OneWayRangeModel "${SRCROOT}${OBSERVATIONMODELSDIR}/tests/unitTestOneWayRangeObservationModel.cpp")
#    setup_custom_test_program(test_OneWayRangeModel "${SRCROOT}${OBSERVATIONMODELSDIR}")
#    target_link_libraries(test_OneWayRangeModel ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#    add_executable(test_OneWayDifferencedRangeModel "${SRCROOT}${OBSERVATIONMODELSDIR}/tests/unitTestOneWayDifferencedRangeModel.cpp")
#    setup_custom_test_program(test_OneWayDifferencedRangeModel "${SRCROOT}${OBSERVATIONMODELSDIR}")
#    target_link_libraries(test_OneWayDifferencedRangeModel ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#    add_executable(test_DopplerModels "${SRCROOT}${OBSERVATIONMODELSDIR}/tests/unitTestDopplerModels.cpp")
#    setup_custom_test_program(test_DopplerModels "${SRCROOT}${OBSERVATIONMODELSDIR}")
#    target_link_libraries(test_DopplerModels ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#    add_executable(test_PositionObservationModel "${SRCROOT}${OBSERVATIONMODELSDIR}/tests/unitTestPositionObservationModel.cpp")
#    setup_custom_test_program(test_PositionObservationModel "${SRCROOT}${OBSERVATIONMODELSDIR}")
#    target_link_libraries(test_PositionObservationModel ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#    add_executable(test_EulerAngleObservationModel "${SRCROOT}${OBSERVATIONMODELSDIR}/tests/unitTestEulerAngleObservationModel.cpp")
#    setup_custom_test_program(test_EulerAngleObservationModel "${SRCROOT}${OBSERVATIONMODELSDIR}")
#    target_link_libraries(test_EulerAngleObservationModel ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#    add_executable(test_VelocityObservationModel "${SRCROOT}${OBSERVATIONMODELSDIR}/tests/unitTestVelocityObservationModel.cpp")
#    setup_custom_test_program(test_VelocityObservationModel "${SRCROOT}${OBSERVATIONMODELSDIR}")
#    target_link_libraries(test_VelocityObservationModel ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#    add_executable(test_ObservationNoiseSimulation "${SRCROOT}${OBSERVATIONMODELSDIR}/tests/unitTestSimulatedObservationNoise.cpp")
#    setup_custom_test_program(test_ObservationNoiseSimulation "${SRCROOT}${OBSERVATIONMODELSDIR}")
#    target_link_libraries(test_ObservationNoiseSimulation ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#    add_executable(test_NWayRangeObservationModel "${SRCROOT}${OBSERVATIONMODELSDIR}/tests/unitTestNWayRangeObservationModel.cpp")
#    setup_custom_test_program(test_NWayRangeObservationModel "${SRCROOT}${OBSERVATIONMODELSDIR}")
#    target_link_libraries(test_NWayRangeObservationModel ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#    add_executable(test_ObservationViabilityCalculators "${SRCROOT}${OBSERVATIONMODELSDIR}/tests/unitTestObservationViabilityCalculators.cpp")
#    setup_custom_test_program(test_ObservationViabilityCalculators "${SRCROOT}${OBSERVATIONMODELSDIR}")
#    target_link_libraries(test_ObservationViabilityCalculators ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#endif( )
