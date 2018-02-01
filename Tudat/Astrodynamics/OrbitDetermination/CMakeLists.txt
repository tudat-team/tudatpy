 #    Copyright (c) 2010-2018, Delft University of Technology
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
  "${SRCROOT}${ORBITDETERMINATIONDIR}/podProcessing.h"
  "${SRCROOT}${ORBITDETERMINATIONDIR}/determinePostFitParameterInfluence.h"
  "${SRCROOT}${ORBITDETERMINATIONDIR}/UnitTests/orbitDeterminationTestCases.h"
)


# Add static libraries.
add_library(tudat_orbit_determination STATIC ${ORBIT_DETERMINATION_SOURCES} ${ORBIT_DETERMINATION_HEADERS})
setup_tudat_library_target(tudat_orbit_determination "${SRCROOT}{ORBITDETERMINATIONDIR}")

if( COMPILE_PROPAGATION_TESTS )

add_executable(test_EstimationFromPositionDoubleDouble "${SRCROOT}${ORBITDETERMINATIONDIR}/UnitTests/unitTestEstimationFromIdealDataDoubleDouble.cpp")
setup_custom_test_program(test_EstimationFromPositionDoubleDouble "${SRCROOT}${ORBITDETERMINATIONDIR}")
target_link_libraries(test_EstimationFromPositionDoubleDouble ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_EstimationInput "${SRCROOT}${ORBITDETERMINATIONDIR}/UnitTests/unitTestEstimationInput.cpp")
setup_custom_test_program(test_EstimationInput "${SRCROOT}${ORBITDETERMINATIONDIR}")
target_link_libraries(test_EstimationInput ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})


add_executable(test_TidalPropertyEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}/UnitTests/unitTestTidalPropertyEstimation.cpp")
setup_custom_test_program(test_TidalPropertyEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}")
target_link_libraries(test_TidalPropertyEstimation ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

if( COMPILE_HIGH_ACCURACY_ESTIMATION_TESTS )

    add_executable(test_EstimationFromPositionDoubleLongDouble "${SRCROOT}${ORBITDETERMINATIONDIR}/UnitTests/unitTestEstimationFromIdealDataDoubleLongDouble.cpp")
    setup_custom_test_program(test_EstimationFromPositionDoubleLongDouble "${SRCROOT}${ORBITDETERMINATIONDIR}")
    target_link_libraries(test_EstimationFromPositionDoubleLongDouble ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

    add_executable(test_EstimationFromPositionTimeDouble "${SRCROOT}${ORBITDETERMINATIONDIR}/UnitTests/unitTestEstimationFromIdealDataTimeDouble.cpp")
    setup_custom_test_program(test_EstimationFromPositionTimeDouble "${SRCROOT}${ORBITDETERMINATIONDIR}")
    target_link_libraries(test_EstimationFromPositionTimeDouble ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

    add_executable(test_EstimationFromPositionTimeLongDouble "${SRCROOT}${ORBITDETERMINATIONDIR}/UnitTests/unitTestEstimationFromIdealDataTimeLongDouble.cpp")
    setup_custom_test_program(test_EstimationFromPositionTimeLongDouble "${SRCROOT}${ORBITDETERMINATIONDIR}")
    target_link_libraries(test_EstimationFromPositionTimeLongDouble ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

    add_executable(test_MultiArcStateEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}/UnitTests/unitTestMultiArcStateEstimation.cpp")
    setup_custom_test_program(test_MultiArcStateEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}")
    target_link_libraries(test_MultiArcStateEstimation ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

    add_executable(test_ParameterInfluenceDetermination "${SRCROOT}${ORBITDETERMINATIONDIR}/UnitTests/unitTestParameterInfluenceDetermination.cpp")
    setup_custom_test_program(test_ParameterInfluenceDetermination "${SRCROOT}${ORBITDETERMINATIONDIR}")
    target_link_libraries(test_ParameterInfluenceDetermination ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

    add_executable(test_BiasEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}/UnitTests/unitTestBiasEstimation.cpp")
    setup_custom_test_program(test_BiasEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}")
    target_link_libraries(test_BiasEstimation ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
endif( )

endif( )


