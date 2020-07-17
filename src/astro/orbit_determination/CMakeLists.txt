 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

 # TODO: Add orbit determination library back in

 add_subdirectory("acceleration_partials")
 add_subdirectory("rotational_dynamics_partials")
 add_subdirectory("estimatable_parameters")
 add_subdirectory("observation_partials")

## Set the source files.
set(orbit_determination_SOURCES
  "stateDerivativePartial.cpp"
  "podInputOutputTypes.cpp"
)

# Set the header files.
set(orbit_determination_HEADERS
  "stateDerivativePartial.h"
  "podInputOutputTypes.h"
)
#
#
# Add library.
TUDAT_ADD_LIBRARY("orbit_determination"
        "${orbit_determination_SOURCES}"
        "${orbit_determination_HEADERS}"
        PUBLIC_LINKS
#        tudat_input_output
#        tudat_geometric_shapes
#        tudat_reference_frames
#        tudat_interpolators
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )

#if( BUILD_WITH_PROPAGATION_TESTS )
#
#add_executable(test_EstimationFromPositionDoubleDouble "${SRCROOT}${ORBITDETERMINATIONDIR}/tests/unitTestEstimationFromIdealDataDoubleDouble.cpp")
#setup_custom_test_program(test_EstimationFromPositionDoubleDouble "${SRCROOT}${ORBITDETERMINATIONDIR}")
#target_link_libraries(test_EstimationFromPositionDoubleDouble ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_EstimationInput "${SRCROOT}${ORBITDETERMINATIONDIR}/tests/unitTestEstimationInput.cpp")
#setup_custom_test_program(test_EstimationInput "${SRCROOT}${ORBITDETERMINATIONDIR}")
#target_link_libraries(test_EstimationInput ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#
#add_executable(test_TidalPropertyEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}/tests/unitTestTidalPropertyEstimation.cpp")
#setup_custom_test_program(test_TidalPropertyEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}")
#target_link_libraries(test_TidalPropertyEstimation ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_MultiArcStateEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}/tests/unitTestMultiArcStateEstimation.cpp")
#setup_custom_test_program(test_MultiArcStateEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}")
#target_link_libraries(test_MultiArcStateEstimation ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_HybridArcStateEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}/tests/unitTestHybridArcStateEstimation.cpp")
#setup_custom_test_program(test_HybridArcStateEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}")
#target_link_libraries(test_HybridArcStateEstimation ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_ParameterInfluenceDetermination "${SRCROOT}${ORBITDETERMINATIONDIR}/tests/unitTestParameterInfluenceDetermination.cpp")
#setup_custom_test_program(test_ParameterInfluenceDetermination "${SRCROOT}${ORBITDETERMINATIONDIR}")
#target_link_libraries(test_ParameterInfluenceDetermination ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_BiasEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}/tests/unitTestBiasEstimation.cpp")
#setup_custom_test_program(test_BiasEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}")
#target_link_libraries(test_BiasEstimation ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_RotationalStateEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}/tests/unitTestRotationalStateEstimation.cpp")
#setup_custom_test_program(test_RotationalStateEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}")
#target_link_libraries(test_RotationalStateEstimation ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_ArcwiseEnvironmentEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}/tests/unitTestArcwiseEnvironmentParameters.cpp")
#setup_custom_test_program(test_ArcwiseEnvironmentEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}")
#target_link_libraries(test_ArcwiseEnvironmentEstimation ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_DesaturationDeltaVsEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}/tests/unitTestDesaturationDeltaVsEstimation.cpp")
#setup_custom_test_program(test_DesaturationDeltaVsEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}")
#target_link_libraries(test_DesaturationDeltaVsEstimation ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_FullPlanetaryRotationModelEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}/tests/unitTestFullPlanetaryRotationModelEstimation.cpp")
#setup_custom_test_program(test_FullPlanetaryRotationModelEstimation "${SRCROOT}${ORBITDETERMINATIONDIR}")
#target_link_libraries(test_FullPlanetaryRotationModelEstimation ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#if( BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
#
#add_executable(test_EstimationFromPositionDoubleLongDouble "${SRCROOT}${ORBITDETERMINATIONDIR}/tests/unitTestEstimationFromIdealDataDoubleLongDouble.cpp")
#setup_custom_test_program(test_EstimationFromPositionDoubleLongDouble "${SRCROOT}${ORBITDETERMINATIONDIR}")
#target_link_libraries(test_EstimationFromPositionDoubleLongDouble ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_EstimationFromPositionTimeDouble "${SRCROOT}${ORBITDETERMINATIONDIR}/tests/unitTestEstimationFromIdealDataTimeDouble.cpp")
#setup_custom_test_program(test_EstimationFromPositionTimeDouble "${SRCROOT}${ORBITDETERMINATIONDIR}")
#target_link_libraries(test_EstimationFromPositionTimeDouble ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_EstimationFromPositionTimeLongDouble "${SRCROOT}${ORBITDETERMINATIONDIR}/tests/unitTestEstimationFromIdealDataTimeLongDouble.cpp")
#setup_custom_test_program(test_EstimationFromPositionTimeLongDouble "${SRCROOT}${ORBITDETERMINATIONDIR}")
#target_link_libraries(test_EstimationFromPositionTimeLongDouble ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#endif( )
#
#endif( )

