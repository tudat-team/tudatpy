 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Add source files.
set(PROPAGATORS_SOURCES
  "${SRCROOT}${PROPAGATORSDIR}/nBodyStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyCowellStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyEnckeStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyGaussKeplerStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyGaussModifiedEquinoctialStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/dynamicsStateDerivativeModel.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyUnifiedStateModelQuaternionsStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyUnifiedStateModelModifiedRodriguesParametersStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyUnifiedStateModelExponentialMapStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/variationalEquations.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/stateTransitionMatrixInterface.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/environmentUpdateTypes.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/singleStateTypeDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/rotationalMotionStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/rotationalMotionQuaternionsStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/rotationalMotionModifiedRodriguesParametersStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/rotationalMotionExponentialMapStateDerivative.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/stateDerivativeCircularRestrictedThreeBodyProblem.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/integrateEquations.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/dynamicsStateDerivativeModel.cpp"
)

# Add header files.
set(PROPAGATORS_HEADERS
  "${SRCROOT}${PROPAGATORSDIR}/centralBodyData.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyCowellStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyEnckeStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyGaussKeplerStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyGaussModifiedEquinoctialStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyUnifiedStateModelQuaternionsStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyUnifiedStateModelModifiedRodriguesParametersStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/nBodyUnifiedStateModelExponentialMapStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/dynamicsStateDerivativeModel.h"
  "${SRCROOT}${PROPAGATORSDIR}/singleStateTypeDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/integrateEquations.h"
  "${SRCROOT}${PROPAGATORSDIR}/bodyMassStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/variationalEquations.h"
  "${SRCROOT}${PROPAGATORSDIR}/stateTransitionMatrixInterface.h"
  "${SRCROOT}${PROPAGATORSDIR}/environmentUpdateTypes.h"
  "${SRCROOT}${PROPAGATORSDIR}/customStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/rotationalMotionStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/rotationalMotionQuaternionsStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/rotationalMotionModifiedRodriguesParametersStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/rotationalMotionExponentialMapStateDerivative.h"
  "${SRCROOT}${PROPAGATORSDIR}/stateDerivativeCircularRestrictedThreeBodyProblem.h"
  "${SRCROOT}${PROPAGATORSDIR}/getZeroProperModeRotationalInitialState.h"
)

# Add static libraries.
add_library(tudat_propagators STATIC ${PROPAGATORS_SOURCES} ${PROPAGATORS_HEADERS})
setup_tudat_library_target(tudat_propagators "${SRCROOT}${PROPAGATORSDIR}")

# Add unit tests.
add_executable(test_CentralBodyData "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestCentralBodyData.cpp")
setup_custom_test_program(test_CentralBodyData "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_CentralBodyData tudat_propagators tudat_basic_mathematics ${Boost_LIBRARIES})

if(USE_CSPICE)

if( BUILD_PROPAGATION_TESTS )

add_executable(test_CowellStateDerivative "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestCowellStateDerivative.cpp")
setup_custom_test_program(test_CowellStateDerivative "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_CowellStateDerivative ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_EnckeStateDerivative "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestEnckeStateDerivative.cpp")
setup_custom_test_program(test_EnckeStateDerivative "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_EnckeStateDerivative ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_GaussStateDerivative "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestGaussStateDerivative.cpp")
setup_custom_test_program(test_GaussStateDerivative "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_GaussStateDerivative ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_UnifiedStateModelStateDerivative "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestUnifiedStateModelStateDerivative.cpp")
setup_custom_test_program(test_UnifiedStateModelStateDerivative "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_UnifiedStateModelStateDerivative ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

if( BUILD_WITH_ESTIMATION_TOOLS )
add_executable(test_SequentialVariationEquationIntegration "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestSequentialVariationalEquationIntegration.cpp")
setup_custom_test_program(test_SequentialVariationEquationIntegration "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_SequentialVariationEquationIntegration ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_VariationalEquations "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestVariationalEquationPropagation.cpp")
setup_custom_test_program(test_VariationalEquations "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_VariationalEquations ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_MultiArcVariationalEquations "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestMultiArcVariationalEquationPropagation.cpp")
setup_custom_test_program(test_MultiArcVariationalEquations "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_MultiArcVariationalEquations ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_HybridArcVariationalEquations "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestHybridArcVariationalEquationPropagation.cpp")
setup_custom_test_program(test_HybridArcVariationalEquations "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_HybridArcVariationalEquations ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})


add_executable(test_DependentVariableOutput "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestDependentVariableOutput.cpp")
setup_custom_test_program(test_DependentVariableOutput "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_DependentVariableOutput ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

endif( )

add_executable(test_EnvironmentModelComputations "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestEnvironmentUpdater.cpp")
setup_custom_test_program(test_EnvironmentModelComputations "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_EnvironmentModelComputations ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_BodyMassPropagation "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestBodyMassPropagation.cpp")
setup_custom_test_program(test_BodyMassPropagation "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_BodyMassPropagation ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_MultiTypeStatePropagation "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestMultiTypeStatePropagation.cpp")
setup_custom_test_program(test_MultiTypeStatePropagation "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_MultiTypeStatePropagation ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_StoppingConditions "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestStoppingConditions.cpp")
setup_custom_test_program(test_StoppingConditions "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_StoppingConditions ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_CustomStatePropagation "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestCustomStatePropagation.cpp")
setup_custom_test_program(test_CustomStatePropagation "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_CustomStatePropagation ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_MultiArcDynamics "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestMultiArcDynamics.cpp")
setup_custom_test_program(test_MultiArcDynamics "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_MultiArcDynamics ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_HybridArcDynamics "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestHybridArcDynamics.cpp")
setup_custom_test_program(test_HybridArcDynamics "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_HybridArcDynamics ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_PropagationTerminationReason "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestPropagationTerminationReason.cpp")
setup_custom_test_program(test_PropagationTerminationReason "${SRCROOT}${PROPAGATORSDIR}/")
target_link_libraries(test_PropagationTerminationReason ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_RotationalDynamicsPropagation "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestRotationalDynamicsPropagator.cpp")
setup_custom_test_program(test_RotationalDynamicsPropagation "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_RotationalDynamicsPropagation ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_CR3BPPropagation "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestCR3BPPropagation.cpp")
setup_custom_test_program(test_CR3BPPropagation "${SRCROOT}${PROPAGATORSDIR}/")
target_link_libraries(test_CR3BPPropagation ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_PropagationTerminationCheckOnFinalStep "${SRCROOT}${PROPAGATORSDIR}/UnitTests/uniTestPropagationTerminationCheckOnFinalStep.cpp")
setup_custom_test_program(test_PropagationTerminationCheckOnFinalStep "${SRCROOT}${PROPAGATORSDIR}/")
target_link_libraries(test_PropagationTerminationCheckOnFinalStep ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_GlobalFrameOrigin "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestGlobalFrameOrigin.cpp")
setup_custom_test_program(test_GlobalFrameOrigin "${SRCROOT}${PROPAGATORSDIR}/")
target_link_libraries(test_GlobalFrameOrigin ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_ForwardBackwardIntegration "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestForwardsBackwardsIntegration.cpp")
setup_custom_test_program(test_ForwardBackwardIntegration "${SRCROOT}${PROPAGATORSDIR}/")
target_link_libraries(test_ForwardBackwardIntegration ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_ExactTermination "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestExactTermination.cpp")
setup_custom_test_program(test_ExactTermination "${SRCROOT}${PROPAGATORSDIR}/")
target_link_libraries(test_ExactTermination ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

endif( )

endif( )

add_executable(test_StateDerivativeRestrictedThreeBodyProblem "${SRCROOT}${PROPAGATORSDIR}/UnitTests/unitTestStateDerivativeCircularRestrictedThreeBodyProblem.cpp")
setup_custom_test_program(test_StateDerivativeRestrictedThreeBodyProblem "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_StateDerivativeRestrictedThreeBodyProblem tudat_mission_segments tudat_root_finders tudat_propagators tudat_numerical_integrators tudat_basic_astrodynamics tudat_input_output ${Boost_LIBRARIES})

