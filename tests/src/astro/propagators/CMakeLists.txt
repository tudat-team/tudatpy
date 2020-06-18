#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# TODO: Add propagators library tests

# Add static libraries.
#add_library(tudat_propagators STATIC ${PROPAGATORS_SOURCES} ${PROPAGATORS_HEADERS})
#setup_tudat_library_target(tudat_propagators "${SRCROOT}${PROPAGATORSDIR}")
#
#
#include(tudatLinkLibraries)
#
## Add unit tests.
#add_executable(test_CentralBodyData "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestCentralBodyData.cpp")
#setup_custom_test_program(test_CentralBodyData "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_CentralBodyData tudat_propagators tudat_basic_mathematics ${Boost_LIBRARIES})
#
#if(TUDAT_BUILD_WITH_SPICE_INTERFACE)
#
#if( BUILD_WITH_PROPAGATION_TESTS )
#
#add_executable(test_CowellStateDerivative "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestCowellStateDerivative.cpp")
#setup_custom_test_program(test_CowellStateDerivative "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_CowellStateDerivative ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_EnckeStateDerivative "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestEnckeStateDerivative.cpp")
#setup_custom_test_program(test_EnckeStateDerivative "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_EnckeStateDerivative ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_GaussStateDerivative "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestGaussStateDerivative.cpp")
#setup_custom_test_program(test_GaussStateDerivative "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_GaussStateDerivative ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_UnifiedStateModelStateDerivative "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestUnifiedStateModelStateDerivative.cpp")
#setup_custom_test_program(test_UnifiedStateModelStateDerivative "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_UnifiedStateModelStateDerivative ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#if(TUDAT_BUILD_WITH_ESTIMATION_TOOLS )
#add_executable(test_SequentialVariationEquationIntegration "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestSequentialVariationalEquationIntegration.cpp")
#setup_custom_test_program(test_SequentialVariationEquationIntegration "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_SequentialVariationEquationIntegration ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_VariationalEquations "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestVariationalEquationPropagation.cpp")
#setup_custom_test_program(test_VariationalEquations "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_VariationalEquations ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_MultiArcVariationalEquations "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestMultiArcVariationalEquationPropagation.cpp")
#setup_custom_test_program(test_MultiArcVariationalEquations "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_MultiArcVariationalEquations ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_HybridArcVariationalEquations "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestHybridArcVariationalEquationPropagation.cpp")
#setup_custom_test_program(test_HybridArcVariationalEquations "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_HybridArcVariationalEquations ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#
#add_executable(test_DependentVariableOutput "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestDependentVariableOutput.cpp")
#setup_custom_test_program(test_DependentVariableOutput "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_DependentVariableOutput ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
#
#endif( )
#
#add_executable(test_EnvironmentModelComputations "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestEnvironmentUpdater.cpp")
#setup_custom_test_program(test_EnvironmentModelComputations "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_EnvironmentModelComputations ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_BodyMassPropagation "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestBodyMassPropagation.cpp")
#setup_custom_test_program(test_BodyMassPropagation "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_BodyMassPropagation ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_MultiTypeStatePropagation "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestMultiTypeStatePropagation.cpp")
#setup_custom_test_program(test_MultiTypeStatePropagation "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_MultiTypeStatePropagation ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_StoppingConditions "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestStoppingConditions.cpp")
#setup_custom_test_program(test_StoppingConditions "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_StoppingConditions ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_CustomStatePropagation "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestCustomStatePropagation.cpp")
#setup_custom_test_program(test_CustomStatePropagation "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_CustomStatePropagation ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_MultiArcDynamics "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestMultiArcDynamics.cpp")
#setup_custom_test_program(test_MultiArcDynamics "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_MultiArcDynamics ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_HybridArcDynamics "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestHybridArcDynamics.cpp")
#setup_custom_test_program(test_HybridArcDynamics "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_HybridArcDynamics ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_PropagationTerminationReason "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestPropagationTerminationReason.cpp")
#setup_custom_test_program(test_PropagationTerminationReason "${SRCROOT}${PROPAGATORSDIR}/")
#target_link_libraries(test_PropagationTerminationReason ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_RotationalDynamicsPropagation "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestRotationalDynamicsPropagator.cpp")
#setup_custom_test_program(test_RotationalDynamicsPropagation "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_RotationalDynamicsPropagation ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_CR3BPPropagation "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestCR3BPPropagation.cpp")
#setup_custom_test_program(test_CR3BPPropagation "${SRCROOT}${PROPAGATORSDIR}/")
#target_link_libraries(test_CR3BPPropagation ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_PropagationTerminationCheckOnFinalStep "${SRCROOT}${PROPAGATORSDIR}/tests/uniTestPropagationTerminationCheckOnFinalStep.cpp")
#setup_custom_test_program(test_PropagationTerminationCheckOnFinalStep "${SRCROOT}${PROPAGATORSDIR}/")
#target_link_libraries(test_PropagationTerminationCheckOnFinalStep ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_GlobalFrameOrigin "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestGlobalFrameOrigin.cpp")
#setup_custom_test_program(test_GlobalFrameOrigin "${SRCROOT}${PROPAGATORSDIR}/")
#target_link_libraries(test_GlobalFrameOrigin ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_ForwardBackwardIntegration "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestForwardsBackwardsIntegration.cpp")
#setup_custom_test_program(test_ForwardBackwardIntegration "${SRCROOT}${PROPAGATORSDIR}/")
#target_link_libraries(test_ForwardBackwardIntegration ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_ExactTermination "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestExactTermination.cpp")
#setup_custom_test_program(test_ExactTermination "${SRCROOT}${PROPAGATORSDIR}/")
#target_link_libraries(test_ExactTermination ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#endif( )
#
#endif( )
#
#add_executable(test_StateDerivativeRestrictedThreeBodyProblem "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestStateDerivativeCircularRestrictedThreeBodyProblem.cpp")
#setup_custom_test_program(test_StateDerivativeRestrictedThreeBodyProblem "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_StateDerivativeRestrictedThreeBodyProblem tudat_mission_segments tudat_root_finders tudat_propagators tudat_numerical_integrators tudat_basic_astrodynamics tudat_input_output ${Boost_LIBRARIES})
#
#add_executable(test_FullPropagationRestrictedThreeBodyProblem "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestFullPropagationCircularRestrictedThreeBodyProblem.cpp")
#setup_custom_test_program(test_FullPropagationRestrictedThreeBodyProblem "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_FullPropagationRestrictedThreeBodyProblem ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_FullPropagationLambertTargeter "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestFullPropagationLambertTargeter.cpp")
#setup_custom_test_program(test_FullPropagationLambertTargeter "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_FullPropagationLambertTargeter tudat_trajectory_design tudat_mission_segments tudat_ephemerides tudat_basic_astrodynamics tudat_basic_mathematics ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_FullPropagationPatchedConicsTrajectory "${SRCROOT}${PROPAGATORSDIR}/tests/unitTestFullPropagationPatchedConicsTrajectory.cpp")
#setup_custom_test_program(test_FullPropagationPatchedConicsTrajectory "${SRCROOT}${PROPAGATORSDIR}")
#target_link_libraries(test_FullPropagationPatchedConicsTrajectory tudat_trajectory_design tudat_mission_segments tudat_ephemerides tudat_basic_astrodynamics tudat_basic_mathematics ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#
