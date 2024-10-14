#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

TUDAT_ADD_TEST_CASE(CentralBodyData
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(CowellStateDerivative
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(EnckeStateDerivative PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(GaussStateDerivative PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(UnifiedStateModelStateDerivative PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(EnvironmentUpdater PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(BodyMassPropagation PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(RigidBodyProperties PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(MultiTypeStatePropagation PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(StoppingConditions PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(CustomStatePropagation PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(MultiArcDynamics PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(HybridArcDynamics PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(PropagationTerminationReason PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(RotationalDynamicsPropagator PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(CR3BPPropagation PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(PropagationTerminationCheckOnFinalStep PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(GlobalFrameOrigin PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(ForwardsBackwardsIntegration PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(ExactTermination PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(NonSequentialPropagation PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(PropagationResultsSaving PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(IntegratorSteps PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(RadiationPressurePropagation PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(StateDerivativeRestrictedThreeBodyProblem PRIVATE_LINKS tudat_mission_segments tudat_root_finders tudat_propagators tudat_numerical_integrators tudat_basic_astrodynamics tudat_input_output)

#TUDAT_ADD_TEST_CASE(FullPropagationRestrictedThreeBodyProblem PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

#TUDAT_ADD_TEST_CASE(FullPropagationLambertTargeter PRIVATE_LINKS tudat_trajectory_design tudat_mission_segments tudat_ephemerides tudat_basic_astrodynamics tudat_basic_mathematics PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

#TUDAT_ADD_TEST_CASE(FullPropagationPatchedConicsTrajectory PRIVATE_LINKS tudat_trajectory_design tudat_mission_segments tudat_ephemerides tudat_basic_astrodynamics tudat_basic_mathematics PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})


if(TUDAT_BUILD_WITH_ESTIMATION_TOOLS )
#TUDAT_ADD_TEST_CASE(SequentialVariationalEquations PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(VariationalEquations PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(MultiArcVariationalEquations PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

#TUDAT_ADD_TEST_CASE(MultiArcMultiBodyVariationalEquations PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(HybridArcVariationalEquations PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(HybridArcStateTransitionMatrixInterface PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

#TUDAT_ADD_TEST_CASE(HybridArcMultiBodyVariationalEquations PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

#TUDAT_ADD_TEST_CASE(HybridArcMultiBodyStateTransitionMatrixInterface PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(DependentVariableOutput PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(DependentVariablesInterface PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(NonSequentialVariationalEquations PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

endif( )
