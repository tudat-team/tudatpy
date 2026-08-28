#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

if (TUDAT_BUILD_WITH_PCH AND TUDAT_BUILD_WITH_TEST_PCH)
    add_library(tudat_propagators_test_case_support
        propagatorsTestCaseSupport.cpp
    )

    target_include_directories(tudat_propagators_test_case_support PUBLIC
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/tests/test_tudat/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    )

    target_include_directories(tudat_propagators_test_case_support SYSTEM PRIVATE
        "${EIGEN3_INCLUDE_DIRS}"
        "${Boost_INCLUDE_DIRS}"
        "${CSpice_INCLUDE_DIRS}"
        "${Sofa_INCLUDE_DIRS}"
        "${TudatResources_INCLUDE_DIRS}"
    )

    target_link_libraries(tudat_propagators_test_case_support PUBLIC
        ${Tudat_PROPAGATION_LIBRARIES}
    )
    if (TUDAT_BUILD_WITH_ESTIMATION_TOOLS)
        target_link_libraries(tudat_propagators_test_case_support PUBLIC
            ${Tudat_ESTIMATION_LIBRARIES}
        )
    endif ()
    target_link_libraries(tudat_propagators_test_case_support PRIVATE
        "${Boost_LIBRARIES}"
    )

    set_target_properties(tudat_propagators_test_case_support PROPERTIES POSITION_INDEPENDENT_CODE OFF)
    target_compile_options(tudat_propagators_test_case_support PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-fPIE>
    )

    target_precompile_headers(tudat_propagators_test_case_support PRIVATE
        <tudat/math/basic/linearAlgebra.h>
        <tudat/math/integrators/createNumericalIntegrator.h>
        <tudat/math/integrators/rungeKuttaCoefficients.h>
        <tudat/simulation/environment_setup/body.h>
        <tudat/simulation/environment_setup/createBodiesFactory.h>
        <tudat/simulation/environment_setup/defaultBodies.h>
        <tudat/simulation/estimation_setup/createNumericalSimulator.h>
        <tudat/simulation/propagation_setup/createAccelerationModels.h>
        <tudat/simulation/propagation_setup/propagationSettings.h>
        <tudat/simulation/propagation_setup/singleArcDynamicsSimulator.h>
        <tudat/astro/basic_astro/accelerationModel.h>
        <tudat/interface/spice/spiceInterface.h>
        <Eigen/Core>
        <limits>
        <map>
        <memory>
        <string>
        <thread>
        <vector>
    )

    set(TUDAT_TEST_PCH_REUSE_TARGET tudat_propagators_test_case_support)
endif ()

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

TUDAT_ADD_TEST_CASE(PropagationTermination PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(RotationalDynamicsPropagator PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(FullTwoBodyCoupledDynamicsPropagation PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(CR3BPPropagation PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(PropagationTerminationCheckOnFinalStep PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(GlobalFrameOrigin PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(ForwardsBackwardsIntegration PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(ExactTermination PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(NonSequentialPropagation PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(PropagationResultsSaving PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(IntegratorSteps PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(StateDerivativeModifier PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

if(TUDAT_BUILD_WITH_ESTIMATION_TOOLS)
    TUDAT_ADD_TEST_CASE(RadiationPressurePropagation PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})
else()
    TUDAT_ADD_TEST_CASE(RadiationPressurePropagation PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})
endif()

TUDAT_ADD_TEST_CASE(RelativisticTimePropagation PRIVATE_LINKS  ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(IntegratedStateInterpolator PRIVATE_LINKS ${Tudat_PROPAGATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(StateDerivativeRestrictedThreeBodyProblem PRIVATE_LINKS tudat_mission_segments tudat_root_finders tudat_propagators tudat_numerical_integrators tudat_basic_astrodynamics tudat_basic_mathematics tudat_statistics tudat_input_output)

if(TUDAT_BUILD_WITH_ESTIMATION_TOOLS )
TUDAT_ADD_TEST_CASE(VariationalEquations PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(FullTwoBodyVariationalEquations PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(MultiArcVariationalEquations PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(MultiArcMultiBodyVariationalEquations PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(HybridArcVariationalEquations PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(HybridArcStateTransitionMatrixInterface PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(HybridArcMultiBodyVariationalEquations PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(HybridArcMultiBodyStateTransitionMatrixInterface PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(DependentVariableOutput PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(DependentVariablesInterface PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(NonSequentialVariationalEquations PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

endif( )
