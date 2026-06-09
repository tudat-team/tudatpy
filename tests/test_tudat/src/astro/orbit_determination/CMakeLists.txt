 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

#set(ACCELERATIONPARTIALSDIR "${ORBITDETERMINATIONDIR}/AccelerationPartials")
#set(TORQUEPARTIALSDIR "${ORBITDETERMINATIONDIR}/RotationalDynamicsPartials")
#set(ESTIMATABLEPARAMETERSDIR "${ORBITDETERMINATIONDIR}/EstimatableParameters")
#set(OBSERVATIONPARTIALSDIR "${ORBITDETERMINATIONDIR}/ObservationPartials")
#set(LIGHTTIMECORRECTIONPARTIALSDIR "${ORBITDETERMINATIONDIR}/LightTimeCorrectionPartials")

add_library(tudat_estimation_test_case_support
    executePlanetaryParameterEstimationTestCase.cpp
)

target_include_directories(tudat_estimation_test_case_support PUBLIC
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/tests/test_tudat/include>
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_include_directories(tudat_estimation_test_case_support SYSTEM PRIVATE
    "${EIGEN3_INCLUDE_DIRS}"
    "${Boost_INCLUDE_DIRS}"
    "${CSpice_INCLUDE_DIRS}"
    "${Sofa_INCLUDE_DIRS}"
    "${TudatResources_INCLUDE_DIRS}"
)

target_link_libraries(tudat_estimation_test_case_support PUBLIC
    ${Tudat_ESTIMATION_LIBRARIES}
)
target_link_libraries(tudat_estimation_test_case_support PRIVATE
    "${Boost_LIBRARIES}"
)

if (TUDAT_BUILD_WITH_PCH AND TUDAT_BUILD_WITH_TEST_PCH)
    set_target_properties(tudat_estimation_test_case_support PROPERTIES POSITION_INDEPENDENT_CODE OFF)
    target_compile_options(tudat_estimation_test_case_support PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-fPIE>
    )

    target_precompile_headers(tudat_estimation_test_case_support PRIVATE
        <tudat/simulation/environment_setup/createBodiesFactory.h>
        <tudat/simulation/environment_setup/defaultBodies.h>
        <tudat/simulation/environment_setup/createGroundStations.h>
        <tudat/simulation/estimation_setup/createEstimatableParametersFactory.h>
        <tudat/simulation/estimation_setup/createObservationModelFactory.h>
        <tudat/simulation/estimation_setup/createNumericalSimulator.h>
        <tudat/simulation/estimation_setup/orbitDeterminationManager.h>
        <tudat/simulation/estimation_setup/podProcessing.h>
        <tudat/simulation/estimation_setup/simulateObservations.h>
        <tudat/astro/propagators/propagateCovariance.h>
        <tudat/interface/spice/spiceInterface.h>
        <Eigen/Core>
        <Eigen/Geometry>
        <limits>
        <map>
        <memory>
        <string>
        <thread>
        <vector>
    )
    set(TUDAT_TEST_PCH_REUSE_TARGET tudat_estimation_test_case_support)
endif ()

# Add subdirectories.
add_subdirectory(acceleration_partials)
add_subdirectory(rotational_dynamics_partials)
add_subdirectory(observation_partials)

set(ORBIT_DETERMINATION_TEST_PRIVATE_LINKS
    ${Tudat_ESTIMATION_LIBRARIES}
    tudat_estimation_test_case_support
)

TUDAT_ADD_TEST_CASE(EstimationInput
    SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/unitTestEstimationInputCovariance.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/unitTestEstimationInputWeightsAndConvergence.cpp
    PRIVATE_LINKS
    ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
    )

TUDAT_ADD_TEST_CASE(EstimationFromIdealDataDoubleDouble
    PRIVATE_LINKS
    ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
    )

TUDAT_ADD_TEST_CASE(TidalPropertyEstimation
    PRIVATE_LINKS
    ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
    )

TUDAT_ADD_TEST_CASE(HybridArcStateEstimation
    PRIVATE_LINKS
    ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
    )

# TUDAT_ADD_TEST_CASE(HybridArcMultiBodyEstimation
#         PRIVATE_LINKS
#         ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
#         )

TUDAT_ADD_TEST_CASE(MultiArcStateEstimation
    PRIVATE_LINKS
    ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
    )

TUDAT_ADD_TEST_CASE(InterArcContinuityConstraint
    PRIVATE_LINKS
    ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
    )

TUDAT_ADD_TEST_CASE(ParameterInfluenceDetermination
    PRIVATE_LINKS
    ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
    )

TUDAT_ADD_TEST_CASE(BiasEstimation
    PRIVATE_LINKS
    ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
    )

TUDAT_ADD_TEST_CASE(TimeBiasEstimation
     PRIVATE_LINKS
     ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
     )

TUDAT_ADD_TEST_CASE(RotationalStateEstimation
    PRIVATE_LINKS
    ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
    )

TUDAT_ADD_TEST_CASE(ArcwiseEnvironmentParameters
    PRIVATE_LINKS
    ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
    )

TUDAT_ADD_TEST_CASE(EstimationArcWiseAtmosphereParameters
     PRIVATE_LINKS
     ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
)

TUDAT_ADD_TEST_CASE(DesaturationDeltaVsEstimation
    PRIVATE_LINKS
    ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
    )

TUDAT_ADD_TEST_CASE(FullPlanetaryRotationModelEstimation
    PRIVATE_LINKS
    ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
    )


TUDAT_ADD_TEST_CASE(NonSequentialStateEstimation
     PRIVATE_LINKS
     ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
     )

 TUDAT_ADD_TEST_CASE(EihEstimation
         PRIVATE_LINKS
         ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
         )
# Example application: estimation of Doppler observations from ODF/DSN data
#TUDAT_ADD_TEST_CASE(EstimationFromDsnData
#    PRIVATE_LINKS
#    ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
#    )

# Example application: estimation of position observables from SPICE data
#TUDAT_ADD_TEST_CASE(EstimationFromSpiceData
#    PRIVATE_LINKS
#    ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
#    )

 TUDAT_ADD_TEST_CASE(ConsiderParameters
         PRIVATE_LINKS
         ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
         )

 TUDAT_ADD_TEST_CASE(FitToSpice
        PRIVATE_LINKS
        ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
 )

 TUDAT_ADD_TEST_CASE(CustomAccelerationVariationalEquations
         PRIVATE_LINKS
         ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
 )


 TUDAT_ADD_TEST_CASE(VehicleReferencePointEstimation
         PRIVATE_LINKS
         ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
 )

 TUDAT_ADD_TEST_CASE(RadiationPressureEstimation
         PRIVATE_LINKS
         ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
 )

 TUDAT_ADD_TEST_CASE(RadiationPressureVariationalEquations
         PRIVATE_LINKS
         ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
 )

 TUDAT_ADD_TEST_CASE(ObservationsProcessing
         PRIVATE_LINKS
         ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
 )

TUDAT_ADD_TEST_CASE(CovariancePropagation
         PRIVATE_LINKS
         ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
 )

 TUDAT_ADD_TEST_CASE(ClockPropertyEstimation
         PRIVATE_LINKS
         ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
         )


 TUDAT_ADD_TEST_CASE(RTGForceVectorEstimation
         PRIVATE_LINKS
         ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
 )


 if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )


    TUDAT_ADD_TEST_CASE(EstimationFromIdealDataTimeDouble
            PRIVATE_LINKS
            ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
    )

    TUDAT_ADD_TEST_CASE(EstimationFromIdealDataDoubleLongDouble
        PRIVATE_LINKS
        ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
        )

    TUDAT_ADD_TEST_CASE(EstimationFromIdealDataTimeLongDouble
        PRIVATE_LINKS
        ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
        )

endif( )

TUDAT_ADD_TEST_CASE(EstimationDragScaling
         PRIVATE_LINKS
         ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
         )

 TUDAT_ADD_TEST_CASE(ParameterPropagatorConsistency
         PRIVATE_LINKS
         ${ORBIT_DETERMINATION_TEST_PRIVATE_LINKS}
 )
