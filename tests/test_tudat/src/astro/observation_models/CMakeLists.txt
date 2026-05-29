#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

if (TUDAT_BUILD_WITH_PCH AND TUDAT_BUILD_WITH_TEST_PCH)
    add_library(tudat_observation_models_test_case_support
        observationModelsTestCaseSupport.cpp
    )

    target_include_directories(tudat_observation_models_test_case_support PUBLIC
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/tests/test_tudat/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    )

    target_include_directories(tudat_observation_models_test_case_support SYSTEM PRIVATE
        "${EIGEN3_INCLUDE_DIRS}"
        "${Boost_INCLUDE_DIRS}"
        "${CSpice_INCLUDE_DIRS}"
        "${Sofa_INCLUDE_DIRS}"
        "${TudatResources_INCLUDE_DIRS}"
    )

    target_link_libraries(tudat_observation_models_test_case_support PUBLIC
        ${Tudat_ESTIMATION_LIBRARIES}
    )
    target_link_libraries(tudat_observation_models_test_case_support PRIVATE
        "${Boost_LIBRARIES}"
    )

    set_target_properties(tudat_observation_models_test_case_support PROPERTIES POSITION_INDEPENDENT_CODE OFF)
    target_compile_options(tudat_observation_models_test_case_support PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-fPIE>
    )

    target_precompile_headers(tudat_observation_models_test_case_support PRIVATE
        <tudat/simulation/environment_setup/body.h>
        <tudat/simulation/environment_setup/createBodiesFactory.h>
        <tudat/simulation/environment_setup/defaultBodies.h>
        <tudat/simulation/estimation_setup/createObservationModelFactory.h>
        <tudat/simulation/estimation_setup/processOdfFile.h>
        <tudat/simulation/estimation_setup/simulateObservations.h>
        <tudat/io/basicInputOutput.h>
        <Eigen/Core>
        <limits>
        <map>
        <memory>
        <string>
        <vector>
    )

    set(TUDAT_TEST_PCH_REUSE_TARGET tudat_observation_models_test_case_support)
endif ()


TUDAT_ADD_TEST_CASE(LightTimeSolution PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(AngularPositionModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(OneWayRangeModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(OneWayDifferencedRangeModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(DopplerModels PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(PositionObservationModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(EulerAngleObservationModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(VelocityObservationModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(SimulatedObservationNoise PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(NWayRangeObservationModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(NWayRangeRateObservationModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(DopplerMeasuredFrequencyObservationModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(RelativeAngularPositionModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(RelativePositionModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(ObservationViabilityCalculators PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(PerArcObservationSimulation PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(TimeBias PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(ClockModels PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(AtmosphereCorrection PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(SolarPlasmaCorrection PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(ObservationRoundingErrors PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(IfmsMex PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(ObservationDependentVariables PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(LightTimeCorrectionComponents PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(ProcessOdfData PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

# TUDAT_ADD_TEST_CASE(DsnOdfGrail PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

#TUDAT_ADD_TEST_CASE(DsnOdfGrailStatistics PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(DsnNWayDopplerObservationModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(DsnNWayRangeObservationModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(DifferencedTimeOfArrivalObservationModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(OneWayDopplerMeasuredFrequencyObservationModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})

TUDAT_ADD_TEST_CASE(DifferencedFrequencyOfArrivalObservationModel PRIVATE_LINKS ${Tudat_ESTIMATION_LIBRARIES})