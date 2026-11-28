#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Add source files.
set(estimation_HEADERS
        simulateObservations.h
        observationOutputSettings.h
        createObservationModel.h
        createObservationViability.h
        createLightTimeCalculator.h
        createAccelerationPartials.h
        createLightTimeCorrection.h
        createAtmosphericLightTimeCorrection.h
        createObservationManager.h
        simulatePseudoObservations.h
        orbitDeterminationManager.h
        estimationInterfacesForwardDeclarations.h
        observationInterfacesForwardDeclarations.h
        orbitDeterminationManagerForwardDeclarations.h
        orbitDeterminationManagerHelpers.h
        orbitDeterminationManagerConstructionImplementation.h
        orbitDeterminationManagerUtilitiesImplementation.h
        orbitDeterminationManagerCovarianceImplementation.h
        orbitDeterminationManagerEstimationImplementation.h
        createLightTimeCorrectionPartials.h
        createEstimatableParameters.h
        createEstimatableParametersInterface.h
        createEstimatableParametersFactory.h
        getAssociatedMultiArcParameter.h
        setInitialStateVectorFromParameterSet.h
        wrapObservationResiduals.h
        createStateTransitionAndSensitivityMatrixInterpolator.h
        createCartesianStatePartials.h
        createTorquePartials.h
        createStateDerivativePartials.h
        estimatableParameterSettings.h
        estimatableParameterSettingsBase.h
        orbitDeterminationTestCaseUtilities.h
        createObservationPartials.h
        podProcessing.h
        determinePostFitParameterInfluence.h
        variationalEquationsSolver.h
        variationalEquationsSolverBase.h
        variationalEquationsSolverBaseImplementation.h
        singleArcVariationalEquationsSolver.h
        singleArcVariationalEquationsSolverImplementation.h
        multiArcVariationalEquationsSolver.h
        multiArcVariationalEquationsSolverImplementation.h
        hybridArcVariationalEquationsSolver.h
        hybridArcVariationalEquationsSolverImplementation.h
        createNumericalSimulator.h
        createEulerAngleObservationPartials.h
        createInverseAprioriCovariance.h
        createDopplerPartials.h
        createNWayRangePartials.h
        observationOutputSettings.h
        observationOutput.h
        observationSimulationSettings.h
        observationsProcessing.h
        observations.h
        createDirectObservationPartials.h
        createPositionPartialScaling.h
        createClockPartials.h
        processOdfFile.h
        processPsfFile.h
        processTrackingTxtFile.h
        trackingDataWeighting.h
        createObservationCollection.h
        interArcStateContinuityConstraintSettings.h
        interArcContinuityConstraint.h
        compressDopplerObservationCollection.h
        )

# Add header files.
set(estimation_SOURCES
        createEulerAngleObservationPartials.cpp
        createTorquePartials.cpp
        variationalEquationsSolver.cpp
        createObservationManager.cpp
        createEstimatableParameters.cpp
        createLightTimeCorrection.cpp
        createAtmosphericLightTimeCorrection.cpp
        createLightTimeCorrectionPartials.cpp
        createCartesianStatePartials.cpp
        createStateDerivativePartials.cpp
        orbitDeterminationTestCases.cpp
        createAccelerationPartials.cpp
        createObservationModel.cpp
        createObservationCollection.cpp
        trackingDataWeighting.cpp
        createInverseAprioriCovariance.cpp
        createObservationViability.cpp
        observationSimulationSettings.cpp
        observationOutputSettings.cpp
        observationOutput.cpp
        simulateObservations.cpp
        createClockPartials.cpp
        processOdfFile.cpp
        processTrackingTxtFile.cpp
        interArcContinuityConstraint.cpp
        interArcStateContinuityConstraintSettings.cpp
        )

if (TUDAT_BUILD_WITH_SERIALIZATION)
    list(APPEND estimation_SOURCES observationOutputSerialization.cpp)
endif ()

if (TUDAT_BUILD_TESTS_WITH_EXPLICIT_INSTANTIATIONS)
    set(estimation_INSTANTIATION_SOURCES
            variationalEquationsSolverDoubleDoubleInstantiation.cpp
            singleArcVariationalEquationsSolverDoubleDoubleInstantiation.cpp
            multiArcVariationalEquationsSolverDoubleDoubleInstantiation.cpp
            hybridArcVariationalEquationsSolverDoubleDoubleInstantiation.cpp
            createObservationManagersBaseDoubleDoubleInstantiation.cpp
            )
    list(APPEND estimation_SOURCES ${estimation_INSTANTIATION_SOURCES})
    # The instantiations largely overlap, so compile them as one translation unit
    set_source_files_properties(${estimation_INSTANTIATION_SOURCES} PROPERTIES UNITY_GROUP double_double_instantiations)
    # These share little with the group above (adding them to it raises its peak memory), so they are compiled on their own
    list(APPEND estimation_SOURCES
            createParametersToEstimateDoubleDoubleInstantiation.cpp
            simulateObservationsDoubleDoubleInstantiation.cpp
            orbitDeterminationManager.cpp
            createNumericalSimulator.cpp
            )
endif ()

# Add library.
TUDAT_ADD_LIBRARY("estimation_setup"
        "${estimation_SOURCES}"
        "${estimation_HEADERS}"
        PUBLIC_SERIALIZATION
#        PRIVATE_LINKS "Boost::boost"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )

# Not with a global unity build (CMAKE_UNITY_BUILD), which already batches all sources of the library
if (TUDAT_BUILD_TESTS_WITH_EXPLICIT_INSTANTIATIONS AND TUDAT_BUILD_WITH_UNITY_GROUPS AND NOT CMAKE_UNITY_BUILD)
    set_target_properties(tudat_estimation_setup PROPERTIES UNITY_BUILD ON UNITY_BUILD_MODE GROUP)
endif ()
