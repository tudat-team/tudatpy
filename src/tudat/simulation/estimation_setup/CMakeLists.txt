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
        createCartesianStatePartials.h
        createTorquePartials.h
        createStateDerivativePartials.h
        estimatableParameterSettings.h
        orbitDeterminationTestCaseUtilities.h
        createObservationPartials.h
        podProcessing.h
        determinePostFitParameterInfluence.h
        variationalEquationsSolver.h
        variationalEquationsSolverBase.h
        singleArcVariationalEquationsSolver.h
        multiArcVariationalEquationsSolver.h
        hybridArcVariationalEquationsSolver.h
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
        processTrackingTxtFile.h
        interArcStateContinuityConstraintSettings.h
        interArcContinuityConstraint.h
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
        orbitDeterminationManager.cpp
        createCartesianStatePartials.cpp
        createStateDerivativePartials.cpp
        orbitDeterminationTestCases.cpp
        createNumericalSimulator.cpp
        createAccelerationPartials.cpp
        createObservationModel.cpp
        createInverseAprioriCovariance.cpp
        createObservationViability.cpp
        observationSimulationSettings.cpp
        observationOutputSettings.cpp
        observationOutput.cpp
        simulateObservations.cpp
        createClockPartials.cpp
        processOdfFile.cpp
        processTrackingTxtFile.cpp
        interArcStateContinuityConstraintSettings.cpp
        )

if (TUDAT_BUILD_EXPLICIT_INSTANTIATIONS)
    list(APPEND estimation_SOURCES
            variationalEquationsSolverDoubleDoubleInstantiation.cpp
            singleArcVariationalEquationsSolverDoubleDoubleInstantiation.cpp
            multiArcVariationalEquationsSolverDoubleDoubleInstantiation.cpp
            hybridArcVariationalEquationsSolverDoubleDoubleInstantiation.cpp
            createObservationManagersBaseDoubleDoubleInstantiation.cpp
            )
endif ()

# Add library.
TUDAT_ADD_LIBRARY("estimation_setup"
        "${estimation_SOURCES}"
        "${estimation_HEADERS}"
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
