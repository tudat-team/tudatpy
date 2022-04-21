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
        createObservationManager.h
        orbitDeterminationManager.h
        createLightTimeCorrectionPartials.h
        createEstimatableParameters.h
        createCartesianStatePartials.h
        createTorquePartials.h
        createStateDerivativePartials.h
        estimatableParameterSettings.h
        orbitDeterminationTestCases.h
        createOneWayRangePartials.h
        createObservationPartials.h
        podProcessing.h
        determinePostFitParameterInfluence.h
        variationalEquationsSolver.h
        createNumericalSimulator.h
        createEulerAngleObservationPartials.h
        createDopplerPartials.h
        createAngularPositionPartials.h
        createDifferencedOneWayRangeRatePartials.h
        createNWayRangePartials.h
        observationOutputSettings.h
        observationOutput.h
        observationSimulationSettings.h
        observations.h
        )

# Add header files.
set(estimation_SOURCES
        createEulerAngleObservationPartials.cpp
        createTorquePartials.cpp
        variationalEquationsSolver.cpp
        createObservationManager.cpp
        createEstimatableParameters.cpp
        createLightTimeCorrection.cpp
        createDifferencedOneWayRangeRatePartials.cpp
        createLightTimeCorrectionPartials.cpp
        orbitDeterminationManager.cpp
        createCartesianStatePartials.cpp
        createStateDerivativePartials.cpp
        orbitDeterminationTestCases.cpp
        createOneWayRangePartials.cpp
        createNumericalSimulator.cpp
        createAccelerationPartials.cpp
        createObservationModel.cpp
        createObservationViability.cpp
        createDopplerPartials.cpp
        createAngularPositionPartials.cpp
        observationSimulationSettings.cpp
        observationOutputSettings.cpp
        observationOutput.cpp
        simulateObservations.cpp
        )

# Add library.
TUDAT_ADD_LIBRARY("estimation_setup"
        "${estimation_SOURCES}"
        "${estimation_HEADERS}"
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
