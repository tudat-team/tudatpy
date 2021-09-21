#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

# Set the source files.
set(observation_models_SOURCES
        "lightTimeSolution.cpp"
        "observableTypes.cpp"
        "linkTypeDefs.cpp"
        "observationViabilityCalculator.cpp"
        "observationManager.cpp"
        "observationSimulator.cpp"
        "observationModel.cpp"
        "corrections/firstOrderRelativisticCorrection.cpp"
        )

# Set the header files.
set(observation_models_HEADERS
        "angularPositionObservationModel.h"
        "lightTimeSolution.h"
        "linkTypeDefs.h"
        "observableTypes.h"
        "observationModel.h"
        "observationManager.h"
        "observationSimulator.h"
        "observationBias.h"
        "nWayRangeObservationModel.h"
        "oneWayRangeObservationModel.h"
        "oneWayDopplerObservationModel.h"
        "oneWayDifferencedRangeRateObservationModel.h"
        "positionObservationModel.h"
        "velocityObservationModel.h"
        "twoWayDopplerObservationModel.h"
        "eulerAngleObservationModel.h"
        "observationViabilityCalculator.h"
        "corrections/lightTimeCorrection.h"
        "corrections/firstOrderRelativisticCorrection.h"
        )

TUDAT_ADD_LIBRARY("observation_models"
        "${observation_models_SOURCES}"
        "${observation_models_HEADERS}"
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
