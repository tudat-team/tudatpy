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
set(propagators_SOURCES
        "nBodyStateDerivative.cpp"
        "nBodyCowellStateDerivative.cpp"
        "nBodyEnckeStateDerivative.cpp"
        "nBodyGaussKeplerStateDerivative.cpp"
        "nBodyGaussModifiedEquinoctialStateDerivative.cpp"
        "dynamicsStateDerivativeModel.cpp"
        "nBodyUnifiedStateModelQuaternionsStateDerivative.cpp"
        "nBodyUnifiedStateModelModifiedRodriguesParametersStateDerivative.cpp"
        "nBodyUnifiedStateModelExponentialMapStateDerivative.cpp"
        "variationalEquations.cpp"
        "stateTransitionMatrixInterface.cpp"
        "environmentUpdateTypes.cpp"
        "singleStateTypeDerivative.cpp"
        "rotationalMotionStateDerivative.cpp"
        "rotationalMotionQuaternionsStateDerivative.cpp"
        "rotationalMotionModifiedRodriguesParametersStateDerivative.cpp"
        "rotationalMotionExponentialMapStateDerivative.cpp"
        "stateDerivativeCircularRestrictedThreeBodyProblem.cpp"
        "integrateEquations.cpp"
        "dynamicsStateDerivativeModel.cpp"
        "propagateCovariance.cpp"
        )

# Add header files.
set(propagators_HEADERS
        "centralBodyData.h"
        "nBodyStateDerivative.h"
        "nBodyCowellStateDerivative.h"
        "nBodyEnckeStateDerivative.h"
        "nBodyGaussKeplerStateDerivative.h"
        "nBodyGaussModifiedEquinoctialStateDerivative.h"
        "nBodyUnifiedStateModelQuaternionsStateDerivative.h"
        "nBodyUnifiedStateModelModifiedRodriguesParametersStateDerivative.h"
        "nBodyUnifiedStateModelExponentialMapStateDerivative.h"
        "dynamicsStateDerivativeModel.h"
        "singleStateTypeDerivative.h"
        "integrateEquations.h"
        "bodyMassStateDerivative.h"
        "variationalEquations.h"
        "stateTransitionMatrixInterface.h"
        "environmentUpdateTypes.h"
        "customStateDerivative.h"
        "rotationalMotionStateDerivative.h"
        "rotationalMotionQuaternionsStateDerivative.h"
        "rotationalMotionModifiedRodriguesParametersStateDerivative.h"
        "rotationalMotionExponentialMapStateDerivative.h"
        "stateDerivativeCircularRestrictedThreeBodyProblem.h"
        "getZeroProperModeRotationalInitialState.h"
        "propagateCovariance.h"
        )

TUDAT_ADD_LIBRARY("propagators"
        "${propagators_SOURCES}"
        "${propagators_HEADERS}"
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
