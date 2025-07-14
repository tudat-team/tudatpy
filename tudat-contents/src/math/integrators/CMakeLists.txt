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
set(numerical_integrators_SOURCES
        "rungeKuttaCoefficients.cpp"
        "createNumericalIntegrator.cpp"
        "bulirschStoerVariableStepsizeIntegrator.cpp"
        "numericalIntegrator.cpp"
        "reinitializableNumericalIntegrator.cpp"
        "rungeKutta4Integrator.cpp"
        "rungeKuttaFixedStepSizeIntegrator.cpp"
        "rungeKuttaVariableStepSizeIntegrator.cpp"
        "adamsBashforthMoultonIntegrator.cpp"
        "bulirschStoerVariableStepsizeIntegrator.cpp"
        )

# Add header files.
set(numerical_integrators_HEADERS
        "adamsBashforthMoultonIntegrator.h"
        "createNumericalIntegrator.h"
        "bulirschStoerVariableStepsizeIntegrator.h"
        "euler.h"
        "numericalIntegrator.h"
        "reinitializableNumericalIntegrator.h"
        "rungeKutta4Integrator.h"
        "rungeKuttaCoefficients.h"
        "rungeKuttaFixedStepSizeIntegrator.h"
        "rungeKuttaVariableStepSizeIntegrator.h"
        )

# Add library.
TUDAT_ADD_LIBRARY("numerical_integrators"
        "${numerical_integrators_SOURCES}"
        "${numerical_integrators_HEADERS}"
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
