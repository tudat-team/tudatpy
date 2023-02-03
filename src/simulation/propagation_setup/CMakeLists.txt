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
set(propagation_HEADERS
        createAccelerationModels.h
        dynamicsSimulator.h
        createTorqueModel.h
        createStateDerivativeModel.h
        createEnvironmentUpdater.h
        propagationOutput.h
        propagationTerminationSettings.h
        torqueSettings.h
        createMassRateModels.h
        propagationTermination.h
        propagationSettings.h
        propagationPrintSettings.h
        propagationProcessingSettings.h
        accelerationSettings.h
        propagationOutputSettings.h
        setNumericallyIntegratedStates.h
        environmentUpdater.h
        dependentVariablesInterface.h
        )

# Add header files.
set(propagation_SOURCES
        createStateDerivativeModel.cpp
        propagationSettings.cpp
        setNumericallyIntegratedStates.cpp
        dynamicsSimulator.cpp
        createTorqueModel.cpp
        propagationOutputSettings.cpp
        createMassRateModels.cpp
        createAccelerationModels.cpp
        createEnvironmentUpdater.cpp
        propagationTermination.cpp
        propagationOutput.cpp
        environmentUpdater.cpp
        dependentVariablesInterface.cpp
        )

# Add library.
TUDAT_ADD_LIBRARY("propagation_setup"
        "${propagation_SOURCES}"
        "${propagation_HEADERS}"
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}" "${CSpice_INCLUDE_DIRS}" "${Sofa_INCLUDE_DIRS}"
        )
