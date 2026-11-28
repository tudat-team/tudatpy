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
        dynamicsSimulatorBase.h
        dynamicsSimulatorBaseImplementation.h
        singleArcDynamicsSimulator.h
        singleArcDynamicsSimulatorImplementation.h
        multiArcDynamicsSimulator.h
        multiArcDynamicsSimulatorImplementation.h
        hybridArcDynamicsSimulator.h
        hybridArcDynamicsSimulatorImplementation.h
        createTorqueModel.h
        createStateDerivativeModel.h
        createEnvironmentUpdater.h
        propagationOutput.h
        propagationTerminationSettings.h
        torqueSettings.h
        createMassRateModels.h
        propagationTermination.h
        propagationTerminationDetails.h
        propagationSettings.h
        propagatorType.h
        getProcessedStateStrings.h
        propagationPrintSettings.h
        propagationProcessingSettings.h
        accelerationSettings.h
        propagationOutputSettings.h
        setNumericallyIntegratedStates.h
        environmentUpdater.h
        dependentVariablesInterface.h
        getInitialStateOfBody.h
        addEmptyTabulatedEphemeris.h
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
        accelerationSettings.cpp
        createEnvironmentUpdater.cpp
        propagationTermination.cpp
        propagationOutput.cpp
        environmentUpdater.cpp
        dependentVariablesInterface.cpp
        )

set(PUBLIC_LINKS_LIST)
if(TUDAT_BUILD_WITH_MCD_INTERFACE)
    list(APPEND PUBLIC_LINKS_LIST tudat_mcd_interface)
endif()

if (TUDAT_BUILD_WITH_SERIALIZATION)
    list(APPEND propagation_SOURCES propagationSettingsSerialization.cpp)
endif ()

if (TUDAT_BUILD_TESTS_WITH_EXPLICIT_INSTANTIATIONS)
    set(propagation_INSTANTIATION_SOURCES
            dynamicsSimulatorDoubleDoubleInstantiation.cpp
            singleArcDynamicsSimulatorDoubleDoubleInstantiation.cpp
            multiArcDynamicsSimulatorDoubleDoubleInstantiation.cpp
            hybridArcDynamicsSimulatorDoubleDoubleInstantiation.cpp
            )
    list(APPEND propagation_SOURCES ${propagation_INSTANTIATION_SOURCES})
    # The instantiations largely overlap, so compile them as one translation unit
    set_source_files_properties(${propagation_INSTANTIATION_SOURCES} PROPERTIES UNITY_GROUP double_double_instantiations)
endif ()

# Add library.
TUDAT_ADD_LIBRARY("propagation_setup"
        "${propagation_SOURCES}"
        "${propagation_HEADERS}"
        PUBLIC_LINKS ${PUBLIC_LINKS_LIST}
        PUBLIC_SERIALIZATION
#        PRIVATE_LINKS "Boost::boost"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}" "${CSpice_INCLUDE_DIRS}" "${Sofa_INCLUDE_DIRS}"
        )

# Not with a global unity build (CMAKE_UNITY_BUILD), which already batches all sources of the library
if (TUDAT_BUILD_TESTS_WITH_EXPLICIT_INSTANTIATIONS AND TUDAT_BUILD_WITH_UNITY_GROUPS AND NOT CMAKE_UNITY_BUILD)
    set_target_properties(tudat_propagation_setup PROPERTIES UNITY_BUILD ON UNITY_BUILD_MODE GROUP)
endif ()
