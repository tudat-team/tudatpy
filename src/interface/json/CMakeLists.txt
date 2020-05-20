#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Set the source files.
set(json_interface_SOURCES
        # support
        "support/keys.cpp"
        "support/deserialization.cpp"
        "support/path.cpp"
        "support/valueAccess.cpp"
        "support/options.cpp"
        # environment
        "environment/spice.cpp"
        "environment/body.cpp"
        "environment/atmosphere.cpp"
        "environment/ephemeris.cpp"
        "environment/aerodynamics.cpp"
        "environment/gravityField.cpp"
        "environment/rotationModel.cpp"
        "environment/shapeModel.cpp"
        "environment/radiationPressure.cpp"
        "environment/gravityFieldVariation.cpp"
        "environment/groundStations.cpp"
        # propagation
        "propagation/variable.cpp"
        "propagation/termination.cpp"
        "propagation/acceleration.cpp"
        "propagation/massRateModel.cpp"
        "propagation/thrust.cpp"
        "propagation/torque.cpp"
        "propagation/export.cpp"
        # math
        "math/interpolation.cpp"
        # Executable source file
        "jsonInterface.cpp"
        )

# Set the header files.
set(json_interface_HEADERS
        # support
        "support/errorHandling.h"
        "support/keys.h"
        "support/deserialization.h"
        "support/path.h"
        "support/utilities.h"
        "support/valueAccess.h"
        "support/valueConversions.h"
        "support/options.h"
        # environment
        "environment/spice.h"
        "environment/body.h"
        "environment/atmosphere.h"
        "environment/ephemeris.h"
        "environment/aerodynamics.h"
        "environment/gravityField.h"
        "environment/rotationModel.h"
        "environment/shapeModel.h"
        "environment/radiationPressure.h"
        "environment/gravityFieldVariation.h"
        "environment/groundStations.h"
        # propagation
        "propagation/variable.h"
        "propagation/state.h"
        "propagation/propagator.h"
        "propagation/termination.h"
        "propagation/referenceFrames.h"
        "propagation/acceleration.h"
        "propagation/thrust.h"
        "propagation/massRateModel.h"
        "propagation/torque.h"
        "propagation/export.h"
        # math
        "math/integrator.h"
        "math/interpolation.h"
        # Unittests
        "tests/unitTestSupport.h"
        # Executable header file
        "jsonInterface.h"
        )

if (TUDAT_BUILD_WITH_ESTIMATION_TOOLS)
    set(json_interface_SOURCES
            ${json_interface_SOURCES}
            "estimation/parameter.cpp"
            "estimation/observation.cpp"
            #"estimation/orbitDetermination.cpp"
            "jsonInterfaceVariational.cpp"
            #"jsonEstimationInterface.cpp"
            )

    # Set the header files.
    set(json_interface_HEADERS
            ${json_interface_HEADERS}
            "estimation/parameter.h"
            "estimation/observation.h"
            #"estimation/orbitDetermination.h"
            "jsonInterfaceVariational.h"
            #"jsonEstimationInterface.h"
            )
endif ()



if (TUDAT_BUILD_WITH_ESTIMATION_TOOLS)
    list(APPEND json_propagation_LIBRARIES
            json_interface_library
            ${Tudat_ESTIMATION_LIBRARIES}
            ${Boost_LIBRARIES})
else ()
    list(APPEND json_propagation_LIBRARIES
            json_interface_library
            ${Tudat_PROPAGATION_LIBRARIES}
            ${Boost_LIBRARIES})
endif ()

# Add library.
TUDAT_ADD_LIBRARY("json_interface"
        "${json_interface_SOURCES}"
        "${json_interface_HEADERS}"
        json_propagation_LIBRARIES
        )

# Add JSON-based application.
TUDAT_ADD_EXECUTABLE(json_interface
        jsonInterfaceExecute.cpp
        tudat_json_interface
        )

TUDAT_ADD_TESTCASE(Acceleration
        tudat_json_interface)

#set(JSON_TESTS_DIR "${SRCROOT}${JSONINTERFACEDIR}/tests")

# Acceleration
#add_executable(test_JsonInterfaceAcceleration "${JSON_TESTS_DIR}/unitTestAcceleration.cpp")
#setup_custom_test_program(test_JsonInterfaceAcceleration "")
#target_link_libraries(test_JsonInterfaceAcceleration ${JSON_PROPAGATION_LIBRARIES})
#
## aerodynamics
#add_executable(test_JsonInterfaceAerodynamics "${JSON_TESTS_DIR}/unitTestAerodynamics.cpp")
#setup_custom_test_program(test_JsonInterfaceAerodynamics "")
#target_link_libraries(test_JsonInterfaceAerodynamics ${JSON_PROPAGATION_LIBRARIES})
#
## Atmosphere
#add_executable(test_JsonInterfaceAtmosphere "${JSON_TESTS_DIR}/unitTestAtmosphere.cpp")
#setup_custom_test_program(test_JsonInterfaceAtmosphere "")
#target_link_libraries(test_JsonInterfaceAtmosphere ${JSON_PROPAGATION_LIBRARIES})
#
## Body
#add_executable(test_JsonInterfaceBody "${JSON_TESTS_DIR}/unitTestBody.cpp")
#setup_custom_test_program(test_JsonInterfaceBody "")
#target_link_libraries(test_JsonInterfaceBody ${JSON_PROPAGATION_LIBRARIES})
#
## Deserialization
#add_executable(test_JsonInterfaceDeserialization "${JSON_TESTS_DIR}/unitTestDeserialization.cpp")
#setup_custom_test_program(test_JsonInterfaceDeserialization "")
#target_link_libraries(test_JsonInterfaceDeserialization ${JSON_PROPAGATION_LIBRARIES})
#
## Ephemeris
#add_executable(test_JsonInterfaceEphemeris "${JSON_TESTS_DIR}/unitTestEphemeris.cpp")
#setup_custom_test_program(test_JsonInterfaceEphemeris "")
#target_link_libraries(test_JsonInterfaceEphemeris ${JSON_PROPAGATION_LIBRARIES})
#
## Export
#add_executable(test_JsonInterfaceExport "${JSON_TESTS_DIR}/unitTestExport.cpp")
#setup_custom_test_program(test_JsonInterfaceExport "")
#target_link_libraries(test_JsonInterfaceExport ${JSON_PROPAGATION_LIBRARIES})
#
## GravityField
#add_executable(test_JsonInterfaceGravityField "${JSON_TESTS_DIR}/unitTestGravityField.cpp")
#setup_custom_test_program(test_JsonInterfaceGravityField "")
#target_link_libraries(test_JsonInterfaceGravityField ${JSON_PROPAGATION_LIBRARIES})
#
## GravityFieldVariation
#add_executable(test_JsonInterfaceGravityFieldVariation "${JSON_TESTS_DIR}/unitTestGravityFieldVariation.cpp")
#setup_custom_test_program(test_JsonInterfaceGravityFieldVariation "")
#target_link_libraries(test_JsonInterfaceGravityFieldVariation ${JSON_PROPAGATION_LIBRARIES})
#
## GroundStation
#add_executable(test_JsonInterfaceGroundStation "${JSON_TESTS_DIR}/unitTestGroundStation.cpp")
#setup_custom_test_program(test_JsonInterfaceGroundStation "")
#target_link_libraries(test_JsonInterfaceGroundStation ${JSON_PROPAGATION_LIBRARIES})
#
## Integrator
#add_executable(test_JsonInterfaceIntegrator "${JSON_TESTS_DIR}/unitTestIntegrator.cpp")
#setup_custom_test_program(test_JsonInterfaceIntegrator "")
#target_link_libraries(test_JsonInterfaceIntegrator ${JSON_PROPAGATION_LIBRARIES})
#
## Interpolation
#add_executable(test_JsonInterfaceInterpolation "${JSON_TESTS_DIR}/unitTestInterpolation.cpp")
#setup_custom_test_program(test_JsonInterfaceInterpolation "")
#target_link_libraries(test_JsonInterfaceInterpolation ${JSON_PROPAGATION_LIBRARIES})
#
## MassRateModel
#add_executable(test_JsonInterfaceMassRateModel "${JSON_TESTS_DIR}/unitTestMassRateModel.cpp")
#setup_custom_test_program(test_JsonInterfaceMassRateModel "")
#target_link_libraries(test_JsonInterfaceMassRateModel ${JSON_PROPAGATION_LIBRARIES})
#
## Propagator
#add_executable(test_JsonInterfacePropagator "${JSON_TESTS_DIR}/unitTestPropagator.cpp")
#setup_custom_test_program(test_JsonInterfacePropagator "")
#target_link_libraries(test_JsonInterfacePropagator ${JSON_PROPAGATION_LIBRARIES})
#
## RadiationPressure
#add_executable(test_JsonInterfaceRadiationPressure "${JSON_TESTS_DIR}/unitTestRadiationPressure.cpp")
#setup_custom_test_program(test_JsonInterfaceRadiationPressure "")
#target_link_libraries(test_JsonInterfaceRadiationPressure ${JSON_PROPAGATION_LIBRARIES})
#
## reference_frames
#add_executable(test_JsonInterfaceReferenceFrames "${JSON_TESTS_DIR}/unitTestReferenceFrames.cpp")
#setup_custom_test_program(test_JsonInterfaceReferenceFrames "")
#target_link_libraries(test_JsonInterfaceReferenceFrames ${JSON_PROPAGATION_LIBRARIES})
#
## RotationModel
#add_executable(test_JsonInterfaceRotationModel "${JSON_TESTS_DIR}/unitTestRotationModel.cpp")
#setup_custom_test_program(test_JsonInterfaceRotationModel "")
#target_link_libraries(test_JsonInterfaceRotationModel ${JSON_PROPAGATION_LIBRARIES})
#
## ShapeModel
#add_executable(test_JsonInterfaceShapeModel "${JSON_TESTS_DIR}/unitTestShapeModel.cpp")
#setup_custom_test_program(test_JsonInterfaceShapeModel "")
#target_link_libraries(test_JsonInterfaceShapeModel ${JSON_PROPAGATION_LIBRARIES})
#
## SimulationSingleSatellite
#add_executable(test_JsonInterfaceSimulationSingleSatellite "${JSON_TESTS_DIR}/unitTestSimulationSingleSatellite.cpp")
#setup_custom_test_program(test_JsonInterfaceSimulationSingleSatellite "")
#target_link_libraries(test_JsonInterfaceSimulationSingleSatellite ${JSON_PROPAGATION_LIBRARIES})
#
## SimulationSinglePerturbedSatellite
#add_executable(test_JsonInterfaceSimulationSinglePerturbedSatellite "${JSON_TESTS_DIR}/unitTestSimulationSinglePerturbedSatellite.cpp")
#setup_custom_test_program(test_JsonInterfaceSimulationSinglePerturbedSatellite "")
#target_link_libraries(test_JsonInterfaceSimulationSinglePerturbedSatellite ${JSON_PROPAGATION_LIBRARIES})
#
## SimulationInnerSolarSystem
#add_executable(test_JsonInterfaceSimulationInnerSolarSystem "${JSON_TESTS_DIR}/unitTestSimulationInnerSolarSystem.cpp")
#setup_custom_test_program(test_JsonInterfaceSimulationInnerSolarSystem "")
#target_link_libraries(test_JsonInterfaceSimulationInnerSolarSystem ${JSON_PROPAGATION_LIBRARIES})
#
## SimulationGalileoConstellation
#add_executable(test_JsonInterfaceSimulationGalileoConstellation "${JSON_TESTS_DIR}/unitTestSimulationGalileoConstellation.cpp")
#setup_custom_test_program(test_JsonInterfaceSimulationGalileoConstellation "")
#target_link_libraries(test_JsonInterfaceSimulationGalileoConstellation ${JSON_PROPAGATION_LIBRARIES})
#
## SimulationThrustAlongVelocityVector
#add_executable(test_JsonInterfaceSimulationThrustAlongVelocityVector "${JSON_TESTS_DIR}/unitTestSimulationThrustAlongVelocityVector.cpp")
#setup_custom_test_program(test_JsonInterfaceSimulationThrustAlongVelocityVector "")
#target_link_libraries(test_JsonInterfaceSimulationThrustAlongVelocityVector ${JSON_PROPAGATION_LIBRARIES})
#
## SimulationThrustAccelerationFromFile
#add_executable(test_JsonInterfaceSimulationThrustAccelerationFromFile "${JSON_TESTS_DIR}/unitTestSimulationThrustAccelerationFromFile.cpp")
#setup_custom_test_program(test_JsonInterfaceSimulationThrustAccelerationFromFile "")
#target_link_libraries(test_JsonInterfaceSimulationThrustAccelerationFromFile ${JSON_PROPAGATION_LIBRARIES})
#
## Spice
#add_executable(test_JsonInterfaceSpice "${JSON_TESTS_DIR}/unitTestSpice.cpp")
#setup_custom_test_program(test_JsonInterfaceSpice "")
#target_link_libraries(test_JsonInterfaceSpice ${JSON_PROPAGATION_LIBRARIES})
#
## State
#add_executable(test_JsonInterfaceState "${JSON_TESTS_DIR}/unitTestState.cpp")
#setup_custom_test_program(test_JsonInterfaceState "")
#target_link_libraries(test_JsonInterfaceState ${JSON_PROPAGATION_LIBRARIES})
#
## Termination
#add_executable(test_JsonInterfaceTermination "${JSON_TESTS_DIR}/unitTestTermination.cpp")
#setup_custom_test_program(test_JsonInterfaceTermination "")
#target_link_libraries(test_JsonInterfaceTermination ${JSON_PROPAGATION_LIBRARIES})
#
## Thrust
#add_executable(test_JsonInterfaceThrust "${JSON_TESTS_DIR}/unitTestThrust.cpp")
#setup_custom_test_program(test_JsonInterfaceThrust "")
#target_link_libraries(test_JsonInterfaceThrust ${JSON_PROPAGATION_LIBRARIES})
#
## Torque
#add_executable(test_JsonInterfaceTorque "${JSON_TESTS_DIR}/unitTestTorque.cpp")
#setup_custom_test_program(test_JsonInterfaceTorque "")
#target_link_libraries(test_JsonInterfaceTorque ${JSON_PROPAGATION_LIBRARIES})
#
## Variable
#add_executable(test_JsonInterfaceVariable "${JSON_TESTS_DIR}/unitTestVariable.cpp")
#setup_custom_test_program(test_JsonInterfaceVariable "")
#target_link_libraries(test_JsonInterfaceVariable ${JSON_PROPAGATION_LIBRARIES})
#
#if (BUILD_WITH_ESTIMATION_TOOLS)
#    add_executable(test_JsonInterfaceObservation "${JSON_TESTS_DIR}/unitTestObservation.cpp")
#    setup_custom_test_program(test_JsonInterfaceObservation "")
#    target_link_libraries(test_JsonInterfaceObservation ${JSON_PROPAGATION_LIBRARIES})
#
#    add_executable(test_JsonInterfaceParameter "${JSON_TESTS_DIR}/unitTestParameter.cpp")
#    setup_custom_test_program(test_JsonInterfaceParameter "")
#    target_link_libraries(test_JsonInterfaceParameter ${JSON_PROPAGATION_LIBRARIES})
#
#    add_executable(test_JsonInterfaceSimulationSingleSatelliteVariational "${JSON_TESTS_DIR}/unitTestSimulationSingleSatelliteVariational.cpp")
#    setup_custom_test_program(test_JsonInterfaceSimulationSingleSatelliteVariational "")
#    target_link_libraries(test_JsonInterfaceSimulationSingleSatelliteVariational ${JSON_PROPAGATION_LIBRARIES})
#
#    #add_executable(test_JsonInterfaceSimulationSingleSatelliteEstimation "${JSON_TESTS_DIR}/unitTestSimulationSingleSatelliteEstimation.cpp")
#    #setup_custom_test_program(test_JsonInterfaceSimulationSingleSatelliteEstimation "")
#    #target_link_libraries(test_JsonInterfaceSimulationSingleSatelliteEstimation ${JSON_PROPAGATION_LIBRARIES})
#endif ()
#

