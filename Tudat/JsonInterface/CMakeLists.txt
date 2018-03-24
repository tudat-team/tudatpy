 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Set the source files.
set(JSONINTERFACE_SOURCES
    "${SRCROOT}${JSONINTERFACEDIR}/Support/keys.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/deserialization.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/path.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/valueAccess.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/spice.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/body.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/atmosphere.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/ephemeris.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/aerodynamics.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/gravityField.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/rotationModel.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/shapeModel.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/radiationPressure.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/gravityFieldVariation.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/variable.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/termination.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/acceleration.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/massRateModel.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/thrust.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/torque.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Mathematics/interpolation.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/export.cpp"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/options.cpp"
)

# Set the header files.
set(JSONINTERFACE_HEADERS
    "${SRCROOT}${JSONINTERFACEDIR}/Support/errorHandling.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/keys.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/deserialization.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/path.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/utilities.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/valueAccess.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/valueConversions.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/spice.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/body.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/atmosphere.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/ephemeris.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/aerodynamics.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/gravityField.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/rotationModel.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/shapeModel.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/radiationPressure.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Environment/gravityFieldVariation.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/variable.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/state.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/propagator.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/termination.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/referenceFrames.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/acceleration.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/thrust.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/massRateModel.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/torque.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Mathematics/integrator.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Mathematics/interpolation.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Propagation/export.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/options.h"
    "${SRCROOT}${JSONINTERFACEDIR}/UnitTests/unitTestSupport.h"
    "${SRCROOT}${JSONINTERFACEDIR}/jsonInterface.h"
)

# Add static libraries.
add_library(json_interface_library STATIC ${JSONINTERFACE_SOURCES} ${JSONINTERFACE_HEADERS})
setup_tudat_library_target(json_interface_library "")

list(APPEND JSON_PROPAGATION_LIBRARIES json_interface_library ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

# Set the global macros for setting up targets.
macro(setup_executable_target target_name)
  set_property(TARGET ${target_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY "${BINROOT}")
  install(TARGETS ${target_name} RUNTIME DESTINATION "${BINROOT}")
endmacro()


# Add JSON-based application.

add_executable(json_interface "${SRCROOT}${JSONINTERFACEDIR}/jsonInterface.cpp")
setup_executable_target(json_interface)
target_link_libraries(json_interface ${JSON_PROPAGATION_LIBRARIES})


# Add unit tests.

set(JSON_TESTS_DIR "${SRCROOT}${JSONINTERFACEDIR}/UnitTests")

# Acceleration
add_executable(test_JsonInterfaceAcceleration "${JSON_TESTS_DIR}/unitTestAcceleration.cpp")
setup_custom_test_program(test_JsonInterfaceAcceleration "")
target_link_libraries(test_JsonInterfaceAcceleration ${JSON_PROPAGATION_LIBRARIES})

# Aerodynamics
add_executable(test_JsonInterfaceAerodynamics "${JSON_TESTS_DIR}/unitTestAerodynamics.cpp")
setup_custom_test_program(test_JsonInterfaceAerodynamics "")
target_link_libraries(test_JsonInterfaceAerodynamics ${JSON_PROPAGATION_LIBRARIES})

# Atmosphere
add_executable(test_JsonInterfaceAtmosphere "${JSON_TESTS_DIR}/unitTestAtmosphere.cpp")
setup_custom_test_program(test_JsonInterfaceAtmosphere "")
target_link_libraries(test_JsonInterfaceAtmosphere ${JSON_PROPAGATION_LIBRARIES})

# Body
add_executable(test_JsonInterfaceBody "${JSON_TESTS_DIR}/unitTestBody.cpp")
setup_custom_test_program(test_JsonInterfaceBody "")
target_link_libraries(test_JsonInterfaceBody ${JSON_PROPAGATION_LIBRARIES})

# Deserialization
add_executable(test_JsonInterfaceDeserialization "${JSON_TESTS_DIR}/unitTestDeserialization.cpp")
setup_custom_test_program(test_JsonInterfaceDeserialization "")
target_link_libraries(test_JsonInterfaceDeserialization ${JSON_PROPAGATION_LIBRARIES})

# Ephemeris
add_executable(test_JsonInterfaceEphemeris "${JSON_TESTS_DIR}/unitTestEphemeris.cpp")
setup_custom_test_program(test_JsonInterfaceEphemeris "")
target_link_libraries(test_JsonInterfaceEphemeris ${JSON_PROPAGATION_LIBRARIES})

# Export
add_executable(test_JsonInterfaceExport "${JSON_TESTS_DIR}/unitTestExport.cpp")
setup_custom_test_program(test_JsonInterfaceExport "")
target_link_libraries(test_JsonInterfaceExport ${JSON_PROPAGATION_LIBRARIES})

# GravityField
add_executable(test_JsonInterfaceGravityField "${JSON_TESTS_DIR}/unitTestGravityField.cpp")
setup_custom_test_program(test_JsonInterfaceGravityField "")
target_link_libraries(test_JsonInterfaceGravityField ${JSON_PROPAGATION_LIBRARIES})

# GravityFieldVariation
add_executable(test_JsonInterfaceGravityFieldVariation "${JSON_TESTS_DIR}/unitTestGravityFieldVariation.cpp")
setup_custom_test_program(test_JsonInterfaceGravityFieldVariation "")
target_link_libraries(test_JsonInterfaceGravityFieldVariation ${JSON_PROPAGATION_LIBRARIES})

# Integrator
add_executable(test_JsonInterfaceIntegrator "${JSON_TESTS_DIR}/unitTestIntegrator.cpp")
setup_custom_test_program(test_JsonInterfaceIntegrator "")
target_link_libraries(test_JsonInterfaceIntegrator ${JSON_PROPAGATION_LIBRARIES})

# Interpolation
add_executable(test_JsonInterfaceInterpolation "${JSON_TESTS_DIR}/unitTestInterpolation.cpp")
setup_custom_test_program(test_JsonInterfaceInterpolation "")
target_link_libraries(test_JsonInterfaceInterpolation ${JSON_PROPAGATION_LIBRARIES})

# MassRateModel
add_executable(test_JsonInterfaceMassRateModel "${JSON_TESTS_DIR}/unitTestMassRateModel.cpp")
setup_custom_test_program(test_JsonInterfaceMassRateModel "")
target_link_libraries(test_JsonInterfaceMassRateModel ${JSON_PROPAGATION_LIBRARIES})

# Propagator
add_executable(test_JsonInterfacePropagator "${JSON_TESTS_DIR}/unitTestPropagator.cpp")
setup_custom_test_program(test_JsonInterfacePropagator "")
target_link_libraries(test_JsonInterfacePropagator ${JSON_PROPAGATION_LIBRARIES})

# RadiationPressure
add_executable(test_JsonInterfaceRadiationPressure "${JSON_TESTS_DIR}/unitTestRadiationPressure.cpp")
setup_custom_test_program(test_JsonInterfaceRadiationPressure "")
target_link_libraries(test_JsonInterfaceRadiationPressure ${JSON_PROPAGATION_LIBRARIES})

# ReferenceFrames
add_executable(test_JsonInterfaceReferenceFrames "${JSON_TESTS_DIR}/unitTestReferenceFrames.cpp")
setup_custom_test_program(test_JsonInterfaceReferenceFrames "")
target_link_libraries(test_JsonInterfaceReferenceFrames ${JSON_PROPAGATION_LIBRARIES})

# RotationModel
add_executable(test_JsonInterfaceRotationModel "${JSON_TESTS_DIR}/unitTestRotationModel.cpp")
setup_custom_test_program(test_JsonInterfaceRotationModel "")
target_link_libraries(test_JsonInterfaceRotationModel ${JSON_PROPAGATION_LIBRARIES})

# ShapeModel
add_executable(test_JsonInterfaceShapeModel "${JSON_TESTS_DIR}/unitTestShapeModel.cpp")
setup_custom_test_program(test_JsonInterfaceShapeModel "")
target_link_libraries(test_JsonInterfaceShapeModel ${JSON_PROPAGATION_LIBRARIES})

# SimulationSingleSatellite
add_executable(test_JsonInterfaceSimulationSingleSatellite "${JSON_TESTS_DIR}/unitTestSimulationSingleSatellite.cpp")
setup_custom_test_program(test_JsonInterfaceSimulationSingleSatellite "")
target_link_libraries(test_JsonInterfaceSimulationSingleSatellite ${JSON_PROPAGATION_LIBRARIES})

# SimulationSinglePerturbedSatellite
add_executable(test_JsonInterfaceSimulationSinglePerturbedSatellite "${JSON_TESTS_DIR}/unitTestSimulationSinglePerturbedSatellite.cpp")
setup_custom_test_program(test_JsonInterfaceSimulationSinglePerturbedSatellite "")
target_link_libraries(test_JsonInterfaceSimulationSinglePerturbedSatellite ${JSON_PROPAGATION_LIBRARIES})

# SimulationInnerSolarSystem
add_executable(test_JsonInterfaceSimulationInnerSolarSystem "${JSON_TESTS_DIR}/unitTestSimulationInnerSolarSystem.cpp")
setup_custom_test_program(test_JsonInterfaceSimulationInnerSolarSystem "")
target_link_libraries(test_JsonInterfaceSimulationInnerSolarSystem ${JSON_PROPAGATION_LIBRARIES})

# SimulationGalileoConstellation
add_executable(test_JsonInterfaceSimulationGalileoConstellation "${JSON_TESTS_DIR}/unitTestSimulationGalileoConstellation.cpp")
setup_custom_test_program(test_JsonInterfaceSimulationGalileoConstellation "")
target_link_libraries(test_JsonInterfaceSimulationGalileoConstellation ${JSON_PROPAGATION_LIBRARIES})

# SimulationThrustAlongVelocityVector
add_executable(test_JsonInterfaceSimulationThrustAlongVelocityVector "${JSON_TESTS_DIR}/unitTestSimulationThrustAlongVelocityVector.cpp")
setup_custom_test_program(test_JsonInterfaceSimulationThrustAlongVelocityVector "")
target_link_libraries(test_JsonInterfaceSimulationThrustAlongVelocityVector ${JSON_PROPAGATION_LIBRARIES})

# SimulationThrustAccelerationFromFile
add_executable(test_JsonInterfaceSimulationThrustAccelerationFromFile "${JSON_TESTS_DIR}/unitTestSimulationThrustAccelerationFromFile.cpp")
setup_custom_test_program(test_JsonInterfaceSimulationThrustAccelerationFromFile "")
target_link_libraries(test_JsonInterfaceSimulationThrustAccelerationFromFile ${JSON_PROPAGATION_LIBRARIES})

# Spice
add_executable(test_JsonInterfaceSpice "${JSON_TESTS_DIR}/unitTestSpice.cpp")
setup_custom_test_program(test_JsonInterfaceSpice "")
target_link_libraries(test_JsonInterfaceSpice ${JSON_PROPAGATION_LIBRARIES})

# State
add_executable(test_JsonInterfaceState "${JSON_TESTS_DIR}/unitTestState.cpp")
setup_custom_test_program(test_JsonInterfaceState "")
target_link_libraries(test_JsonInterfaceState ${JSON_PROPAGATION_LIBRARIES})

# Termination
add_executable(test_JsonInterfaceTermination "${JSON_TESTS_DIR}/unitTestTermination.cpp")
setup_custom_test_program(test_JsonInterfaceTermination "")
target_link_libraries(test_JsonInterfaceTermination ${JSON_PROPAGATION_LIBRARIES})

# Thrust
add_executable(test_JsonInterfaceThrust "${JSON_TESTS_DIR}/unitTestThrust.cpp")
setup_custom_test_program(test_JsonInterfaceThrust "")
target_link_libraries(test_JsonInterfaceThrust ${JSON_PROPAGATION_LIBRARIES})

# Torque
add_executable(test_JsonInterfaceTorque "${JSON_TESTS_DIR}/unitTestTorque.cpp")
setup_custom_test_program(test_JsonInterfaceTorque "")
target_link_libraries(test_JsonInterfaceTorque ${JSON_PROPAGATION_LIBRARIES})

# Variable
add_executable(test_JsonInterfaceVariable "${JSON_TESTS_DIR}/unitTestVariable.cpp")
setup_custom_test_program(test_JsonInterfaceVariable "")
target_link_libraries(test_JsonInterfaceVariable ${JSON_PROPAGATION_LIBRARIES})



