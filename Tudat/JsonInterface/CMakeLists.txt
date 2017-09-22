 #    Copyright (c) 2010-2017, Delft University of Technology
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
endmacro(setup_custom_test_program)


# Add JSON-based application.

add_executable(json_interface "${SRCROOT}${JSONINTERFACEDIR}/jsonInterface.cpp")
setup_executable_target(json_interface)
target_link_libraries(json_interface ${JSON_PROPAGATION_LIBRARIES})

set(JSON_TESTS_DIR "${SRCROOT}${JSONINTERFACEDIR}/UnitTests")

# Add unit tests.

# valueAccess
add_executable(test_json_valueAccess "${JSON_TESTS_DIR}/valueAccess.cpp")
setup_custom_test_program(test_json_valueAccess "")
target_link_libraries(test_json_valueAccess ${JSON_PROPAGATION_LIBRARIES})

# # modularFiles
# add_executable(test_json_modularFiles "${JSON_TESTS_DIR}/modularFiles.cpp")
# setup_custom_test_program(test_json_modularFiles "")
# target_link_libraries(test_json_modularFiles ${JSON_PROPAGATION_LIBRARIES})
#
# # mergeableFiles
# add_executable(test_json_mergeableFiles "${JSON_TESTS_DIR}/mergeableFiles.cpp")
# setup_custom_test_program(test_json_mergeableFiles "")
# target_link_libraries(test_json_mergeableFiles ${JSON_PROPAGATION_LIBRARIES})

# aerodynamics
add_executable(test_json_aerodynamics "${JSON_TESTS_DIR}/aerodynamics.cpp")
setup_custom_test_program(test_json_aerodynamics "")
target_link_libraries(test_json_aerodynamics ${JSON_PROPAGATION_LIBRARIES})

# atmosphere
add_executable(test_json_atmosphere "${JSON_TESTS_DIR}/atmosphere.cpp")
setup_custom_test_program(test_json_atmosphere "")
target_link_libraries(test_json_atmosphere ${JSON_PROPAGATION_LIBRARIES})

# ephemeris
add_executable(test_json_ephemeris "${JSON_TESTS_DIR}/ephemeris.cpp")
setup_custom_test_program(test_json_ephemeris "")
target_link_libraries(test_json_ephemeris ${JSON_PROPAGATION_LIBRARIES})

# gravityField
add_executable(test_json_gravityField "${JSON_TESTS_DIR}/gravityField.cpp")
setup_custom_test_program(test_json_gravityField "")
target_link_libraries(test_json_gravityField ${JSON_PROPAGATION_LIBRARIES})

# gravityFieldVariation
add_executable(test_json_gravityFieldVariation "${JSON_TESTS_DIR}/gravityFieldVariation.cpp")
setup_custom_test_program(test_json_gravityFieldVariation "")
target_link_libraries(test_json_gravityFieldVariation ${JSON_PROPAGATION_LIBRARIES})

# radiationPressure
add_executable(test_json_radiationPressure "${JSON_TESTS_DIR}/radiationPressure.cpp")
setup_custom_test_program(test_json_radiationPressure "")
target_link_libraries(test_json_radiationPressure ${JSON_PROPAGATION_LIBRARIES})

# rotationModel
add_executable(test_json_rotationModel "${JSON_TESTS_DIR}/rotationModel.cpp")
setup_custom_test_program(test_json_rotationModel "")
target_link_libraries(test_json_rotationModel ${JSON_PROPAGATION_LIBRARIES})

# shapeModel
add_executable(test_json_shapeModel "${JSON_TESTS_DIR}/shapeModel.cpp")
setup_custom_test_program(test_json_shapeModel "")
target_link_libraries(test_json_shapeModel ${JSON_PROPAGATION_LIBRARIES})

# simulationSingleSatellite
add_executable(test_json_simulationSingleSatellite "${JSON_TESTS_DIR}/simulationSingleSatellite.cpp")
setup_custom_test_program(test_json_simulationSingleSatellite "")
target_link_libraries(test_json_simulationSingleSatellite ${JSON_PROPAGATION_LIBRARIES})

# simulationSinglePerturbedSatellite
add_executable(test_json_simulationSinglePerturbedSatellite "${JSON_TESTS_DIR}/simulationSinglePerturbedSatellite.cpp")
setup_custom_test_program(test_json_simulationSinglePerturbedSatellite "")
target_link_libraries(test_json_simulationSinglePerturbedSatellite ${JSON_PROPAGATION_LIBRARIES})

# simulationInnerSolarSystem
add_executable(test_json_simulationInnerSolarSystem "${JSON_TESTS_DIR}/simulationInnerSolarSystem.cpp")
setup_custom_test_program(test_json_simulationInnerSolarSystem "")
target_link_libraries(test_json_simulationInnerSolarSystem ${JSON_PROPAGATION_LIBRARIES})

# simulationGalileoConstellation
add_executable(test_json_simulationGalileoConstellation "${JSON_TESTS_DIR}/simulationGalileoConstellation.cpp")
setup_custom_test_program(test_json_simulationGalileoConstellation "")
target_link_libraries(test_json_simulationGalileoConstellation ${JSON_PROPAGATION_LIBRARIES})

# simulationThrustAlongVelocityVector
add_executable(test_json_simulationThrustAlongVelocityVector "${JSON_TESTS_DIR}/simulationThrustAlongVelocityVector.cpp")
setup_custom_test_program(test_json_simulationThrustAlongVelocityVector "")
target_link_libraries(test_json_simulationThrustAlongVelocityVector ${JSON_PROPAGATION_LIBRARIES})

# simulationThrustAccelerationFromFile
add_executable(test_json_simulationThrustAccelerationFromFile "${JSON_TESTS_DIR}/simulationThrustAccelerationFromFile.cpp")
setup_custom_test_program(test_json_simulationThrustAccelerationFromFile "")
target_link_libraries(test_json_simulationThrustAccelerationFromFile ${JSON_PROPAGATION_LIBRARIES})






