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
    "${SRCROOT}${JSONINTERFACEDIR}/Support/modular.cpp"
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
    "${SRCROOT}${JSONINTERFACEDIR}/Support/modular.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/path.h"
    "${SRCROOT}${JSONINTERFACEDIR}/Support/units.h"
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
    "${SRCROOT}${JSONINTERFACEDIR}/simulation.h"
)

# Add static libraries.
add_library(tudat_json_interface STATIC ${JSONINTERFACE_SOURCES} ${JSONINTERFACE_HEADERS})
setup_tudat_library_target(tudat_json_interface "${SRCROOT}${JSONINTERFACEDIR}")


# Add unit tests.

# aerodynamics
add_executable(test_json_aerodynamics "${SRCROOT}${JSONINTERFACEDIR}/UnitTests/aerodynamics.cpp")
setup_custom_test_program(test_json_aerodynamics "${SRCROOT}${JSONINTERFACEDIR}")
target_link_libraries(test_json_aerodynamics tudat_json_interface ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

# atmosphere
add_executable(test_json_atmosphere "${SRCROOT}${JSONINTERFACEDIR}/UnitTests/atmosphere.cpp")
setup_custom_test_program(test_json_atmosphere "${SRCROOT}${JSONINTERFACEDIR}")
target_link_libraries(test_json_atmosphere tudat_json_interface ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

# ephemeris
add_executable(test_json_ephemeris "${SRCROOT}${JSONINTERFACEDIR}/UnitTests/ephemeris.cpp")
setup_custom_test_program(test_json_ephemeris "${SRCROOT}${JSONINTERFACEDIR}")
target_link_libraries(test_json_ephemeris tudat_json_interface ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

# gravityField
add_executable(test_json_gravityField "${SRCROOT}${JSONINTERFACEDIR}/UnitTests/gravityField.cpp")
setup_custom_test_program(test_json_gravityField "${SRCROOT}${JSONINTERFACEDIR}")
target_link_libraries(test_json_gravityField tudat_json_interface ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

# gravityFieldVariation
add_executable(test_json_gravityFieldVariation "${SRCROOT}${JSONINTERFACEDIR}/UnitTests/gravityFieldVariation.cpp")
setup_custom_test_program(test_json_gravityFieldVariation "${SRCROOT}${JSONINTERFACEDIR}")
target_link_libraries(test_json_gravityFieldVariation tudat_json_interface ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

# radiationPressure
add_executable(test_json_radiationPressure "${SRCROOT}${JSONINTERFACEDIR}/UnitTests/radiationPressure.cpp")
setup_custom_test_program(test_json_radiationPressure "${SRCROOT}${JSONINTERFACEDIR}")
target_link_libraries(test_json_radiationPressure tudat_json_interface ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

# rotationModel
add_executable(test_json_rotationModel "${SRCROOT}${JSONINTERFACEDIR}/UnitTests/rotationModel.cpp")
setup_custom_test_program(test_json_rotationModel "${SRCROOT}${JSONINTERFACEDIR}")
target_link_libraries(test_json_rotationModel tudat_json_interface ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

# shapeModel
add_executable(test_json_shapeModel "${SRCROOT}${JSONINTERFACEDIR}/UnitTests/shapeModel.cpp")
setup_custom_test_program(test_json_shapeModel "${SRCROOT}${JSONINTERFACEDIR}")
target_link_libraries(test_json_shapeModel tudat_json_interface ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

# simulationSingleSatellite
add_executable(test_json_simulationSingleSatellite "${SRCROOT}${JSONINTERFACEDIR}/UnitTests/simulationSingleSatellite.cpp")
setup_custom_test_program(test_json_simulationSingleSatellite "${SRCROOT}${JSONINTERFACEDIR}")
target_link_libraries(test_json_simulationSingleSatellite tudat_json_interface ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

# simulationSinglePerturbedSatellite
add_executable(test_json_simulationSinglePerturbedSatellite "${SRCROOT}${JSONINTERFACEDIR}/UnitTests/simulationSinglePerturbedSatellite.cpp")
setup_custom_test_program(test_json_simulationSinglePerturbedSatellite "${SRCROOT}${JSONINTERFACEDIR}")
target_link_libraries(test_json_simulationSinglePerturbedSatellite tudat_json_interface ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

# simulationInnerSolarSystem
add_executable(test_json_simulationInnerSolarSystem "${SRCROOT}${JSONINTERFACEDIR}/UnitTests/simulationInnerSolarSystem.cpp")
setup_custom_test_program(test_json_simulationInnerSolarSystem "${SRCROOT}${JSONINTERFACEDIR}")
target_link_libraries(test_json_simulationInnerSolarSystem tudat_json_interface ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

# simulationGalileoConstellation
add_executable(test_json_simulationGalileoConstellation "${SRCROOT}${JSONINTERFACEDIR}/UnitTests/simulationGalileoConstellation.cpp")
setup_custom_test_program(test_json_simulationGalileoConstellation "${SRCROOT}${JSONINTERFACEDIR}")
target_link_libraries(test_json_simulationGalileoConstellation tudat_json_interface ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

# simulationThrustAlongVelocityVector
add_executable(test_json_simulationThrustAlongVelocityVector "${SRCROOT}${JSONINTERFACEDIR}/UnitTests/simulationThrustAlongVelocityVector.cpp")
setup_custom_test_program(test_json_simulationThrustAlongVelocityVector "${SRCROOT}${JSONINTERFACEDIR}")
target_link_libraries(test_json_simulationThrustAlongVelocityVector tudat_json_interface ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

# # simulationThrustAccelerationFromFile
# add_executable(test_json_simulationThrustAccelerationFromFile "${SRCROOT}${JSONINTERFACEDIR}/UnitTests/simulationThrustAccelerationFromFile.cpp")
# setup_custom_test_program(test_json_simulationThrustAccelerationFromFile "${SRCROOT}${JSONINTERFACEDIR}")
# target_link_libraries(test_json_simulationThrustAccelerationFromFile tudat_json_interface ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})






