#    Copyright (c) 2010-2016, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#


# Add source files.
set(SIMULATION_SETUP_SOURCES
  "${SRCROOT}${SIMULATIONSETUPDIR}/createAccelerationModels.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createAtmosphereModel.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createEphemeris.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createFlightConditions.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createGravityFieldVariations.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createGravityField.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createRotationModel.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createBodies.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createBodyShapeModel.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createRadiationPressureInterface.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/defaultBodies.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/body.cpp"
)

# Add header files.
set(SIMULATION_SETUP_HEADERS 
  "${SRCROOT}${SIMULATIONSETUPDIR}/accelerationSettings.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createAccelerationModels.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createAtmosphereModel.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createEphemeris.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createFlightConditions.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createGravityFieldVariations.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createGravityField.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createRotationModel.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createBodies.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createBodyShapeModel.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createRadiationPressureInterface.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/defaultBodies.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/body.h"
)

# Add static libraries.
add_library(tudat_simulation_setup STATIC ${SIMULATION_SETUP_SOURCES} ${SIMULATION_SETUP_HEADERS} )
setup_tudat_library_target(tudat_simulation_setup "${SRCROOT}${SIMULATIONSETUPDIR}")

list(APPEND TUDAT_CSPICE_LIBRARIES "")
if(USE_CSPICE)
  list(APPEND TUDAT_CSPICE_LIBRARIES tudat_spice_interface cspice)
endif()

# Add unit tests.
add_executable(test_EnvironmentSetup "${SRCROOT}${SIMULATIONSETUPDIR}/UnitTests/unitTestEnvironmentModelSetup.cpp")
setup_custom_test_program(test_EnvironmentSetup "${SRCROOT}${SIMULATIONSETUPDIR}/")
target_link_libraries(test_EnvironmentSetup tudat_simulation_setup tudat_gravitation tudat_ephemerides
                      tudat_reference_frames tudat_electro_magnetism tudat_aerodynamics tudat_geometric_shapes
                      tudat_input_output tudat_basic_astrodynamics
                      tudat_basic_mathematics ${TUDAT_CSPICE_LIBRARIES} ${Boost_LIBRARIES})

if(USE_CSPICE)
add_executable(test_AccelerationModelSetup "${SRCROOT}${SIMULATIONSETUPDIR}/UnitTests/unitTestAccelerationModelSetup.cpp")
setup_custom_test_program(test_AccelerationModelSetup "${SRCROOT}${SIMULATIONSETUPDIR}/")
target_link_libraries(test_AccelerationModelSetup tudat_simulation_setup tudat_gravitation tudat_ephemerides
                      tudat_reference_frames tudat_electro_magnetism tudat_aerodynamics tudat_geometric_shapes
                      tudat_input_output tudat_basic_astrodynamics
                      tudat_basic_mathematics ${TUDAT_CSPICE_LIBRARIES} ${Boost_LIBRARIES})
endif()
