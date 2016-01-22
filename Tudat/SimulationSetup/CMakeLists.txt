 #    Copyright (c) 2010-2012 Delft University of Technology.
 #
 #    This software is protected by national and international copyright.
 #    Any unauthorized use, reproduction or modification is unlawful and
 #    will be prosecuted. Commercial and non-private application of the
 #    software in any form is strictly prohibited unless otherwise granted
 #    by the authors.
 #
 #    The code is provided without any warranty; without even the implied
 #    warranty of merchantibility or fitness for a particular purpose.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      110820    S.M. Persson      File created.
 #      111025    K. Kumar          Adapted file to work with Revision 194.
 #      111026    K. Kumar          Adapted file so all headers show in project tree in Qt Creator.
 #

# Add source files.
set(SIMULATION_SETUP_SOURCES
  "${SRCROOT}${SIMULATIONSETUPDIR}/accelerationModelTypes.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createAccelerationModels.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createAtmosphereModel.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createEphemeris.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createGravityField.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createRotationModel.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createBodies.cpp"
  "${SRCROOT}${SIMULATIONSETUPDIR}/defaultBodies.cpp"
)

# Add header files.
set(SIMULATION_SETUP_HEADERS 
  "${SRCROOT}${SIMULATIONSETUPDIR}/accelerationModelTypes.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createAccelerationModels.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createAtmosphereModel.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createEphemeris.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createGravityField.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createRotationModel.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/createBodies.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/defaultBodies.h"
  "${SRCROOT}${SIMULATIONSETUPDIR}/body.h"
)

# Add static libraries.
add_library(tudat_simulation_setup STATIC ${SIMULATION_SETUP_SOURCES} ${SIMULATION_SETUP_HEADERS} )
setup_tudat_library_target(tudat_simulation_setup "${SRCROOT}${SIMULATIONSETUPDIR}")

# Add unit tests.
add_executable(test_EnvironmentSetup "${SRCROOT}${SIMULATIONSETUPDIR}/UnitTests/unitTestEnvironmentModelSetup.cpp")
setup_custom_test_program(test_EnvironmentSetup "${SRCROOT}${SIMULATIONSETUPDIR}/")
target_link_libraries(test_EnvironmentSetup tudat_simulation_setup tudat_gravitation tudat_ephemerides tudat_reference_frames tudat_aerodynamics tudat_spice_interface tudat_input_output tudat_basic_astrodynamics tudat_basic_mathematics cspice ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_AccelerationModelSetup "${SRCROOT}${SIMULATIONSETUPDIR}/UnitTests/unitTestAccelerationModelSetup.cpp")
setup_custom_test_program(test_AccelerationModelSetup "${SRCROOT}${SIMULATIONSETUPDIR}/")
target_link_libraries(test_AccelerationModelSetup tudat_simulation_setup tudat_gravitation tudat_ephemerides tudat_reference_frames tudat_aerodynamics tudat_spice_interface tudat_input_output tudat_basic_astrodynamics tudat_basic_mathematics cspice ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
