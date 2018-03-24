#    Copyright (c) 2010-2018, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

set(ESTIMATIONSETUPSDIR "${SIMULATIONSETUPDIR}/EstimationSetup")

# Add source files.
file(GLOB_RECURSE SIMULATION_SETUP_HEADERS ${SRCROOT}/SimulationSetup ABSOLUTE ${CODEROOT} *.h)
file(GLOB_RECURSE SIMULATION_SETUP_SOURCES ${SRCROOT}/SimulationSetup ABSOLUTE ${CODEROOT} *.cpp)


# Add static libraries.
add_library(tudat_simulation_setup STATIC ${SIMULATION_SETUP_SOURCES} ${SIMULATION_SETUP_HEADERS} )
setup_tudat_library_target(tudat_simulation_setup "${SRCROOT}${SIMULATIONSETUPDIR}")

# Add unit tests.
if(USE_CSPICE)
    add_executable(test_EnvironmentCreation "${SRCROOT}${SIMULATIONSETUPDIR}/UnitTests/unitTestEnvironmentModelSetup.cpp")
    setup_custom_test_program(test_EnvironmentCreation "${SRCROOT}${SIMULATIONSETUPDIR}/")
    target_link_libraries(test_EnvironmentCreation ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

    add_executable(test_AccelerationModelCreation "${SRCROOT}${SIMULATIONSETUPDIR}/UnitTests/unitTestAccelerationModelSetup.cpp")
    setup_custom_test_program(test_AccelerationModelCreation "${SRCROOT}${SIMULATIONSETUPDIR}/")
    target_link_libraries(test_AccelerationModelCreation ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
endif()
