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
 #      111123    B. Tong Minh      Added unit test for State.
 #

# Add source files.
set(STATES_SOURCES
  "${SRCROOT}${STATESDIR}/writeStateHistoryToFile.cpp"
)

# Add header files.
set(STATES_HEADERS 
  "${SRCROOT}${STATESDIR}/cartesianElements.h"
  "${SRCROOT}${STATESDIR}/cartesianPositionElements.h"
  "${SRCROOT}${STATESDIR}/cartesianVelocityElements.h"
  "${SRCROOT}${STATESDIR}/keplerianElements.h"
  "${SRCROOT}${STATESDIR}/state.h"
  "${SRCROOT}${STATESDIR}/writeStateHistoryToFile.h"
)

# Add static libraries.
add_library(tudat_states STATIC ${STATES_SOURCES} ${STATES_HEADERS})
setup_tudat_library_target(tudat_states "${SRCROOT}${STATESDIR}")

# Add unit tests.
add_executable(test_CartesianElements "${SRCROOT}${STATESDIR}/UnitTests/unitTestCartesianElements.cpp")
setup_custom_test_program(test_CartesianElements "${SRCROOT}${STATESDIR}")
target_link_libraries(test_CartesianElements tudat_states ${Boost_LIBRARIES})

add_executable(test_KeplerianElements "${SRCROOT}${STATESDIR}/UnitTests/unitTestKeplerianElements.cpp")
setup_custom_test_program(test_KeplerianElements "${SRCROOT}${STATESDIR}")
target_link_libraries(test_KeplerianElements tudat_states ${Boost_LIBRARIES})

add_executable(test_State "${SRCROOT}${STATESDIR}/UnitTests/unitTestState.cpp")
setup_custom_test_program(test_State "${SRCROOT}${STATESDIR}")
target_link_libraries(test_State tudat_states ${Boost_LIBRARIES})
