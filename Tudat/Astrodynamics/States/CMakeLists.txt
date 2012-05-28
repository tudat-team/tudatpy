 #    Copyright (c) 2010-2012, Delft University of Technology
 #    All rights reserved.
 #
 #    Redistribution and use in source and binary forms, with or without modification, are
 #    permitted provided that the following conditions are met:
 #      - Redistributions of source code must retain the above copyright notice, this list of
 #        conditions and the following disclaimer.
 #      - Redistributions in binary form must reproduce the above copyright notice, this list of
 #        conditions and the following disclaimer in the documentation and/or other materials
 #        provided with the distribution.
 #      - Neither the name of the Delft University of Technology nor the names of its contributors
 #        may be used to endorse or promote products derived from this software without specific
 #        prior written permission.
 #
 #    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 #    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 #    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 #    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 #    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 #    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 #    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 #    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 #    OF THE POSSIBILITY OF SUCH DAMAGE.
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
