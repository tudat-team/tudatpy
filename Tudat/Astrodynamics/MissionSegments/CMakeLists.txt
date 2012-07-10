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
 #      120531    P. Musegaas       Deleted obsolete escapePhase, capturePhase and deepSpaceManeuvre
 #                                  headers. Removed unitTestDeepSpaceManeuvre.
 #

# Set the source files.
set(MISSIONSEGMENTS_SOURCES
  "${SRCROOT}${MISSIONSEGMENTSDIR}/escapeAndCapture.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/gravityAssist.cpp"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertTargeter.cpp"
)

# Set the header files.
set(MISSIONSEGMENTS_HEADERS 
  "${SRCROOT}${MISSIONSEGMENTSDIR}/escapeAndCapture.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/gravityAssist.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/lambertTargeter.h"
  "${SRCROOT}${MISSIONSEGMENTSDIR}/trajectoryDesignMethod.h"
)

# Add static libraries.
add_library(tudat_mission_segments STATIC ${MISSIONSEGMENTS_SOURCES} ${MISSIONSEGMENTS_HEADERS})
setup_tudat_library_target(tudat_mission_segments "${SRCROOT}${MISSIONSEGMENTSDIR}")

# Add unit tests.
add_executable(test_EscapeAndCapture "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests/unitTestEscapeAndCapture.cpp")
setup_custom_test_program(test_EscapeAndCapture "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_EscapeAndCapture tudat_bodies tudat_mission_segments tudat_gravitation tudat_basic_astrodynamics tudat_root_finding_methods tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_GravityAssist "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests/unitTestGravityAssist.cpp")
setup_custom_test_program(test_GravityAssist "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_GravityAssist tudat_bodies tudat_gravitation tudat_mission_segments tudat_basic_astrodynamics tudat_root_finding_methods tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_LambertTargeter "${SRCROOT}${MISSIONSEGMENTSDIR}/UnitTests/unitTestLambertTargeter.cpp")
setup_custom_test_program(test_LambertTargeter "${SRCROOT}${MISSIONSEGMENTSDIR}")
target_link_libraries(test_LambertTargeter tudat_bodies tudat_gravitation tudat_mission_segments tudat_root_finding_methods tudat_basic_astrodynamics tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
