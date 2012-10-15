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
 #	120203    K. Kumar	    Adapted for new BasicAstrodynamics package.
 #      120216    K. Kumar          Removed propagator files.
 #
 #    References
 #
 #    Notes
 #

# Set the source files.
set(BASICASTRODYNAMICS_SOURCES
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/convertMeanAnomalyToEccentricAnomaly.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/convertMeanAnomalyToHyperbolicEccentricAnomaly.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/clohessyWiltshirePropagator.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/keplerPropagator.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/missionGeometry.cpp"
)

# Set the header files.
set(BASICASTRODYNAMICS_HEADERS 
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/accelerationModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/convertMeanAnomalyToEccentricAnomaly.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/convertMeanAnomalyToHyperbolicEccentricAnomaly.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/clohessyWiltshirePropagator.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/forceModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/generalizedForceModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/keplerPropagator.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/missionGeometry.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/momentDueToForceModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/pureMomentModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/testAccelerationModels.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/testBody.h"
)

# Add static libraries.
add_library(tudat_basic_astrodynamics STATIC ${BASICASTRODYNAMICS_SOURCES} ${BASICASTRODYNAMICS_HEADERS})
setup_tudat_library_target(tudat_basic_astrodynamics "${SRCROOT}${BASICASTRODYNAMICSDIR}")

# Add unit tests.
add_executable(test_MeanToEccenticAnomalyConversion "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestConvertMeanAnomalyToEccentricAnomaly.cpp")
setup_custom_test_program(test_MeanToEccenticAnomalyConversion "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_MeanToEccenticAnomalyConversion tudat_basic_astrodynamics tudat_root_finding_methods ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_KeplerPropagator "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestKeplerPropagator.cpp")
setup_custom_test_program(test_KeplerPropagator "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_KeplerPropagator tudat_input_output tudat_bodies tudat_gravitation tudat_basic_astrodynamics tudat_states tudat_root_finding_methods ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_MomentDueToForce "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestMomentDueToForce.cpp")
setup_custom_test_program(test_MomentDueToForce "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_MomentDueToForce tudat_basic_astrodynamics ${Boost_LIBRARIES})

add_executable(test_AccelerationModel "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestAccelerationModel.cpp")
setup_custom_test_program(test_AccelerationModel "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_AccelerationModel tudat_basic_astrodynamics ${Boost_LIBRARIES})

add_executable(test_ClohessyWiltshirePropagator "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestClohessyWiltshirePropagator.cpp")
setup_custom_test_program(test_ClohessyWiltshirePropagator "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_ClohessyWiltshirePropagator tudat_basic_astrodynamics ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_MissionGeometry "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestMissionGeometry.cpp")
setup_custom_test_program(test_MissionGeometry "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_MissionGeometry tudat_basic_astrodynamics ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
