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
 #      120606    T. Secretin       Adapted to new Ephemeris folder.
 #
 #    References
 #
 #    Notes
 #

# Set the source files.
set(EPHEMERIDES_SOURCES
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositionsBase.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositions.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositionsCircularCoplanar.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/cartesianStateExtractor.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/keplerStateExtractor.cpp"
)

# Set the header files.
set(EPHEMERIDES_HEADERS
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositionsBase.h"
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositions.h"
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositionsCircularCoplanar.h"
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositionsDataContainer.h"
  "${SRCROOT}${EPHEMERIDESDIR}/ephemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/cartesianStateExtractor.h"
  "${SRCROOT}${EPHEMERIDESDIR}/keplerStateExtractor.h"
)

# Add static libraries.
add_library(tudat_ephemerides STATIC ${EPHEMERIDES_SOURCES} ${EPHEMERIDES_HEADERS})
setup_tudat_library_target(tudat_ephemerides "${SRCROOT}${EPHEMERIDESDIR}")

# Add unit tests.
add_executable(test_ApproximatePlanetPositions "${SRCROOT}${EPHEMERIDESDIR}/UnitTests/unitTestApproximatePlanetPositions.cpp")
setup_custom_test_program(test_ApproximatePlanetPositions "${SRCROOT}${EPHEMERIDESDIR}")
target_link_libraries(test_ApproximatePlanetPositions tudat_ephemerides tudat_gravitation tudat_basic_astrodynamics tudat_input_output tudat_root_finders  ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_CartesianStateExtractor "${SRCROOT}${EPHEMERIDESDIR}/UnitTests/unitTestCartesianStateExtractor.cpp")
setup_custom_test_program(test_CartesianStateExtractor "${SRCROOT}${EPHEMERIDESDIR}")
target_link_libraries(test_CartesianStateExtractor tudat_input_output tudat_ephemerides ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_KeplerStateExtractor "${SRCROOT}${EPHEMERIDESDIR}/UnitTests/unitTestKeplerStateExtractor.cpp")
setup_custom_test_program(test_KeplerStateExtractor "${SRCROOT}${EPHEMERIDESDIR}")
target_link_libraries(test_KeplerStateExtractor tudat_ephemerides tudat_input_output tudat_basic_astrodynamics tudat_root_finders ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
