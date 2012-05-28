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
 #      120209    B. Tong Minh      Merge Ephemeris package
 #

# Set Ephemeris path.
set(EPHEMERISDIR "${BODIESDIR}/Ephemeris")

# Set the source files.
set(BODIES_SOURCES
  "${SRCROOT}${BODIESDIR}/planet.cpp"
  "${SRCROOT}${EPHEMERISDIR}/approximatePlanetPositionsBase.cpp"
  "${SRCROOT}${EPHEMERISDIR}/approximatePlanetPositions.cpp"
  "${SRCROOT}${EPHEMERISDIR}/approximatePlanetPositionsCircularCoplanar.cpp"
)

# Set the header files.
set(BODIES_HEADERS 
  "${SRCROOT}${BODIESDIR}/body.h"
  "${SRCROOT}${BODIESDIR}/celestialBody.h"
  "${SRCROOT}${BODIESDIR}/planet.h"
  "${SRCROOT}${BODIESDIR}/vehicle.h"
  "${SRCROOT}${BODIESDIR}/vehicleExternalModel.h"
  "${SRCROOT}${EPHEMERISDIR}/approximatePlanetPositionsBase.h"
  "${SRCROOT}${EPHEMERISDIR}/approximatePlanetPositions.h"
  "${SRCROOT}${EPHEMERISDIR}/approximatePlanetPositionsCircularCoplanar.h"
  "${SRCROOT}${EPHEMERISDIR}/approximatePlanetPositionsDataContainer.h"
  "${SRCROOT}${EPHEMERISDIR}/ephemeris.h"
)

# Add static libraries.
add_library(tudat_bodies STATIC ${BODIES_SOURCES} ${BODIES_HEADERS})
setup_tudat_library_target(tudat_bodies "${SRCROOT}${BODIESDIR}")

# Add unit tests.
add_executable(test_ApproximatePlanetPositions "${SRCROOT}${EPHEMERISDIR}/UnitTests/unitTestApproximatePlanetPositions.cpp")
setup_custom_test_program(test_ApproximatePlanetPositions "${SRCROOT}${EPHEMERISDIR}")
target_link_libraries(test_ApproximatePlanetPositions tudat_states tudat_bodies tudat_gravitation tudat_basic_astrodynamics tudat_input_output tudat_root_finding_methods  ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
