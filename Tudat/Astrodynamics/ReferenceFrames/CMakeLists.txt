 #    Copyright (c) 2010-2015, Delft University of Technology
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
 #
 #    References
 #
 #    Notes
 #

# Add source files.
set(REFERENCEFRAMES_SOURCES
  "${SRCROOT}${REFERENCEFRAMESDIR}/referenceFrameTransformations.cpp"
  "${SRCROOT}${REFERENCEFRAMESDIR}/apparentAccelerationModel.cpp"
  "${SRCROOT}${REFERENCEFRAMESDIR}/aerodynamicAngleCalculator.cpp"
)

# Add header files.
set(REFERENCEFRAMES_HEADERS 
  "${SRCROOT}${REFERENCEFRAMESDIR}/referenceFrameTransformations.h"
  "${SRCROOT}${REFERENCEFRAMESDIR}/apparentAccelerationModel.h"
  "${SRCROOT}${REFERENCEFRAMESDIR}/aerodynamicAngleCalculator.h"
)

# Add static libraries.
add_library(tudat_reference_frames STATIC ${REFERENCEFRAMES_SOURCES} ${REFERENCEFRAMES_HEADERS})
setup_tudat_library_target(tudat_reference_frames "${SRCROOT}${REFERENCEFRAMESDIR}")

# Add unit tests.
add_executable(test_ReferenceFrameTransformations "${SRCROOT}${REFERENCEFRAMESDIR}/UnitTests/unitTestReferenceFrameTransformations.cpp")
setup_custom_test_program(test_ReferenceFrameTransformations "${SRCROOT}${REFERENCEFRAMESDIR}")
target_link_libraries(test_ReferenceFrameTransformations tudat_reference_frames ${Boost_LIBRARIES})

add_executable(test_ApparentAccelerationModel "${SRCROOT}${REFERENCEFRAMESDIR}/UnitTests/unitTestApparentAccelerationModel.cpp")
setup_custom_test_program(test_ApparentAccelerationModel "${SRCROOT}${REFERENCEFRAMESDIR}")
target_link_libraries(test_ApparentAccelerationModel tudat_reference_frames ${Boost_LIBRARIES})

add_executable(test_AerodynamicAngleCalculator "${SRCROOT}${REFERENCEFRAMESDIR}/UnitTests/unitTestAerodynamicAngleCalculator.cpp")
setup_custom_test_program(test_AerodynamicAngleCalculator "${SRCROOT}${REFERENCEFRAMESDIR}")
target_link_libraries(test_AerodynamicAngleCalculator tudat_reference_frames tudat_basic_mathematics ${Boost_LIBRARIES})
