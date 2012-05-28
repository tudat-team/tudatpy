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
 #	120202    K. Kumar	    Adapted for new Interpolators sub-directory.
 #

# Add source files.
set(INTERPOLATORS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/cubicSplineInterpolation.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/linearInterpolation.cpp"
)

# Add header files.
set(INTERPOLATORS_HEADERS 
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/cubicSplineInterpolation.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/linearInterpolation.h"
)

# Add static libraries.
add_library(tudat_interpolators STATIC ${INTERPOLATORS_SOURCES})
setup_tudat_library_target(tudat_interpolators "${SRCROOT}${MATHEMATICSDIR}/Interpolators")

# Add unit tests.
add_executable(test_CubicSplineInterpolation "${SRCROOT}${MATHEMATICSDIR}/Interpolators/unitTestCubicSplineInterpolation.cpp")
setup_custom_test_program(test_CubicSplineInterpolation "${SRCROOT}${MATHEMATICSDIR}")
target_link_libraries(test_CubicSplineInterpolation tudat_interpolators tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_LinearInterpolation "${SRCROOT}${MATHEMATICSDIR}/Interpolators/unitTestLinearInterpolation.cpp")
setup_custom_test_program(test_LinearInterpolation "${SRCROOT}${MATHEMATICSDIR}")
target_link_libraries(test_LinearInterpolation tudat_interpolators tudat_basic_mathematics ${Boost_LIBRARIES})
