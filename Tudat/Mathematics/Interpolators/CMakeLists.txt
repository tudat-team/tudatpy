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
 #	120202    K. Kumar	    Adapted for new Interpolators sub-directory.
 #      120716    D. Dirkx          Updated with new interpolator architecture.
 #
 #    References
 #
 #    Notes
 #

# Add source files.
set(INTERPOLATORS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/linearInterpolator.cpp"
)

# Add header files.
set(INTERPOLATORS_HEADERS
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/cubicSplineInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/linearInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/lagrangeInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/interpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/lookupScheme.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/oneDimensionalInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/multiLinearInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/createInterpolator.h"
)

# Add static libraries.
add_library(tudat_interpolators STATIC ${INTERPOLATORS_SOURCES} ${INTERPOLATORS_HEADERS})
setup_tudat_library_target(tudat_interpolators "${SRCROOT}${MATHEMATICSDIR}/Interpolators")

# Add unit tests.
add_executable(test_CubicSplineInterpolator "${SRCROOT}${MATHEMATICSDIR}/Interpolators/UnitTests/unitTestCubicSplineInterpolator.cpp")
setup_custom_test_program(test_CubicSplineInterpolator "${SRCROOT}${MATHEMATICSDIR}")
target_link_libraries(test_CubicSplineInterpolator tudat_input_output tudat_interpolators tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_LinearInterpolator "${SRCROOT}${MATHEMATICSDIR}/Interpolators/UnitTests/unitTestLinearInterpolator.cpp")
setup_custom_test_program(test_LinearInterpolator "${SRCROOT}${MATHEMATICSDIR}")
target_link_libraries(test_LinearInterpolator tudat_input_output tudat_interpolators tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_MultiLinearInterpolator "${SRCROOT}${MATHEMATICSDIR}/Interpolators/UnitTests/unitTestMultiLinearInterpolator.cpp")
setup_custom_test_program(test_MultiLinearInterpolator "${SRCROOT}${MATHEMATICSDIR}")
target_link_libraries(test_MultiLinearInterpolator tudat_input_output tudat_interpolators tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_LagrangeInterpolator "${SRCROOT}${MATHEMATICSDIR}/Interpolators/UnitTests/unitTestLagrangeInterpolators.cpp")
setup_custom_test_program(test_LagrangeInterpolator "${SRCROOT}${MATHEMATICSDIR}")
target_link_libraries(test_LagrangeInterpolator tudat_input_output tudat_interpolators tudat_basic_mathematics ${Boost_LIBRARIES})


