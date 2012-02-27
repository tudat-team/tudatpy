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
setup_headers("${INTERPOLATORS_HEADERS}" "${MATHEMATICSDIR}/Interpolators")

# Add unit tests.
add_executable(test_CubicSplineInterpolation "${SRCROOT}${MATHEMATICSDIR}/Interpolators/unitTestCubicSplineInterpolation.cpp")
setup_custom_test_program(test_CubicSplineInterpolation "${SRCROOT}${MATHEMATICSDIR}")
target_link_libraries(test_CubicSplineInterpolation tudat_interpolators tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_LinearInterpolation "${SRCROOT}${MATHEMATICSDIR}/Interpolators/unitTestLinearInterpolation.cpp")
setup_custom_test_program(test_LinearInterpolation "${SRCROOT}${MATHEMATICSDIR}")
target_link_libraries(test_LinearInterpolation tudat_interpolators tudat_basic_mathematics ${Boost_LIBRARIES})
