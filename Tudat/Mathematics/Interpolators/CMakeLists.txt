 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Add source files.
set(INTERPOLATORS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/cubicSplineInterpolator.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/linearInterpolator.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/lagrangeInterpolator.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/interpolator.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/multiLinearInterpolator.cpp"
)

# Add header files.
set(INTERPOLATORS_HEADERS
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/cubicSplineInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/hermiteCubicSplineInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/linearInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/lagrangeInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/interpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/lookupScheme.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/multiDimensionalInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/oneDimensionalInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/multiLinearInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/piecewiseConstantInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/jumpDataLinearInterpolator.h"
  "${SRCROOT}${MATHEMATICSDIR}/Interpolators/createInterpolator.h"
)

# Add static libraries.
add_library(tudat_interpolators STATIC ${INTERPOLATORS_SOURCES} ${INTERPOLATORS_HEADERS})
setup_tudat_library_target(tudat_interpolators "${SRCROOT}${MATHEMATICSDIR}/Interpolators")

# Add unit tests.
add_executable(test_CubicSplineInterpolator "${SRCROOT}${MATHEMATICSDIR}/Interpolators/UnitTests/unitTestCubicSplineInterpolator.cpp")
setup_custom_test_program(test_CubicSplineInterpolator "${SRCROOT}${MATHEMATICSDIR}")
target_link_libraries(test_CubicSplineInterpolator tudat_input_output tudat_interpolators tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_HermiteCubicSplineInterpolator "${SRCROOT}${MATHEMATICSDIR}/Interpolators/UnitTests/unitTestHermiteCubicSplineInterpolator.cpp")
setup_custom_test_program(test_HermiteCubicSplineInterpolator "${SRCROOT}${MATHEMATICSDIR}")
target_link_libraries(test_HermiteCubicSplineInterpolator tudat_input_output tudat_interpolators tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_LinearInterpolator "${SRCROOT}${MATHEMATICSDIR}/Interpolators/UnitTests/unitTestLinearInterpolator.cpp")
setup_custom_test_program(test_LinearInterpolator "${SRCROOT}${MATHEMATICSDIR}")
target_link_libraries(test_LinearInterpolator tudat_input_output tudat_interpolators tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_MultiLinearInterpolator "${SRCROOT}${MATHEMATICSDIR}/Interpolators/UnitTests/unitTestMultiLinearInterpolator.cpp")
setup_custom_test_program(test_MultiLinearInterpolator "${SRCROOT}${MATHEMATICSDIR}")
target_link_libraries(test_MultiLinearInterpolator tudat_input_output tudat_interpolators tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_LagrangeInterpolator "${SRCROOT}${MATHEMATICSDIR}/Interpolators/UnitTests/unitTestLagrangeInterpolators.cpp")
setup_custom_test_program(test_LagrangeInterpolator "${SRCROOT}${MATHEMATICSDIR}")
target_link_libraries(test_LagrangeInterpolator tudat_input_output tudat_interpolators tudat_basic_mathematics ${Boost_LIBRARIES})


