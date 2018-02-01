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
set(BASICMATHEMATICS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/coordinateConversions.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/legendrePolynomials.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/nearestNeighbourSearch.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/numericalDerivative.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/sphericalHarmonics.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/rotationAboutArbitraryAxis.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/basicMathematicsFunctions.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/coordinateConversions.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/linearAlgebra.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/leastSquaresEstimation.cpp"
)

# Add header files.
set(BASICMATHEMATICS_HEADERS
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/basicFunction.h"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/convergenceException.h"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/coordinateConversions.h"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/function.h"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/functionProxy.h"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/legendrePolynomials.h"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/linearAlgebra.h"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/nearestNeighbourSearch.h"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/numericalDerivative.h"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/sphericalHarmonics.h"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/rotationAboutArbitraryAxis.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/basicMathematicsFunctions.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/coordinateConversions.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/linearAlgebra.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/mathematicalConstants.h"
  "${SRCROOT}${BASICMATHEMATICSDIR}/leastSquaresEstimation.h"
)

# Add static libraries.
add_library(tudat_basic_mathematics STATIC ${BASICMATHEMATICS_SOURCES} ${BASICMATHEMATICS_HEADERS})
setup_tudat_library_target(tudat_basic_mathematics "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")

# Add unit tests.
add_executable(test_MathematicalConstants "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/UnitTests/unitTestMathematicalConstants.cpp")
setup_custom_test_program(test_MathematicalConstants "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")
target_link_libraries(test_MathematicalConstants tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_LinearAlgebra "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/UnitTests/unitTestLinearAlgebra.cpp")
setup_custom_test_program(test_LinearAlgebra "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")
target_link_libraries(test_LinearAlgebra tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_CoordinateConversions "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/UnitTests/unitTestCoordinateConversions.cpp")
setup_custom_test_program(test_CoordinateConversions "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")
target_link_libraries(test_CoordinateConversions tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_NearestNeighbourSearch "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/UnitTests/unitTestNearestNeighbourSearch.cpp")
setup_custom_test_program(test_NearestNeighbourSearch "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")
target_link_libraries(test_NearestNeighbourSearch tudat_input_output tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_NumericalDerivative "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/UnitTests/unitTestNumericalDerivative.cpp")
setup_custom_test_program(test_NumericalDerivative "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")
target_link_libraries(test_NumericalDerivative tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_LegendrePolynomials "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/UnitTests/unitTestLegendrePolynomials.cpp")
setup_custom_test_program(test_LegendrePolynomials "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")
target_link_libraries(test_LegendrePolynomials tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_SphericalHarmonics "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/UnitTests/unitTestSphericalHarmonics.cpp")
setup_custom_test_program(test_SphericalHarmonics "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")
target_link_libraries(test_SphericalHarmonics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_RotationAboutArbitraryAxis "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/UnitTests/unitTestRotationAboutArbitraryAxis.cpp")
setup_custom_test_program(test_RotationAboutArbitraryAxis "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")
target_link_libraries(test_RotationAboutArbitraryAxis tudat_basic_mathematics ${Boost_LIBRARIES})
