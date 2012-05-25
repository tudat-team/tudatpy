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
 #      111025    K. Kumar          Adapted file to work with Revision 194.
 #      111026    K. Kumar          Adapted file so all headers show in project tree in Qt Creator.
 #      120909    B. Tong Minh      Moved to BasicMathematics
 #

# Add source files.
set(BASICMATHEMATICS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/coordinateConversions.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/nearestNeighbourSearch.cpp"
)

# Add header files.
set(BASICMATHEMATICS_HEADERS 
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/coordinateConversions.h"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/nearestNeighbourSearch.h"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/numericalDerivative.h"
)

# Add static libraries.
add_library(tudat_basic_mathematics STATIC ${BASICMATHEMATICS_SOURCES} ${BASICMATHEMATICS_HEADERS})
setup_tudat_library_target(tudat_basic_mathematics "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")

# Add unit tests.
add_executable(test_CoordinateConversions "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/UnitTests/unitTestCoordinateConversions.cpp")
setup_custom_test_program(test_CoordinateConversions "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")
target_link_libraries(test_CoordinateConversions tudat_basic_mathematics ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_NearestNeighbourSearch "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/UnitTests/unitTestNearestNeighbourSearch.cpp")
setup_custom_test_program(test_NearestNeighbourSearch "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")
target_link_libraries(test_NearestNeighbourSearch tudat_input_output tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_NumericalDerivative "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/UnitTests/unitTestNumericalDerivative.cpp")
setup_custom_test_program(test_NumericalDerivative "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")
target_link_libraries(test_NumericalDerivative ${Boost_LIBRARIES})
