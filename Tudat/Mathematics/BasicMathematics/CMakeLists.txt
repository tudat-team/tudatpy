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
 #      120909    B. Tong Minh      Moved to BasicMathematics
 #      120920    E. Dekens         Added 'legendre' and 'sphericalHarmonics' files.
 #      140212    S. Hirsh          Added 'rotation about arbitrary axis' files.
 #
 #    References
 #
 #    Notes
 #

# Add source files.
set(BASICMATHEMATICS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/coordinateConversions.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/legendrePolynomials.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/nearestNeighbourSearch.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/sphericalHarmonics.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/rotationAboutArbitraryAxis.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/basicMathematicsFunctions.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/coordinateConversions.cpp"
  "${SRCROOT}${BASICMATHEMATICSDIR}/linearAlgebra.cpp"
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
target_link_libraries(test_NumericalDerivative ${Boost_LIBRARIES})

add_executable(test_LegendrePolynomials "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/UnitTests/unitTestLegendrePolynomials.cpp")
setup_custom_test_program(test_LegendrePolynomials "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")
target_link_libraries(test_LegendrePolynomials tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_SphericalHarmonics "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/UnitTests/unitTestSphericalHarmonics.cpp")
setup_custom_test_program(test_SphericalHarmonics "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")
target_link_libraries(test_SphericalHarmonics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_RotationAboutArbitraryAxis "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics/UnitTests/unitTestRotationAboutArbitraryAxis.cpp")
setup_custom_test_program(test_RotationAboutArbitraryAxis "${SRCROOT}${MATHEMATICSDIR}/BasicMathematics")
target_link_libraries(test_RotationAboutArbitraryAxis tudat_basic_mathematics ${Boost_LIBRARIES})
