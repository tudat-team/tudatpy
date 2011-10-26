# \file CMakeLists.txt
 #    This is the CMake file for the Mathematics sub-directory.
 #
 #    Path              : /Mathematics/
 #    Version           : 1
 #    Check status      : Checked
 #
 #    Author            : S.M. Persson
 #    Affiliation       : McGill University, Department of Mechanical Engineering
 #    E-mail address    : mikael.s.persson@gmail.com
 #
 #    Author            : K. Kumar
 #    Affiliation       : Delft University of Technology
 #    E-mail address    : K.Kumar@tudelft.nl
 #
 #    Author            : S.M. Persson
 #    Affiliation       : McGill University, Department of Mechanical Engineering
 #    E-mail address    : mikael.s.persson@gmail.com
 #
 #    Date created      : 20 August, 2011
 #    Last modified     : 25 October, 2011
 #
 #    References
 #
 #    Notes
 #
 #    Copyright (c) 2010-2011 Delft University of Technology.
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
 #/

# Add source files.
set(MATHEMATICS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/basicMathematicsFunctions.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/cubicSplineInterpolation.cpp"

  "${SRCROOT}${MATHEMATICSDIR}/unitConversions.h"
)

# Add header files.
set(MATHEMATICS_HEADERS 
  "${MATHEMATICSDIR}/basicMathematicsFunctions.h"
  "${MATHEMATICSDIR}/cubicSplineInterpolation.h"
)

# Add sub-directories.
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/GeometricShapes")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/LinearAlgebra")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/RandomNumberGenerators")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/RootFindingMethods")
add_subdirectory("${SRCROOT}${MATHEMATICSDIR}/Statistics")

# Get target properties for static libraries.
get_target_property(GEOMETRICSHAPESSOURCES tudat_geometric_shapes SOURCES)
get_target_property(LINEARALGEBRASOURCES tudat_linear_algebra SOURCES)
get_target_property(NUMERICALINTEGRATORSSOURCES tudat_numerical_integrators SOURCES)
get_target_property(RANDOMNUMBERGENERATORSSOURCES tudat_random_number_generators SOURCES)
get_target_property(ROOTFINDINGMETHODSSOURCES tudat_root_finding_methods SOURCES)
get_target_property(STATISTICSSOURCES tudat_statistics SOURCES)

# Add static libraries.
add_library(tudat_mathematics STATIC ${MATHEMATICS_SOURCES} ${GEOMETRICSHAPESSOURCES} ${LINEARALGEBRASOURCES} ${NUMERICALINTEGRATORSSOURCES} ${RANDOMNUMBERGENERATORSSOURCES} ${ROOTFINDINGMETHODSSOURCES} ${STATISTICSSOURCES})
setup_tudat_library_target(tudat_mathematics "${SRCROOT}${MATHEMATICSDIR}")
setup_headers("${MATHEMATICS_HEADERS}" "${MATHEMATICSDIR}")

# Add unit tests.
add_executable(test_BasicMathematicsFunctions "${SRCROOT}${MATHEMATICSDIR}/unitTestBasicMathematicsFunctions.cpp")
setup_custom_test_program(test_BasicMathematicsFunctions "${SRCROOT}${MATHEMATICSDIR}")
target_link_libraries(test_BasicMathematicsFunctions tudat_astrodynamics tudat_input tudat_output tudat_mathematics tudat_basics)

add_executable(test_CubicSplineInterpolation "${SRCROOT}${MATHEMATICSDIR}/unitTestCubicSplineInterpolation.cpp")
setup_custom_test_program(test_CubicSplineInterpolation "${SRCROOT}${MATHEMATICSDIR}")
target_link_libraries(test_CubicSplineInterpolation tudat_astrodynamics tudat_input tudat_output tudat_mathematics tudat_basics)

add_executable(test_UnitConversions "${SRCROOT}${MATHEMATICSDIR}/unitTestUnitConversions.cpp")
setup_custom_test_program(test_UnitConversions "${SRCROOT}${MATHEMATICSDIR}")
target_link_libraries(test_UnitConversions tudat_input tudat_output tudat_mathematics tudat_astrodynamics tudat_basics)

# End of file.
