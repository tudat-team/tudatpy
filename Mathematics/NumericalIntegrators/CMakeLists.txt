# \file CMakeLists.txt
 #    This is the CMake file for the NumericalIntegrators sub-directory.
 #
 #    Path              : /Mathematics/NumericalIntegrators/
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
set(NUMERICALINTEGRATORS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/euler.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/integrator.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/rungeKutta4thOrderFixedStepsize.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/singleStepIntegrationMethods.cpp"

  "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/stateDerivativeBase.h"
)

# Add header files.
set(NUMERICALINTEGRATORS_HEADERS 
  "${MATHEMATICSDIR}/NumericalIntegrators/euler.h"
  "${MATHEMATICSDIR}/NumericalIntegrators/integrator.h"
  "${MATHEMATICSDIR}/NumericalIntegrators/rungeKutta4thOrderFixedStepsize.h"
  "${MATHEMATICSDIR}/NumericalIntegrators/singleStepIntegrationMethods.h"
)

# Add static libraries.
add_library(tudat_numerical_integrators STATIC ${NUMERICALINTEGRATORS_SOURCES})
setup_tudat_library_target(tudat_numerical_integrators "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
setup_headers("${NUMERICALINTEGRATORS_HEADERS}" "${MATHEMATICSDIR}/NumericalIntegrators")

# Add unit tests.
add_executable(test_EulerIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/unitTestEulerIntegrator.cpp")
setup_custom_test_program(test_EulerIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_EulerIntegrator tudat_input tudat_output tudat_mathematics tudat_astrodynamics tudat_basics)

add_executable(test_RungeKutta4thOrderFixedStepsizeIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators/unitTestRungeKutta4thOrderFixedStepsizeIntegrator.cpp")
setup_custom_test_program(test_RungeKutta4thOrderFixedStepsizeIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalIntegrators")
target_link_libraries(test_RungeKutta4thOrderFixedStepsizeIntegrator tudat_input tudat_output tudat_mathematics tudat_astrodynamics tudat_basics)

# End of file.
