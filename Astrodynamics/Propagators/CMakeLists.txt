# \file CMakeLists.txt
 #    This is the CMake file for the Propagators sub-directory.
 #
 #    Path              : /Astrodynamics/Propagators/
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
set(PROPAGATORS_SOURCES
  "${SRCROOT}${PROPAGATORSDIR}/cartesianStateNumericalPropagator.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/keplerPropagator.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/numericalPropagator.cpp"
  "${SRCROOT}${PROPAGATORSDIR}/seriesPropagator.cpp"

  "${SRCROOT}${PROPAGATORSDIR}/propagator.h"
  "${SRCROOT}${PROPAGATORSDIR}/propagatorDataContainer.h"
)

# Add header files.
set(PROPAGATORS_HEADERS 
  "${PROPAGATORSDIR}/cartesianStateNumericalPropagator.h"
  "${PROPAGATORSDIR}/keplerPropagator.h"
  "${PROPAGATORSDIR}/numericalPropagator.h"
  "${PROPAGATORSDIR}/seriesPropagator.h"
)

# Add static libraries.
add_library(tudat_propagators STATIC ${PROPAGATORS_SOURCES})
setup_tudat_library_target(tudat_propagators "${SRCROOT}${PROPAGATORSDIR}")
setup_headers("${PROPAGATORS_HEADERS}" "${PROPAGATORSDIR}")

# Add unit tests.
add_executable(test_KeplerPropagator "${SRCROOT}${PROPAGATORSDIR}/unitTestKeplerPropagator.cpp")
setup_custom_test_program(test_KeplerPropagator "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_KeplerPropagator tudat_astrodynamics tudat_input tudat_output tudat_mathematics tudat_basics)

add_executable(test_NumericalPropagator "${SRCROOT}${PROPAGATORSDIR}/unitTestNumericalPropagator.cpp")
setup_custom_test_program(test_NumericalPropagator "${SRCROOT}${PROPAGATORSDIR}")
target_link_libraries(test_NumericalPropagator tudat_astrodynamics tudat_input tudat_output tudat_mathematics tudat_basics)

# End of file.
