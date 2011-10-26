# \file CMakeLists.txt
 #    This is the CMake file for the RandomNumberGenerators sub-directory.
 #
 #    Path              : /Mathematics/RandomNumberGenerators/
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
set(RANDOMNUMBERGENERATORS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/RandomNumberGenerators/exponentialRandomNumberGenerator.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/RandomNumberGenerators/normalRandomNumberGenerator.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/RandomNumberGenerators/uniformRandomNumberGenerator.cpp"
)

# Add header files.
set(RANDOMNUMBERGENERATORS_HEADERS 
  "${MATHEMATICSDIR}/RandomNumberGenerators/exponentialRandomNumberGenerator.h"
  "${MATHEMATICSDIR}/RandomNumberGenerators/normalRandomNumberGenerator.h"
  "${MATHEMATICSDIR}/RandomNumberGenerators/randomNumberGenerator.h"
  "${MATHEMATICSDIR}/RandomNumberGenerators/uniformRandomNumberGenerator.h"
)

# Add static libraries.
add_library(tudat_random_number_generators STATIC ${RANDOMNUMBERGENERATORS_SOURCES})
setup_tudat_library_target(tudat_random_number_generators "${SRCROOT}${MATHEMATICSDIR}/RandomNumberGenerators")
setup_headers("${RANDOMNUMBERGENERATORS_HEADERS}" "${MATHEMATICSDIR}/RandomNumberGenerators")

# Add unit tests.
add_executable(test_RandomNumberGenerator "${SRCROOT}${MATHEMATICSDIR}/RandomNumberGenerators/unitTestRandomNumberGenerator.cpp")
setup_custom_test_program(test_RandomNumberGenerator "${SRCROOT}${MATHEMATICSDIR}/RandomNumberGenerators")
target_link_libraries(test_RandomNumberGenerator tudat_input tudat_output tudat_mathematics tudat_astrodynamics tudat_basics)

# End of file.
