# \file CMakeLists.txt
 #    This is the CMake file for the EnvironmentModels sub-directory.
 #
 #    Path              : /Astrodynamics/EnvironmentModels/
 #    Version           : 3
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
 #    Last modified     : 26 October, 2011
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
 #      111026    K. Kumar          Adapted file so all headers show in project tree in Qt Creator.
 #/

# Set the source files.
set(ENVIRONMENTMODELS_SOURCES
  "${SRCROOT}${ENVIRONMENTMODELSDIR}/centralGravityField.cpp"
  "${SRCROOT}${ENVIRONMENTMODELSDIR}/exponentialAtmosphere.cpp"
  "${SRCROOT}${ENVIRONMENTMODELSDIR}/sphericalHarmonicsGravityField.cpp"
  "${SRCROOT}${ENVIRONMENTMODELSDIR}/tabulatedAtmosphere.cpp"
)

# Set the header files.
set(ENVIRONMENTMODELS_HEADERS 
  "${SRCROOT}${ENVIRONMENTMODELSDIR}/atmosphereModel.h"
  "${SRCROOT}${ENVIRONMENTMODELSDIR}/centralGravityField.h"
  "${SRCROOT}${ENVIRONMENTMODELSDIR}/environmentModel.h"
  "${SRCROOT}${ENVIRONMENTMODELSDIR}/exponentialAtmosphere.h"
  "${SRCROOT}${ENVIRONMENTMODELSDIR}/gravityFieldModel.h"
  "${SRCROOT}${ENVIRONMENTMODELSDIR}/sphericalHarmonicsGravityField.h"
  "${SRCROOT}${ENVIRONMENTMODELSDIR}/tabulatedAtmosphere.h"
)

# Add static libraries
add_library(tudat_environment_models STATIC ${ENVIRONMENTMODELS_SOURCES} ${ENVIRONMENTMODELS_HEADERS})
setup_tudat_library_target(tudat_environment_models "${SRCROOT}${ENVIRONMENTMODELSDIR}")
setup_headers("${ENVIRONMENTMODELS_HEADERS}" "${ENVIRONMENTMODELSDIR}")

# Add unit tests.
add_executable(test_ExponentialAtmosphere "${SRCROOT}${ENVIRONMENTMODELSDIR}/unitTestExponentialAtmosphere.cpp")
setup_custom_test_program(test_ExponentialAtmosphere "${SRCROOT}${ENVIRONMENTMODELSDIR}")
target_link_libraries(test_ExponentialAtmosphere tudat_input tudat_output tudat_mathematics tudat_astrodynamics tudat_basics)

add_executable(test_TabulatedAtmosphere "${SRCROOT}${ENVIRONMENTMODELSDIR}/unitTestTabulatedAtmosphere.cpp")
setup_custom_test_program(test_TabulatedAtmosphere "${SRCROOT}${ENVIRONMENTMODELSDIR}")
target_link_libraries(test_TabulatedAtmosphere tudat_input tudat_output tudat_mathematics tudat_astrodynamics tudat_basics tudat_output tudat_input ${Boost_LIBRARIES})

add_executable(test_SphericalHarmonicsGravityField "${SRCROOT}${ENVIRONMENTMODELSDIR}/unitTestSphericalHarmonicsGravityField.cpp")
setup_custom_test_program(test_SphericalHarmonicsGravityField "${SRCROOT}${ENVIRONMENTMODELSDIR}")
target_link_libraries(test_SphericalHarmonicsGravityField tudat_input tudat_output tudat_mathematics tudat_astrodynamics tudat_basics ${Boost_LIBRARIES})

# End of file.
