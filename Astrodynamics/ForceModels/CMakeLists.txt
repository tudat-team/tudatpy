# \file CMakeLists.txt
 #    This is the CMake file for the ForceModels sub-directory.
 #
 #    Path              : /Astrodynamics/ForceModels/
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
set(FORCEMODELS_SOURCES
  "${SRCROOT}${FORCEMODELSDIR}/aerodynamicCoefficientGenerator.cpp"
  "${SRCROOT}${FORCEMODELSDIR}/aerodynamics.cpp"
  "${SRCROOT}${FORCEMODELSDIR}/hypersonicLocalInclinationAnalysis.cpp"
)

# Set the header files.
set(FORCEMODELS_HEADERS 
  "${SRCROOT}${FORCEMODELSDIR}/aerodynamicCoefficientGenerator.h"
  "${SRCROOT}${FORCEMODELSDIR}/aerodynamicCoefficientInterface.h"
  "${SRCROOT}${FORCEMODELSDIR}/aerodynamicForce.h"
  "${SRCROOT}${FORCEMODELSDIR}/aerodynamics.h"
  "${SRCROOT}${FORCEMODELSDIR}/forceModel.h"
  "${SRCROOT}${FORCEMODELSDIR}/gravitationalForceModel.h"
  "${SRCROOT}${FORCEMODELSDIR}/hypersonicLocalInclinationAnalysis.h"
)

# Add static libraries.
add_library(tudat_force_models STATIC ${FORCEMODELS_SOURCES} ${FORCEMODELS_HEADERS})
setup_tudat_library_target(tudat_force_models "${SRCROOT}${FORCEMODELSDIR}")
setup_headers("${FORCEMODELS_HEADERS}" "${FORCEMODELSDIR}")

# Add unit tests.
add_executable(test_AerodynamicMomentAndAerodynamicForce "${SRCROOT}${FORCEMODELSDIR}/unitTestAerodynamicMomentAndAerodynamicForce.cpp")
setup_custom_test_program(test_AerodynamicMomentAndAerodynamicForce "${SRCROOT}${FORCEMODELSDIR}")
target_link_libraries(test_AerodynamicMomentAndAerodynamicForce tudat_input tudat_output tudat_mathematics tudat_astrodynamics tudat_basics tudat_force_models tudat_moment_models)

add_executable(test_AerodynamicsNamespace "${SRCROOT}${FORCEMODELSDIR}/unitTestAerodynamicsNamespace.cpp")
setup_custom_test_program(test_AerodynamicsNamespace "${SRCROOT}${FORCEMODELSDIR}")
target_link_libraries(test_AerodynamicsNamespace tudat_input tudat_output tudat_mathematics tudat_astrodynamics tudat_basics tudat_force_models)

add_executable(test_CoefficientGenerator "${SRCROOT}${FORCEMODELSDIR}/unitTestCoefficientGenerator.cpp")
setup_custom_test_program(test_CoefficientGenerator "${SRCROOT}${FORCEMODELSDIR}")
target_link_libraries(test_CoefficientGenerator tudat_input tudat_output tudat_astrodynamics tudat_mathematics tudat_astrodynamics tudat_basics tudat_force_models ${Boost_LIBRARIES})

# End of file.
