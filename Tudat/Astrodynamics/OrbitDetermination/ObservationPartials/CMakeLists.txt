# \file CMakeLists.txt
 #    This is the CMake file for the MissionSegments sub-directory.
 #
 #    Path              : /Astrodynamics/MissionSegments/
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
set(OBSERVATION_PARTIALS_SOURCES
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/rotationMatrixPartial.cpp"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/createPositionPartials.cpp"
)

# Set the header files.
set(OBSERVATION_PARTIALS_HEADERS
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/rotationMatrixPartial.h"
  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/createPositionPartials.h"
)

# Add static libraries.
add_library(tudat_observation_partials STATIC ${OBSERVATION_PARTIALS_SOURCES} ${OBSERVATION_PARTIALS_HEADERS})
setup_tudat_library_target(tudat_observation_partials "${SRCROOT}{OBSERVATIONPARTIALSDIR}")

add_executable(test_RotationMatrixPartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}/UnitTests/unitTestRotationMatrixPartials.cpp")
setup_custom_test_program(test_RotationMatrixPartials "${SRCROOT}${OBSERVATIONPARTIALSDIR}")
target_link_libraries(test_RotationMatrixPartials ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})
