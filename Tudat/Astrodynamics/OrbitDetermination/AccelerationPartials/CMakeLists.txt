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
set(ACCELERATION_PARTIALS_SOURCES
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/accelerationPartial.cpp"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/centralGravityAccelerationPartial.cpp"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/numericalAccelerationPartial.cpp"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/radiationPressureAccelerationPartial.cpp"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/sphericalHarmonicPartialFunctions.cpp"
)

# Set the header files.
set(ACCELERATION_PARTIALS_HEADERS
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/accelerationPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/thirdBodyGravityPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/centralGravityAccelerationPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/numericalAccelerationPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/radiationPressureAccelerationPartial.h"
  "${SRCROOT}${ACCELERATIONPARTIALSDIR}/sphericalHarmonicPartialFunctions.h"
)

# Add static libraries.
add_library(tudat_acceleration_partials STATIC ${ACCELERATION_PARTIALS_SOURCES} ${ACCELERATION_PARTIALS_HEADERS})
setup_tudat_library_target(tudat_acceleration_partials "${SRCROOT}{ACCELERATIONPARTIALSDIR}")

if(USE_CSPICE)

# Add unit tests
add_executable(test_AccelerationPartials "${SRCROOT}${ACCELERATIONPARTIALSDIR}/UnitTests/unitTestAccelerationPartials.cpp")
setup_custom_test_program(test_AccelerationPartials "${SRCROOT}${ACCELERATIONPARTIALSDIR}")
target_link_libraries(test_AccelerationPartials ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

# Add unit tests
add_executable(test_SphericalHarmonicPartials "${SRCROOT}${ACCELERATIONPARTIALSDIR}/UnitTests/unitTestSphericalHarmonicPartials.cpp")
setup_custom_test_program(test_SphericalHarmonicPartials "${SRCROOT}${ACCELERATIONPARTIALSDIR}")
target_link_libraries(test_SphericalHarmonicPartials ${TUDAT_ESTIMATION_LIBRARIES} ${Boost_LIBRARIES})

endif( )


