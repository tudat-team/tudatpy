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
 #      120209    B. Tong Minh      Merge Ephemeris package
 #
# Set the source files.
set(EARTH_ORIENTATION_SOURCES
  "${SRCROOT}${EARTHORIENTATIONDIR}/earthOrientationCalculator.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/terrestrialTimeScaleConverter.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/eopReader.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/polarMotionCalculator.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/precessionNutationCalculator.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/readAmplitudeAndDoodsonNumber.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/UnitTests/sofaEarthOrientationCookbookExamples.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/shortPeriodEarthOrientationCorrectionCalculator.cpp"
)

# Set the header files.
set(EARTH_ORIENTATION_HEADERS
  "${SRCROOT}${EARTHORIENTATIONDIR}/earthOrientationCalculator.h"
  "${SRCROOT}${EARTHORIENTATIONDIR}/terrestrialTimeScaleConverter.h"
  "${SRCROOT}${EARTHORIENTATIONDIR}/eopReader.h"
  "${SRCROOT}${EARTHORIENTATIONDIR}/polarMotionCalculator.h"
  "${SRCROOT}${EARTHORIENTATIONDIR}/precessionNutationCalculator.h"
  "${SRCROOT}${EARTHORIENTATIONDIR}/readAmplitudeAndDoodsonNumber.h"
  "${SRCROOT}${EARTHORIENTATIONDIR}/shortPeriodEarthOrientationCorrectionCalculator.h"
)

# Add static libraries.
add_library(tudat_earth_orientation STATIC ${EARTH_ORIENTATION_SOURCES} ${EARTH_ORIENTATION_HEADERS})
setup_tudat_library_target(tudat_earth_orientation "${SRCROOT}${EARTHORIENTATIONDIR}")

#add_executable(test_EarthOrientation "${SRCROOT}${EARTHORIENTATIONDIR}/earthOrientationTest.cpp")
#setup_custom_test_program(test_EarthOrientation "${SRCROOT}${ORBITDETERMINATIONDIR}")
#target_link_libraries(test_EarthOrientation tudat_earth_orientation tudat_spice_interface tudat_sofa_interface tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_input_output sofa cspice ${Boost_LIBRARIES})

add_executable(test_EarthOrientationCalculator "${SRCROOT}${EARTHORIENTATIONDIR}/UnitTests/unitTestEarthOrientationCalculator.cpp")
setup_custom_test_program(test_EarthOrientationCalculator "${SRCROOT}${ORBITDETERMINATIONDIR}")
target_link_libraries(test_EarthOrientationCalculator tudat_earth_orientation tudat_spice_interface tudat_sofa_interface tudat_interpolators tudat_basic_astrodynamics tudat_basic_mathematics tudat_input_output tudat_input_output sofa cspice ${Boost_LIBRARIES})

add_executable(test_EopReader "${SRCROOT}${EARTHORIENTATIONDIR}/UnitTests/unitTestEopReader.cpp")
setup_custom_test_program(test_EopReader "${SRCROOT}${ORBITDETERMINATIONDIR}")
target_link_libraries(test_EopReader tudat_earth_orientation tudat_sofa_interface tudat_interpolators tudat_basic_astrodynamics tudat_basic_mathematics tudat_input_output sofa ${Boost_LIBRARIES})

add_executable(test_PolarMotion "${SRCROOT}${EARTHORIENTATIONDIR}/UnitTests/unitTestPolarMotionCalculator.cpp")
setup_custom_test_program(test_PolarMotion "${SRCROOT}${ORBITDETERMINATIONDIR}")
target_link_libraries(test_PolarMotion tudat_earth_orientation tudat_sofa_interface tudat_interpolators tudat_basic_astrodynamics tudat_basic_mathematics tudat_input_output sofa ${Boost_LIBRARIES})

add_executable(test_TimeScaleConverter "${SRCROOT}${EARTHORIENTATIONDIR}/UnitTests/unitTestTimeScaleConverter.cpp")
setup_custom_test_program(test_TimeScaleConverter "${SRCROOT}${EARTHORIENTATIONDIR}")
target_link_libraries(test_TimeScaleConverter tudat_earth_orientation tudat_sofa_interface tudat_interpolators tudat_basic_astrodynamics tudat_basic_mathematics tudat_input_output sofa ${Boost_LIBRARIES})

add_executable(test_ShortPeriodPolarMotion "${SRCROOT}${EARTHORIENTATIONDIR}/UnitTests/unitTestShortPeriodPolarMotion.cpp")
setup_custom_test_program(test_ShortPeriodPolarMotion "${SRCROOT}${EARTHORIENTATIONDIR}")
target_link_libraries(test_ShortPeriodPolarMotion tudat_earth_orientation tudat_sofa_interface tudat_basic_astrodynamics tudat_basic_astrodynamics tudat_basic_mathematics tudat_input_output sofa ${Boost_LIBRARIES})

add_executable(test_ShortPeriodUT1Correction "${SRCROOT}${EARTHORIENTATIONDIR}/UnitTests/unitTestShortPeriodUT1Corrections.cpp")
setup_custom_test_program(test_ShortPeriodUT1Correction "${SRCROOT}${EARTHORIENTATIONDIR}")
target_link_libraries(test_ShortPeriodUT1Correction tudat_earth_orientation tudat_sofa_interface tudat_basic_astrodynamics tudat_basic_astrodynamics tudat_basic_mathematics tudat_input_output sofa ${Boost_LIBRARIES})

#add_executable(test_EarthOrientationInterpolator "${SRCROOT}${EARTHORIENTATIONDIR}/UnitTests/testEarthOrientationInterpolator.cpp")
#setup_application_executable_target(test_EarthOrientationInterpolator "${SRCROOT}${EARTHORIENTATIONDIR}")
#target_link_libraries(test_EarthOrientationInterpolator tudat_earth_orientation tudat_sofa_interface tudat_interpolators tudat_basic_astrodynamics tudat_basic_mathematics tudat_input_output tudat_input_output sofa ${Boost_LIBRARIES})



