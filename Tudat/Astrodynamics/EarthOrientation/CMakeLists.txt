#    Copyright (c) 2010-2018, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.

# Set the source files.
set(EARTH_ORIENTATION_SOURCES
  "${SRCROOT}${EARTHORIENTATIONDIR}/earthOrientationCalculator.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/terrestrialTimeScaleConverter.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/eopReader.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/polarMotionCalculator.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/precessionNutationCalculator.cpp"
  "${SRCROOT}${EARTHORIENTATIONDIR}/readAmplitudeAndArgumentMultipliers.cpp"
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
  "${SRCROOT}${EARTHORIENTATIONDIR}/readAmplitudeAndArgumentMultipliers.h"
  "${SRCROOT}${EARTHORIENTATIONDIR}/shortPeriodEarthOrientationCorrectionCalculator.h"
  "${SRCROOT}${EARTHORIENTATIONDIR}/UnitTests/sofaEarthOrientationCookbookExamples.h"
)

# Add static libraries.
add_library(tudat_earth_orientation STATIC ${EARTH_ORIENTATION_SOURCES} ${EARTH_ORIENTATION_HEADERS})
setup_tudat_library_target(tudat_earth_orientation "${SRCROOT}${EARTHORIENTATIONDIR}")

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

add_executable(test_ShortPeriodEopCorrections "${SRCROOT}${EARTHORIENTATIONDIR}/UnitTests/unitTestShortPeriodEopCorrections.cpp")
setup_custom_test_program(test_ShortPeriodEopCorrections "${SRCROOT}${EARTHORIENTATIONDIR}")
target_link_libraries(test_ShortPeriodEopCorrections tudat_earth_orientation tudat_sofa_interface tudat_basic_astrodynamics tudat_basic_astrodynamics tudat_basic_mathematics tudat_input_output sofa ${Boost_LIBRARIES})


