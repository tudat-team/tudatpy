 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Set the source files.
set(EPHEMERIDES_SOURCES
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositionsBase.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositions.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositionsCircularCoplanar.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/ephemeris.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/rotationalEphemeris.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/cartesianStateExtractor.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/keplerStateExtractor.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/keplerEphemeris.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/rotationalEphemeris.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/simpleRotationalEphemeris.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/tabulatedEphemeris.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/frameManager.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/compositeEphemeris.cpp"
  "${SRCROOT}${EPHEMERIDESDIR}/tabulatedRotationalEphemeris.cpp"
)

# Set the header files.
set(EPHEMERIDES_HEADERS
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositionsBase.h"
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositions.h"
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositionsCircularCoplanar.h"
  "${SRCROOT}${EPHEMERIDESDIR}/approximatePlanetPositionsDataContainer.h"
  "${SRCROOT}${EPHEMERIDESDIR}/ephemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/constantEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/cartesianStateExtractor.h"
  "${SRCROOT}${EPHEMERIDESDIR}/keplerStateExtractor.h"
  "${SRCROOT}${EPHEMERIDESDIR}/keplerEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/rotationalEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/simpleRotationalEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/tabulatedEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/frameManager.h"
  "${SRCROOT}${EPHEMERIDESDIR}/itrsToGcrsRotationModel.h"
  "${SRCROOT}${EPHEMERIDESDIR}/compositeEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/constantEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/constantRotationalEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/multiArcEphemeris.h"
  "${SRCROOT}${EPHEMERIDESDIR}/tabulatedRotationalEphemeris.h"
)

# Add static libraries.
add_library(tudat_ephemerides STATIC ${EPHEMERIDES_SOURCES} ${EPHEMERIDES_HEADERS})
setup_tudat_library_target(tudat_ephemerides "${SRCROOT}${EPHEMERIDESDIR}")

# Add unit tests.
add_executable(test_ApproximatePlanetPositions "${SRCROOT}${EPHEMERIDESDIR}/UnitTests/unitTestApproximatePlanetPositions.cpp")
setup_custom_test_program(test_ApproximatePlanetPositions "${SRCROOT}${EPHEMERIDESDIR}")
target_link_libraries(test_ApproximatePlanetPositions tudat_ephemerides tudat_gravitation tudat_basic_astrodynamics tudat_basic_mathematics tudat_input_output tudat_root_finders  ${Boost_LIBRARIES})

add_executable(test_TabulatedEphemeris "${SRCROOT}${EPHEMERIDESDIR}/UnitTests/unitTestTabulatedEphemeris.cpp")
setup_custom_test_program(test_TabulatedEphemeris "${SRCROOT}${EPHEMERIDESDIR}")
target_link_libraries(test_TabulatedEphemeris tudat_ephemerides tudat_interpolators tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_CartesianStateExtractor "${SRCROOT}${EPHEMERIDESDIR}/UnitTests/unitTestCartesianStateExtractor.cpp")
setup_custom_test_program(test_CartesianStateExtractor "${SRCROOT}${EPHEMERIDESDIR}")
target_link_libraries(test_CartesianStateExtractor tudat_input_output tudat_ephemerides ${Boost_LIBRARIES})

add_executable(test_KeplerStateExtractor "${SRCROOT}${EPHEMERIDESDIR}/UnitTests/unitTestKeplerStateExtractor.cpp")
setup_custom_test_program(test_KeplerStateExtractor "${SRCROOT}${EPHEMERIDESDIR}")
target_link_libraries(test_KeplerStateExtractor tudat_ephemerides tudat_input_output tudat_basic_astrodynamics tudat_basic_mathematics tudat_root_finders ${Boost_LIBRARIES})

add_executable(test_RotationalEphemeris "${SRCROOT}${EPHEMERIDESDIR}/UnitTests/unitTestRotationalEphemeris.cpp")
setup_custom_test_program(test_RotationalEphemeris "${SRCROOT}${EPHEMERIDESDIR}")
target_link_libraries(test_RotationalEphemeris tudat_ephemerides tudat_reference_frames tudat_input_output tudat_basic_astrodynamics tudat_basic_mathematics tudat_spice_interface ${TUDAT_EXTERNAL_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_SimpleRotationalEphemeris "${SRCROOT}${EPHEMERIDESDIR}/UnitTests/unitTestSimpleRotationalEphemeris.cpp")
setup_custom_test_program(test_SimpleRotationalEphemeris "${SRCROOT}${EPHEMERIDESDIR}")
target_link_libraries(test_SimpleRotationalEphemeris tudat_ephemerides tudat_reference_frames tudat_input_output tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

if(USE_CSPICE)
add_executable(test_FrameManager "${SRCROOT}${EPHEMERIDESDIR}/UnitTests/unitTestFrameManager.cpp")
setup_custom_test_program(test_FrameManager "${SRCROOT}${EPHEMERIDESDIR}")
target_link_libraries(test_FrameManager tudat_ephemerides tudat_reference_frames ${TUDAT_EXTERNAL_INTERFACE_LIBRARIES} tudat_input_output tudat_basic_astrodynamics tudat_basic_mathematics ${TUDAT_EXTERNAL_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_CompositeEphemeris "${SRCROOT}${EPHEMERIDESDIR}/UnitTests/unitTestCompositeEphemeris.cpp")
setup_custom_test_program(test_CompositeEphemeris "${SRCROOT}${EPHEMERIDESDIR}")
target_link_libraries(test_CompositeEphemeris ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_TabulatedRotationalEphemeris "${SRCROOT}${EPHEMERIDESDIR}/UnitTests/unitTestTabulatedRotationalEphemeris.cpp")
setup_custom_test_program(test_TabulatedRotationalEphemeris "${SRCROOT}${EPHEMERIDESDIR}")
target_link_libraries(test_TabulatedRotationalEphemeris tudat_ephemerides tudat_reference_frames tudat_input_output tudat_basic_astrodynamics tudat_basic_mathematics tudat_spice_interface ${TUDAT_EXTERNAL_LIBRARIES} ${Boost_LIBRARIES})

endif()

add_executable(test_KeplerEphemeris "${SRCROOT}${EPHEMERIDESDIR}/UnitTests/unitTestKeplerEphemeris.cpp")
setup_custom_test_program(test_KeplerEphemeris "${SRCROOT}${EPHEMERIDESDIR}")
target_link_libraries(test_KeplerEphemeris tudat_ephemerides tudat_reference_frames tudat_input_output tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

if(USE_SOFA)
add_executable(test_GcrsToItrsRotation "${SRCROOT}${EPHEMERIDESDIR}/UnitTests/unitTestItrsToGcrsRotationModel.cpp")
setup_custom_test_program(test_GcrsToItrsRotation "${SRCROOT}${EPHEMERIDESDIR}")
target_link_libraries(test_GcrsToItrsRotation tudat_ephemerides tudat_spice_interface tudat_reference_frames tudat_earth_orientation tudat_sofa_interface tudat_basic_astrodynamics tudat_basic_mathematics tudat_input_output cspice sofa ${Boost_LIBRARIES})
endif( )
