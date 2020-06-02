#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Set the source files.
set(ephemerides_SOURCES
        "approximatePlanetPositionsBase.cpp"
        "approximatePlanetPositions.cpp"
        "approximatePlanetPositionsCircularCoplanar.cpp"
        "ephemeris.cpp"
        "rotationalEphemeris.cpp"
        "cartesianStateExtractor.cpp"
        "keplerStateExtractor.cpp"
        "keplerEphemeris.cpp"
        "rotationalEphemeris.cpp"
        "simpleRotationalEphemeris.cpp"
        "tabulatedEphemeris.cpp"
        "frameManager.cpp"
        "compositeEphemeris.cpp"
        "tabulatedRotationalEphemeris.cpp"
        "synchronousRotationalEphemeris.cpp"
        "fullPlanetaryRotationModel.cpp"
        )

# Set the header files.
set(ephemerides_HEADERS
        "approximatePlanetPositionsBase.h"
        "approximatePlanetPositions.h"
        "approximatePlanetPositionsCircularCoplanar.h"
        "approximatePlanetPositionsDataContainer.h"
        "ephemeris.h"
        "constantEphemeris.h"
        "cartesianStateExtractor.h"
        "keplerStateExtractor.h"
        "keplerEphemeris.h"
        "rotationalEphemeris.h"
        "constantRotationalEphemeris.h"
        "simpleRotationalEphemeris.h"
        "tabulatedEphemeris.h"
        "frameManager.h"
        "itrsToGcrsRotationModel.h"
        "compositeEphemeris.h"
        "constantEphemeris.h"
        "constantRotationalEphemeris.h"
        "multiArcEphemeris.h"
        "tabulatedRotationalEphemeris.h"
        "fullPlanetaryRotationModel.h"
        "synchronousRotationalEphemeris.h"
        )


TUDAT_ADD_LIBRARY("ephemerides"
        "${ephemerides_SOURCES}"
        "${ephemerides_HEADERS}"
        PUBLIC_LINKS tudat_spice_interface tudat_reference_frames tudat_interpolators
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}" "${CSpice_INCLUDE_DIRS}"
        )

# Add tests.
TUDAT_ADD_TEST_CASE(ApproximatePlanetPositions
        PRIVATE_LINKS
        tudat_ephemerides
        tudat_gravitation
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        tudat_input_output
        tudat_root_finders
        )

TUDAT_ADD_TEST_CASE(TabulatedEphemeris
        PRIVATE_LINKS
        tudat_ephemerides
        tudat_gravitation
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        tudat_input_output
        tudat_interpolators
        )

TUDAT_ADD_TEST_CASE(CartesianStateExtractor
        PRIVATE_LINKS
        tudat_ephemerides
        tudat_gravitation
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        tudat_input_output
        tudat_interpolators
        )

## Add unit tests.
#add_executable(test_ApproximatePlanetPositions "${SRCROOT}${EPHEMERIDESDIR}/tests/unitTestApproximatePlanetPositions.cpp")
#setup_custom_test_program(test_ApproximatePlanetPositions "${SRCROOT}${EPHEMERIDESDIR}")
#target_link_libraries(test_ApproximatePlanetPositions tudat_ephemerides tudat_gravitation tudat_basic_astrodynamics tudat_basic_mathematics tudat_input_output tudat_root_finders  ${Boost_LIBRARIES})
#
#add_executable(test_TabulatedEphemeris "${SRCROOT}${EPHEMERIDESDIR}/tests/unitTestTabulatedEphemeris.cpp")
#setup_custom_test_program(test_TabulatedEphemeris "${SRCROOT}${EPHEMERIDESDIR}")
#target_link_libraries(test_TabulatedEphemeris tudat_ephemerides tudat_interpolators tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})
#
#add_executable(test_CartesianStateExtractor "${SRCROOT}${EPHEMERIDESDIR}/tests/unitTestCartesianStateExtractor.cpp")
#setup_custom_test_program(test_CartesianStateExtractor "${SRCROOT}${EPHEMERIDESDIR}")
#target_link_libraries(test_CartesianStateExtractor tudat_input_output tudat_ephemerides ${Boost_LIBRARIES})
#
#add_executable(test_KeplerStateExtractor "${SRCROOT}${EPHEMERIDESDIR}/tests/unitTestKeplerStateExtractor.cpp")
#setup_custom_test_program(test_KeplerStateExtractor "${SRCROOT}${EPHEMERIDESDIR}")
#target_link_libraries(test_KeplerStateExtractor tudat_ephemerides tudat_input_output tudat_basic_astrodynamics tudat_basic_mathematics tudat_root_finders ${Boost_LIBRARIES})
#
#add_executable(test_RotationalEphemeris "${SRCROOT}${EPHEMERIDESDIR}/tests/unitTestRotationalEphemeris.cpp")
#setup_custom_test_program(test_RotationalEphemeris "${SRCROOT}${EPHEMERIDESDIR}")
#target_link_libraries(test_RotationalEphemeris tudat_ephemerides tudat_reference_frames tudat_input_output tudat_basic_astrodynamics tudat_basic_mathematics tudat_spice_interface ${TUDAT_EXTERNAL_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_SimpleRotationalEphemeris "${SRCROOT}${EPHEMERIDESDIR}/tests/unitTestSimpleRotationalEphemeris.cpp")
#setup_custom_test_program(test_SimpleRotationalEphemeris "${SRCROOT}${EPHEMERIDESDIR}")
#target_link_libraries(test_SimpleRotationalEphemeris tudat_ephemerides tudat_reference_frames tudat_input_output tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})
#
#if(TUDAT_BUILD_WITH_SPICE_INTERFACE)
#add_executable(test_FrameManager "${SRCROOT}${EPHEMERIDESDIR}/tests/unitTestFrameManager.cpp")
#setup_custom_test_program(test_FrameManager "${SRCROOT}${EPHEMERIDESDIR}")
#target_link_libraries(test_FrameManager tudat_ephemerides tudat_reference_frames ${TUDAT_EXTERNAL_INTERFACE_LIBRARIES} tudat_input_output tudat_basic_astrodynamics tudat_basic_mathematics ${TUDAT_EXTERNAL_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_CompositeEphemeris "${SRCROOT}${EPHEMERIDESDIR}/tests/unitTestCompositeEphemeris.cpp")
#setup_custom_test_program(test_CompositeEphemeris "${SRCROOT}${EPHEMERIDESDIR}")
#target_link_libraries(test_CompositeEphemeris ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_TabulatedRotationalEphemeris "${SRCROOT}${EPHEMERIDESDIR}/tests/unitTestTabulatedRotationalEphemeris.cpp")
#setup_custom_test_program(test_TabulatedRotationalEphemeris "${SRCROOT}${EPHEMERIDESDIR}")
#target_link_libraries(test_TabulatedRotationalEphemeris tudat_ephemerides tudat_reference_frames tudat_input_output tudat_basic_astrodynamics tudat_basic_mathematics tudat_spice_interface ${TUDAT_EXTERNAL_LIBRARIES} ${Boost_LIBRARIES})
#
#add_executable(test_SynchronousRotationalEphemeris "${SRCROOT}${EPHEMERIDESDIR}/tests/unitTestSynchronousRotationalEphemeris.cpp")
#setup_custom_test_program(test_SynchronousRotationalEphemeris "${SRCROOT}${EPHEMERIDESDIR}")
#target_link_libraries(test_SynchronousRotationalEphemeris ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#endif()
#
#add_executable(test_KeplerEphemeris "${SRCROOT}${EPHEMERIDESDIR}/tests/unitTestKeplerEphemeris.cpp")
#setup_custom_test_program(test_KeplerEphemeris "${SRCROOT}${EPHEMERIDESDIR}")
#target_link_libraries(test_KeplerEphemeris tudat_ephemerides tudat_reference_frames tudat_input_output tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})
#
#if(USE_SOFA)
#add_executable(test_GcrsToItrsRotation "${SRCROOT}${EPHEMERIDESDIR}/tests/unitTestItrsToGcrsRotationModel.cpp")
#setup_custom_test_program(test_GcrsToItrsRotation "${SRCROOT}${EPHEMERIDESDIR}")
#target_link_libraries(test_GcrsToItrsRotation tudat_ephemerides tudat_spice_interface tudat_reference_frames tudat_earth_orientation tudat_sofa_interface tudat_basic_astrodynamics tudat_basic_mathematics tudat_input_output sofa ${Boost_LIBRARIES})
#endif( )
#
#add_executable(test_PlanetaryRotationModel "${SRCROOT}${EPHEMERIDESDIR}/tests/unitTestPlanetaryRotationModel.cpp")
#setup_custom_test_program(test_PlanetaryRotationModel "${SRCROOT}${EPHEMERIDESDIR}")
#target_link_libraries(test_PlanetaryRotationModel ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
