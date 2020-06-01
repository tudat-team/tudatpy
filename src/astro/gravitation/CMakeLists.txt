#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#
#    Notes
#
#

# Set the source files.
set(gravitation_SOURCES
        "basicSolidBodyTideGravityFieldVariations.cpp"
        "gravityFieldVariations.cpp"
        "centralGravityModel.cpp"
        "centralJ2GravityModel.cpp"
        "centralJ2J3GravityModel.cpp"
        "centralJ2J3J4GravityModel.cpp"
        "gravityFieldModel.cpp"
        "jacobiEnergy.cpp"
        "librationPoint.cpp"
        "sphericalHarmonicsGravityModel.cpp"
        "sphericalHarmonicsGravityField.cpp"
        "thirdBodyPerturbation.cpp"
        "timeDependentSphericalHarmonicsGravityField.cpp"
        "unitConversionsCircularRestrictedThreeBodyProblem.cpp"
        "tests/planetTestData.cpp"
        "unitConversionsCircularRestrictedThreeBodyProblem.cpp"
        "triAxialEllipsoidGravity.cpp"
        "tabulatedGravityFieldVariations.cpp"
        "mutualSphericalHarmonicGravityModel.cpp"
        "secondDegreeGravitationalTorque.cpp"
        "directTidalDissipationAcceleration.cpp"
        )

# Set the header files.
set(gravitation_HEADERS
        "basicSolidBodyTideGravityFieldVariations.h"
        "gravityFieldVariations.h"
        "centralGravityModel.h"
        "centralJ2GravityModel.h"
        "centralJ2J3GravityModel.h"
        "centralJ2J3J4GravityModel.h"
        "gravityFieldModel.h"
        "jacobiEnergy.h"
        "librationPoint.h"
        "sphericalHarmonicsGravityModel.h"
        "sphericalHarmonicsGravityModelBase.h"
        "sphericalHarmonicsGravityField.h"
        "thirdBodyPerturbation.h"
        "timeDependentSphericalHarmonicsGravityField.h"
        "unitConversionsCircularRestrictedThreeBodyProblem.h"
        "tests/planetTestData.h"
        "triAxialEllipsoidGravity.h"
        "tabulatedGravityFieldVariations.h"
        "mutualSphericalHarmonicGravityModel.h"
        "secondDegreeGravitationalTorque.h"
        "directTidalDissipationAcceleration.h"
        "sphericalHarmonicGravitationalTorque.h"
        )

TUDAT_ADD_LIBRARY("gravitation"
        "${gravitation_SOURCES}"
        "${gravitation_HEADERS}"
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )

# Add static libraries.
#add_library(tudat_gravitation STATIC ${GRAVITATION_SOURCES} ${GRAVITATION_HEADERS})
#setup_tudat_library_target(tudat_gravitation "${SRCROOT}${GRAVITATIONDIR}")
#
## Add unit tests.
#add_executable(test_SphericalHarmonicsGravityField "${SRCROOT}${GRAVITATIONDIR}/tests/unitTestSphericalHarmonicsGravityField.cpp")
#setup_custom_test_program(test_SphericalHarmonicsGravityField "${SRCROOT}${GRAVITATIONDIR}")
#target_link_libraries(test_SphericalHarmonicsGravityField tudat_gravitation tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})
#
#add_executable(test_GravitationalForce "${SRCROOT}${GRAVITATIONDIR}/tests/unitTestGravitationalForce.cpp")
#setup_custom_test_program(test_GravitationalForce "${SRCROOT}${GRAVITATIONDIR}")
#target_link_libraries(test_GravitationalForce tudat_gravitation tudat_basic_astrodynamics ${Boost_LIBRARIES})
#
#add_executable(test_CentralAndZonalGravityModel "${SRCROOT}${GRAVITATIONDIR}/tests/unitTestCentralAndZonalGravityModel.cpp")
#setup_custom_test_program(test_CentralAndZonalGravityModel "${SRCROOT}${GRAVITATIONDIR}")
#target_link_libraries(test_CentralAndZonalGravityModel tudat_gravitation tudat_basic_astrodynamics tudat_basic_astrodynamics tudat_numerical_integrators ${Boost_LIBRARIES})
#
#add_executable(test_UnitConversionsCircularRestrictedThreeBodyProblem "${SRCROOT}${GRAVITATIONDIR}/tests/unitTestUnitConversionsCircularRestrictedThreeBodyProblem.cpp")
#setup_custom_test_program(test_UnitConversionsCircularRestrictedThreeBodyProblem "${SRCROOT}${GRAVITATIONDIR}")
#target_link_libraries(test_UnitConversionsCircularRestrictedThreeBodyProblem tudat_gravitation ${Boost_LIBRARIES})
#
#add_executable(test_JacobiEnergy "${SRCROOT}${GRAVITATIONDIR}/tests/unitTestJacobiEnergy.cpp")
#setup_custom_test_program(test_JacobiEnergy "${SRCROOT}${GRAVITATIONDIR}")
#target_link_libraries(test_JacobiEnergy tudat_gravitation tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_input_output ${Boost_LIBRARIES})
#
#add_executable(test_LibrationPoints "${SRCROOT}${GRAVITATIONDIR}/tests/unitTestLibrationPoints.cpp")
#setup_custom_test_program(test_LibrationPoints "${SRCROOT}${GRAVITATIONDIR}")
#target_link_libraries(test_LibrationPoints tudat_ephemerides tudat_root_finders tudat_gravitation tudat_input_output ${Boost_LIBRARIES})
#
#add_executable(test_SphericalHarmonicsGravityModel "${SRCROOT}${GRAVITATIONDIR}/tests/unitTestSphericalHarmonicsGravityModel.cpp")
#setup_custom_test_program(test_SphericalHarmonicsGravityModel "${SRCROOT}${GRAVITATIONDIR}")
#target_link_libraries(test_SphericalHarmonicsGravityModel tudat_gravitation tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})
#
#add_executable(test_ThirdBodyPerturbation "${SRCROOT}${GRAVITATIONDIR}/tests/unitTestThirdBodyPerturbation.cpp")
#setup_custom_test_program(test_ThirdBodyPerturbation "${SRCROOT}${GRAVITATIONDIR}")
#target_link_libraries(test_ThirdBodyPerturbation tudat_gravitation tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})
#
#add_executable(test_TriAxialEllipsoidGravity "${SRCROOT}${GRAVITATIONDIR}/tests/unitTestTriAxialEllipsoidGravity.cpp")
#setup_custom_test_program(test_TriAxialEllipsoidGravity "${SRCROOT}${GRAVITATIONDIR}")
#target_link_libraries(test_TriAxialEllipsoidGravity tudat_gravitation tudat_basic_mathematics ${Boost_LIBRARIES})
#
#add_executable(test_DirectTidalDissipationAcceleration "${SRCROOT}${GRAVITATIONDIR}/tests/unitTestDirectTidalDissipation.cpp")
#setup_custom_test_program(test_DirectTidalDissipationAcceleration "${SRCROOT}${GRAVITATIONDIR}")
#target_link_libraries(test_DirectTidalDissipationAcceleration ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#
#if (TUDAT_BUILD_WITH_SPICE_INTERFACE)
#    add_executable(test_GravityFieldVariations "${SRCROOT}${GRAVITATIONDIR}/tests/unitTestGravityFieldVariations.cpp")
#    setup_custom_test_program(test_GravityFieldVariations "${SRCROOT}${GRAVITATIONDIR}")
#    target_link_libraries(test_GravityFieldVariations tudat_gravitation tudat_interpolators tudat_basic_mathematics tudat_spice_interface ${SPICE_LIBRARIES} ${Boost_LIBRARIES})
#
#    add_executable(test_MutualSphericalHarmonicsGravityModel "${SRCROOT}${GRAVITATIONDIR}/tests/unitTestMutualSphericalHarmonicAcceleration.cpp")
#    setup_custom_test_program(test_MutualSphericalHarmonicsGravityModel "${SRCROOT}${GRAVITATIONDIR}")
#    target_link_libraries(test_MutualSphericalHarmonicsGravityModel ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#    add_executable(test_GravitationalTorques "${SRCROOT}${GRAVITATIONDIR}/tests/unitTestGravitationalTorques.cpp")
#    setup_custom_test_program(test_GravitationalTorques "${SRCROOT}${GRAVITATIONDIR}")
#    target_link_libraries(test_GravitationalTorques ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
#
#endif ()
