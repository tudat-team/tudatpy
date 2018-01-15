 #    Copyright (c) 2010-2018, Delft University of Technology
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
set(GRAVITATION_SOURCES
  "${SRCROOT}${GRAVITATIONDIR}/basicSolidBodyTideGravityFieldVariations.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/gravityFieldVariations.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/centralGravityModel.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/centralJ2GravityModel.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/centralJ2J3GravityModel.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/centralJ2J3J4GravityModel.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/gravityFieldModel.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/jacobiEnergy.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/librationPoint.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/sphericalHarmonicsGravityModel.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/sphericalHarmonicsGravityField.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/thirdBodyPerturbation.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/timeDependentSphericalHarmonicsGravityField.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/unitConversionsCircularRestrictedThreeBodyProblem.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/UnitTests/planetTestData.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/unitConversionsCircularRestrictedThreeBodyProblem.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/triAxialEllipsoidGravity.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/tabulatedGravityFieldVariations.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/mutualSphericalHarmonicGravityModel.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/secondDegreeGravitationalTorque.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/directTidalDissipationAcceleration.cpp"
)

# Set the header files.
set(GRAVITATION_HEADERS
  "${SRCROOT}${GRAVITATIONDIR}/basicSolidBodyTideGravityFieldVariations.h"
  "${SRCROOT}${GRAVITATIONDIR}/gravityFieldVariations.h"
  "${SRCROOT}${GRAVITATIONDIR}/centralGravityModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/centralJ2GravityModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/centralJ2J3GravityModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/centralJ2J3J4GravityModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/gravityFieldModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/jacobiEnergy.h"
  "${SRCROOT}${GRAVITATIONDIR}/librationPoint.h"
  "${SRCROOT}${GRAVITATIONDIR}/sphericalHarmonicsGravityModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/sphericalHarmonicsGravityModelBase.h"
  "${SRCROOT}${GRAVITATIONDIR}/sphericalHarmonicsGravityField.h"
  "${SRCROOT}${GRAVITATIONDIR}/thirdBodyPerturbation.h"
  "${SRCROOT}${GRAVITATIONDIR}/timeDependentSphericalHarmonicsGravityField.h"
  "${SRCROOT}${GRAVITATIONDIR}/unitConversionsCircularRestrictedThreeBodyProblem.h"
  "${SRCROOT}${GRAVITATIONDIR}/UnitTests/planetTestData.h"
  "${SRCROOT}${GRAVITATIONDIR}/triAxialEllipsoidGravity.h"
  "${SRCROOT}${GRAVITATIONDIR}/tabulatedGravityFieldVariations.h"
  "${SRCROOT}${GRAVITATIONDIR}/mutualSphericalHarmonicGravityModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/secondDegreeGravitationalTorque.h"
  "${SRCROOT}${GRAVITATIONDIR}/directTidalDissipationAcceleration.h"
)

# Add static libraries.
add_library(tudat_gravitation STATIC ${GRAVITATION_SOURCES} ${GRAVITATION_HEADERS})
setup_tudat_library_target(tudat_gravitation "${SRCROOT}${GRAVITATIONDIR}")

# Add unit tests.
add_executable(test_SphericalHarmonicsGravityField "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestSphericalHarmonicsGravityField.cpp")
setup_custom_test_program(test_SphericalHarmonicsGravityField "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_SphericalHarmonicsGravityField tudat_gravitation tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_GravitationalForce "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestGravitationalForce.cpp")
setup_custom_test_program(test_GravitationalForce "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_GravitationalForce tudat_gravitation ${Boost_LIBRARIES})

add_executable(test_CentralAndZonalGravityModel "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestCentralAndZonalGravityModel.cpp")
setup_custom_test_program(test_CentralAndZonalGravityModel "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_CentralAndZonalGravityModel tudat_gravitation tudat_numerical_integrators ${Boost_LIBRARIES})

add_executable(test_UnitConversionsCircularRestrictedThreeBodyProblem "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestUnitConversionsCircularRestrictedThreeBodyProblem.cpp")
setup_custom_test_program(test_UnitConversionsCircularRestrictedThreeBodyProblem "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_UnitConversionsCircularRestrictedThreeBodyProblem tudat_gravitation ${Boost_LIBRARIES})

add_executable(test_JacobiEnergy "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestJacobiEnergy.cpp")
setup_custom_test_program(test_JacobiEnergy "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_JacobiEnergy tudat_gravitation tudat_mission_segments tudat_root_finders tudat_basic_astrodynamics tudat_input_output ${Boost_LIBRARIES})

add_executable(test_LibrationPoints "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestLibrationPoints.cpp")
setup_custom_test_program(test_LibrationPoints "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_LibrationPoints tudat_ephemerides tudat_root_finders tudat_gravitation tudat_basic_astrodynamics tudat_input_output ${Boost_LIBRARIES})

add_executable(test_SphericalHarmonicsGravityModel "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestSphericalHarmonicsGravityModel.cpp")
setup_custom_test_program(test_SphericalHarmonicsGravityModel "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_SphericalHarmonicsGravityModel tudat_gravitation tudat_basic_mathematics ${Boost_LIBRARIES} )

add_executable(test_ThirdBodyPerturbation "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestThirdBodyPerturbation.cpp")
setup_custom_test_program(test_ThirdBodyPerturbation "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_ThirdBodyPerturbation tudat_gravitation tudat_basic_mathematics ${Boost_LIBRARIES} )

add_executable(test_TriAxialEllipsoidGravity "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestTriAxialEllipsoidGravity.cpp")
setup_custom_test_program(test_TriAxialEllipsoidGravity "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_TriAxialEllipsoidGravity tudat_gravitation tudat_basic_mathematics ${Boost_LIBRARIES} )

add_executable(test_DirectTidalDissipationAcceleration "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestDirectTidalDissipation.cpp")
setup_custom_test_program(test_DirectTidalDissipationAcceleration "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_DirectTidalDissipationAcceleration ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

if(USE_CSPICE)
add_executable(test_GravityFieldVariations "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestGravityFieldVariations.cpp")
setup_custom_test_program(test_GravityFieldVariations "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_GravityFieldVariations tudat_gravitation tudat_basic_mathematics tudat_spice_interface ${SPICE_LIBRARIES} ${Boost_LIBRARIES} )

add_executable(test_MutualSphericalHarmonicsGravityModel "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestMutualSphericalHarmonicAcceleration.cpp")
setup_custom_test_program(test_MutualSphericalHarmonicsGravityModel "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_MutualSphericalHarmonicsGravityModel ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
endif()
