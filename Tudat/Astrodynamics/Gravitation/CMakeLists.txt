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
 #      120321    K. Kumar          Added header, source and unit test files for state-derivative
 #                                  of CRTBP.
 #

# Set the source files.
set(GRAVITATION_SOURCES
  "${SRCROOT}${GRAVITATIONDIR}/centralGravityField.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/gravitationalAccelerationModel.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/gravitationalForceModel.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/jacobiEnergy.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/librationPoint.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/sphericalHarmonicsGravityField.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/stateDerivativeCircularRestrictedThreeBodyProblem.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/unitConversionsCircularRestrictedThreeBodyProblem.cpp"
)

# Set the header files.
set(GRAVITATION_HEADERS 
  "${SRCROOT}${GRAVITATIONDIR}/gravityFieldModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/centralGravityField.h"
  "${SRCROOT}${GRAVITATIONDIR}/gravitationalAccelerationModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/gravitationalForceModel.h"
  "${SRCROOT}${GRAVITATIONDIR}/jacobiEnergy.h"
  "${SRCROOT}${GRAVITATIONDIR}/librationPoint.h"
  "${SRCROOT}${GRAVITATIONDIR}/sphericalHarmonicsGravityField.h"
  "${SRCROOT}${GRAVITATIONDIR}/stateDerivativeCircularRestrictedThreeBodyProblem.h"
  "${SRCROOT}${GRAVITATIONDIR}/unitConversionsCircularRestrictedThreeBodyProblem.h"
)

# Add static libraries.
add_library(tudat_gravitation STATIC ${GRAVITATION_SOURCES} ${GRAVITATION_HEADERS})
setup_tudat_library_target(tudat_gravitation "${SRCROOT}${GRAVITATIONDIR}")

# Add unit tests.
add_executable(test_SphericalHarmonicsGravityField "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestSphericalHarmonicsGravityField.cpp")
setup_custom_test_program(test_SphericalHarmonicsGravityField "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_SphericalHarmonicsGravityField tudat_gravitation tudat_bodies ${Boost_LIBRARIES})

add_executable(test_GravitationalForce "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestGravitationalForce.cpp")
setup_custom_test_program(test_GravitationalForce "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_GravitationalForce tudat_gravitation tudat_bodies ${Boost_LIBRARIES})

add_executable(test_GravitationalAcceleration "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestGravitationalAcceleration.cpp")
setup_custom_test_program(test_GravitationalAcceleration "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_GravitationalAcceleration tudat_gravitation tudat_bodies ${Boost_LIBRARIES})

add_executable(test_UnitConversionsCircularRestrictedThreeBodyProblem "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestUnitConversionsCircularRestrictedThreeBodyProblem.cpp")
setup_custom_test_program(test_UnitConversionsCircularRestrictedThreeBodyProblem "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_UnitConversionsCircularRestrictedThreeBodyProblem tudat_gravitation ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_JacobiEnergy "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestJacobiEnergy.cpp")
setup_custom_test_program(test_JacobiEnergy "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_JacobiEnergy tudat_gravitation tudat_bodies tudat_mission_segments tudat_root_finding_methods tudat_basic_astrodynamics tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_LibrationPoints "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestLibrationPoints.cpp")
setup_custom_test_program(test_LibrationPoints "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_LibrationPoints tudat_bodies tudat_mission_segments tudat_root_finding_methods tudat_gravitation tudat_basic_astrodynamics tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

add_executable(test_StateDerivativeRestrictedThreeBodyProblem "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestStateDerivativeCircularRestrictedThreeBodyProblem.cpp")
setup_custom_test_program(test_StateDerivativeRestrictedThreeBodyProblem "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_StateDerivativeRestrictedThreeBodyProblem tudat_bodies tudat_mission_segments tudat_root_finding_methods tudat_gravitation tudat_basic_astrodynamics tudat_input_output ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
