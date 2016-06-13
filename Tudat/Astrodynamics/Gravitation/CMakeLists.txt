 #    Copyright (c) 2010-2015, Delft University of Technology
 #    All rights reserved.
 #
 #    Redistribution and use in source and binary forms, with or without modification, are
 #    permitted provided that the following conditions are met:
 #      - Redistributions of source code must retain the above copyright notice, this list of
 #        conditions and the following disclaimer.
 #      - Redistributions in binary form must reproduce the above copyright notice, this list of
 #        conditions and the following disclaimer in the documentation and/or other materials
 #        provided with the distribution.
 #      - Neither the name of the Delft University of Technology nor the names of its contributors
 #        may be used to endorse or promote products derived from this software without specific
 #        prior written permission.
 #
 #    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 #    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 #    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 #    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 #    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 #    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 #    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 #    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 #    OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      110820    S.M. Persson      File created.
 #      111025    K. Kumar          Adapted file to work with Revision 194.
 #      111026    K. Kumar          Adapted file so all headers show in project tree in Qt Creator.
 #      120321    K. Kumar          Added header, source and unit test files for state-derivative
 #                                  of CRTBP.
 #      121017    E. Dekens         Added sphericalHarmonicsGravitationalAccelerationModel files.
 #      150205                      Last modified
 #    References
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
  "${SRCROOT}${GRAVITATIONDIR}/stateDerivativeCircularRestrictedThreeBodyProblem.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/thirdBodyPerturbation.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/timeDependentSphericalHarmonicsGravityField.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/unitConversionsCircularRestrictedThreeBodyProblem.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/UnitTests/planetTestData.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/unitConversionsCircularRestrictedThreeBodyProblem.cpp"
  "${SRCROOT}${GRAVITATIONDIR}/tabulatedGravityFieldVariations.cpp"
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
  "${SRCROOT}${GRAVITATIONDIR}/stateDerivativeCircularRestrictedThreeBodyProblem.h"
  "${SRCROOT}${GRAVITATIONDIR}/thirdBodyPerturbation.h"
  "${SRCROOT}${GRAVITATIONDIR}/timeDependentSphericalHarmonicsGravityField.h"
  "${SRCROOT}${GRAVITATIONDIR}/unitConversionsCircularRestrictedThreeBodyProblem.h"
  "${SRCROOT}${GRAVITATIONDIR}/UnitTests/planetTestData.h"
  "${SRCROOT}${GRAVITATIONDIR}/tabulatedGravityFieldVariations.h"
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

add_executable(test_StateDerivativeRestrictedThreeBodyProblem "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestStateDerivativeCircularRestrictedThreeBodyProblem.cpp")
setup_custom_test_program(test_StateDerivativeRestrictedThreeBodyProblem "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_StateDerivativeRestrictedThreeBodyProblem tudat_mission_segments tudat_root_finders tudat_gravitation tudat_basic_astrodynamics tudat_input_output ${Boost_LIBRARIES})

add_executable(test_SphericalHarmonicsGravityModel "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestSphericalHarmonicsGravityModel.cpp")
setup_custom_test_program(test_SphericalHarmonicsGravityModel "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_SphericalHarmonicsGravityModel tudat_gravitation tudat_basic_mathematics ${Boost_LIBRARIES} )

add_executable(test_ThirdBodyPerturbation "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestThirdBodyPerturbation.cpp")
setup_custom_test_program(test_ThirdBodyPerturbation "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_ThirdBodyPerturbation tudat_gravitation tudat_basic_mathematics ${Boost_LIBRARIES} )

if(USE_CSPICE)
add_executable(test_GravityFieldVariations "${SRCROOT}${GRAVITATIONDIR}/UnitTests/unitTestGravityFieldVariations.cpp")
setup_custom_test_program(test_GravityFieldVariations "${SRCROOT}${GRAVITATIONDIR}")
target_link_libraries(test_GravityFieldVariations tudat_gravitation tudat_basic_mathematics tudat_spice_interface ${SPICE_LIBRARIES} ${Boost_LIBRARIES} )
endif()
