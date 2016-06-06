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
 #	120203    K. Kumar	    Adapted for new BasicAstrodynamics package.
 #      120216    K. Kumar          Removed propagator files.
 #      120822    P. Musegaas       Added unit test for mean to hyperbolic eccentric anomaly
 #                                  conversion.
 #      130218    D. Dirkx          Added unit test for Julian date conversions.
 #      130214    R.C.A. Boon       Added modified equinoctial element files (.h/.cpp, unit test).
 #      130301    D. Dirkx          Added unit test for geodetic coordinate conversions.
 #      160406    M. Van den Broeck Added Unified State Model element files (.h/.cpp, unit test).
 #
 #    References
 #
 #    Notes
 #

# Set the source files.
set(BASICASTRODYNAMICS_SOURCES
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/accelerationModelTypes.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/clohessyWiltshirePropagator.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/geodeticCoordinateConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/missionGeometry.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/modifiedEquinoctialElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/timeConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/astrodynamicsFunctions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/physicalConstants.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/bodyShapeModel.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/unifiedStateModelElementConversions.cpp"
)

# Set the header files.
set(BASICASTRODYNAMICS_HEADERS
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/accelerationModelTypes.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/accelerationModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/celestialBodyConstants.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/convertMeanToEccentricAnomalies.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/clohessyWiltshirePropagator.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/geodeticCoordinateConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/keplerPropagator.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/missionGeometry.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/modifiedEquinoctialElementConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/stateVectorIndices.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/timeConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/testAccelerationModels.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/testBody.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/keplerPropagatorTestData.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/astrodynamicsFunctions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/orbitalElementConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/physicalConstants.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/unitConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/bodyShapeModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/oblateSpheroidBodyShapeModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/sphericalBodyShapeModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/unifiedStateModelElementConversions.h"
)

# Add static libraries.
add_library(tudat_basic_astrodynamics STATIC ${BASICASTRODYNAMICS_SOURCES} ${BASICASTRODYNAMICS_HEADERS})
setup_tudat_library_target(tudat_basic_astrodynamics "${SRCROOT}${BASICASTRODYNAMICSDIR}")

# Add unit tests.
set(BASICASTRODYNAMICS_UNITTESTS
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestBasicAstrodynamics.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestAstrodynamicsFunctions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestOrbitalElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestPhysicalConstants.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestUnitConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestUnifiedStateModelElementConversions.cpp"
)

add_executable(test_AstrodynamicsFunctions "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestAstrodynamicsFunctions.cpp")
setup_custom_test_program(test_AstrodynamicsFunctions "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_AstrodynamicsFunctions tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_OrbitalElementConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestOrbitalElementConversions.cpp")
setup_custom_test_program(test_OrbitalElementConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_OrbitalElementConversions tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_PhysicalConstants "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestPhysicalConstants.cpp")
setup_custom_test_program(test_PhysicalConstants "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_PhysicalConstants tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_UnitConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestUnitConversions.cpp")
setup_custom_test_program(test_UnitConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_UnitConversions tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_MeanToEccentricAnomalyConversion "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestConvertMeanAnomalyToEccentricAnomaly.cpp")
setup_custom_test_program(test_MeanToEccentricAnomalyConversion "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_MeanToEccentricAnomalyConversion tudat_basic_astrodynamics tudat_basic_mathematics tudat_root_finders ${Boost_LIBRARIES})

add_executable(test_MeanToHyperbolicEccentricAnomalyConversion "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestConvertMeanAnomalyToHyperbolicEccentricAnomaly.cpp")
setup_custom_test_program(test_MeanToHyperbolicEccentricAnomalyConversion "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_MeanToHyperbolicEccentricAnomalyConversion tudat_basic_astrodynamics tudat_basic_mathematics tudat_root_finders ${Boost_LIBRARIES})

add_executable(test_KeplerPropagator "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestKeplerPropagator.cpp")
setup_custom_test_program(test_KeplerPropagator "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_KeplerPropagator tudat_input_output tudat_gravitation tudat_basic_astrodynamics tudat_basic_mathematics tudat_root_finders ${Boost_LIBRARIES})

add_executable(test_AccelerationModel "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestAccelerationModel.cpp")
setup_custom_test_program(test_AccelerationModel "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_AccelerationModel tudat_basic_astrodynamics ${Boost_LIBRARIES})

add_executable(test_ClohessyWiltshirePropagator "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestClohessyWiltshirePropagator.cpp")
setup_custom_test_program(test_ClohessyWiltshirePropagator "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_ClohessyWiltshirePropagator tudat_basic_astrodynamics ${Boost_LIBRARIES})

add_executable(test_MissionGeometry "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestMissionGeometry.cpp")
setup_custom_test_program(test_MissionGeometry "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_MissionGeometry tudat_basic_astrodynamics ${Boost_LIBRARIES})

add_executable(test_TimeConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestTimeConversions.cpp")
setup_custom_test_program(test_TimeConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_TimeConversions tudat_basic_astrodynamics ${Boost_LIBRARIES})

add_executable(test_CelestialBodyConstants "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestCelestialBodyConstants.cpp")
setup_custom_test_program(test_CelestialBodyConstants "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_CelestialBodyConstants tudat_basic_astrodynamics ${Boost_LIBRARIES})

add_executable(test_ModifiedEquinoctialElementConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestModifiedEquinoctialElementConversions.cpp")
setup_custom_test_program(test_ModifiedEquinoctialElementConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_ModifiedEquinoctialElementConversions tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_GeodeticCoordinateConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestGeodeticCoordinateConversions.cpp")
setup_custom_test_program(test_GeodeticCoordinateConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_GeodeticCoordinateConversions tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_BodyShapeModels "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestBodyShapeModels.cpp")
setup_custom_test_program(test_BodyShapeModels "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_BodyShapeModels tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_UnifiedStateModelElementConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestUnifiedStateModelElementConversions.cpp")
setup_custom_test_program(test_UnifiedStateModelElementConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_UnifiedStateModelElementConversions tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})
