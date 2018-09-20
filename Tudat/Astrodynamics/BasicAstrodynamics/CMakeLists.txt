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
set(BASICASTRODYNAMICS_SOURCES
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/accelerationModelTypes.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/accelerationModel.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/attitudeElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/clohessyWiltshirePropagator.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/geodeticCoordinateConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/missionGeometry.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/modifiedEquinoctialElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/timeConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/astrodynamicsFunctions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/physicalConstants.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/bodyShapeModel.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/sphericalStateConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/unifiedStateModelQuaternionElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/unifiedStateModelModifiedRodriguesParameterElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/unifiedStateModelExponentialMapElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/stateRepresentationConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/torqueModelTypes.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/torqueModel.cpp"
)

# Set the header files.
set(BASICASTRODYNAMICS_HEADERS
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/accelerationModelTypes.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/accelerationModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/attitudeElementConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/celestialBodyConstants.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/convertMeanToEccentricAnomalies.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/clohessyWiltshirePropagator.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/customTorque.h"
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
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/massRateModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/sphericalStateConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/unifiedStateModelQuaternionElementConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/unifiedStateModelModifiedRodriguesParameterElementConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/unifiedStateModelExponentialMapElementConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/stateRepresentationConversions.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/empiricalAcceleration.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/torqueModel.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/torqueModelTypes.h"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/dissipativeTorqueModel.h"
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
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestUnifiedStateModelQuaternionsElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestUnifiedStateModelModifiedRodriguesParametersElementConversions.cpp"
  "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestUnifiedStateModelExponentialMapElementConversions.cpp"
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
setup_custom_test_program(test_MeanToEccentricAnomalyConversion "${SRCRthrustFunctionsOOT}${BASICASTRODYNAMICSDIR}")
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

add_executable(test_StateElementConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestStateConversions.cpp")
setup_custom_test_program(test_StateElementConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_StateElementConversions tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_SphericalOrbitalStateConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestSphericalOrbitStateConversions.cpp")
setup_custom_test_program(test_SphericalOrbitalStateConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_SphericalOrbitalStateConversions tudat_basic_astrodynamics tudat_reference_frames tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_BodyShapeModels "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestBodyShapeModels.cpp")
setup_custom_test_program(test_BodyShapeModels "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_BodyShapeModels tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_UnifiedStateModelQuaternionElementConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestUnifiedStateModelQuaternionElementConversions.cpp")
setup_custom_test_program(test_UnifiedStateModelQuaternionElementConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_UnifiedStateModelQuaternionElementConversions tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_UnifiedStateModelModifiedRodriguesParameterElementConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestUnifiedStateModelModifiedRodriguesParameterElementConversions.cpp")
setup_custom_test_program(test_UnifiedStateModelModifiedRodriguesParameterElementConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_UnifiedStateModelModifiedRodriguesParameterElementConversions tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_UnifiedStateModelExponentialMapElementConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestUnifiedStateModelExponentialMapElementConversions.cpp")
setup_custom_test_program(test_UnifiedStateModelExponentialMapElementConversions "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_UnifiedStateModelExponentialMapElementConversions tudat_basic_astrodynamics tudat_basic_mathematics ${Boost_LIBRARIES})

add_executable(test_EmpiricalAccelerations "${SRCROOT}${BASICASTRODYNAMICSDIR}/UnitTests/unitTestEmpiricalAcceleration.cpp")
setup_custom_test_program(test_EmpiricalAccelerations "${SRCROOT}${BASICASTRODYNAMICSDIR}")
target_link_libraries(test_EmpiricalAccelerations ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})
