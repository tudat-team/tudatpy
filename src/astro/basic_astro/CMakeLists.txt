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
set(basic_astrodynamics_SOURCES
        "accelerationModelTypes.cpp"
        "accelerationModel.cpp"
        "attitudeElementConversions.cpp"
        "clohessyWiltshirePropagator.cpp"
        "geodeticCoordinateConversions.cpp"
        "missionGeometry.cpp"
        "modifiedEquinoctialElementConversions.cpp"
        "timeConversions.cpp"
        "astrodynamicsFunctions.cpp"
        "physicalConstants.cpp"
        "bodyShapeModel.cpp"
        "sphericalStateConversions.cpp"
        "unifiedStateModelQuaternionElementConversions.cpp"
        "unifiedStateModelModifiedRodriguesParameterElementConversions.cpp"
        "unifiedStateModelExponentialMapElementConversions.cpp"
        "stateRepresentationConversions.cpp"
        "torqueModelTypes.cpp"
        "torqueModel.cpp"
        )

# Set the header files.
set(basic_astrodynamics_HEADERS
        "accelerationModelTypes.h"
        "accelerationModel.h"
        "attitudeElementConversions.h"
        "celestialBodyConstants.h"
        "convertMeanToEccentricAnomalies.h"
        "clohessyWiltshirePropagator.h"
        "customTorque.h"
        "geodeticCoordinateConversions.h"
        "keplerPropagator.h"
        "missionGeometry.h"
        "modifiedEquinoctialElementConversions.h"
        "stateVectorIndices.h"
        "timeConversions.h"
        "tests/testAccelerationModels.h"
        "tests/testBody.h"
        "tests/keplerPropagatorTestData.h"
        "astrodynamicsFunctions.h"
        "orbitalElementConversions.h"
        "physicalConstants.h"
        "unitConversions.h"
        "bodyShapeModel.h"
        "oblateSpheroidBodyShapeModel.h"
        "sphericalBodyShapeModel.h"
        "massRateModel.h"
        "sphericalStateConversions.h"
        "unifiedStateModelQuaternionElementConversions.h"
        "unifiedStateModelModifiedRodriguesParameterElementConversions.h"
        "unifiedStateModelExponentialMapElementConversions.h"
        "stateRepresentationConversions.h"
        "empiricalAcceleration.h"
        "torqueModel.h"
        "torqueModelTypes.h"
        "dissipativeTorqueModel.h"
        )


TUDAT_ADD_LIBRARY("basic_astrodynamics"
        "${basic_astrodynamics_SOURCES}"
        "${basic_astrodynamics_HEADERS}"
        PUBLIC_LINKS tudat_basic_mathematics
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )

# Add tests.
TUDAT_ADD_TEST_CASE(AstrodynamicsFunctions
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(OrbitalElementConversions
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(UnitConversions
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(PhysicalConstants
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(ConvertMeanAnomalyToEccentricAnomaly
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        tudat_root_finders
        )

TUDAT_ADD_TEST_CASE(ConvertMeanAnomalyToHyperbolicEccentricAnomaly
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        tudat_root_finders
        )

TUDAT_ADD_TEST_CASE(KeplerPropagator
        PRIVATE_LINKS
        tudat_input_output
        tudat_gravitation
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        tudat_root_finders
        )

TUDAT_ADD_TEST_CASE(AccelerationModel
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(ClohessyWiltshirePropagator
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(MissionGeometry
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(TimeConversions
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(CelestialBodyConstants
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(ModifiedEquinoctialElementConversions
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(GeodeticCoordinateConversions
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(StateConversions
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(SphericalOrbitStateConversions
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        tudat_reference_frames
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(BodyShapeModels
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(UnifiedStateModelQuaternionElementConversions
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(UnifiedStateModelMRPElementConversions
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(UnifiedStateModelEMElementConversions
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(EmpiricalAcceleration
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        tudat_electromagnetism
        ${CSpice_LIBRARIES}
        ${Boost_LIBRARIES}
        ${Sofa_LIBRARIES}
        )

