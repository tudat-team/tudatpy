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
        # For .so
        tudat_basic_mathematics
        )

# Add tests.
TUDAT_ADD_TESTCASE(AstrodynamicsFunctions
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(OrbitalElementConversions
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(UnitConversions
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(PhysicalConstants
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(ConvertMeanAnomalyToEccentricAnomaly
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        tudat_root_finders
        )

TUDAT_ADD_TESTCASE(ConvertMeanAnomalyToHyperbolicEccentricAnomaly
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        tudat_root_finders
        )

TUDAT_ADD_TESTCASE(KeplerPropagator
        tudat_input_output
        tudat_gravitation
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        tudat_root_finders
        )

TUDAT_ADD_TESTCASE(AccelerationModel
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TESTCASE(ClohessyWiltshirePropagator
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TESTCASE(MissionGeometry
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TESTCASE(TimeConversions
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TESTCASE(CelestialBodyConstants
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TESTCASE(ModifiedEquinoctialElementConversions
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(GeodeticCoordinateConversions
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(StateConversions
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(SphericalOrbitStateConversions
        tudat_basic_astrodynamics
        tudat_reference_frames
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(BodyShapeModels
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(UnifiedStateModelQuaternionElementConversions
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(UnifiedStateModelMRPElementConversions
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(UnifiedStateModelEMElementConversions
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(EmpiricalAcceleration
        ${Tudat_PROPAGATION_LIBRARIES}
        tudat_electromagnetism
        ${CSpice_LIBRARIES}
        ${Boost_LIBRARIES}
        ${Sofa_LIBRARIES}
        )

