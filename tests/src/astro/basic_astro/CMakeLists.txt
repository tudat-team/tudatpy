#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

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
        tudat_input_output
        )

TUDAT_ADD_TEST_CASE(ConvertMeanAnomalyToHyperbolicEccentricAnomaly
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        tudat_root_finders
        tudat_input_output
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

TUDAT_ADD_TEST_CASE(CustomAccelerationAndTorqueModel
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
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
        )

TUDAT_ADD_TEST_CASE(PolyhedronFunctions
        PRIVATE_LINKS
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

