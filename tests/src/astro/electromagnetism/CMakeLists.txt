#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

TUDAT_ADD_TEST_CASE(CannonBallRadiationPressureAccelerationAndForce
        PRIVATE_LINKS
        tudat_electromagnetism
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(LorentzStaticMagneticAccelerationAndForce
        PRIVATE_LINKS
        tudat_electromagnetism
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(RadiationPressureInterface
        PRIVATE_LINKS
        tudat_electromagnetism
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(PanelledRadiationPressure
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(SolarSailAccelerationAndForce
        PRIVATE_LINKS
        tudat_electromagnetism
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(RadiationSourceModel
        PRIVATE_LINKS
        tudat_electromagnetism
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(LuminosityModel
        PRIVATE_LINKS
        tudat_electromagnetism
        )

TUDAT_ADD_TEST_CASE(SourcePanelRadiosityModel
        PRIVATE_LINKS
        tudat_electromagnetism
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(RadiationPressureTargetModel
        PRIVATE_LINKS
        tudat_electromagnetism
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(RadiationPressureAcceleration
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(ReflectionLaw
        PRIVATE_LINKS
        tudat_electromagnetism
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(OccultationModel
        PRIVATE_LINKS
        tudat_electromagnetism
        tudat_basic_mathematics
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(SurfacePropertyDistribution
        PRIVATE_LINKS
        tudat_electromagnetism
        tudat_basic_mathematics
        tudat_basic_astrodynamics
        tudat_spice_interface
        )
