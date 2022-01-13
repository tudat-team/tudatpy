#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

TUDAT_ADD_TEST_CASE(ApproximatePlanetPositions
        PRIVATE_LINKS
        tudat_ephemerides
        tudat_gravitation
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        tudat_input_output
        tudat_root_finders
        )

TUDAT_ADD_TEST_CASE(TabulatedEphemeris
        PRIVATE_LINKS
        tudat_ephemerides
        tudat_gravitation
        tudat_basic_astrodynamics
        tudat_input_output
        tudat_interpolators
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(CartesianStateExtractor
        PRIVATE_LINKS
        tudat_ephemerides
        tudat_gravitation
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        tudat_input_output
        tudat_interpolators
        )


TUDAT_ADD_TEST_CASE(KeplerStateExtractor
        PRIVATE_LINKS
        tudat_ephemerides
        tudat_input_output
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        tudat_root_finders
        )

TUDAT_ADD_TEST_CASE(RotationalEphemeris
        PRIVATE_LINKS
        tudat_ephemerides
        tudat_reference_frames
        tudat_input_output
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        tudat_spice_interface
        )

TUDAT_ADD_TEST_CASE(SimpleRotationalEphemeris
        PRIVATE_LINKS
        tudat_ephemerides
        tudat_reference_frames
        tudat_input_output
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(KeplerEphemeris
        PRIVATE_LINKS
        tudat_ephemerides
        tudat_reference_frames
        tudat_input_output
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(PlanetaryRotationModel
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(TabulatedRotationalEphemeris
        PRIVATE_LINKS
        tudat_ephemerides
        tudat_reference_frames
        tudat_input_output
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        tudat_spice_interface
        )

TUDAT_ADD_TEST_CASE(CompositeEphemeris
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(FrameManager
        PRIVATE_LINKS
        tudat_ephemerides
        tudat_reference_frames
        tudat_spice_interface
        )

TUDAT_ADD_TEST_CASE(SynchronousRotationalEphemeris
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(TwoLineElementsEphemeris
        PRIVATE_LINKS
       ${Tudat_PROPAGATION_LIBRARIES}
        )

if(TUDAT_BUILD_WITH_SOFA_INTERFACE)

    TUDAT_ADD_TEST_CASE(ItrsToGcrsRotationModel
            PRIVATE_LINKS
            tudat_ephemerides
            tudat_spice_interface
            tudat_reference_frames
            tudat_earth_orientation
            tudat_sofa_interface
            tudat_basic_astrodynamics
            tudat_basic_mathematics
            tudat_input_output
            )

endif( )
