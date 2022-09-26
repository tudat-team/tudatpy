#    Copyright (c) 2010-2019, Delft University of Technology
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

TUDAT_ADD_TEST_CASE(SphericalHarmonicsGravityField
        PRIVATE_LINKS
        tudat_gravitation
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(CentralGravityModel
        PRIVATE_LINKS
        tudat_gravitation
        tudat_basic_astrodynamics
        tudat_basic_mathematics
        )

#TUDAT_ADD_TEST_CASE(JacobiEnergy
#        PRIVATE_LINKS
#        tudat_gravitation
#        tudat_mission_segments
#        tudat_root_finders
#        tudat_basic_astrodynamics
#        tudat_input_output
#        )
#
#TUDAT_ADD_TEST_CASE(UnitConversionsCircularRestrictedThreeBodyProblem
#        PRIVATE_LINKS
#        tudat_gravitation
#        )
#
#TUDAT_ADD_TEST_CASE(LibrationPoints
#    PRIVATE_LINKS
#    tudat_ephemerides tudat_root_finders tudat_gravitation tudat_input_output )

TUDAT_ADD_TEST_CASE(SphericalHarmonicsGravityModel
    PRIVATE_LINKS
    tudat_gravitation tudat_basic_astrodynamics tudat_basic_mathematics)

#TUDAT_ADD_TEST_CASE(ThirdBodyPerturbation
#    PRIVATE_LINKS
#    tudat_gravitation tudat_basic_astrodynamics tudat_basic_mathematics)
#
#TUDAT_ADD_TEST_CASE(TriAxialEllipsoidGravity
#    PRIVATE_LINKS
#    tudat_gravitation tudat_basic_mathematics)
#
#TUDAT_ADD_TEST_CASE(DirectTidalDissipationAcceleration
#    PRIVATE_LINKS
#    ${Tudat_PROPAGATION_LIBRARIES})
#
##TUDAT_ADD_TEST_CASE(CentralAndZonalGravityModel
##        PRIVATE_LINKS
##        tudat_gravitation
##        tudat_basic_astrodynamics
##        tudat_basic_astrodynamics
##        tudat_numerical_integrators
##        SOURCES
##        "planetTestData.cpp"
##        )
#
#TUDAT_ADD_TEST_CASE(GravityFieldVariations
#        PRIVATE_LINKS
#        ${Tudat_PROPAGATION_LIBRARIES}
#        )
#
#TUDAT_ADD_TEST_CASE(MutualSphericalHarmonicAcceleration
#        PRIVATE_LINKS
#        ${Tudat_PROPAGATION_LIBRARIES}
#        )
#
#TUDAT_ADD_TEST_CASE(GravitationalTorques
#        PRIVATE_LINKS
#        ${Tudat_PROPAGATION_LIBRARIES}
#        )

TUDAT_ADD_TEST_CASE(PolyhedronGravityField
        PRIVATE_LINKS
        tudat_gravitation
        tudat_basic_astrodynamics
        tudat_basic_mathematics)

TUDAT_ADD_TEST_CASE(PolyhedronGravityModel
        PRIVATE_LINKS
        tudat_gravitation
        tudat_basic_astrodynamics
        tudat_basic_mathematics)