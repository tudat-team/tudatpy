#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

TUDAT_ADD_TEST_CASE(AerodynamicMomentAndAerodynamicForce
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        Tudat::tudat_aerodynamics
        )

TUDAT_ADD_TEST_CASE(AerodynamicsNamespace
        PRIVATE_LINKS
        Tudat::tudat_aerodynamics
        Tudat::tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(CoefficientGenerator
        PRIVATE_LINKS
        Tudat::tudat_aerodynamics
        Tudat::tudat_geometric_shapes
        Tudat::tudat_interpolators
        Tudat::tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(ExponentialAtmosphere
        PRIVATE_LINKS
        Tudat::tudat_aerodynamics
        Tudat::tudat_basic_mathematics
        )

TUDAT_ADD_TEST_CASE(CustomConstantTemperatureAtmosphere
        PRIVATE_LINKS
        Tudat::tudat_aerodynamics
        Tudat::tudat_basic_mathematics
        )

# TUDAT_ADD_TEST_CASE(TabulatedAtmosphere
#         PRIVATE_LINKS
#         Tudat::tudat_aerodynamics
#         Tudat::tudat_input_output
#         Tudat::tudat_interpolators
#         Tudat::tudat_basic_mathematics
#         Tudat::tudat_basic_astrodynamics
#         )

TUDAT_ADD_TEST_CASE(HeatTransfer
        PRIVATE_LINKS
        Tudat::tudat_aerodynamics
        Tudat::tudat_root_finders
        Tudat::tudat_basic_mathematics
        Tudat::tudat_basic_astrodynamics
        )

TUDAT_ADD_TEST_CASE(ControlSurfaceIncrements
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(AerodynamicCoefficientsFromFile
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_TEST_CASE(WindModel
        PRIVATE_LINKS
        ${Tudat_PROPAGATION_LIBRARIES}
        )

if (TUDAT_BUILD_WITH_NRLMSISE00)
    TUDAT_ADD_TEST_CASE(NRLMSISE00Atmosphere
            PRIVATE_LINKS
            ${Tudat_PROPAGATION_LIBRARIES}
            )
endif ()
