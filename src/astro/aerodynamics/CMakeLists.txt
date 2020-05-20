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
set(aerodynamics_SOURCES
        "aerodynamicAcceleration.cpp"
        "aerodynamicCoefficientGenerator.cpp"
        "aerodynamicTorque.cpp"
        "aerodynamicForce.cpp"
        "aerodynamics.cpp"
        "customConstantTemperatureAtmosphere.cpp"
        "exponentialAtmosphere.cpp"
        "hypersonicLocalInclinationAnalysis.cpp"
        "tabulatedAtmosphere.cpp"
        "flightConditions.cpp"
        "trimOrientation.cpp"
        "equilibriumWallTemperature.cpp"
        )

# Set the header files.
set(aerodynamics_HEADERS
        "aerodynamicAcceleration.h"
        "aerodynamicTorque.h"
        "aerodynamicCoefficientGenerator.h"
        "aerodynamicCoefficientInterface.h"
        "aerodynamicAcceleration.h"
        "aerodynamicForce.h"
        "aerodynamics.h"
        "atmosphereModel.h"
        "customConstantTemperatureAtmosphere.h"
        "exponentialAtmosphere.h"
        "hypersonicLocalInclinationAnalysis.h"
        "tabulatedAtmosphere.h"
        "standardAtmosphere.h"
        "customAerodynamicCoefficientInterface.h"
        "controlSurfaceAerodynamicCoefficientInterface.h"
        "flightConditions.h"
        "tests/testApolloCapsuleCoefficients.h"
        "flightConditions.cpp"
        "aerodynamicGuidance.h"
        "equilibriumWallTemperature.h"
        "windModel.h"
        )

if (TUDAT_WITH_NRLMSISE00)
    set(aerodynamics_SOURCES
            "${aerodynamics_SOURCES}"
            "nrlmsise00Atmosphere.cpp"
            "nrlmsise00InputFunctions.cpp"
            )
    set(aerodynamics_HEADERS
            "${aerodynamics_HEADERS}"
            "nrlmsise00Atmosphere.h"
            "nrlmsise00InputFunctions.h"
            )
endif ()

# Add library.
TUDAT_ADD_LIBRARY("aerodynamics"
        "${aerodynamics_SOURCES}"
        "${aerodynamics_HEADERS}"
        # For SHARED
        tudat_input_output
        tudat_geometric_shapes
        tudat_reference_frames
        tudat_interpolators
        )

# Add test data.
#TUDAT_ADD_TEST_DATA(tests/aerodynamics_test_data)

# Add tests.
TUDAT_ADD_TESTCASE(AerodynamicMomentAndAerodynamicForce
        ${Tudat_PROPAGATION_LIBRARIES}
        tudat_aerodynamics
        )

TUDAT_ADD_TESTCASE(AerodynamicsNamespace
        tudat_aerodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(CoefficientGenerator
        tudat_aerodynamics
        tudat_geometric_shapes
        tudat_interpolators
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(ExponentialAtmosphere
        tudat_aerodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(CustomConstantTemperatureAtmosphere
        tudat_aerodynamics
        tudat_basic_mathematics
        )

TUDAT_ADD_TESTCASE(TabulatedAtmosphere
        tudat_aerodynamics
        tudat_input_output
        tudat_interpolators
        tudat_basic_mathematics
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TESTCASE(HeatTransfer
        tudat_aerodynamics
        tudat_root_finders
        tudat_basic_mathematics
        tudat_basic_astrodynamics
        )

TUDAT_ADD_TESTCASE(ControlSurfaceIncrements
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_TESTCASE(AerodynamicCoefficientsFromFile
        ${Tudat_PROPAGATION_LIBRARIES}
        )

TUDAT_ADD_TESTCASE(WindModel
        ${Tudat_PROPAGATION_LIBRARIES}
        )

if (TUDAT_BUILD_WITH_NRLMSISE00)
    TUDAT_ADD_TESTCASE(NRLMSISE00Atmosphere
            tudat_aerodynamics
            tudat_interpolators
            tudat_basic_mathematics
            nrlmsise00
            tudat_input_output
            tudat_basic_astrodynamics)
endif ()
