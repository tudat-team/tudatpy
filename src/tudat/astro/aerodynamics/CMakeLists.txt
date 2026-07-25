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
        "aerodynamicCoefficientInterface.cpp"
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
        "gasSurfaceInteractionModel.cpp"
        "marsDtmAtmosphereModel.cpp"
        "comaModel.cpp"
        "comaWindModel.cpp"
        "nrlmsise00Atmosphere.cpp"
        "nrlmsise00InputFunctions.cpp"
)

# Set the header files.
set(aerodynamics_HEADERS
        "aerodynamicAcceleration.h"
        "aerodynamicTorque.h"
        "aerodynamicCoefficientGenerator.h"
        "aerodynamicCoefficientInterface.h"
        "panelledAerodynamicCoefficientInterface.h"
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
        "aerodynamicGuidance.h"
        "equilibriumWallTemperature.h"
        "windModel.h"
        "comaWindModel.h"
        "gasSurfaceInteractionModel.h"
        "marsDtmAtmosphereModel.h"
        "comaModel.h"
        "nrlmsise00Atmosphere.h"
        "nrlmsise00InputFunctions.h"
)

# Conditionally add MCD sources and headers
if(TUDAT_BUILD_WITH_MCD_INTERFACE)
    list(APPEND aerodynamics_SOURCES "mcdAtmosphereModel.cpp")
    list(APPEND aerodynamics_HEADERS "mcdAtmosphereModel.h")
endif()

set(PUBLIC_LINKS_LIST)
if(TUDAT_BUILD_WITH_MCD_INTERFACE)
    list(APPEND PUBLIC_LINKS_LIST tudat_mcd_interface)
endif()

TUDAT_ADD_LIBRARY("aerodynamics"
            "${aerodynamics_SOURCES}"
            "${aerodynamics_HEADERS}"
            PRIVATE_LINKS
            "${NRLMSISE00_LIBRARIES}"
            PUBLIC_LINKS
            ${PUBLIC_LINKS_LIST}
            #        tudat_input_output
            #        tudat_geometric_shapes
            #        tudat_reference_frames
            #        tudat_interpolators
            #        PRIVATE_LINKS "${Boost_LIBRARIES}"
            )
