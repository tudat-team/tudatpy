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
        "gasSurfaceInteractionModel.h"
        "marsDtmAtmosphereModel.h"
        "nrlmsise00Atmosphere.h"
        "nrlmsise00InputFunctions.h"
)

# Conditionally add MCD sources and headers
if(TUDAT_BUILD_WITH_MCD)
    list(APPEND aerodynamics_SOURCES "mcdAtmosphereModel.cpp")
    list(APPEND aerodynamics_HEADERS "mcdAtmosphereModel.h")
endif()

# Create the library
TUDAT_ADD_LIBRARY("aerodynamics"
    "${aerodynamics_SOURCES}"
    "${aerodynamics_HEADERS}"
    PRIVATE_LINKS
    "${NRLMSISE00_LIBRARIES}"
)

# Link MCD if enabled - this should be the ONLY place we link MCD to aerodynamics
if(TUDAT_BUILD_WITH_MCD)
    # Link the MCD Fortran library
    target_link_libraries(tudat_aerodynamics PRIVATE mcd_fortran)
    
    # Include MCD directories for building tudat_aerodynamics only
    # These are not propagated to consumers of tudat_aerodynamics
    target_include_directories(tudat_aerodynamics PRIVATE
        ${CMAKE_BINARY_DIR}/mcd_modules
        ${CMAKE_SOURCE_DIR}/third_parties/mcd
    )
    
        # Get the MCD_DATA_PATH from mcd_fortran target
    get_target_property(MCD_DATA_PATH_DEF mcd_fortran INTERFACE_COMPILE_DEFINITIONS)
    
    # Propagate MCD_DATA_PATH to tudat_aerodynamics
    # The path is already set by mcd_fortran's INTERFACE_COMPILE_DEFINITIONS
    # but we explicitly add it here for clarity
    set(MCD_DATA_DIR "${CMAKE_SOURCE_DIR}/third_parties/mcd/data")
    target_compile_definitions(tudat_aerodynamics PRIVATE
        MCD_DATA_PATH="${MCD_DATA_DIR}/"
    )
    
    message(STATUS "MCD library linked to tudat_aerodynamics")
    message(STATUS "MCD data path for tudat_aerodynamics: ${MCD_DATA_DIR}/")
endif()
