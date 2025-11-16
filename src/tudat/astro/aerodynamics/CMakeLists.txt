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
)

# Conditionally add MCD atmosphere model (follows NRLMSISE pattern)
if(TUDAT_BUILD_WITH_MCD)
    list(APPEND aerodynamics_SOURCES "mcdAtmosphereModel.cpp")
endif()

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
)
# Conditionally add MCD header
if(TUDAT_BUILD_WITH_MCD)
    list(APPEND aerodynamics_HEADERS "mcdAtmosphereModel.h")
endif()

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
    TUDAT_ADD_LIBRARY("aerodynamics"
            "${aerodynamics_SOURCES}"
            "${aerodynamics_HEADERS}"
            PRIVATE_LINKS
            "${NRLMSISE00_LIBRARIES}"
            #        tudat_input_output
            #        tudat_geometric_shapes
            #        tudat_reference_frames
            #        tudat_interpolators
            #        PRIVATE_LINKS "${Boost_LIBRARIES}"
            )

            # Add MCD-specific linking if enabled

if(TUDAT_BUILD_WITH_MCD)
    # Link NetCDF libraries
    target_link_libraries(tudat_aerodynamics PUBLIC 
        NetCDF::NetCDF_Fortran
        NetCDF::NetCDF_C
    )
    
    # Include NetCDF headers
    target_include_directories(tudat_aerodynamics PUBLIC 
        ${NetCDF_INCLUDE_DIRS}
    )
    
    # Link MCD library if it exists as separate target
    if(TARGET mcd_library)
        target_link_libraries(tudat_aerodynamics PUBLIC mcd_library)
    endif()
    
    message(STATUS "MCD support enabled for tudat_aerodynamics")
else()
    # Add compile definition to disable MCD
    target_compile_definitions(tudat_aerodynamics PUBLIC TUDAT_BUILD_WITH_MCD=0)
endif()
