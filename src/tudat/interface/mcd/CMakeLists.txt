#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

#Set the source files.
set(mcd_interface_SOURCES
        "marsClimateDatabaseClimateModel.cpp"
        )

# Set the header files.
set(mcd_interface_HEADERS
        "marsClimateDatabaseClimateModel.h"
        )

# Add static libraries.
TUDAT_ADD_LIBRARY("mcd_interface"
        "${mcd_interface_SOURCES}"
        "${mcd_interface_HEADERS}"
        PUBLIC_LINKS
        tudat_basic_astrodynamics
        PRIVATE_INCLUDES
        ${CMAKE_BINARY_DIR}/mcd_modules
        ${CMAKE_SOURCE_DIR}/third_parties/mcd
        )

# Expose MCD headers to dependent targets (e.g., tudat_propagation_setup)
# Use generator expressions to handle build vs install scenarios
target_include_directories(tudat_mcd_interface 
        SYSTEM PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/mcd_modules>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/third_parties/mcd>
        )

target_link_libraries(tudat_mcd_interface PRIVATE mcd_fortran)

set(MCD_DATA_DIR "${CMAKE_SOURCE_DIR}/third_parties/mcd/data")
target_compile_definitions(tudat_mcd_interface PRIVATE
    MCD_DATA_PATH="${MCD_DATA_DIR}/"
)

message(STATUS "MCD fortran library linked to mcd_interface")
message(STATUS "MCD data path for mcd_interface: ${MCD_DATA_DIR}/")