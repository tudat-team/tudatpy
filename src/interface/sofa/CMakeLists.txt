#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# http://www.iausofa.org/current_C.html#Downloads

#Set the source files.
set(sofa_interface_SOURCES
        "sofaTimeConversions.cpp"
        "fundamentalArguments.cpp"
        "earthOrientation.cpp"
        )

# Set the header files.
set(sofa_interface_HEADERS
        "sofaTimeConversions.h"
        "fundamentalArguments.h"
        "earthOrientation.h"
        )

# Add static libraries.
TUDAT_ADD_LIBRARY("sofa_interface"
        "${sofa_interface_SOURCES}"
        "${sofa_interface_HEADERS}"
#        PUBLIC_LINKS
#        tudat_basic_astrodynamics
        PRIVATE_LINKS "${Sofa_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}" "${Sofa_INCLUDE_DIRS}"
        )
