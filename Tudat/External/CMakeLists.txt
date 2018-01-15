#    Copyright (c) 2010-2018, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Define the main sub-directories.

if(USE_SOFA)
    # Set SPICE interface directory.
    set(SOFAINTERFACEDIR "${EXTERNALDIR}/SofaInterface")

    # Add subdirectories.
    add_subdirectory("${SRCROOT}${SOFAINTERFACEDIR}")

    # Get target properties for static libraries.
    get_target_property(SOFAINTERFACESOURCES tudat_sofa_interface SOURCES)

endif( )


if(USE_CSPICE)
    # Set SPICE interface directory.
    set(SPICEINTERFACEDIR "${EXTERNALDIR}/SpiceInterface")

    # Add subdirectories.
    add_subdirectory("${SRCROOT}${SPICEINTERFACEDIR}")

    # Get target properties for static libraries.
    get_target_property(SPICEINTERFACESOURCES tudat_spice_interface SOURCES)

endif( )
