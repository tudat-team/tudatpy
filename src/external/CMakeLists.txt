#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

## Define the main sub-directories.
#if (TUDAT_USE_SOFA)
#    # Set SPICE interface directory.
#    set(SOFAINTERFACEDIR "${EXTERNALDIR}/SofaInterface")
#
#    # Add subdirectories.
#    add_subdirectory("${SRCROOT}${SOFAINTERFACEDIR}")
#
#    # Get target properties for static libraries.
#    get_target_property(SOFAINTERFACESOURCES tudat_sofa_interface SOURCES)
#
#endif ()
#
#
#if (TUDAT_USE_SPICE)
#    # Set SPICE interface directory.
#    set(SPICEINTERFACEDIR "${EXTERNALDIR}/SpiceInterface")
#
#    # Add subdirectories.
#    add_subdirectory("${SRCROOT}${SPICEINTERFACEDIR}")
#
#    # Get target properties for static libraries.
#    get_target_property(SPICEINTERFACESOURCES tudat_spice_interface SOURCES)
#endif ()

TUDAT_ADD_DATA_DIR("atmosphere_tables")
TUDAT_ADD_DATA_DIR("ephemeris_data")
TUDAT_ADD_DATA_DIR("gravity_models")
TUDAT_ADD_DATA_DIR("space_weather_data")