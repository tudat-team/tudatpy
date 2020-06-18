#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

if (TUDAT_BUILD_WITH_JSON_INTERFACE)
    add_subdirectory(json)
endif ()

if (TUDAT_BUILD_WITH_SPICE_INTERFACE)
    add_subdirectory(spice)
endif ()

if (TUDAT_BUILD_WITH_SOFA_INTERFACE)
    add_subdirectory(sofa)
endif ()