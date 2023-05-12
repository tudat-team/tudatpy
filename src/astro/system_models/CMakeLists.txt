#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Set the header files.
set(system_models_HEADERS
        "engineModel.h"
        "vehicleSystems.h"
        "timingSystem.h"
        )

set(system_models_SOURCES
        "dummySystemModel.cpp"
        "timingSystem.cpp"
        )

TUDAT_ADD_LIBRARY("system_models"
        "${system_models_SOURCES}"
        "${system_models_HEADERS}"
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
