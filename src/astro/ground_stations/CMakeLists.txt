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
set(ground_stations_SOURCES
        "groundStation.cpp"
        "groundStationState.cpp"
        "pointingAnglesCalculator.cpp"        
        "basicTidalBodyDeformation.cpp"
        "iers2010SolidTidalBodyDeformation.cpp"
        "transmittingFrequencies.cpp"
        "oceanTideEarthDeformation.cpp"
        "poleTideDeformation.cpp"
)

# Set the header files.
set(ground_stations_HEADERS
        "groundStation.h"
        "groundStationState.h"
        "pointingAnglesCalculator.h"        
        "basicTidalBodyDeformation.h"
        "bodyDeformationModel.h"
        "iers2010SolidTidalBodyDeformation.h"
        "transmittingFrequencies.h"
        "oceanTideEarthDeformation.h"
        "poleTideDeformation.h"
)

TUDAT_ADD_LIBRARY("ground_stations"
        "${ground_stations_SOURCES}"
        "${ground_stations_HEADERS}"
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
