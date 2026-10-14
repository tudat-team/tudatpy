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
        "bodyDeformationModel.cpp"
        "stationLoadingDeformation.cpp"
        "earthDeformationTime.cpp"
        "groundStation.cpp"
        "groundStationState.cpp"
        "stationHistoryTrajectory.cpp"
        "stationHistoryTime.cpp"
        "transformedStationHistoryTrajectory.cpp"
        "stationReferencePoint.cpp"
        "pointingAnglesCalculator.cpp"        
        "basicTidalBodyDeformation.cpp"
        "iers2010SolidTidalBodyDeformation.cpp"
        "transmittingFrequencies.cpp"
        "oceanTideEarthDeformation.cpp"
        "poleTideDeformation.cpp"
        "oceanPoleTideDeformation.cpp"
        "atmosphericTideLoading.cpp"
        "atmosphericTideDeformation.cpp"
        "meteorologicalConditions.cpp"
)

# Set the header files.
set(ground_stations_HEADERS
        "groundStation.h"
        "groundStationState.h"
        "stationHistoryTrajectory.h"
        "stationHistoryMotionModel.h"
        "transformedStationHistoryTrajectory.h"
        "stationReferencePoint.h"
        "pointingAnglesCalculator.h"        
        "basicTidalBodyDeformation.h"
        "bodyDeformationModel.h"
        "stationLoadingDeformation.h"
        "earthDeformationTime.h"
        "iers2010SolidTidalBodyDeformation.h"
        "transmittingFrequencies.h"
        "oceanTideEarthDeformation.h"
        "poleTideDeformation.h"
        "oceanPoleTideDeformation.h"
        "atmosphericTideLoading.h"
        "atmosphericTideDeformation.h"
        "meteorologicalConditions.h"
)

TUDAT_ADD_LIBRARY("ground_stations"
        "${ground_stations_SOURCES}"
        "${ground_stations_HEADERS}"
        PUBLIC_LINKS tudat_environment_data tudat_reference_frames
#        PRIVATE_LINKS "Boost::boost"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
