#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# # Add subdirectories.
add_subdirectory(shape_based)

set(low_thrust_trajectories_SOURCES
    "lowThrustLeg.cpp"
    )

set(low_thrust_trajectories_HEADERS
    "lowThrustLeg.h"
    )


# Add library.
TUDAT_ADD_LIBRARY("low_thrust_trajectories"
        "${low_thrust_trajectories_SOURCES}"
        "${low_thrust_trajectories_HEADERS}"
        )
