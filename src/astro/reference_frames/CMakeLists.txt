#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Add source files.
set(reference_frames_SOURCES
        "referenceFrameTransformations.cpp"
        "apparentAccelerationModel.cpp"
        "aerodynamicAngleCalculator.cpp"
        )

# Add header files.
set(reference_frames_HEADERS
        "referenceFrameTransformations.h"
        "apparentAccelerationModel.h"
        "aerodynamicAngleCalculator.h"
        "dependentOrientationCalculator.h"
        )

TUDAT_ADD_LIBRARY("reference_frames"
        "${reference_frames_SOURCES}"
        "${reference_frames_HEADERS}"
        PUBLIC_LINKS tudat_basic_astrodynamics
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
