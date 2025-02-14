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
set(relativity_SOURCES
        "relativisticLightTimeCorrection.cpp"
        "relativisticAccelerationCorrection.cpp"
        "relativisticTimeConversion.cpp"
        "metric.cpp"
        "einsteinInfeldHoffmannEquations.cpp"
        )

# Set the header files.
set(relativity_HEADERS
        "relativisticLightTimeCorrection.h"
        "relativisticAccelerationCorrection.h"
        "relativisticTimeConversion.h"
        "metric.h"
        "einsteinInfeldHoffmannEquations.h"
        "einsteinInfeldHoffmannAcceleration.h"
        )

TUDAT_ADD_LIBRARY("relativity"
        "${relativity_SOURCES}"
        "${relativity_HEADERS}"
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
