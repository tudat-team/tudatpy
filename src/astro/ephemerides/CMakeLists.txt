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
set(ephemerides_SOURCES
        "approximatePlanetPositionsBase.cpp"
        "approximatePlanetPositions.cpp"
        "approximatePlanetPositionsCircularCoplanar.cpp"
        "ephemeris.cpp"
        "rotationalEphemeris.cpp"
        "cartesianStateExtractor.cpp"
        "keplerStateExtractor.cpp"
        "keplerEphemeris.cpp"
        "rotationalEphemeris.cpp"
        "simpleRotationalEphemeris.cpp"
        "tabulatedEphemeris.cpp"
        "frameManager.cpp"
        "compositeEphemeris.cpp"
        "tabulatedRotationalEphemeris.cpp"
        "synchronousRotationalEphemeris.cpp"
        "fullPlanetaryRotationModel.cpp"
        "tleEphemeris.cpp"
        "aeordynamicAngleRotationalEphemeris.cpp"
        "directionBasedRotationalEphemeris.cpp"
        )

# Set the header files.
set(ephemerides_HEADERS
        "approximatePlanetPositionsBase.h"
        "approximatePlanetPositions.h"
        "approximatePlanetPositionsCircularCoplanar.h"
        "approximatePlanetPositionsDataContainer.h"
        "ephemeris.h"
        "constantEphemeris.h"
        "cartesianStateExtractor.h"
        "keplerStateExtractor.h"
        "keplerEphemeris.h"
        "rotationalEphemeris.h"
        "constantRotationalEphemeris.h"
        "simpleRotationalEphemeris.h"
        "tabulatedEphemeris.h"
        "frameManager.h"
        "itrsToGcrsRotationModel.h"
        "compositeEphemeris.h"
        "constantEphemeris.h"
        "constantRotationalEphemeris.h"
        "multiArcEphemeris.h"
        "tabulatedRotationalEphemeris.h"
        "fullPlanetaryRotationModel.h"
        "synchronousRotationalEphemeris.h"
        "tleEphemeris.h"
        "aeordynamicAngleRotationalEphemeris.h"
        "directionBasedRotationalEphemeris.h"
        )

TUDAT_ADD_LIBRARY("ephemerides"
        "${ephemerides_SOURCES}"
        "${ephemerides_HEADERS}"
#        PUBLIC_LINKS
#        tudat_spice_interface
#        tudat_reference_frames
#        tudat_interpolators
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}" "${CSpice_INCLUDE_DIRS}"
        )
