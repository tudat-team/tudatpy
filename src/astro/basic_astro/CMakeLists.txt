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
set(basic_astrodynamics_SOURCES
        "accelerationModelTypes.cpp"
        "accelerationModel.cpp"
        "attitudeElementConversions.cpp"
        "clohessyWiltshirePropagator.cpp"
        "geodeticCoordinateConversions.cpp"
        "missionGeometry.cpp"
        "modifiedEquinoctialElementConversions.cpp"
        "timeConversions.cpp"
        "astrodynamicsFunctions.cpp"
        "physicalConstants.cpp"
        "polyhedronFunctions.cpp"
        "bodyShapeModel.cpp"
        "polyhedronBodyShapeModel.cpp"
        "sphericalStateConversions.cpp"
        "unifiedStateModelQuaternionElementConversions.cpp"
        "unifiedStateModelModifiedRodriguesParameterElementConversions.cpp"
        "unifiedStateModelExponentialMapElementConversions.cpp"
        "stateRepresentationConversions.cpp"
        "torqueModelTypes.cpp"
        "torqueModel.cpp"
        )

# Set the header files.
set(basic_astrodynamics_HEADERS
        "accelerationModelTypes.h"
        "accelerationModel.h"
        "customAccelerationModel.h"
        "attitudeElementConversions.h"
        "celestialBodyConstants.h"
        "convertMeanToEccentricAnomalies.h"
        "clohessyWiltshirePropagator.h"
        "customTorque.h"
        "geodeticCoordinateConversions.h"
        "keplerPropagator.h"
        "missionGeometry.h"
        "modifiedEquinoctialElementConversions.h"
        "stateVectorIndices.h"
        "timeConversions.h"
        "tests/testAccelerationModels.h"
        "tests/testBody.h"
        "tests/keplerPropagatorTestData.h"
        "astrodynamicsFunctions.h"
        "orbitalElementConversions.h"
        "physicalConstants.h"
        "polyhedronFunctions.h"
        "unitConversions.h"
        "bodyShapeModel.h"
        "oblateSpheroidBodyShapeModel.h"
        "sphericalBodyShapeModel.h"
        "polyhedronBodyShapeModel.h"
        "hybridBodyShapeModel.h"
        "massRateModel.h"
        "sphericalStateConversions.h"
        "unifiedStateModelQuaternionElementConversions.h"
        "unifiedStateModelModifiedRodriguesParameterElementConversions.h"
        "unifiedStateModelExponentialMapElementConversions.h"
        "stateRepresentationConversions.h"
        "empiricalAcceleration.h"
        "torqueModel.h"
        "torqueModelTypes.h"
        "dissipativeTorqueModel.h"
        )

TUDAT_ADD_LIBRARY("basic_astrodynamics"
        "${basic_astrodynamics_SOURCES}"
        "${basic_astrodynamics_HEADERS}")
