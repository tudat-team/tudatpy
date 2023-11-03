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
set(environment_HEADERS
        createAerodynamicCoefficientInterface.h
        createBodyShapeModel.h
        createEphemeris.h
        createGravityField.h
        createAerodynamicControlSurfaces.h
        createRotationModel.h
        defaultBodies.h
        createBodies.h
        createGroundStations.h
        createBodyDeformationModel.h
        body.h
        createRadiationPressureInterface.h
        createGravityFieldVariations.h
        createAtmosphereModel.h
        createFlightConditions.h
        createRadiationSourceModel.h
        createRadiationPressureTargetModel.h
        createOccultationModel.h
        createSurfacePropertyDistribution.h
        createSystemModel.h
        createThrustModelGuidance.h
        thrustSettings.h
        )

# Add header files.
set(environment_SOURCES
        createRadiationPressureInterface.cpp
        createEphemeris.cpp
        createRotationModel.cpp
        createGravityFieldVariations.cpp
        createBodies.cpp
        createGravityField.cpp
        createBodyShapeModel.cpp
        createBodyDeformationModel.cpp
        createFlightConditions.cpp
        createGroundStations.cpp
        createAerodynamicCoefficientInterface.cpp
        createAerodynamicControlSurfaces.cpp
        defaultBodies.cpp
        body.cpp
        createAtmosphereModel.cpp
        createSystemModel.cpp
        createThrustModelGuidance.cpp
        thrustSettings.cpp
        createRadiationSourceModel.cpp
        createRadiationPressureTargetModel.cpp
        createOccultationModel.cpp
        createSurfacePropertyDistribution.cpp
        )

# Add library.
TUDAT_ADD_LIBRARY("environment_setup"
        "${environment_SOURCES}"
        "${environment_HEADERS}"
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}" ${CSpice_INCLUDE_DIRS} ${Sofa_INCLUDE_DIRS}
        )
