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
        defaultGroundStationSettings.h
        createBodies.h
        createBodiesFactory.h
        createBodiesSettings.h
        createGroundStations.h
        createBodyDeformationModel.h
        baseStateInterface.h
        rigidBodyProperties.h
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
        createClimateModel.h
        createCameras.h
        createRelativisticTimeConverter.h
        createMetric.h
        )

# Add header files.
set(environment_SOURCES
        createRadiationPressureInterface.cpp
        createEphemeris.cpp
        createRotationModel.cpp
        createRotationModelEarthOrientation.cpp
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
        defaultGroundStationSettings.cpp
        baseStateInterface.cpp
        rigidBodyProperties.cpp
        body.cpp
        createAtmosphereModel.cpp
        createSystemModel.cpp
        createThrustModelGuidance.cpp
        thrustSettings.cpp
        createRadiationSourceModel.cpp
        createRadiationPressureTargetModel.cpp
        createOccultationModel.cpp
        createSurfacePropertyDistribution.cpp
        createClimateModel.cpp
        createCameras.cpp
        createRelativisticTimeConverter.cpp
        createMetric.cpp
        )

if (TUDAT_BUILD_WITH_SERIALIZATION)
    list(APPEND environment_SOURCES gravityFieldVariationSettingsSerialization.cpp)
endif ()

# Add library.
TUDAT_ADD_LIBRARY("environment_setup"
        "${environment_SOURCES}"
        "${environment_HEADERS}"
        PUBLIC_SERIALIZATION
#        PRIVATE_LINKS "Boost::boost"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}" ${CSpice_INCLUDE_DIRS} ${Sofa_INCLUDE_DIRS}
        )

if(TUDAT_BUILD_WITH_MCD_INTERFACE)
    target_link_libraries(tudat_environment_setup PRIVATE
        mcd_fortran
        tudat_mcd_interface
        ${NetCDF_Fortran_LIBRARIES}
        ${NetCDF_C_LIBRARIES}
    )
    target_include_directories(tudat_environment_setup SYSTEM PRIVATE 
        "${NetCDF_INCLUDE_DIRS}"
        "${CMAKE_SOURCE_DIR}/third_parties/mcd")
endif()
