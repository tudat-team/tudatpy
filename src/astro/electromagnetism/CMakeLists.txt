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
set(electromagnetism_HEADERS
        "lorentzStaticMagneticForce.h"
        "lorentzStaticMagneticAcceleration.h"
        "radiationPressureInterface.h"
        "basicElectroMagnetism.h"
        "radiationSourceModel.h"
        "luminosityModel.h"
        "sourcePanelRadiosityModel.h"
        "radiationPressureTargetModel.h"
        "radiationPressureAcceleration.h"
        "reflectionLaw.h"
        "occultationModel.h"
        "surfacePropertyDistribution.h"
        "yarkovskyAcceleration.h"
        )

# Set the source files.
set(electromagnetism_SOURCES
        "lorentzStaticMagneticForce.cpp"
        "lorentzStaticMagneticAcceleration.cpp"
        "radiationSourceModel.cpp"
        "luminosityModel.cpp"
        "sourcePanelRadiosityModel.cpp"
        "radiationPressureTargetModel.cpp"
        "radiationPressureAcceleration.cpp"
        "reflectionLaw.cpp"
        "occultationModel.cpp"
        "surfacePropertyDistribution.cpp"
        "yarkovskyAcceleration.cpp"
        )

# Add tests.
TUDAT_ADD_LIBRARY("electromagnetism"
        "${electromagnetism_SOURCES}"
        "${electromagnetism_HEADERS}"
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
