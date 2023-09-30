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
        "cannonBallRadiationPressureAcceleration.h"
        "cannonBallRadiationPressureForce.h"
        "lorentzStaticMagneticForce.h"
        "lorentzStaticMagneticAcceleration.h"
        "radiationPressureInterface.h"
        "basicElectroMagnetism.h"
        "panelledRadiationPressure.h"
        "solarSailAcceleration.h"
        "solarSailForce.h"
        "yarkovskyAcceleration.h"
        )

# Set the source files.
set(electromagnetism_SOURCES
        "cannonBallRadiationPressureAcceleration.cpp"
        "cannonBallRadiationPressureForce.cpp"
        "lorentzStaticMagneticForce.cpp"
        "lorentzStaticMagneticAcceleration.cpp"
        "radiationPressureInterface.cpp"
        "panelledRadiationPressure.cpp"
        "solarSailAcceleration.cpp"
        "solarSailForce.cpp"
        "yarkovskyAcceleration.cpp"
        )

# Add tests.
TUDAT_ADD_LIBRARY("electromagnetism"
        "${electromagnetism_SOURCES}"
        "${electromagnetism_HEADERS}"
#        PRIVATE_LINKS "${Boost_LIBRARIES}"
#        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )
