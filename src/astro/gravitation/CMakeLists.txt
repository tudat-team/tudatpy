#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#
#    Notes
#
#

# Set the source files.
set(gravitation_SOURCES
        "basicSolidBodyTideGravityFieldVariations.cpp"
        "gravityFieldVariations.cpp"
        "centralGravityModel.cpp"
#        "centralJ2GravityModel.cpp"
#        "centralJ2J3GravityModel.cpp"
#        "centralJ2J3J4GravityModel.cpp"
        "gravityFieldModel.cpp"
        "jacobiEnergy.cpp"
        "librationPoint.cpp"
        "sphericalHarmonicsGravityModel.cpp"
        "sphericalHarmonicsGravityField.cpp"
        "thirdBodyPerturbation.cpp"
        "timeDependentSphericalHarmonicsGravityField.cpp"
        "unitConversionsCircularRestrictedThreeBodyProblem.cpp"
        "unitConversionsCircularRestrictedThreeBodyProblem.cpp"
        "triAxialEllipsoidGravity.cpp"
        "tabulatedGravityFieldVariations.cpp"
        "mutualSphericalHarmonicGravityModel.cpp"
        "secondDegreeGravitationalTorque.cpp"
        "directTidalDissipationAcceleration.cpp"
        "periodicGravityFieldVariations.cpp"
        "polynomialGravityFieldVariations.cpp"
        "polyhedronGravityField.cpp"
        "polyhedronGravityModel.cpp"
        "ringGravityField.cpp"
        "ringGravityModel.cpp"
        )

# Set the header files.
set(gravitation_HEADERS
        "basicSolidBodyTideGravityFieldVariations.h"
        "gravityFieldVariations.h"
        "centralGravityModel.h"
#        "centralJ2GravityModel.h"
#        "centralJ2J3GravityModel.h"
#        "centralJ2J3J4GravityModel.h"
        "gravityFieldModel.h"
        "jacobiEnergy.h"
        "librationPoint.h"
        "sphericalHarmonicsGravityModel.h"
        "sphericalHarmonicsGravityModelBase.h"
        "sphericalHarmonicsGravityField.h"
        "thirdBodyPerturbation.h"
        "timeDependentSphericalHarmonicsGravityField.h"
        "unitConversionsCircularRestrictedThreeBodyProblem.h"
        "triAxialEllipsoidGravity.h"
        "tabulatedGravityFieldVariations.h"
        "mutualSphericalHarmonicGravityModel.h"
        "secondDegreeGravitationalTorque.h"
        "directTidalDissipationAcceleration.h"
        "sphericalHarmonicGravitationalTorque.h"
        "periodicGravityFieldVariations.h"
        "polynomialGravityFieldVariations.h"
        "polyhedronGravityField.h"
        "polyhedronGravityModel.h"
        "ringGravityField.h"
        "ringGravityModel.h"
        )

TUDAT_ADD_LIBRARY("gravitation"
        "${gravitation_SOURCES}"
        "${gravitation_HEADERS}")
