 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

# Set the source files.
set(acceleration_partials_SOURCES
  "aerodynamicAccelerationPartial.cpp"
  "centralGravityAccelerationPartial.cpp"
  "numericalAccelerationPartial.cpp"
  "relativisticAccelerationPartial.cpp"
  "radiationPressureAccelerationPartial.cpp"
  "sphericalHarmonicPartialFunctions.cpp"
  "sphericalHarmonicAccelerationPartial.cpp"
  "mutualSphericalHarmonicGravityPartial.cpp"
  "empiricalAccelerationPartial.cpp"
  "tidalLoveNumberPartialInterface.cpp"
  "directTidalDissipationAccelerationPartial.cpp"
  "panelledRadiationPressureAccelerationPartial.cpp"
  "thrustAccelerationPartial.cpp"
  "polyhedronAccelerationPartial.cpp"
  "ringAccelerationPartial.cpp"
  "einsteinInfeldHoffmannPartials.cpp"
  "yarkovskyAccelerationPartial.cpp"
  "fullRadiationPressureAccelerationPartial.cpp"
)

# Set the header files.
set(acceleration_partials_HEADERS
  "accelerationPartial.h"
  "aerodynamicAccelerationPartial.h"
  "thirdBodyGravityPartial.h"
  "centralGravityAccelerationPartial.h"
  "numericalAccelerationPartial.h"
  "relativisticAccelerationPartial.h"
  "radiationPressureAccelerationPartial.h"
  "sphericalHarmonicPartialFunctions.h"
  "sphericalHarmonicAccelerationPartial.h"
  "mutualSphericalHarmonicGravityPartial.h"
  "empiricalAccelerationPartial.h"
  "tidalLoveNumberPartialInterface.h"
  "directTidalDissipationAccelerationPartial.h"
  "panelledRadiationPressureAccelerationPartial.h"
  "thrustAccelerationPartial.h"
  "polyhedronAccelerationPartial.h"
  "ringAccelerationPartial.h"
  "einsteinInfeldHoffmannPartials.h"
  "yarkovskyAccelerationPartial.h"
  "fullRadiationPressureAccelerationPartial.h"
)

TUDAT_ADD_LIBRARY("acceleration_partials"
        "${acceleration_partials_SOURCES}"
        "${acceleration_partials_HEADERS}"
        )

