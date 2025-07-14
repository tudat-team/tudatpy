 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

# Set the source files.
set(observation_partials_SOURCES
  "angularPositionPartial.cpp"
  "rotationMatrixPartial.cpp"
  "positionPartials.cpp"
  "differencedObservationPartial.cpp"
  "oneWayRangePartial.cpp"
  "twoWayDopplerPartial.cpp"
  "nWayRangePartial.cpp"
  "oneWayDopplerPartial.cpp"
  "relativeAngularPositionPartial.cpp"
  "observationPartial.cpp"
  "lightTimeCorrectionPartial.cpp"
  "firstOrderRelativisticPartial.cpp"
  "clockParameterPartials.cpp"
)

# Set the header files.
set(observation_partials_HEADERS
  "angularPositionPartial.h"
  "rotationMatrixPartial.h"
  "positionPartials.h"
  "oneWayRangePartial.h"
  "nWayRangePartial.h"
  "differencedObservationPartial.h"
  "directObservationPartial.h"
  "oneWayDopplerPartial.h"
  "twoWayDopplerPartial.h"
  "observationPartial.h"
  "eulerAngleObservablePartials.h"
  "relativeAngularPositionPartial.h"
  "lightTimeCorrectionPartial.h"
  "firstOrderRelativisticPartial.h"
  "clockParameterPartials.h"
        )

TUDAT_ADD_LIBRARY("observation_partials"
        "${observation_partials_SOURCES}"
        "${observation_partials_HEADERS}"
        )

