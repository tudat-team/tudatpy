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
  "oneWayRangePartial.cpp"
  "twoWayDopplerPartial.cpp"
  "nWayRangePartial.cpp"
  "differencedOneWayRangeRatePartial.cpp"
  "oneWayDopplerPartial.cpp"
  "observationPartial.cpp"
  "lightTimeCorrectionPartial.cpp"
  "firstOrderRelativisticPartial.cpp"
)

# Set the header files.
set(observation_partials_HEADERS
  "angularPositionPartial.h"
  "rotationMatrixPartial.h"
  "positionPartials.h"
  "oneWayRangePartial.h"
  "nWayRangePartial.h"
  "differencedOneWayRangeRatePartial.h"
  "oneWayDopplerPartial.h"
  "twoWayDopplerPartial.h"
  "observationPartial.h"
  "eulerAngleObservablePartials.h"
  "lightTimeCorrectionPartial.h"
  "firstOrderRelativisticPartial.h"
)

## Set the header files.
#set(observation_partials_TEST_HEADERS
#  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/tests/observationPartialTestFunctions.h"
#  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/tests/numericalObservationPartial.h"
#)

#set(observation_partials_TEST_SOURCES
#  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/tests/observationPartialTestFunctions.cpp"
#  "${SRCROOT}${OBSERVATIONPARTIALSDIR}/tests/numericalObservationPartial.cpp"
#)

TUDAT_ADD_LIBRARY("observation_partials"
        "${observation_partials_SOURCES}"
        "${observation_partials_HEADERS}"
#        PUBLIC_LINKS
#        tudat_input_output
#        tudat_geometric_shapes
#        tudat_reference_frames
#        tudat_interpolators
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )


#add_library(tudat_observation_partials_test_functions STATIC ${observation_partials_TEST_SOURCES} ${observation_partials_TEST_HEADERS})
#setup_tudat_library_target(tudat_observation_partials "${SRCROOT}{OBSERVATIONPARTIALSDIR}/UnitTests")
