 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

 add_subdirectory("acceleration_partials")
 add_subdirectory("rotational_dynamics_partials")
 add_subdirectory("estimatable_parameters")
 add_subdirectory("observation_partials")

## Set the source files.
set(orbit_determination_SOURCES
  "massDerivativePartial.cpp"
  "stateDerivativePartial.cpp"
  "podInputOutputTypes.cpp"
)

# Set the header files.
set(orbit_determination_HEADERS
  "massDerivativePartial.h"
  "stateDerivativePartial.h"
  "podInputOutputTypes.h"
)
#
#
# Add library.
TUDAT_ADD_LIBRARY("orbit_determination"
        "${orbit_determination_SOURCES}"
        "${orbit_determination_HEADERS}"
        )
