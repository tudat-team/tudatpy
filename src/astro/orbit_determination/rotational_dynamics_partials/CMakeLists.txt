 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.

# Set the source files.
set(rotational_dynamics_partials_SOURCES
  "inertialTorquePartial.cpp"
  "constantTorquePartial.cpp"
  "secondDegreeGravitationalTorquePartial.cpp"
  "sphericalHarmonicGravitationalTorquePartial.cpp"
)

# Set the header files.
set(rotational_dynamics_partials_HEADERS
  "torquePartial.h"
  "inertiaTensorPartial.h"
  "inertialTorquePartial.h"
  "constantTorquePartial.h"
  "secondDegreeGravitationalTorquePartial.h"
  "sphericalHarmonicGravitationalTorquePartial.h"
)


TUDAT_ADD_LIBRARY("torque_partials"
        "${rotational_dynamics_partials_SOURCES}"
        "${rotational_dynamics_partials_HEADERS}"
        PUBLIC_LINKS
#        tudat_input_output
#        tudat_geometric_shapes
#        tudat_reference_frames
#        tudat_interpolators
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )

