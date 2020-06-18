#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Set the source files.
set(json_interface_SOURCES
        # support
        "support/keys.cpp"
        "support/deserialization.cpp"
        "support/path.cpp"
        "support/valueAccess.cpp"
        "support/options.cpp"
        # environment
        "environment/spice.cpp"
        "environment/body.cpp"
        "environment/atmosphere.cpp"
        "environment/ephemeris.cpp"
        "environment/aerodynamics.cpp"
        "environment/gravityField.cpp"
        "environment/rotationModel.cpp"
        "environment/shapeModel.cpp"
        "environment/radiationPressure.cpp"
        "environment/gravityFieldVariation.cpp"
        "environment/groundStations.cpp"
        # propagation
        "propagation/variable.cpp"
        "propagation/termination.cpp"
        "propagation/acceleration.cpp"
        "propagation/massRateModel.cpp"
        "propagation/thrust.cpp"
        "propagation/torque.cpp"
        "propagation/export.cpp"
        # math
        "math/interpolation.cpp"
        # Executable source file
        "jsonInterface.cpp"
        )

# Set the header files.
set(json_interface_HEADERS
        # support
        "support/errorHandling.h"
        "support/keys.h"
        "support/deserialization.h"
        "support/path.h"
        "support/utilities.h"
        "support/valueAccess.h"
        "support/valueConversions.h"
        "support/options.h"
        # environment
        "environment/spice.h"
        "environment/body.h"
        "environment/atmosphere.h"
        "environment/ephemeris.h"
        "environment/aerodynamics.h"
        "environment/gravityField.h"
        "environment/rotationModel.h"
        "environment/shapeModel.h"
        "environment/radiationPressure.h"
        "environment/gravityFieldVariation.h"
        "environment/groundStations.h"
        # propagation
        "propagation/variable.h"
        "propagation/state.h"
        "propagation/propagator.h"
        "propagation/termination.h"
        "propagation/referenceFrames.h"
        "propagation/acceleration.h"
        "propagation/thrust.h"
        "propagation/massRateModel.h"
        "propagation/torque.h"
        "propagation/export.h"
        # math
        "math/integrator.h"
        "math/interpolation.h"
        # Unittests
        "tests/unitTestSupport.h"
        # Executable header file
        "jsonInterface.h"
        )

if (TUDAT_BUILD_WITH_ESTIMATION_TOOLS)
    set(json_interface_SOURCES
            ${json_interface_SOURCES}
            "estimation/parameter.cpp"
            "estimation/observation.cpp"
            #"estimation/orbitDetermination.cpp"
            "jsonInterfaceVariational.cpp"
            #"jsonEstimationInterface.cpp"
            )

    # Set the header files.
    set(json_interface_HEADERS
            ${json_interface_HEADERS}
            "estimation/parameter.h"
            "estimation/observation.h"
            #"estimation/orbitDetermination.h"
            "jsonInterfaceVariational.h"
            #"jsonEstimationInterface.h"
            )
endif ()


if (TUDAT_BUILD_WITH_ESTIMATION_TOOLS)
    list(APPEND json_propagation_LIBRARIES
            json_interface_library
            ${Tudat_ESTIMATION_LIBRARIES}
            ${Boost_LIBRARIES})
else ()
    list(APPEND json_propagation_LIBRARIES
            json_interface_library
            ${Tudat_PROPAGATION_LIBRARIES}
            ${Boost_LIBRARIES})
endif ()

# Add library.
TUDAT_ADD_LIBRARY("json_interface"
        "${json_interface_SOURCES}"
        "${json_interface_HEADERS}"
        PUBLIC_INTERFACE
        json_propagation_LIBRARIES
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )

# Add JSON-based application.
TUDAT_ADD_EXECUTABLE(json_interface
        jsonInterfaceExecute.cpp
        tudat_json_interface
        )
