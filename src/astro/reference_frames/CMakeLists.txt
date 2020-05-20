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
set(reference_frames_SOURCES
        "referenceFrameTransformations.cpp"
        "apparentAccelerationModel.cpp"
        "aerodynamicAngleCalculator.cpp"
        )

# Add header files.
set(reference_frames_HEADERS
        "referenceFrameTransformations.h"
        "apparentAccelerationModel.h"
        "aerodynamicAngleCalculator.h"
        "dependentOrientationCalculator.h"
        )

TUDAT_ADD_LIBRARY("reference_frames"
        "${reference_frames_SOURCES}"
        "${reference_frames_HEADERS}"
        # For SHARED
        tudat_basic_astrodynamics
        )

# Add static libraries.
#add_library(tudat_reference_frames STATIC ${REFERENCEFRAMES_SOURCES} ${REFERENCEFRAMES_HEADERS})
#setup_tudat_library_target(tudat_reference_frames "${SRCROOT}${REFERENCEFRAMESDIR}")

## Add unit tests.
#add_executable(test_ReferenceFrameTransformations "${SRCROOT}${REFERENCEFRAMESDIR}/tests/unitTestReferenceFrameTransformations.cpp")
#setup_custom_test_program(test_ReferenceFrameTransformations "${SRCROOT}${REFERENCEFRAMESDIR}")
#target_link_libraries(test_ReferenceFrameTransformations tudat_reference_frames ${Boost_LIBRARIES})
#
#add_executable(test_ApparentAccelerationModel "${SRCROOT}${REFERENCEFRAMESDIR}/tests/unitTestApparentAccelerationModel.cpp")
#setup_custom_test_program(test_ApparentAccelerationModel "${SRCROOT}${REFERENCEFRAMESDIR}")
#target_link_libraries(test_ApparentAccelerationModel tudat_reference_frames tudat_basic_astrodynamics ${Boost_LIBRARIES})
#
#add_executable(test_AerodynamicAngleCalculator "${SRCROOT}${REFERENCEFRAMESDIR}/tests/unitTestAerodynamicAngleCalculator.cpp")
#setup_custom_test_program(test_AerodynamicAngleCalculator "${SRCROOT}${REFERENCEFRAMESDIR}")
#target_link_libraries(test_AerodynamicAngleCalculator tudat_environment_setup tudat_propagators tudat_aerodynamics tudat_geometric_shapes tudat_gravitation tudat_electromagnetism tudat_propulsion tudat_ephemerides tudat_numerical_integrators tudat_reference_frames tudat_basic_astrodynamics tudat_input_output tudat_basic_mathematics tudat_spice_interface ${CSpice_LIBRARIES} ${Boost_LIBRARIES})
#
