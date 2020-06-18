#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# TODO: Add relativity library tests

# Add static libraries.
#add_library(tudat_relativity STATIC ${RELATIVITY_SOURCES} ${RELATIVITY_HEADERS})
#setup_tudat_library_target(tudat_relativity "${SRCROOT}${RELATIVITYDIR}")

#if( TUDAT_BUILD_WITH_ESTIMATION_TOOLS )
#
## Add unit tests
#add_executable(test_ShapiroTimeDelay "${SRCROOT}${RELATIVITYDIR}/tests/unitTestShapiroTimeDelay.cpp")
#setup_custom_test_program(test_ShapiroTimeDelay "${SRCROOT}${RELATIVITYDIR}")
#target_link_libraries(test_ShapiroTimeDelay tudat_relativity tudat_gravitation tudat_observation_models tudat_ephemerides ${Boost_LIBRARIES})
#
#endif( )
#
#if( TUDAT_BUILD_WITH_ESTIMATION_TOOLS )
#
#add_executable(test_RelativisticAccelerationCorrection "${SRCROOT}${RELATIVITYDIR}/tests/unitTestRelativisticAccelerationCorrection.cpp")
#setup_custom_test_program(test_RelativisticAccelerationCorrection "${SRCROOT}${RELATIVITYDIR}")
#target_link_libraries(test_RelativisticAccelerationCorrection ${Tudat_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

#endif( )
