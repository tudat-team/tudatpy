 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Set the source files.
set(RELATIVITY_SOURCES
  "${SRCROOT}${RELATIVITYDIR}/relativisticLightTimeCorrection.cpp"
  "${SRCROOT}${RELATIVITYDIR}/relativisticAccelerationCorrection.cpp"
  "${SRCROOT}${RELATIVITYDIR}/relativisticTimeConversion.cpp"
  "${SRCROOT}${RELATIVITYDIR}/metric.cpp"
)

# Set the header files.
set(RELATIVITY_HEADERS 
  "${SRCROOT}${RELATIVITYDIR}/relativisticLightTimeCorrection.h"
  "${SRCROOT}${RELATIVITYDIR}/relativisticAccelerationCorrection.h"
  "${SRCROOT}${RELATIVITYDIR}/relativisticTimeConversion.h"
  "${SRCROOT}${RELATIVITYDIR}/metric.h"
)

# Add static libraries.
add_library(tudat_relativity STATIC ${RELATIVITY_SOURCES} ${RELATIVITY_HEADERS})
setup_tudat_library_target(tudat_relativity "${SRCROOT}${RELATIVITYDIR}")


# Add unit tests
add_executable(test_ShapiroTimeDelay "${SRCROOT}${RELATIVITYDIR}/UnitTests/unitTestShapiroTimeDelay.cpp")
setup_custom_test_program(test_ShapiroTimeDelay "${SRCROOT}${RELATIVITYDIR}")
target_link_libraries(test_ShapiroTimeDelay tudat_relativity tudat_gravitation tudat_observation_models tudat_ephemerides ${Boost_LIBRARIES})

if( COMPILE_PROPAGATION_TESTS )

add_executable(test_RelativisticAccelerationCorrection "${SRCROOT}${RELATIVITYDIR}/UnitTests/unitTestRelativisticAccelerationCorrection.cpp")
setup_custom_test_program(test_RelativisticAccelerationCorrection "${SRCROOT}${RELATIVITYDIR}")
target_link_libraries(test_RelativisticAccelerationCorrection ${TUDAT_PROPAGATION_LIBRARIES} ${Boost_LIBRARIES})

endif( )
