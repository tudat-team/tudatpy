 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Add header files.
set(FILTERS_HEADERS
#  "${SRCROOT}${MATHEMATICSDIR}/Filters/extendedKalmanFilter.h"
  "${SRCROOT}${MATHEMATICSDIR}/Filters/kalmanFilter.h"
  "${SRCROOT}${MATHEMATICSDIR}/Filters/linearKalmanFilter.h"
#  "${SRCROOT}${MATHEMATICSDIR}/Filters/unscentedKalmanFilter.h"
)

set(FILTERS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/Filters/linearKalmanFilter.cpp"
)


# Add unit test files.
set(FILTERS_TESTS
  "${SRCROOT}${MATHEMATICSDIR}/Filters/UnitTests/unitTestLinearKalmanFilter.cpp"
)

# Add static libraries.
add_library(tudat_filters STATIC ${FILTERS_SOURCES} ${FILTERS_HEADERS})
setup_tudat_library_target(tudat_filters "${SRCROOT}${MATHEMATICSDIR}/Filters")
set_target_properties(tudat_filters PROPERTIES LINKER_LANGUAGE CXX)

# Add unit tests.
add_executable(test_Filters ${FILTERS_TESTS})
setup_custom_test_program(test_Filters "${SRCROOT}${MATHEMATICSDIR}/Filters")
target_link_libraries(test_Filters tudat_filters tudat_numerical_integrators tudat_statistics tudat_input_output ${Boost_LIBRARIES})
