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
  "${SRCROOT}${MATHEMATICSDIR}/Filters/createFilter.h"
  "${SRCROOT}${MATHEMATICSDIR}/Filters/extendedKalmanFilter.h"
  "${SRCROOT}${MATHEMATICSDIR}/Filters/filter.h"
  "${SRCROOT}${MATHEMATICSDIR}/Filters/kalmanFilter.h"
  "${SRCROOT}${MATHEMATICSDIR}/Filters/linearKalmanFilter.h"
  "${SRCROOT}${MATHEMATICSDIR}/Filters/unscentedKalmanFilter.h"
  "${SRCROOT}${MATHEMATICSDIR}/Filters/UnitTests/controlClass.h"
)

set(FILTERS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/Filters/filter.cpp"
)

# Add static libraries.
add_library(tudat_filters STATIC ${FILTERS_SOURCES} ${FILTERS_HEADERS})
setup_tudat_library_target(tudat_filters "${SRCROOT}${MATHEMATICSDIR}/Filters")
set_target_properties(tudat_filters PROPERTIES LINKER_LANGUAGE CXX)

# Add unit tests.
add_executable(test_LinearKalmanFilter "${SRCROOT}${MATHEMATICSDIR}/Filters/UnitTests/unitTestLinearKalmanFilter.cpp")
setup_custom_test_program(test_LinearKalmanFilter "${SRCROOT}${MATHEMATICSDIR}/Filters")
target_link_libraries(test_LinearKalmanFilter tudat_filters tudat_numerical_integrators tudat_statistics tudat_basics tudat_basic_mathematics
    tudat_input_output ${Boost_LIBRARIES})

add_executable(test_ExtendedKalmanFilter "${SRCROOT}${MATHEMATICSDIR}/Filters/UnitTests/unitTestExtendedKalmanFilter.cpp")
setup_custom_test_program(test_ExtendedKalmanFilter "${SRCROOT}${MATHEMATICSDIR}/Filters")
target_link_libraries(test_ExtendedKalmanFilter tudat_filters tudat_numerical_integrators tudat_statistics tudat_basics tudat_basic_mathematics
    tudat_input_output ${Boost_LIBRARIES})

add_executable(test_UnscentedKalmanFilter "${SRCROOT}${MATHEMATICSDIR}/Filters/UnitTests/unitTestUnscentedKalmanFilter.cpp")
setup_custom_test_program(test_UnscentedKalmanFilter "${SRCROOT}${MATHEMATICSDIR}/Filters")
target_link_libraries(test_UnscentedKalmanFilter tudat_filters tudat_numerical_integrators tudat_statistics tudat_basics tudat_basic_mathematics
    tudat_input_output ${Boost_LIBRARIES})
