#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Add header files.
set(filters_HEADERS
        "createFilter.h"
        "extendedKalmanFilter.h"
        "filter.h"
        "kalmanFilter.h"
        "linearKalmanFilter.h"
        "unscentedKalmanFilter.h"
        "tests/controlClass.h"
        )

set(filters_SOURCES
        "filter.cpp"
        )

# Add library.
TUDAT_ADD_LIBRARY("filters"
        "${filters_SOURCES}"
        "${filters_HEADERS}"
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )

# Add test library interface.
add_library(filters_test_library INTERFACE)
target_link_libraries(filters_test_library
        INTERFACE
        tudat_filters
        tudat_numerical_integrators
        tudat_statistics
        tudat_basics
        tudat_basic_mathematics
        tudat_input_output)

# Add tests.
TUDAT_ADD_TEST_CASE(LinearKalmanFilter PRIVATE_LINKSfilters_test_library)
TUDAT_ADD_TEST_CASE(ExtendedKalmanFilter PRIVATE_LINKS filters_test_library)
TUDAT_ADD_TEST_CASE(UnscentedKalmanFilter PRIVATE_LINKS filters_test_library)
