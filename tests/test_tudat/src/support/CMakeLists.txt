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
set(test_support_SOURCES
	"numericalObservationPartial.cpp"
	"observationPartialTestFunctions.cpp"
        )

# Set the header files.
set(test_support_HEADERS
	"numericalIntegratorTests.h"
	"numericalObservationPartial.h"
	"observationPartialTestFunctions.h"
	"sumLmkTestHelpers.h"
        )

TUDAT_ADD_TEST_LIBRARY("support"
        "${test_support_SOURCES}"
        "${test_support_HEADERS}")

# Public Tudat model headers include archive-independent Cereal declarations,
# also when building without the serialization implementations.
target_link_libraries(tudat_test_support PUBLIC cereal::cereal)

