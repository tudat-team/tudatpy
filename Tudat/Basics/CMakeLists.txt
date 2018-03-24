 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Add source files.
set(BASICSDIR_SOURCES
  "${SRCROOT}${BASICSDIR}/utilities.cpp"
)

# Add header files.
set(BASICSDIR_HEADERS 
  "${SRCROOT}${BASICSDIR}/utilities.h"
  "${SRCROOT}${BASICSDIR}/testMacros.h"
  "${SRCROOT}${BASICSDIR}/utilityMacros.h"
  "${SRCROOT}${BASICSDIR}/timeType.h"
  "${SRCROOT}${BASICSDIR}/basicTypedefs.h"
)

# Add unit test files.
set(BASICSDIR_UNITTESTS
  "${SRCROOT}${BASICSDIR}/UnitTests/unitTestBasics.cpp"
  "${SRCROOT}${BASICSDIR}/UnitTests/unitTestTestMacros.cpp"
)

# Add static libraries.
add_library(tudat_basics STATIC ${BASICSDIR_SOURCES} ${BASICSDIR_HEADERS})
setup_tudat_library_target(tudat_basics "${SRCROOT}${BASICSDIR}")

add_executable(test_TimeTypes "${SRCROOT}${BASICSDIR}/UnitTests/unitTestTimeTypes.cpp")
setup_custom_test_program(test_TimeTypes "${SRCROOT}${BASICSDIR}")
target_link_libraries(test_TimeTypes ${Boost_LIBRARIES})

