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
set(ROOTFINDERS_HEADERS 
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/halleyRootFinder.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/newtonRaphson.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/rootFinder.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/secantRootFinder.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/bisection.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/terminationConditions.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/createRootFinder.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/testFunction.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/testFunction1.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/testFunction2.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/testFunction3.h"
  #"${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/testFunctionWithLargeRootDifference.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/testFunctionWithZeroRoot.h"
)

set(ROOTFINDERS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/createRootFinder.cpp"
)


# Add unit test files.
set(ROOTFINDERS_TESTS
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/unitTestBisection.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/unitTestHalleyRootFinder.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/unitTestNewtonRaphson.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/unitTestRootFinders.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/unitTestSecantRootFinder.cpp"
)

# Add static libraries.
add_library(tudat_root_finders STATIC ${ROOTFINDERS_SOURCES} ${ROOTFINDERS_HEADERS})
setup_tudat_library_target(tudat_root_finders "${SRCROOT}${MATHEMATICSDIR}/RootFinders")
set_target_properties(tudat_root_finders PROPERTIES LINKER_LANGUAGE CXX)

# Add unit tests.
add_executable(test_RootFinders ${ROOTFINDERS_TESTS})
setup_custom_test_program(test_RootFinders "${SRCROOT}${MATHEMATICSDIR}/RootFinders")
target_link_libraries(test_RootFinders tudat_root_finders ${Boost_LIBRARIES})
