 #    Copyright (c) 2010-2014, Delft University of Technology
 #    All rights reserved.
 #
 #    Redistribution and use in source and binary forms, with or without modification, are
 #    permitted provided that the following conditions are met:
 #      - Redistributions of source code must retain the above copyright notice, this list of
 #        conditions and the following disclaimer.
 #      - Redistributions in binary form must reproduce the above copyright notice, this list of
 #        conditions and the following disclaimer in the documentation and/or other materials
 #        provided with the distribution.
 #      - Neither the name of the Delft University of Technology nor the names of its contributors
 #        may be used to endorse or promote products derived from this software without specific
 #        prior written permission.
 #
 #    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 #    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 #    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 #    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 #    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 #    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 #    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 #    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 #    OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 #    Changelog
 #      YYMMDD    Author            Comment
 #      110820    S.M. Persson      File created.
 #      111025    K. Kumar          Adapted file to work with Revision 194.
 #      111026    K. Kumar          Adapted file so all headers show in project tree in Qt Creator.
 #      120726    S. Billemont      Added files regarding new root-finder methods and unit tests.
 #      120813    P. Musegaas       Added Newton Raphson unit test executable, added additional test
 #                                  functions.
 #	120920    K. Kumar          Changed sub-directory name to RootFinders.
 #
 #    References
 #
 #    Notes
 #

# Add source files.
set(ROOTFINDERS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/void.cpp"
)

# Add header files.
set(ROOTFINDERS_HEADERS 
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/newtonRaphson.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/rootFinder.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/terminationConditions.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/testFunction.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/testFunction1.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/testFunction2.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/testFunction3.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/testFunctionWithLargeRootDifference.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/testFunctionWithZeroRoot.h"
)

# Add unit test files.
set(ROOTFINDERS_TESTS
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/unitTestNewtonRaphson.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/RootFinders/UnitTests/unitTestRootFinders.cpp"
)

# Add static libraries.
add_library(tudat_root_finders STATIC ${ROOTFINDERS_SOURCES} ${ROOTFINDERS_HEADERS})
setup_tudat_library_target(tudat_root_finders "${SRCROOT}${MATHEMATICSDIR}/RootFinders")

# Add unit tests.
add_executable(test_RootFinders ${ROOTFINDERS_TESTS})
setup_custom_test_program(test_RootFinders "${SRCROOT}${MATHEMATICSDIR}/RootFinders")
target_link_libraries(test_RootFinders tudat_root_finders ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})

# Add unit tests.
add_executable(test_NewtonRaphson ${ROOTFINDERS_TESTS})
setup_custom_test_program(test_NewtonRaphson "${SRCROOT}${MATHEMATICSDIR}/RootFinders")
target_link_libraries(test_NewtonRaphson tudat_root_finders ${TUDAT_CORE_LIBRARIES} ${Boost_LIBRARIES})
