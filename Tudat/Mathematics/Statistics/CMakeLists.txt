 #    Copyright (c) 2010-2012, Delft University of Technology
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
 #

# Add source files.
set(STATISTICS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/basicStatistics.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/simpleLinearRegression.cpp"
)

# Add header files.
set(STATISTICS_HEADERS 
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/basicStatistics.h"
  "${SRCROOT}${MATHEMATICSDIR}/Statistics/simpleLinearRegression.h"
)

# Add static libraries.
add_library(tudat_statistics STATIC ${STATISTICS_SOURCES} ${STATISTICS_HEADERS})
setup_tudat_library_target(tudat_statistics "${SRCROOT}${MATHEMATICSDIR}/Statistics")

# Add unit tests.
add_executable(test_SimpleLinearRegression "${SRCROOT}${MATHEMATICSDIR}/Statistics/UnitTests/unitTestSimpleLinearRegression.cpp")
setup_custom_test_program(test_SimpleLinearRegression "${SRCROOT}${MATHEMATICSDIR}/Statistics")
target_link_libraries(test_SimpleLinearRegression tudat_statistics ${Boost_LIBRARIES})

add_executable(test_BasicStatistics "${SRCROOT}${MATHEMATICSDIR}/Statistics/UnitTests/unitTestBasicStatistics.cpp")
setup_custom_test_program(test_BasicStatistics "${SRCROOT}${MATHEMATICSDIR}/Statistics")
target_link_libraries(test_BasicStatistics tudat_statistics ${Boost_LIBRARIES})
