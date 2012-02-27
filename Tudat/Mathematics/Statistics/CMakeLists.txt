 #    Copyright (c) 2010-2012 Delft University of Technology.
 #
 #    This software is protected by national and international copyright.
 #    Any unauthorized use, reproduction or modification is unlawful and
 #    will be prosecuted. Commercial and non-private application of the
 #    software in any form is strictly prohibited unless otherwise granted
 #    by the authors.
 #
 #    The code is provided without any warranty; without even the implied
 #    warranty of merchantibility or fitness for a particular purpose.
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
setup_headers("${STATISTICS_HEADERS}" "${MATHEMATICSDIR}/Statistics")

# Add unit tests.
add_executable(test_SimpleLinearRegression "${SRCROOT}${MATHEMATICSDIR}/Statistics/unitTestSimpleLinearRegression.cpp")
setup_custom_test_program(test_SimpleLinearRegression "${SRCROOT}${MATHEMATICSDIR}/Statistics")
target_link_libraries(test_SimpleLinearRegression tudat_statistics)

add_executable(test_BasicStatistics "${SRCROOT}${MATHEMATICSDIR}/Statistics/UnitTests/unitTestBasicStatistics.cpp")
setup_custom_test_program(test_BasicStatistics "${SRCROOT}${MATHEMATICSDIR}/Statistics")
target_link_libraries(test_BasicStatistics tudat_statistics ${Boost_LIBRARIES})
