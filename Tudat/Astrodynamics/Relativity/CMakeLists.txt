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
 #      120321    K. Kumar          Added header, source and unit test files for state-derivative
 #                                  of CRTBP.
 #

# Set the source files.
set(RELATIVITY_SOURCES
  "${SRCROOT}${RELATIVITYDIR}/relativisticLightTimeCorrection.cpp"
)

# Set the header files.
set(RELATIVITY_HEADERS 
  "${SRCROOT}${RELATIVITYDIR}/relativisticLightTimeCorrection.h"
)

# Add static libraries.
add_library(tudat_relativity STATIC ${RELATIVITY_SOURCES} ${RELATIVITY_HEADERS})
setup_tudat_library_target(tudat_relativity "${SRCROOT}${RELATIVITYDIR}")


# Add unit tests
add_executable(test_ShapiroTimeDelay "${SRCROOT}${RELATIVITYDIR}/UnitTests/unitTestShapiroTimeDelay.cpp")
setup_custom_test_program(test_ShapiroTimeDelay "${SRCROOT}${RELATIVITYDIR}")
target_link_libraries(test_ShapiroTimeDelay tudat_relativity tudat_gravitation tudat_observation_models tudat_ephemerides ${Boost_LIBRARIES})
