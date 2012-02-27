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
set(REFERENCEFRAMES_SOURCES
  "${SRCROOT}${REFERENCEFRAMESDIR}/referenceFrameTransformations.cpp"
)

# Add header files.
set(REFERENCEFRAMES_HEADERS 
  "${SRCROOT}${REFERENCEFRAMESDIR}/referenceFrameTransformations.h"
)

# Add static libraries.
add_library(tudat_reference_frames STATIC ${REFERENCEFRAMES_SOURCES} ${REFERENCEFRAMES_HEADERS})
setup_tudat_library_target(tudat_reference_frames "${SRCROOT}${REFERENCEFRAMESDIR}")
setup_headers("${REFERENCEFRAMES_HEADERS}" "${REFERENCEFRAMESDIR}")

# Add unit tests.
add_executable(test_ReferenceFrameTransformations "${SRCROOT}${REFERENCEFRAMESDIR}/unitTestReferenceFrameTransformations.cpp")
setup_custom_test_program(test_ReferenceFrameTransformations "${SRCROOT}${REFERENCEFRAMESDIR}")
target_link_libraries(test_ReferenceFrameTransformations tudat_reference_frames ${TUDAT_CORE_LIBRARIES})
