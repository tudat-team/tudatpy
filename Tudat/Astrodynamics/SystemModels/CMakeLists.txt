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
set(SYSTEMMODELS_SOURCES
  "${SRCROOT}${SYSTEMMODELSDIR}/engineModel.cpp"
)

# Set the header files.
set(SYSTEMMODELS_HEADERS 
  "${SRCROOT}${SYSTEMMODELSDIR}/engineModel.h"
  "${SRCROOT}${SYSTEMMODELSDIR}/vehicleSystems.h"
)

# Add static libraries.
add_library(tudat_system_models STATIC ${SYSTEMMODELS_SOURCES} ${SYSTEMMODELS_HEADERS})
setup_tudat_library_target(tudat_system_models "${SRCROOT}${SYSTEMMODELSDIR}")


