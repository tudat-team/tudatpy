# \file CMakeLists.txt
 #    This is the CMake file for the RootFindingMethods sub-directory.
 #
 #    Path              : /Mathematics/RootFindingMethods/
 #    Version           : 3
 #    Check status      : Checked
 #
 #    Author            : S.M. Persson
 #    Affiliation       : McGill University, Department of Mechanical Engineering
 #    E-mail address    : mikael.s.persson@gmail.com
 #
 #    Author            : K. Kumar
 #    Affiliation       : Delft University of Technology
 #    E-mail address    : K.Kumar@tudelft.nl
 #
 #    Author            : S.M. Persson
 #    Affiliation       : McGill University, Department of Mechanical Engineering
 #    E-mail address    : mikael.s.persson@gmail.com
 #
 #    Date created      : 20 August, 2011
 #    Last modified     : 26 October, 2011
 #
 #    References
 #
 #    Notes
 #
 #    Copyright (c) 2010-2011 Delft University of Technology.
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
 #/

# Add source files.
set(ROOTFINDINGMETHODS_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/RootFindingMethods/newtonRaphson.cpp"
)

# Add header files.
set(ROOTFINDINGMETHODS_HEADERS 
  "${SRCROOT}${MATHEMATICSDIR}/RootFindingMethods/newtonRaphson.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFindingMethods/newtonRaphsonAdaptor.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFindingMethods/newtonRaphsonBase.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFindingMethods/rootFinder.h"
  "${SRCROOT}${MATHEMATICSDIR}/RootFindingMethods/rootFinderBase.h"
)

# Add static libraries.
add_library(tudat_root_finding_methods STATIC ${ROOTFINDINGMETHODS_SOURCES} ${ROOTFINDINGMETHODS_HEADERS})
setup_tudat_library_target(tudat_root_finding_methods "${SRCROOT}${MATHEMATICSDIR}/RootFindingMethods")
setup_headers("${ROOTFINDINGMETHODS_HEADERS}" "${MATHEMATICSDIR}/RootFindingMethods")

# Add unit tests.
add_executable(test_NewtonRaphson "${SRCROOT}${MATHEMATICSDIR}/RootFindingMethods/unitTestNewtonRaphson.cpp")
setup_custom_test_program(test_NewtonRaphson "${SRCROOT}${MATHEMATICSDIR}/RootFindingMethods")
target_link_libraries(test_NewtonRaphson tudat_input tudat_output tudat_mathematics tudat_astrodynamics tudat_basics)

# End of file.
