# \file CMakeLists.txt
 #    This is the CMake file for the Basics sub-directory.
 #
 #    Path              : /Basics/
 #    Version           : 1
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
 #    Last modified     : 25 October, 2011
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
 #/

# Add source files.
set(BASICS_SOURCES
  "${SRCROOT}${BASICSDIR}/basicFunctions.cpp"
)

# Add header files.
set(BASICS_HEADERS 
  "${BASICSDIR}/basicFunctions.h"
)

# Add static libraries.
add_library(tudat_basics STATIC ${BASICS_SOURCES})
setup_tudat_library_target(tudat_basics "${SRCROOT}${BASICSDIR}")
setup_headers("${BASICS_HEADERS}" "${BASICSDIR}")

# End of file.
