# \file CMakeLists.txt
 #    This is the CMake file for the Bodies sub-directory.
 #
 #    Path              : /Astrodynamics/Bodies/
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

# Set the source files.
set(BODIES_SOURCES
  "${SRCROOT}${BODIESDIR}/planet.cpp"

  "${SRCROOT}${BODIESDIR}/body.h"
  "${SRCROOT}${BODIESDIR}/celestialBody.h"
  "${SRCROOT}${BODIESDIR}/vehicle.h"
  "${SRCROOT}${BODIESDIR}/vehicleExternalModel.h"
)

# Set the header files.
set(BODIES_HEADERS 
  "${BODIESDIR}/planet.h"
)

# Add static libraries.
add_library(tudat_bodies STATIC ${BODIES_SOURCES})
setup_tudat_library_target(tudat_bodies "${SRCROOT}${BODIESDIR}")
setup_headers("${BODIES_HEADERS}" "${BODIESDIR}")

# End of file.
