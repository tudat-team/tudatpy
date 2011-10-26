# \file CMakeLists.txt
 #    This is the CMake file for the MomentModels sub-directory.
 #
 #    Path              : /Astrodynamics/MomentModels/
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
set(MOMENTMODELS_SOURCES
  "${SRCROOT}${MOMENTMODELSDIR}/aerodynamicMoment.cpp"

  "${SRCROOT}${MOMENTMODELSDIR}/momentModel.h"
)

# Add header files.
set(MOMENTMODELS_HEADERS 
  "${MOMENTMODELSDIR}/aerodynamicMoment.h"
)

# Add static libraries.
add_library(tudat_moment_models STATIC ${MOMENTMODELS_SOURCES})
setup_tudat_library_target(tudat_moment_models "${SRCROOT}${MOMENTMODELSDIR}")
setup_headers("${MOMENTMODELS_HEADERS}" "${MOMENTMODELSDIR}")

# End of file.
