# \file CMakeLists.txt
 #    This is the CMake file for the GeometricShapes sub-directory.
 #
 #    Path              : /Mathematics/GeometricShapes/
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
set(GEOMETRICSHAPES_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/capsule.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/compositeSurfaceGeometry.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/conicalFrustum.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/lawgsPartGeometry.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/quadrilateralMeshedSurfaceGeometry.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/singleSurfaceGeometry.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/sphereSegment.cpp"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/torus.cpp"

  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/geometricShape.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/surfaceGeometry.h"
)

# Add header files.
set(GEOMETRICSHAPES_HEADERS 
  "${MATHEMATICSDIR}/GeometricShapes/capsule.h"
  "${MATHEMATICSDIR}/GeometricShapes/compositeSurfaceGeometry.h"
  "${MATHEMATICSDIR}/GeometricShapes/conicalFrustum.h"
  "${MATHEMATICSDIR}/GeometricShapes/lawgsPartGeometry.h"
  "${MATHEMATICSDIR}/GeometricShapes/quadrilateralMeshedSurfaceGeometry.h"
  "${MATHEMATICSDIR}/GeometricShapes/singleSurfaceGeometry.h"
  "${MATHEMATICSDIR}/GeometricShapes/sphereSegment.h"
  "${MATHEMATICSDIR}/GeometricShapes/torus.h"
)

# Add static libraries.
add_library(tudat_geometric_shapes STATIC ${GEOMETRICSHAPES_SOURCES})
setup_tudat_library_target(tudat_geometric_shapes "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes")
setup_headers("${GEOMETRICSHAPES_HEADERS}" "${MATHEMATICSDIR}/GeometricShapes")

# Add unit tests.
add_executable(test_LawgsSurfaceGeometry "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/unitTestLawgsSurfaceGeometry.cpp")
setup_custom_test_program(test_LawgsSurfaceGeometry "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes")
target_link_libraries(test_LawgsSurfaceGeometry tudat_input tudat_output tudat_mathematics tudat_astrodynamics tudat_basics)

# End of file.
