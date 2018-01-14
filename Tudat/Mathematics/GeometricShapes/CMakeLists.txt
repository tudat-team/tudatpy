 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

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
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/geometricShapesToFile.cpp"
)

# Add header files.
set(GEOMETRICSHAPES_HEADERS 
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/capsule.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/compositeSurfaceGeometry.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/conicalFrustum.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/lawgsPartGeometry.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/quadrilateralMeshedSurfaceGeometry.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/singleSurfaceGeometry.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/sphereSegment.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/surfaceGeometry.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/torus.h"
  "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/geometricShapesToFile.h"
)

# Add static libraries.
add_library(tudat_geometric_shapes STATIC ${GEOMETRICSHAPES_SOURCES})
setup_tudat_library_target(tudat_geometric_shapes "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes")

# Add unit tests.
add_executable(test_LawgsSurfaceGeometry "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes/UnitTests/unitTestLawgsSurfaceGeometry.cpp")
setup_custom_test_program(test_LawgsSurfaceGeometry "${SRCROOT}${MATHEMATICSDIR}/GeometricShapes")
target_link_libraries(test_LawgsSurfaceGeometry tudat_geometric_shapes tudat_basic_mathematics ${Boost_LIBRARIES})
