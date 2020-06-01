#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Add source files.
set(geometric_SOURCES
        "capsule.cpp"
        "compositeSurfaceGeometry.cpp"
        "conicalFrustum.cpp"
        "lawgsPartGeometry.cpp"
        "quadrilateralMeshedSurfaceGeometry.cpp"
        "singleSurfaceGeometry.cpp"
        "sphereSegment.cpp"
        "torus.cpp"
        "geometricShapesToFile.cpp"
        )

# Add header files.
set(geometric_HEADERS
        "capsule.h"
        "compositeSurfaceGeometry.h"
        "conicalFrustum.h"
        "lawgsPartGeometry.h"
        "quadrilateralMeshedSurfaceGeometry.h"
        "singleSurfaceGeometry.h"
        "sphereSegment.h"
        "surfaceGeometry.h"
        "torus.h"
        "geometricShapesToFile.h"
        )

# Add library.
TUDAT_ADD_LIBRARY("geometric_shapes"
        "${geometric_SOURCES}"
        "${geometric_HEADERS}"
        PRIVATE_LINKS "${Boost_LIBRARIES}"
        PRIVATE_INCLUDES "${EIGEN3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}"
        )

# Add tests.
TUDAT_ADD_TEST_CASE(LawgsSurfaceGeometry PRIVATE_LINKS tudat_geometric_shapes tudat_basic_mathematics)
