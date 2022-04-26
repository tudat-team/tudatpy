 #    Copyright (c) 2010-2019, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
# #
#
## Set the source files.
set(SHAPEBASEDMETHODS_SOURCES
  "baseFunctionsHodographicShaping.cpp"
  "createBaseFunctionHodographicShaping.cpp"
  "compositeFunctionHodographicShaping.cpp"
  "hodographicShapingLeg.cpp"
  "baseFunctionsSphericalShaping.cpp"
  "compositeFunctionSphericalShaping.cpp"
  "sphericalShapingLeg.cpp"
  "getRecommendedBaseFunctionsHodographicShaping.cpp"
  )

set(shape_based_methods_HEADERS
  "baseFunctionsHodographicShaping.h"
  "createBaseFunctionHodographicShaping.h"
  "compositeFunctionHodographicShaping.h"
  "hodographicShapingLeg.h"
  "baseFunctionsSphericalShaping.h"
  "compositeFunctionSphericalShaping.h"
  "sphericalShapingLeg.h"
  "getRecommendedBaseFunctionsHodographicShaping.h"
  )

if( TUDAT_BUILD_WITH_PAGMO )
    set(SHAPEBASEDMETHODS_SOURCES
        ${SHAPEBASEDMETHODS_SOURCES}
        "hodographicShapingOptimisationSetup.cpp"
        )

    set(shape_based_methods_HEADERS
        ${shape_based_methods_HEADERS}
        "hodographicShapingOptimisationSetup.cpp"
        )
endif()


# Add library.
TUDAT_ADD_LIBRARY("shape_based_methods"
        "${SHAPEBASEDMETHODS_SOURCES}"
        "${shape_based_methods_HEADERS}"
        )
