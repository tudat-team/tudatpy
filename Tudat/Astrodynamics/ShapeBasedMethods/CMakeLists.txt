 #    Copyright (c) 2010-2018, Delft University of Technology
 #    All rigths reserved
 #
 #    This file is part of the Tudat. Redistribution and use in source and
 #    binary forms, with or without modification, are permitted exclusively
 #    under the terms of the Modified BSD license. You should have received
 #    a copy of the license with this file. If not, please or visit:
 #    http://tudat.tudelft.nl/LICENSE.
 #

# Set the source files.
set(SHAPEBASEDMETHODS_SOURCES
  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsHodographicShaping.cpp"
  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionHodographicShaping.cpp"
  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/hodographicShaping.cpp"
  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsSphericalShaping.cpp"
  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionSphericalShaping.cpp"
  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/sphericalShaping.cpp"
)

# Set the header files.
set(SHAPEBASEDMETHODS_HEADERS
  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsHodographicShaping.h"
  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionHodographicShaping.h"
  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/hodographicShaping.h"
  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/baseFunctionsSphericalShaping.h"
  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/compositeFunctionSphericalShaping.h"
  "${SRCROOT}${SHAPEBASEDMETHODSDIR}/sphericalShaping.h"
)

# Add static libraries.
add_library(tudat_shape_based_methods STATIC ${SHAPEBASEDMETHODS_SOURCES} ${SHAPEBASEDMETHODS_HEADERS})
setup_tudat_library_target(tudat_shape_based_methods "${SRCROOT}${SHAPEBASEDMETHODSDIR}")

include(tudatLinkLibraries)

## Add unit tests.
add_executable(test_HodographicShaping "${SRCROOT}${SHAPEBASEDMETHODSDIR}/UnitTests/unitTestHodographicShaping.cpp")
setup_custom_test_program(test_HodographicShaping "${SRCROOT}${SHAPEBASEDMETHODSDIR}")
target_link_libraries(test_HodographicShaping tudat_shape_based_methods tudat_numerical_quadrature ${TUDAT_PROPAGATION_LIBRARIES}  ${Boost_LIBRARIES})

add_executable(test_SphericalShaping "${SRCROOT}${SHAPEBASEDMETHODSDIR}/UnitTests/unitTestSphericalShaping.cpp")
setup_custom_test_program(test_SphericalShaping "${SRCROOT}${SHAPEBASEDMETHODSDIR}")
target_link_libraries(test_SphericalShaping tudat_shape_based_methods tudat_numerical_quadrature ${TUDAT_PROPAGATION_LIBRARIES}  ${Boost_LIBRARIES})

