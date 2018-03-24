#   Copyright (c) 2010-2018, Delft University of Technology
#   All rigths reserved
#
#   This file is part of the Tudat. Redistribution and use in source and
#   binary forms, with or without modification, are permitted exclusively
#   under the terms of the Modified BSD license. You should have received
#   a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# Add header files.
set(NUMERICAL_QUADRATURE_HEADERS
  "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature/numericalQuadrature.h"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature/trapezoidQuadrature.h"
  "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature/gaussianQuadrature.h"
)

set(NUMERICAL_QUADRATURE_SOURCES
  "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature/dummyNumericalQuadrature.cpp"
)


# Add static libraries.
add_library(tudat_numerical_quadrature STATIC ${NUMERICAL_QUADRATURE_SOURCES} ${NUMERICAL_QUADRATURE_HEADERS} )
setup_tudat_library_target(tudat_numerical_quadrature "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature")
set_target_properties(tudat_numerical_quadrature PROPERTIES LINKER_LANGUAGE CXX)

# Add unit tests.
add_executable(test_TrapezoidalIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature/UnitTests/unitTestTrapezoidQuadrature.cpp")
setup_custom_test_program(test_TrapezoidalIntegrator "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature")
target_link_libraries(test_TrapezoidalIntegrator tudat_numerical_quadrature ${Boost_LIBRARIES})

add_executable(test_GaussianQuadrature "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature/UnitTests/unitTestGaussianQuadrature.cpp")
setup_custom_test_program(test_GaussianQuadrature "${SRCROOT}${MATHEMATICSDIR}/NumericalQuadrature")
target_link_libraries(test_GaussianQuadrature tudat_input_output tudat_numerical_quadrature ${Boost_LIBRARIES})

