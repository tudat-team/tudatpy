#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#

# TODO: Add shape-based low-thrust tests
## Add unit tests.
add_executable(test_HodographicShaping "${SRCROOT}${SHAPEBASEDMETHODSDIR}/tests/unitTestHodographicShaping.cpp")
setup_custom_test_program(test_HodographicShaping "${SRCROOT}${SHAPEBASEDMETHODSDIR}")
target_link_libraries(test_HodographicShaping tudat_shape_based_methods tudat_low_thrust_trajectories tudat_numerical_quadrature ${Tudat_PROPAGATION_LIBRARIES}  ${Boost_LIBRARIES})

add_executable(test_SphericalShaping "${SRCROOT}${SHAPEBASEDMETHODSDIR}/tests/unitTestSphericalShaping.cpp")
setup_custom_test_program(test_SphericalShaping "${SRCROOT}${SHAPEBASEDMETHODSDIR}")
target_link_libraries(test_SphericalShaping tudat_shape_based_methods tudat_low_thrust_trajectories tudat_numerical_quadrature ${Tudat_PROPAGATION_LIBRARIES}  ${Boost_LIBRARIES})
