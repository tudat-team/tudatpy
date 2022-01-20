#    Copyright (c) 2010-2019, Delft University of Technology
#    All rigths reserved
#
#    This file is part of the Tudat. Redistribution and use in source and
#    binary forms, with or without modification, are permitted exclusively
#    under the terms of the Modified BSD license. You should have received
#    a copy of the license with this file. If not, please or visit:
#    http://tudat.tudelft.nl/LICENSE.
#
#

# # Add subdirectories.
add_subdirectory(shape_based)

#### Add unit tests.
#if( TUDAT_WITH_PAGMO )
#    TUDAT_ADD_TEST_CASE(SimsFlanagan PRIVATE_LINKS tudat_low_thrust_trajectories tudat_shape_based_methods tudat_numerical_quadrature pagmo ${Tudat_PROPAGATION_LIBRARIES}  ${Boost_LIBRARIES})
#endif( )
