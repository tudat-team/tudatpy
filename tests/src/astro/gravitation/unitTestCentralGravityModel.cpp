/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "tudat/astro/gravitation/centralGravityModel.h"
#include "tudat/basics/testMacros.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"

#include <cmath>
#include <limits>

#include <Eigen/Core>

namespace tudat
{
namespace unit_tests
{

//! Test the functionality of the central gravity model class.
BOOST_AUTO_TEST_SUITE( test_central_gravity_model )

//! Test the computation of the gravitational potential and acceleration, using the wrapper class.
BOOST_AUTO_TEST_CASE( testWrapperClassPotentialAndAcceleration )
{
    // Short-cuts.
    using namespace gravitation;

    // Define gravitational parameter of Earth [m^3 s^-2]. The value is obtained from the Earth
    // Gravitational Model 2008 as described by Mathworks [2012].
    const double gravitationalParameter = 3.986004418e14;

    // Define arbitrary Cartesian position [m].
    const Eigen::Vector3d position( 7.0e6, 8.0e6, 9.0e6 );

    // Declare spherical harmonics gravitational acceleration class object.
    std::shared_ptr< CentralGravitationalAccelerationModel3d > earthGravity
            = std::make_shared< CentralGravitationalAccelerationModel3d >(
                [ & ]( Eigen::Vector3d& input ){ input = position; }, gravitationalParameter );
    earthGravity->resetUpdatePotential( true );
    earthGravity->updateMembers( );

    // Generate expected result for the potential
    double expectedPotential = gravitationalParameter / position.norm( );
    // Check if expected result matches computed result.
    BOOST_CHECK_EQUAL( expectedPotential, earthGravity->getCurrentPotential( ) );

    // Determine the expected result for the gradient of the potential of myPlanet
    const Eigen::Vector3d expectedAcceleration = - gravitationalParameter
            / std::pow( position.norm( ), 3.0 ) * position;
    // Check if expected result matches computed result.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedAcceleration, earthGravity->getAcceleration( ),
                                       std::numeric_limits< double >::epsilon( ) );

}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace tudat
} // namespace unit_tests
