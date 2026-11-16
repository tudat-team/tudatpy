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

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>

#include "tudat/astro/gravitation/centralGravityModel.h"
#include "tudat/basics/testMacros.h"
#include "tudat/astro/gravitation/gravityFieldModel.h"
#include "tudat/astro/gravitation/gravitationalParameterRateAcceleration.h"

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
    std::shared_ptr< CentralGravitationalAccelerationModel3d > earthGravity = std::make_shared< CentralGravitationalAccelerationModel3d >(
            [ & ]( Eigen::Vector3d& input ) { input = position; }, gravitationalParameter );
    earthGravity->resetUpdatePotential( true );
    earthGravity->updateMembers( );

    // Generate expected result for the potential
    double expectedPotential = gravitationalParameter / position.norm( );
    // Check if expected result matches computed result.
    BOOST_CHECK_EQUAL( expectedPotential, earthGravity->getCurrentPotential( ) );

    // Determine the expected result for the gradient of the potential of myPlanet
    const Eigen::Vector3d expectedAcceleration = -gravitationalParameter / std::pow( position.norm( ), 3.0 ) * position;
    // Check if expected result matches computed result.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedAcceleration, earthGravity->getAcceleration( ), std::numeric_limits< double >::epsilon( ) );
}

//! Test the gravitational parameter rate acceleration against the difference of point-mass accelerations with mu( t ) and mu( t0 ).
BOOST_AUTO_TEST_CASE( testGravitationalParameterRateAcceleration )
{
    using namespace gravitation;

    const double gravitationalParameter = 1.32712440018e20;
    const double relativeRate = -1.0e-10;
    const double referenceEpoch = 1.0e8;
    const double currentTime = 2.0e8;
    const Eigen::Vector3d acceleratedBodyPosition( 1.2e11, -0.7e11, 0.3e11 );
    const Eigen::Vector3d exertingBodyPosition( 1.0e9, 2.0e9, -3.0e9 );

    GravitationalParameterRateAcceleration rateAcceleration( [ & ]( ) { return acceleratedBodyPosition; },
                                                             [ & ]( ) { return exertingBodyPosition; },
                                                             [ & ]( ) { return gravitationalParameter; },
                                                             relativeRate,
                                                             referenceEpoch );
    rateAcceleration.updateMembers( currentTime );

    const double currentGravitationalParameter = gravitationalParameter * ( 1.0 + relativeRate * ( currentTime - referenceEpoch ) );
    const Eigen::Vector3d expectedAcceleration =
            computeGravitationalAcceleration( acceleratedBodyPosition, currentGravitationalParameter, exertingBodyPosition ) -
            computeGravitationalAcceleration( acceleratedBodyPosition, gravitationalParameter, exertingBodyPosition );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedAcceleration, rateAcceleration.getAcceleration( ), 1.0e-12 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
