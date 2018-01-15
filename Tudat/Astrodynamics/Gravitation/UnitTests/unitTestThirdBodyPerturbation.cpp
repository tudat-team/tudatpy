/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include <limits>
#include <vector>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/thirdBodyPerturbation.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_third_body_perturbation )


//! Test if perturbational acceleration is computed correctly.
BOOST_AUTO_TEST_CASE( testComputationOfThirdBodyPerturbation )
{
    // This verifies whether the implementation has been done correctly, using simplified test
    // cases.

    // Define tolerance level for differences
    double tolerance =  std::numeric_limits< double >::epsilon( );

    // Define test cases.

    // Declare vectors for later loop-wise check.
    std::vector< Eigen::Vector3d > positionOfBodyForTestCase;
    std::vector< double > gravitationalParameterOfPerturbingBodyForTestCase;
    std::vector< Eigen::Vector3d > expectedAccelerationForTestCase;

    // Central body at origin assumed.
    // Perturbed body at 0 index.
    Eigen::Vector3d position( 4.0, 0.0, 0.0 ); // Arbitrary body.
    positionOfBodyForTestCase.push_back( position );

    gravitationalParameterOfPerturbingBodyForTestCase.push_back( 1.0 ); // Technically not a
                           // perturbing body, but for corresponding indices still added. Also,
                           // parameter of one to simplify manual computation.

    // Set expected acceleration.
    Eigen::Vector3d acceleration( 0.0, 0.0, 0.0 ); // Should be infinity, but
                                                   // shouldn't be necessary anyways.
    expectedAccelerationForTestCase.push_back( acceleration );

    // Inner planet perturbing body.
    position << 3.0, 0.0, 0.0; // Arbitrary body.
    positionOfBodyForTestCase.push_back( position );

    gravitationalParameterOfPerturbingBodyForTestCase.push_back( 1.0 ); // A parameter of one to
                                                                    // simplify manual computation.

    // Set expected acceleration.
    acceleration << -1.1111111111111111111111111111111, 0.0, 0.0; // Manually computed.
    expectedAccelerationForTestCase.push_back( acceleration );

    // Outer planet perturbing body.
    position << 5.0, 0.0, 0.0; // Arbitrary body
    positionOfBodyForTestCase.push_back( position );

    gravitationalParameterOfPerturbingBodyForTestCase.push_back( 1.0 ); // A parameter of one to
                                                                    // simplify manual computation.

    // Set expected acceleration.
    acceleration << 0.96, 0.0, 0.0; // Manually computed.
    expectedAccelerationForTestCase.push_back( acceleration );

    // Inner planet perturbing body (second position).
    position << 0.0, 3.0, 0.0; // Arbitrary body.
    positionOfBodyForTestCase.push_back( position );

    gravitationalParameterOfPerturbingBodyForTestCase.push_back( 1.0 ); // A parameter of one to
                                                                    // simplify manual computation.

    // Set expected acceleration.
    acceleration << -0.032, -0.08711111111111111111111111111111, 0.0; // Manually computed.
    expectedAccelerationForTestCase.push_back( acceleration );

    // Inner planet perturbing body (third position).
    position << -3.0, 0.0, 0.0; // Arbitrary body.
    positionOfBodyForTestCase.push_back( position );

    gravitationalParameterOfPerturbingBodyForTestCase.push_back( 1.0 ); // A parameter of one to
                                                                     // simplify manual computation.

    // Set expected acceleration.
    acceleration << 0.09070294784580498866213151927438, 0.0, 0.0; // Manually computed.
    expectedAccelerationForTestCase.push_back( acceleration );

    // Heavier outer planet perturbing body.
    position << 5.0, 0.0, 0.0; // Arbitrary body.
    positionOfBodyForTestCase.push_back( position );

    gravitationalParameterOfPerturbingBodyForTestCase.push_back( 2.0 ); // A parameter of two to
                                                                    // simplify manual computation.

    // Set expected acceleration.
    acceleration << 1.92, 0.0, 0.0; // Manually computed.
    expectedAccelerationForTestCase.push_back( acceleration );

    // Examine test cases.
    for ( unsigned int i = 1; i < positionOfBodyForTestCase.size( ); i++ )
    {

        // Create central gravity acceleration objects.
        gravitation::CentralGravitationalAccelerationModel3dPointer directAccelerationModel =
                boost::make_shared< gravitation::CentralGravitationalAccelerationModel3d >(
                    boost::lambda::constant( positionOfBodyForTestCase[ 0 ] ),
                    gravitationalParameterOfPerturbingBodyForTestCase[ i ],
                    boost::lambda::constant( positionOfBodyForTestCase[ i ] ) );

        gravitation::CentralGravitationalAccelerationModel3dPointer
                centralBodyAccelerationModel =
                boost::make_shared< gravitation::CentralGravitationalAccelerationModel3d >(
                    boost::lambda::constant( Eigen::Vector3d::Zero( ) ),
                    gravitationalParameterOfPerturbingBodyForTestCase[ i ],
                    boost::lambda::constant( positionOfBodyForTestCase[ i ] ) );

        // Create third body gravity acceleration objects.
        boost::shared_ptr<
                gravitation::ThirdBodyCentralGravityAcceleration > thirdBodyAcceleration =
                    boost::make_shared< gravitation::ThirdBodyCentralGravityAcceleration >(
                        directAccelerationModel, centralBodyAccelerationModel, "centralBody" );
        thirdBodyAcceleration->updateMembers( );

        // Compute perturbational acceleration for parameters given.
        acceleration = gravitation::computeThirdBodyPerturbingAcceleration(
                    gravitationalParameterOfPerturbingBodyForTestCase[ i ],
                    positionOfBodyForTestCase[ i ],
                    positionOfBodyForTestCase[ 0 ] );

        // Check values with reference data.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedAccelerationForTestCase[ i ], acceleration,
                                           tolerance );

        // Compute from manual sum.
        acceleration = gravitation::computeGravitationalAcceleration(
                    positionOfBodyForTestCase[ 0 ],
                gravitationalParameterOfPerturbingBodyForTestCase[ i ] ,
                positionOfBodyForTestCase[ i ] ) -
                gravitation::computeGravitationalAcceleration(
                    Eigen::Vector3d::Zero( ),
                    gravitationalParameterOfPerturbingBodyForTestCase[ i ],
                    positionOfBodyForTestCase[ i ] );

        // Check values with reference data.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedAccelerationForTestCase[ i ], acceleration,
                                           tolerance );

        acceleration = thirdBodyAcceleration->getAcceleration( );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedAccelerationForTestCase[ i ], acceleration,
                                           tolerance );
    }
}

//! Test if perturbational acceleration is physically correct.
BOOST_AUTO_TEST_CASE( testRealisticThirdBodyPerturbation )
{
    // This is mainly using more realistic orders of magnitude for values. Also, this verifies that
    // the routine handles correctly all three components.

    // Define tolerance.
    // This tolerance is set to higher than machine precision, because the unit tests otherwise
    // fail. This is most likely the case because the accelerations being tested are much smalled
    // than 1.0, which is what the values are normalized to when using a relative measure to
    // check if the computed values match the expected values.
    double tolerance = 1.0e-14;

    // Define the affected position (almost geostationary).
    Eigen::Vector3d realisticTestPosition( -40000000.0, 9000000.0, -9500000.0 );

    // Define the perturbing position (almost like the Moon).
    Eigen::Vector3d realisticPerturberPosition ( 25000000.0, -380000000.0, -55000000.0 );
    double realisticGravitationalParameterOfPerturbingBody = 4900.0e9;

    // Set the expected acceleration values (computed manually).
    Eigen::Vector3d expectedRealisticAcceleration ( 2.9394638802044552864534126789629e-6,
                                                    2.2253948636344714462843917227786e-6,
                                                    1.1680066249314399399805985890762e-6 );

    // Compute the acceleration.
    Eigen::Vector3d realisticComputedAcceleration =
            gravitation::computeThirdBodyPerturbingAcceleration(
                realisticGravitationalParameterOfPerturbingBody,
                realisticPerturberPosition,
                realisticTestPosition );

    // Compare differences.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedRealisticAcceleration, realisticComputedAcceleration,
                                       tolerance );

    // Compute from manual sum.
    realisticComputedAcceleration = gravitation::computeGravitationalAcceleration(
                realisticTestPosition,
                realisticGravitationalParameterOfPerturbingBody,
                realisticPerturberPosition ) -
            gravitation::computeGravitationalAcceleration(
                Eigen::Vector3d::Zero( ),
                realisticGravitationalParameterOfPerturbingBody,
                realisticPerturberPosition );

    // Compare differences.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedRealisticAcceleration,
                                       realisticComputedAcceleration,
                                       tolerance );

    // Redo test in barycentric frame to check whether position of central body is handled
    // correctly.
    Eigen::Vector3d barycentricEarthPosition = Eigen::Vector3d( 100E9, -110E9, 3E9 );
    realisticPerturberPosition += barycentricEarthPosition;
    realisticTestPosition += barycentricEarthPosition;

    // Compute the acceleration.
    realisticComputedAcceleration =
            gravitation::computeThirdBodyPerturbingAcceleration(
                realisticGravitationalParameterOfPerturbingBody,
                realisticPerturberPosition,
                realisticTestPosition,
                barycentricEarthPosition );

    // Compare differences.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedRealisticAcceleration,
                                       realisticComputedAcceleration,
                                       tolerance );

    // Compute from manual sum.
    realisticComputedAcceleration = gravitation::computeGravitationalAcceleration(
                realisticTestPosition,
                realisticGravitationalParameterOfPerturbingBody,
                realisticPerturberPosition ) -
            gravitation::computeGravitationalAcceleration(
                barycentricEarthPosition,
                realisticGravitationalParameterOfPerturbingBody,
                realisticPerturberPosition );

    // Compare differences.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedRealisticAcceleration,
                                       realisticComputedAcceleration,
                                       tolerance );

    // Create central gravity acceleration objects.
    gravitation::CentralGravitationalAccelerationModel3dPointer directAccelerationModel =
            boost::make_shared< gravitation::CentralGravitationalAccelerationModel3d >(
                boost::lambda::constant( realisticTestPosition ),
                realisticGravitationalParameterOfPerturbingBody,
                boost::lambda::constant( realisticPerturberPosition ) );

    gravitation::CentralGravitationalAccelerationModel3dPointer
            centralBodyAccelerationModel =
            boost::make_shared< gravitation::CentralGravitationalAccelerationModel3d >(
                boost::lambda::constant( barycentricEarthPosition ),
                realisticGravitationalParameterOfPerturbingBody,
                boost::lambda::constant( realisticPerturberPosition ) );

    // Create third body gravity acceleration objects.
    boost::shared_ptr<
            gravitation::ThirdBodyCentralGravityAcceleration > thirdBodyAcceleration =
                boost::make_shared< gravitation::ThirdBodyCentralGravityAcceleration >(
                    directAccelerationModel, centralBodyAccelerationModel, "centralBody" );
    thirdBodyAcceleration->updateMembers( );

    realisticComputedAcceleration = thirdBodyAcceleration->getAcceleration( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedRealisticAcceleration,
                                       realisticComputedAcceleration,
                                       tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
