/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license.
 */

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>

#include "tudat/astro/gravitation/integratedGravityFieldVariations.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/legendrePolynomials.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createGravityFieldVariations.h"
#include "tudat/simulation/propagation_setup/setNumericallyIntegratedStates.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_integrated_gravity_field_variations )

BOOST_AUTO_TEST_CASE( testCurrentAndInterpolatedLifecycle )
{
    using namespace gravitation;

    IntegratedGravityFieldVariations variation;
    const Eigen::Vector5d first = ( Eigen::Vector5d( ) << 1.0e-5, 2.0e-5, 3.0e-5, 4.0e-5, 5.0e-5 ).finished( );
    const Eigen::Vector5d last = 3.0 * first;
    const Eigen::Vector5d current = -2.0 * first;
    const Eigen::Vector5d currentDerivative = 0.25 * first;

    variation.setCoefficientCorrectionHistory( { { 0.0, first }, { 10.0, last } } );
    variation.setCurrentCoefficientCorrections( current );
    variation.setCurrentCoefficientCorrectionDerivative( currentDerivative );
    variation.setIsBodyInPropagation( true );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( variation.getCoefficientCorrections( 5.0 ), current, 5.0e-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( variation.getCoefficientCorrectionDerivative( 5.0 ), currentDerivative, 5.0e-15 );

    const auto correctionsInPropagation = variation.calculateSphericalHarmonicsCorrections( 5.0 );
    BOOST_CHECK_CLOSE_FRACTION( correctionsInPropagation.first( 0, 0 ),
                                current( 0 ) / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 ),
                                5.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( correctionsInPropagation.second( 0, 2 ),
                                current( 4 ) / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 ),
                                5.0e-15 );

    variation.setIsBodyInPropagation( false );
    const Eigen::Vector5d interpolated = 0.5 * ( first + last );
    const Eigen::Vector5d interpolatedDerivative = ( last - first ) / 10.0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( variation.getCoefficientCorrections( 5.0 ), interpolated, 5.0e-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( variation.getCoefficientCorrectionDerivative( 5.0 ), interpolatedDerivative, 5.0e-15 );
}

BOOST_AUTO_TEST_CASE( testMultiBodyHistoryAndRigidBodyRefresh )
{
    using namespace gravitation;
    using namespace propagators;
    using namespace simulation_setup;

    SystemOfBodies bodies;
    bodies.createEmptyBody( "A" );
    bodies.createEmptyBody( "B" );

    Eigen::MatrixXd nominalCosine = Eigen::MatrixXd::Zero( 3, 3 );
    Eigen::MatrixXd nominalSine = Eigen::MatrixXd::Zero( 3, 3 );
    nominalCosine( 0, 0 ) = 1.0;
    nominalCosine( 2, 0 ) = -1.0e-3;
    nominalCosine( 2, 2 ) = 2.0e-4;

    bodies.at( "A" )->setGravityFieldModel(
            std::make_shared< SphericalHarmonicsGravityField >( 4.0e5, 2.0e3, nominalCosine, nominalSine, "AFixed", 0.4 ) );
    bodies.at( "B" )->setGravityFieldModel(
            std::make_shared< SphericalHarmonicsGravityField >( 5.0e5, 3.0e3, nominalCosine, nominalSine, "BFixed" ) );

    const std::shared_ptr< IntegratedGravityFieldVariations > variationA = ensureIntegratedGravityFieldVariation( bodies.at( "A" ), "A" );
    const std::shared_ptr< IntegratedGravityFieldVariations > variationB = ensureIntegratedGravityFieldVariation( bodies.at( "B" ), "B" );
    BOOST_CHECK( bodies.at( "A" )->getMassProperties( )->isInertiaTensorAvailable( ) );
    BOOST_CHECK( !bodies.at( "A" )->getMassProperties( )->isInertiaTensorDerivativeAvailable( ) );
    BOOST_CHECK( !bodies.at( "B" )->getMassProperties( )->isInertiaTensorAvailable( ) );

    const Eigen::Vector5d a0 = ( Eigen::Vector5d( ) << 1.0e-5, 2.0e-5, 3.0e-5, 4.0e-5, 5.0e-5 ).finished( );
    const Eigen::Vector5d a1 = 2.0 * a0;
    const Eigen::Vector5d b0 = -3.0 * a0;
    const Eigen::Vector5d b1 = 5.0 * a0;
    std::map< double, Eigen::VectorXd > history;
    history[ 0.0 ] = Eigen::VectorXd::Zero( 12 );
    history[ 10.0 ] = Eigen::VectorXd::Zero( 12 );
    history[ 0.0 ].segment( 1, 5 ) = a0;
    history[ 0.0 ].segment( 6, 5 ) = b0;
    history[ 10.0 ].segment( 1, 5 ) = a1;
    history[ 10.0 ].segment( 6, 5 ) = b1;

    resetIntegratedBodyGravity( bodies, history, { "A", "B" }, { 1, 10 }, nullptr );
    const std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > fieldA =
            std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >( bodies.at( "A" )->getGravityFieldModel( ) );
    const std::shared_ptr< TimeDependentSphericalHarmonicsGravityField > fieldB =
            std::dynamic_pointer_cast< TimeDependentSphericalHarmonicsGravityField >( bodies.at( "B" )->getGravityFieldModel( ) );
    BOOST_REQUIRE( fieldA != nullptr );
    BOOST_REQUIRE( fieldB != nullptr );
    BOOST_CHECK_CLOSE_FRACTION( fieldA->getCosineCoefficients( )( 2, 0 ),
                                nominalCosine( 2, 0 ) + a1( 0 ) / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 ),
                                5.0e-15 );
    BOOST_CHECK( bodies.at( "A" )->getMassProperties( )->isInertiaTensorDerivativeAvailable( ) );

    bodies.at( "A" )->setIsBodyInPropagation( false );
    bodies.at( "B" )->setIsBodyInPropagation( false );
    bodies.at( "A" )->updateCurrentGravityField( 5.0 );
    bodies.at( "B" )->updateCurrentGravityField( 5.0 );

    const Eigen::Vector5d expectedA = 0.5 * ( a0 + a1 );
    const Eigen::Vector5d expectedB = 0.5 * ( b0 + b1 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( variationA->getCoefficientCorrections( 5.0 ), expectedA, 5.0e-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( variationB->getCoefficientCorrections( 5.0 ), expectedB, 5.0e-15 );

    BOOST_CHECK_CLOSE_FRACTION(
            fieldA->getCosineCoefficients( )( 2, 0 ),
            nominalCosine( 2, 0 ) + expectedA( 0 ) / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 ),
            5.0e-15 );
    BOOST_CHECK_CLOSE_FRACTION( fieldB->getSineCoefficients( )( 2, 2 ),
                                expectedB( 4 ) / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 ),
                                5.0e-15 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( bodies.at( "A" )->getBodyInertiaTensor( ), fieldA->getInertiaTensor( ), 5.0e-15 );
    const Eigen::Vector5d expectedDerivativeA = ( a1 - a0 ) / 10.0;
    const Eigen::Matrix3d expectedInertiaDerivative = computeDerivativeInertiaTensor( expectedDerivativeA( 0 ),
                                                                                      expectedDerivativeA( 1 ),
                                                                                      expectedDerivativeA( 2 ),
                                                                                      expectedDerivativeA( 3 ),
                                                                                      expectedDerivativeA( 4 ),
                                                                                      bodies.at( "A" )->getBodyMass( ),
                                                                                      fieldA->getReferenceRadius( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( bodies.at( "A" )->getBodyInertiaTensorDerivative( ), expectedInertiaDerivative, 5.0e-15 );
    BOOST_CHECK_THROW( bodies.at( "B" )->getBodyInertiaTensor( ), std::runtime_error );

    // Switching back into propagation selects the live state instead of the installed history.
    bodies.at( "A" )->setIsBodyInPropagation( true );
    bodies.at( "A" )->setCurrentPropagatedGravityFieldVariation( -a0, 5.0 );
    const Eigen::Vector5d livePropagationState = -a0;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( variationA->getCoefficientCorrections( 5.0 ), livePropagationState, 5.0e-15 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
