/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <memory>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/basics/testMacros.h"
#include "tudat/io/trackingData.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/simulation/estimation_setup/createObservationCollection.h"

namespace tudat
{
namespace unit_tests
{

using namespace observation_models;
using namespace simulation_setup;

BOOST_AUTO_TEST_SUITE( test_astrometry_time_biases_and_uncertainties )

//! Compute right ascension and declination of a relative position
Eigen::Vector2d computeRightAscensionAndDeclination( const Eigen::Vector3d& relativePosition )
{
    return Eigen::Vector2d( std::atan2( relativePosition( 1 ), relativePosition( 0 ) ),
                            std::asin( relativePosition( 2 ) / relativePosition.norm( ) ) );
}

//! Check that observation time biases and uncertainties are applied to angular position observations and weights, using the rates
//! of right ascension and declination w.r.t. J2000 (for an ECLIPJ2000 global frame), for per-observation and set-level weights
BOOST_AUTO_TEST_CASE( testAngularPositionTimeBiasesAndUncertainties )
{
    Eigen::Vector6d observerState, observedBodyState;
    observerState << 1.4E11, 0.3E11, -0.1E11, -5.0E3, 2.9E4, 1.0E3;
    observedBodyState << 1.6E11, -0.2E11, 0.2E11, 1.0E4, 2.2E4, -3.0E3;

    BodyListSettings bodySettings( "SSB", "ECLIPJ2000" );
    bodySettings.addSettings( "Observer" );
    bodySettings.at( "Observer" )->ephemerisSettings = constantEphemerisSettings( observerState, "SSB", "ECLIPJ2000" );
    bodySettings.addSettings( "Asteroid" );
    bodySettings.at( "Asteroid" )->ephemerisSettings = constantEphemerisSettings( observedBodyState, "SSB", "ECLIPJ2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Rates w.r.t. J2000 (constant, since the states of both bodies are constant), from central differences for linear relative motion
    const Eigen::Matrix3d rotationToJ2000 = reference_frames::getJ2000toECLIPJ2000TransformationMatrix( ).transpose( );
    const Eigen::Vector3d relativePosition = rotationToJ2000 * ( observedBodyState - observerState ).segment( 0, 3 );
    const Eigen::Vector3d relativeVelocity = rotationToJ2000 * ( observedBodyState - observerState ).segment( 3, 3 );
    const double timeStep = 10.0;
    const Eigen::Vector2d angularRates = ( computeRightAscensionAndDeclination( relativePosition + timeStep * relativeVelocity ) -
                                           computeRightAscensionAndDeclination( relativePosition - timeStep * relativeVelocity ) ) /
            ( 2.0 * timeStep );

    const std::vector< double > timeUncertainties = { 2.0, 0.5 };
    const std::vector< double > timeBiases = { 3.0, -1.5 };

    // Positional covariance of the two observations, uncorrelated (per-observation weights) or correlated (set-level weights)
    const double positionalVariance = 1.0E-12;
    const Eigen::Matrix4d uncorrelatedCovariance = positionalVariance * Eigen::Matrix4d::Identity( );
    Eigen::Matrix4d correlatedCovariance = uncorrelatedCovariance;
    correlatedCovariance.block( 0, 2, 2, 2 ) = 0.3 * positionalVariance * Eigen::Matrix2d::Identity( );
    correlatedCovariance.block( 2, 0, 2, 2 ) = 0.3 * positionalVariance * Eigen::Matrix2d::Identity( );

    for( const bool useSetLevelWeights : { false, true } )
    {
        const data::PlainLinkDefinition linkEnds = { { { "Asteroid", "" }, "transmitter" }, { { "Observer", "" }, "receiver" } };
        std::shared_ptr< data::TrackingData< double, double > > trackingData =
                std::make_shared< data::TrackingData< double, double > >( "AngularPosition",
                                                                          linkEnds,
                                                                          std::vector< Eigen::VectorXd >( 2, Eigen::Vector2d::Zero( ) ),
                                                                          std::vector< double >{ 0.0, 100.0 },
                                                                          "receiver" );
        trackingData->setObservationWeightSettings(
                useSetLevelWeights ? ObservationWeightSettings::setBlock( correlatedCovariance.inverse( ) )
                                   : ObservationWeightSettings::diagonalPerObservation(
                                             std::vector< Eigen::VectorXd >( 2, Eigen::Vector2d::Constant( 1.0 / positionalVariance ) ) ) );
        trackingData->setObservationTimeUncertainties( timeUncertainties );
        trackingData->setObservationTimeBiases( timeBiases );

        const ObservationVectorData< double, double > vectorData =
                createObservationDatasetFromTrackingData< double, double >( { trackingData }, bodies )->createObservationVectorData( );

        Eigen::Matrix4d expectedCovariance = useSetLevelWeights ? correlatedCovariance : uncorrelatedCovariance;
        for( unsigned int i = 0; i < 2; i++ )
        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    vectorData.getObservationVector( ).segment( 2 * i, 2 ), Eigen::Vector2d( -timeBiases.at( i ) * angularRates ), 1.0E-8 );
            expectedCovariance.block( 2 * i, 2 * i, 2, 2 ) +=
                    timeUncertainties.at( i ) * timeUncertainties.at( i ) * angularRates * angularRates.transpose( );
        }
        const Eigen::MatrixXd computedCovariance = vectorData.getSparseWeightMatrix( ).toDense( ).inverse( );
        BOOST_CHECK_SMALL( ( computedCovariance - expectedCovariance ).norm( ) / expectedCovariance.norm( ), 1.0E-8 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
