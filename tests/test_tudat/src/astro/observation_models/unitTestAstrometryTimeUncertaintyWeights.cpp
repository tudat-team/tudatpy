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

#include "tudat/io/trackingData.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/simulation/estimation_setup/trackingDataWeighting.h"

namespace tudat
{
namespace unit_tests
{

using namespace observation_models;
using namespace simulation_setup;

BOOST_AUTO_TEST_SUITE( test_astrometry_time_uncertainty_weights )

//! Compute right ascension and declination of a relative position
Eigen::Vector2d computeRightAscensionAndDeclination( const Eigen::Vector3d& relativePosition )
{
    return Eigen::Vector2d( std::atan2( relativePosition( 1 ), relativePosition( 0 ) ),
                            std::asin( relativePosition( 2 ) / relativePosition.norm( ) ) );
}

//! Check analytical right ascension and declination rates against central differences for linear relative motion
BOOST_AUTO_TEST_CASE( testRightAscensionAndDeclinationRates )
{
    Eigen::Vector6d relativeState;
    relativeState << 1.2E11, -0.7E11, 0.4E11, 1.5E4, 2.0E4, -0.8E4;

    const double timeStep = 10.0;
    const Eigen::Vector2d numericalRates =
            ( computeRightAscensionAndDeclination( relativeState.segment( 0, 3 ) + timeStep * relativeState.segment( 3, 3 ) ) -
              computeRightAscensionAndDeclination( relativeState.segment( 0, 3 ) - timeStep * relativeState.segment( 3, 3 ) ) ) /
            ( 2.0 * timeStep );
    const Eigen::Vector2d analyticalRates = computeRightAscensionAndDeclinationRates( relativeState );

    for( unsigned int i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( analyticalRates( i ), numericalRates( i ), 1.0E-6 );
    }
}

//! Check that the observation epoch uncertainty is added to the covariance of the angular position observations
BOOST_AUTO_TEST_CASE( testAngularPositionWeightsWithTimeUncertainty )
{
    Eigen::Vector6d observerState, observedBodyState;
    observerState << 1.4E11, 0.3E11, -0.1E11, -5.0E3, 2.9E4, 1.0E3;
    observedBodyState << 1.6E11, -0.2E11, 0.2E11, 1.0E4, 2.2E4, -3.0E3;

    BodyListSettings bodySettings( "SSB", "J2000" );
    bodySettings.addSettings( "Observer" );
    bodySettings.at( "Observer" )->ephemerisSettings = constantEphemerisSettings( observerState, "SSB", "J2000" );
    bodySettings.addSettings( "Asteroid" );
    bodySettings.at( "Asteroid" )->ephemerisSettings = constantEphemerisSettings( observedBodyState, "SSB", "J2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    LinkEnds linkEnds;
    linkEnds[ receiver ] = LinkEndId( "Observer", "" );
    linkEnds[ transmitter ] = LinkEndId( "Asteroid", "" );

    const std::vector< double > epochs = { 0.0, 100.0 };
    std::vector< Eigen::VectorXd > observations( 2, Eigen::Vector2d::Zero( ) );
    std::shared_ptr< data::TrackingData< double, double > > trackingData = std::make_shared< data::TrackingData< double, double > >(
            "angular_position", data::PlainLinkDefinition( ), observations, epochs, "receiver" );

    const double positionalUncertainty = 1.0E-6;
    trackingData->setObservationWeightSettings( ObservationWeightSettings::diagonalPerObservation(
            std::vector< Eigen::VectorXd >( 2, Eigen::Vector2d::Constant( 1.0 / ( positionalUncertainty * positionalUncertainty ) ) ) ) );
    trackingData->setAstrometryTimeUncertainties( { 2.0, 0.5 } );

    const ObservationWeightSettings weightSettings =
            getAngularPositionWeightsWithTimeUncertainty< double, double >( trackingData, bodies, linkEnds, epochs );
    BOOST_REQUIRE( weightSettings.type_ == ObservationWeightSettings::WeightsBlockType::block_per_observation );

    const Eigen::Vector2d angularRates = computeRightAscensionAndDeclinationRates( observedBodyState - observerState );
    for( unsigned int i = 0; i < 2; i++ )
    {
        const double timeUncertainty = trackingData->getAstrometryTimeUncertainties( ).at( i );
        const Eigen::Matrix2d expectedCovariance = positionalUncertainty * positionalUncertainty * Eigen::Matrix2d::Identity( ) +
                timeUncertainty * timeUncertainty * angularRates * angularRates.transpose( );
        const Eigen::Matrix2d computedCovariance = weightSettings.weightBlocks_.at( i ).inverse( );
        BOOST_CHECK_SMALL( ( computedCovariance - expectedCovariance ).norm( ) / expectedCovariance.norm( ), 1.0E-12 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
