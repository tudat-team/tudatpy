/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <stdexcept>
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/observation_models/corrections/laserTargetCorrection.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrection.h"

namespace tudat
{
namespace unit_tests
{
using namespace observation_models;

BOOST_AUTO_TEST_SUITE( test_laser_target_signature )

//! A positive published CoM correction shortens each computed surface leg and shifts the reflection epoch.
BOOST_AUTO_TEST_CASE( scalar_center_of_mass_and_coupled_light_time )
{

    // Use a fixed station and target, so each exact optical path is the nominal distance minus CoM.
    const double distance = 12000000.0;

    const double centerOfMass = 0.251;

    const auto signature = std::make_shared< LaserTargetSignature >( centerOfMass );

    const auto correction = std::make_shared< LaserTargetLightTimeCorrection >( signature, nullptr );

    const Eigen::Vector6d targetState = ( Eigen::Vector6d( ) << distance, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( );

    const auto station = std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ) );

    const auto target = std::make_shared< ephemerides::ConstantEphemeris >( targetState );

    const std::vector< std::shared_ptr< LightTimeCorrection > > corrections = { correction };

    const std::vector< std::shared_ptr< LightTimeCalculator<> > > legs = {
        std::make_shared< LightTimeCalculator<> >( station, target, corrections ),
        std::make_shared< LightTimeCalculator<> >( target, station, corrections )
    };

    FullLinkLightTimeCalculator<> calculator( legs );

    std::vector< double > times;

    std::vector< Eigen::Vector6d > states;

    const double roundTrip = calculator.calculateLightTimeWithLinkEndsStates( 0.0, transmitter, times, states );

    // The optical correction contributes once to each leg, including the uplink reflection time.
    BOOST_CHECK_CLOSE_FRACTION( roundTrip * physical_constants::SPEED_OF_LIGHT, 2.0 * ( distance - centerOfMass ), 2.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( times[ 1 ], ( distance - centerOfMass ) / physical_constants::SPEED_OF_LIGHT, 2.0E-15 );
    BOOST_CHECK( requiresMultiLegIterations( laser_target_signature ) );
    BOOST_CHECK( !calculator.doCorrectionsNeedFrequency( ) );

    const auto response = signature->calculateResponse( Eigen::Vector3d::UnitX( ), Eigen::Vector3d::UnitX( ), 0.0 );

    // A spherical scalar correction has no return-distribution width.
    BOOST_CHECK_EQUAL( response.getMeanRangeCorrection( ), -centerOfMass );
    BOOST_CHECK_EQUAL( response.rangeSpread, 0.0 );
}

//! A symmetric pair has zero centroid shift and finite width; common directional weights change the centroid.
BOOST_AUTO_TEST_CASE( weighted_return_centroid_and_width )
{

    // Equal returns at opposite optical points isolate the difference between a centroid and a width.
    const LaserTargetSignature symmetric( { Eigen::Vector3d( -0.5, 0.0, 0.0 ), Eigen::Vector3d( 0.5, 0.0, 0.0 ) }, { 1.0, 1.0 } );

    const auto pair = symmetric.calculateResponse( Eigen::Vector3d( 12.0, 0.0, 0.0 ), Eigen::Vector3d( 12.0, 0.0, 0.0 ), 0.0 );

    // Collinear ranges are exactly 12 +/- 0.5 m, so width cannot masquerade as a positive range bias.
    BOOST_CHECK_SMALL( pair.getMeanRangeCorrection( ), 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( pair.rangeSpread, 0.5, 2.0E-15 );

    const LaserTargetSignature directional( { Eigen::Vector3d( 3.0, 4.0, 0.0 ), Eigen::Vector3d( 0.0, 0.0, 1.0 ) },
                                            { 1.0, 1.0 },
                                            []( double, const Eigen::Vector3d& up, const Eigen::Vector3d& down ) {
                                                return std::vector< double >{ 1.0, up.dot( down ) < 0.0 ? 3.0 : 1.0 };
                                            } );

    const auto weighted = directional.calculateResponse( Eigen::Vector3d( 0.0, 0.0, 12.0 ), Eigen::Vector3d( 0.0, 0.0, -12.0 ), 0.0 );

    // The 5-12-13 triangle gives point changes (1,1) and (-1,1), with shared weights 1:3.
    BOOST_CHECK_CLOSE_FRACTION( weighted.meanPathCorrections.x( ), -0.5, 2.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( weighted.meanPathCorrections.y( ), 1.0, 2.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( weighted.getMeanRangeCorrection( ), 0.25, 2.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( weighted.rangeSpread, std::sqrt( 3.0 ) / 4.0, 2.0E-15 );
}

//! Small optical offsets survive astronomical ranges and are rotated relative to the solver's target endpoint.
BOOST_AUTO_TEST_CASE( stable_finite_paths_and_target_orientation )
{

    // A micrometre on-axis offset remains measurable when the nominal range is a trillion metres.
    const LaserTargetSignature microscopic( { Eigen::Vector3d( 1.0E-6, 0.0, 0.0 ) }, { 1.0 } );

    const auto distant = microscopic.calculateResponse( Eigen::Vector3d( 1.0E12, 0.0, 0.0 ), Eigen::Vector3d( 1.0E12, 0.0, 0.0 ), 0.0 );

    // Direct subtraction of these two norms would lose the micrometre signal completely.
    BOOST_CHECK_CLOSE_FRACTION( distant.getMeanRangeCorrection( ), -1.0E-6, 2.0E-15 );

    const auto signature = std::make_shared< LaserTargetSignature >( std::vector< Eigen::Vector3d >{ Eigen::Vector3d::UnitX( ) },
                                                                     std::vector< double >{ 1.0 } );

    LaserTargetLightTimeCorrection correction( signature, []( double ) {
        return Eigen::Quaterniond( Eigen::AngleAxisd( -0.5 * std::acos( -1.0 ), Eigen::Vector3d::UnitZ( ) ) );
    } );

    Eigen::Vector6d target = Eigen::Vector6d::Zero( );

    target.head< 3 >( ).setConstant( 1.0E11 );

    Eigen::Vector6d observer = target;

    observer.y( ) += 12.0;

    const std::vector< Eigen::Vector6d > states = { observer, target, target, observer };

    const std::vector< double > times = { 0.0, 1.0, 1.0, 2.0 };

    // Body-fixed +x points toward inertial +y; each retarded leg is shorter by one metre.
    BOOST_CHECK_CLOSE_FRACTION(
            correction.calculateLightTimeCorrectionWithMultiLegLinkEndStates( states, times, 0 ) * physical_constants::SPEED_OF_LIGHT,
            -1.0,
            2.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION(
            correction.calculateLightTimeCorrectionWithMultiLegLinkEndStates( states, times, 2 ) * physical_constants::SPEED_OF_LIGHT,
            -1.0,
            2.0E-15 );
}

//! Reject combining a scalar CoM correction with an array response or applying either to just one leg.
BOOST_AUTO_TEST_CASE( target_response_applied_once_per_leg )
{

    // Alternative target conventions must not enter the same physical light path together.
    const auto scalar = laserTargetLightTimeCorrection( "Target", std::make_shared< LaserTargetSignature >( 0.251 ) );

    const auto array = laserTargetLightTimeCorrection(
            "Target",
            std::make_shared< LaserTargetSignature >( std::vector< Eigen::Vector3d >{ Eigen::Vector3d::Zero( ) },
                                                      std::vector< double >{ 1.0 } ) );

    // One shared response on both legs is valid; mixed conventions would double count the physical target.
    BOOST_CHECK_NO_THROW( validateLaserTargetLightTimeCorrections( { { scalar }, { scalar } } ) );
    BOOST_CHECK_THROW( validateLaserTargetLightTimeCorrections( { { scalar, array }, { scalar, array } } ), std::invalid_argument );
    BOOST_CHECK_THROW( validateLaserTargetLightTimeCorrections( { { scalar }, {} } ), std::invalid_argument );
}

BOOST_AUTO_TEST_SUITE_END( )
}  // namespace unit_tests
}  // namespace tudat
