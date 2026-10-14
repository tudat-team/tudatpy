#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>

#include "tudat/astro/basic_astro/oblateSpheroidBodyShapeModel.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/astro/ground_stations/earthDeformationTime.h"
#include "tudat/astro/ground_stations/oceanTideEarthDeformation.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "tudat/simulation/environment_setup/createBodyDeformationModel.h"

namespace tudat
{
namespace unit_tests
{
namespace
{
using basic_astrodynamics::BodyDeformationModel;
using ground_stations::GroundStationState;

struct FunctionDeformation : BodyDeformationModel {
    explicit FunctionDeformation( const std::function< Eigen::Vector3d( double ) >& function ): function( function ) {}

    Eigen::Vector3d calculateDisplacement( const double time, const Eigen::Vector3d& ) override
    {
        lastTime = time;
        return function( time );
    }

    Eigen::Vector3d calculateDisplacement( const double time, const std::shared_ptr< GroundStationState > ) override
    {
        return calculateDisplacement( time, Eigen::Vector3d::Zero( ) );
    }

    std::function< Eigen::Vector3d( double ) > function;
    double lastTime = std::numeric_limits< double >::quiet_NaN( );
};

struct AnalyticStateDeformation : FunctionDeformation {
    AnalyticStateDeformation( ): FunctionDeformation( []( double ) { return Eigen::Vector3d::Zero( ).eval( ); } ) {}

    Eigen::Vector6d calculateDisplacementState( const double, const std::shared_ptr< GroundStationState > ) override
    {
        Eigen::Vector6d result;
        result << 0.1, 0.2, 0.3, 0.01, 0.02, 0.03;
        return result;
    }
};

struct CoverageError : std::runtime_error {
    CoverageError( ): std::runtime_error( "Unavailable derivative sample" ) {}
};

std::shared_ptr< basic_astrodynamics::Iers2010EarthDeformation > makeIersModel(
        const std::function< Eigen::Vector6d( const double ) > customArguments = {} )
{
    Eigen::Vector6d moon = Eigen::Vector6d::Zero( ), sun = Eigen::Vector6d::Zero( );
    moon.head< 3 >( ) << -179996231.920342, -312468450.131567, -169288918.592160;
    sun.head< 3 >( ) << 137859926952.015, 54228127881.435, 23509422341.696;
    auto rotation = std::make_shared< ephemerides::SimpleRotationalEphemeris >( Eigen::Quaterniond::Identity( ), 7.292115e-5, 0.0 );
    return simulation_setup::createDefaultEarthIers2010DeformationModel(
            std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ) ),
            std::make_shared< ephemerides::ConstantEphemeris >( moon ),
            std::make_shared< ephemerides::ConstantEphemeris >( sun ),
            rotation,
            []( ) { return 1.0; },
            []( ) { return 0.0123000371; },
            []( ) { return 332946.0482; },
            customArguments );
}

input_output::OceanLoadingData makeOceanData( )
{
    input_output::OceanLoadingData data;
    input_output::OceanLoadingStation station;
    station.stationIdentifier = "VELOCITY_TEST";
    station.amplitudesMetres.setConstant( 0.001 );
    data.stations.emplace( station.stationIdentifier, station );
    return data;
}

std::shared_ptr< GroundStationState > stationWithModels( std::vector< std::shared_ptr< BodyDeformationModel > >& models,
                                                         const bool throwWhenUnavailable = true )
{
    auto motion = std::make_shared< ground_stations::BodyDeformationStationMotionModel >(
            [ &models ]( ) -> std::vector< std::shared_ptr< BodyDeformationModel > >& { return models; }, throwWhenUnavailable );
    return std::make_shared< GroundStationState >(
            Eigen::Vector3d( 4075578.385, 931852.890, 4801570.154 ), coordinate_conversions::cartesian_position, nullptr, motion );
}
}  // namespace

BOOST_AUTO_TEST_CASE( analytic_constant_linear_and_sinusoidal_rates )
{
    const Eigen::Vector3d position( 6.4e6, 1.0e6, 2.0e6 );
    FunctionDeformation constant( []( double ) { return Eigen::Vector3d( 0.1, -0.2, 0.3 ); } );
    BOOST_CHECK_SMALL( constant.calculateDisplacementVelocity( 100.0, position ).norm( ), 1.0e-15 );
    const Eigen::Vector3d slope( 1.0e-4, -3.0e-4, 7.0e-4 );
    FunctionDeformation linear( [ = ]( double t ) { return ( Eigen::Vector3d( 0.1, 0.2, 0.3 ) + t * slope ).eval( ); } );
    BOOST_CHECK_SMALL( ( linear.calculateDisplacementVelocity( 100.0, position ) - slope ).norm( ), 1.0e-16 );
    FunctionDeformation sinusoid( []( double t ) {
        return Eigen::Vector3d( 0.3 * std::sin( 0.01 * t ), 0.1 * std::cos( 0.02 * t ), 0.2 * std::sin( 0.03 * t ) );
    } );
    const double epoch = 43.0;
    const Eigen::Vector3d expected( 0.003 * std::cos( 0.01 * epoch ), -0.002 * std::sin( 0.02 * epoch ), 0.006 * std::cos( 0.03 * epoch ) );
    std::array< double, 3 > errors;
    for( std::size_t i = 0; i < errors.size( ); ++i )
    {
        sinusoid.setDisplacementDerivativeTimeStep( 4.0 / std::pow( 2.0, i ) );
        errors[ i ] = ( sinusoid.calculateDisplacementVelocity( epoch, position ) - expected ).norm( );
        BOOST_CHECK_EQUAL( sinusoid.lastTime, epoch );
    }
    BOOST_CHECK_GT( errors[ 0 ] / errors[ 1 ], 15.0 );
    BOOST_CHECK_GT( errors[ 1 ] / errors[ 2 ], 15.0 );
    BOOST_CHECK_SMALL( errors[ 2 ], 2.0e-10 );
    std::cout << "Analytic sinusoid errors h=4,2,1 s: " << errors[ 0 ] << ' ' << errors[ 1 ] << ' ' << errors[ 2 ] << '\n';
}

BOOST_AUTO_TEST_CASE( derivative_step_and_nonfinite_validation )
{
    FunctionDeformation model( []( double t ) { return Eigen::Vector3d( t, t, t ); } );
    for( double step : { 0.0, -1.0, std::numeric_limits< double >::infinity( ), std::numeric_limits< double >::quiet_NaN( ) } )
    {
        BOOST_CHECK_THROW( model.setDisplacementDerivativeTimeStep( step ), std::invalid_argument );
    }
    BOOST_CHECK_EQUAL( model.getDisplacementDerivativeTimeStep( ), 1.0 );
    BOOST_CHECK_THROW( model.calculateDisplacementVelocity( 1.0e20, Eigen::Vector3d::Zero( ) ), std::invalid_argument );
    BOOST_CHECK_THROW( model.calculateDisplacementVelocity( std::numeric_limits< double >::infinity( ), Eigen::Vector3d::Zero( ) ),
                       std::invalid_argument );
    FunctionDeformation invalid(
            []( double ) { return Eigen::Vector3d::Constant( std::numeric_limits< double >::quiet_NaN( ) ).eval( ); } );
    BOOST_CHECK_THROW( invalid.calculateDisplacementVelocity( 0.0, Eigen::Vector3d::Zero( ) ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( nonuniform_epoch_stencil_is_rejected_near_binade_boundary )
{
    const double boundaryEpoch = std::ldexp( 1.0, 53 ) - 1.0;
    const Eigen::Vector3d slope( 0.001, -0.002, 0.003 );
    FunctionDeformation model( [ = ]( double t ) { return ( ( t - boundaryEpoch ) * slope ).eval( ); } );
    model.setDisplacementDerivativeTimeStep( 3.0 );
    // All four points are distinct, but the upper outer point is only 5 seconds
    // from the centre instead of 6. Uniform weights would bias this linear rate
    // by 1/36. Reject it before evaluating any displacement samples.
    BOOST_CHECK_THROW( model.calculateDisplacementVelocity( boundaryEpoch, Eigen::Vector3d::Zero( ) ), std::invalid_argument );
    BOOST_CHECK( std::isnan( model.lastTime ) );
    model.setDisplacementDerivativeTimeStep( 4.0 );
    BOOST_CHECK_SMALL( ( model.calculateDisplacementVelocity( boundaryEpoch + 1.0, Eigen::Vector3d::Zero( ) ) - slope ).norm( ), 1.0e-17 );

    const double contemporaryEpoch = 334368000.125;
    FunctionDeformation contemporary( [ = ]( double t ) { return ( ( t - contemporaryEpoch ) * slope ).eval( ); } );
    BOOST_CHECK_SMALL( ( contemporary.calculateDisplacementVelocity( contemporaryEpoch, Eigen::Vector3d::Zero( ) ) - slope ).norm( ),
                       1.0e-17 );
}

BOOST_AUTO_TEST_CASE( errors_restore_central_cache_and_preserve_exception_type )
{
    FunctionDeformation bounded( []( double t ) {
        if( std::fabs( t ) > 0.5 )
        {
            throw CoverageError( );
        }
        return Eigen::Vector3d( t, t, t );
    } );
    BOOST_CHECK_THROW( bounded.calculateDisplacementVelocity( 0.0, Eigen::Vector3d::Zero( ) ), CoverageError );
    BOOST_CHECK_EQUAL( bounded.lastTime, 0.0 );

    bool failedSample = false;
    bool attemptedRestore = false;
    FunctionDeformation failedRestore( [ & ]( double t ) {
        if( t != 0.0 )
        {
            failedSample = true;
            throw CoverageError( );
        }
        if( failedSample )
        {
            attemptedRestore = true;
            throw std::logic_error( "Restore failure" );
        }
        return Eigen::Vector3d::Zero( ).eval( );
    } );
    BOOST_CHECK_THROW( failedRestore.calculateDisplacementVelocity( 0.0, Eigen::Vector3d::Zero( ) ), CoverageError );
    BOOST_CHECK( attemptedRestore );
}

BOOST_AUTO_TEST_CASE( optional_missing_model_does_not_hide_rate_failure )
{
    auto unavailable = std::make_shared< FunctionDeformation >( []( double ) -> Eigen::Vector3d { throw CoverageError( ); } );
    std::vector< std::shared_ptr< BodyDeformationModel > > models = { unavailable };
    auto station = stationWithModels( models, false );
    BOOST_CHECK_SMALL( station->getCartesianStateInTime( 0.0, "" ).tail< 3 >( ).norm( ), 1.0e-15 );
    models[ 0 ] = std::make_shared< FunctionDeformation >( []( double t ) {
        if( t != 0.0 )
        {
            throw CoverageError( );
        }
        return Eigen::Vector3d::Constant( 0.1 ).eval( );
    } );
    BOOST_CHECK_THROW( station->getCartesianStateInTime( 0.0, "" ), CoverageError );
    models[ 0 ] = unavailable;
    station = stationWithModels( models, true );
    BOOST_CHECK_THROW( station->getCartesianStateInTime( 0.0, "" ), CoverageError );
}

BOOST_AUTO_TEST_CASE( default_station_pipeline_and_state_override )
{
    const Eigen::Vector3d slope( 1.0e-4, 2.0e-4, -3.0e-4 );
    auto linear = std::make_shared< FunctionDeformation >( [ = ]( double t ) { return ( t * slope ).eval( ); } );
    std::vector< std::shared_ptr< BodyDeformationModel > > models = { linear };
    auto station = stationWithModels( models );
    const auto state = station->getCartesianStateInTime( 12.0, "" );
    BOOST_CHECK_SMALL( ( state.tail< 3 >( ) - slope ).norm( ), 1.0e-17 );
    BOOST_CHECK_SMALL( ( state.head< 3 >( ) - station->getNominalCartesianPosition( ) - 12.0 * slope ).norm( ), 1.0e-9 );
    struct AnalyticRateDeformation : FunctionDeformation {
        AnalyticRateDeformation( ): FunctionDeformation( []( double ) { return Eigen::Vector3d::Zero( ).eval( ); } ) {}
        Eigen::Vector3d calculateDisplacementVelocity( const double, const std::shared_ptr< GroundStationState > ) override
        {
            return Eigen::Vector3d( 0.04, 0.05, 0.06 );
        }
    };
    models[ 0 ] = std::make_shared< AnalyticRateDeformation >( );
    BOOST_CHECK_SMALL( ( station->getCartesianStateInTime( 12.0, "" ).tail< 3 >( ) - Eigen::Vector3d( 0.04, 0.05, 0.06 ) ).norm( ),
                       1.0e-17 );
    models[ 0 ] = std::make_shared< AnalyticStateDeformation >( );
    const auto overridden = station->getCartesianStateInTime( 12.0, "" );
    BOOST_CHECK_SMALL( ( overridden.tail< 3 >( ) - Eigen::Vector3d( 0.01, 0.02, 0.03 ) ).norm( ), 1.0e-17 );
}

BOOST_AUTO_TEST_CASE( oblate_nominal_station_basis_is_preserved )
{
    struct LocalDeformation : BodyDeformationModel {
        Eigen::Vector3d calculateDisplacement( const double, const Eigen::Vector3d& ) override
        {
            throw std::logic_error( "Full nominal station state required" );
        }
        Eigen::Vector3d calculateDisplacement( const double t, const std::shared_ptr< GroundStationState > station ) override
        {
            const auto enu = station->getEnuGeocentricUnitVectors( );
            return ( enu[ 0 ] * ( 0.1 + 0.001 * t ) + enu[ 1 ] * ( -0.2 + 0.002 * t ) + enu[ 2 ] * 0.3 ).eval( );
        }
    };
    auto shape = std::make_shared< basic_astrodynamics::OblateSpheroidBodyShapeModel >( 6378137.0, 1.0 / 298.257223563 );
    auto station =
            std::make_shared< GroundStationState >( Eigen::Vector3d( 100.0, 0.8, -1.2 ), coordinate_conversions::geodetic_position, shape );
    const auto xyz = station->getNominalCartesianPosition( );
    const double phi = std::atan2( xyz.z( ), std::hypot( xyz.x( ), xyz.y( ) ) );
    const double lambda = std::atan2( xyz.y( ), xyz.x( ) );
    const Eigen::Vector3d east( -std::sin( lambda ), std::cos( lambda ), 0.0 );
    const Eigen::Vector3d north( -std::sin( phi ) * std::cos( lambda ), -std::sin( phi ) * std::sin( lambda ), std::cos( phi ) );
    BOOST_CHECK_GT( std::fabs( phi - 0.8 ), 0.001 );
    LocalDeformation model;
    BOOST_CHECK_SMALL( ( model.calculateDisplacementVelocity( 123.0, station ) - 0.001 * east - 0.002 * north ).norm( ), 1.0e-16 );
}

BOOST_AUTO_TEST_CASE( production_iers_velocity_step_convergence )
{
    Eigen::Vector6d moon = Eigen::Vector6d::Zero( ), sun = Eigen::Vector6d::Zero( );
    moon.head< 3 >( ) << -179996231.920342, -312468450.131567, -169288918.592160;
    sun.head< 3 >( ) << 137859926952.015, 54228127881.435, 23509422341.696;
    auto rotation = std::make_shared< ephemerides::SimpleRotationalEphemeris >( Eigen::Quaterniond::Identity( ), 7.292115e-5, 0.0 );
    auto model = simulation_setup::createDefaultEarthIers2010DeformationModel(
            std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ) ),
            std::make_shared< ephemerides::ConstantEphemeris >( moon ),
            std::make_shared< ephemerides::ConstantEphemeris >( sun ),
            rotation,
            []( ) { return 1.0; },
            []( ) { return 0.0123000371; },
            []( ) { return 332946.0482; } );
    auto station = std::make_shared< GroundStationState >( Eigen::Vector3d( 4075578.385, 931852.890, 4801570.154 ) );
    const double time = 334324800.0 + 43200.0;
    const Eigen::Vector3d reference = model->calculateDisplacementVelocity( time, station );
    std::array< double, 3 > errors;
    for( std::size_t i = 0; i < errors.size( ); ++i )
    {
        model->setDisplacementDerivativeTimeStep( 1800.0 / std::pow( 2.0, i ) );
        errors[ i ] = ( model->calculateDisplacementVelocity( time, station ) - reference ).norm( );
    }
    BOOST_CHECK_GT( reference.norm( ), 1.0e-7 );
    BOOST_CHECK_GT( errors[ 0 ] / errors[ 1 ], 10.0 );
    BOOST_CHECK_GT( errors[ 1 ] / errors[ 2 ], 10.0 );
    model->setDisplacementDerivativeTimeStep( 0.5 );
    const double fineDifference = ( model->calculateDisplacementVelocity( time, station ) - reference ).norm( );
    BOOST_CHECK_SMALL( fineDifference, 2.0e-10 );
    std::cout << "Production IERS velocity " << reference.transpose( ) << " m/s; errors h=1800,900,450 s: " << errors[ 0 ] << ' '
              << errors[ 1 ] << ' ' << errors[ 2 ] << "; h=1 vs0.5 s: " << fineDifference << '\n';
}

BOOST_AUTO_TEST_CASE( earth_utc_interval_policy_includes_boundaries_and_allows_historical_drift )
{
    using basic_astrodynamics::validateEarthDeformationUtcInterval;
    for( const auto& date : { std::array< int, 3 >{ 2012, 7, 1 }, std::array< int, 3 >{ 2017, 1, 1 } } )
    {
        const double utc = basic_astrodynamics::timeFromDecomposedDateTime< double >( date[ 0 ], date[ 1 ], date[ 2 ], 0, 0, 0.0 );
        const double leap = sofa_interface::convertUTCtoTT( utc );
        // Ambiguity strictly inside the inserted second, even without straddling its end.
        BOOST_CHECK_THROW( validateEarthDeformationUtcInterval( leap - 0.8, leap - 0.2 ), std::domain_error );
        // No ambiguous five-point sample here: the interior midnight jump is detected separately.
        BOOST_CHECK_THROW( validateEarthDeformationUtcInterval( leap - 16.0, leap + 24.0 ), std::domain_error );
        BOOST_CHECK_THROW( validateEarthDeformationUtcInterval( leap, leap + 40.0 ), std::domain_error );
        BOOST_CHECK_THROW( validateEarthDeformationUtcInterval( leap - 40.0, leap ), std::domain_error );
    }
    for( const auto& date : { std::array< int, 3 >{ 2012, 7, 2 }, std::array< int, 3 >{ 1965, 6, 16 } } )
    {
        const double utc = basic_astrodynamics::timeFromDecomposedDateTime< double >( date[ 0 ], date[ 1 ], date[ 2 ], 0, 0, 0.0 );
        const double midnight = sofa_interface::convertUTCtoTT( utc );
        BOOST_CHECK_NO_THROW( validateEarthDeformationUtcInterval( midnight - 2.0, midnight + 2.0 ) );
    }
    const double historicalJump =
            sofa_interface::convertUTCtoTT( basic_astrodynamics::timeFromDecomposedDateTime< double >( 1965, 7, 1, 0, 0, 0.0 ) );
    BOOST_CHECK_THROW( validateEarthDeformationUtcInterval( historicalJump - 16.0, historicalJump + 24.0 ), std::domain_error );
    BOOST_CHECK_THROW( validateEarthDeformationUtcInterval( std::numeric_limits< double >::quiet_NaN( ), 0.0 ), std::invalid_argument );
    BOOST_CHECK_THROW( validateEarthDeformationUtcInterval( -std::numeric_limits< double >::max( ), std::numeric_limits< double >::max( ) ),
                       std::invalid_argument );
    BOOST_CHECK_THROW( validateEarthDeformationUtcInterval( -1.0e20, -1.0e20 ), std::exception );
}

BOOST_AUTO_TEST_CASE( earth_model_hooks_reject_utc_limitations_and_preserve_custom_callbacks )
{
    auto station = std::make_shared< GroundStationState >( Eigen::Vector3d( 4075578.385, 931852.890, 4801570.154 ) );
    station->setSiteId( "VELOCITY_TEST" );
    auto iers = makeIersModel( );
    basic_astrodynamics::OceanTideEarthDeformation ocean( makeOceanData( ) );
    for( const auto& date : { std::array< int, 3 >{ 2012, 7, 1 }, std::array< int, 3 >{ 2017, 1, 1 } } )
    {
        const double utc = basic_astrodynamics::timeFromDecomposedDateTime< double >( date[ 0 ], date[ 1 ], date[ 2 ], 0, 0, 0.0 );
        const double leap = sofa_interface::convertUTCtoTT( utc );
        BOOST_CHECK_THROW( iers->calculateDisplacementVelocity( leap, station ), std::domain_error );
        BOOST_CHECK_THROW( iers->calculateDisplacementVelocity( leap, station->getNominalCartesianPosition( ) ), std::domain_error );
        BOOST_CHECK_THROW( ocean.calculateDisplacementVelocity( leap, station ), std::domain_error );
        const auto customArguments = []( const double ) { return Eigen::Vector6d::Zero( ).eval( ); };
        auto customIers = makeIersModel( customArguments );
        basic_astrodynamics::OceanTideEarthDeformation customOcean( makeOceanData( ), customArguments );
        BOOST_CHECK( customIers->calculateDisplacementVelocity( leap, station ).allFinite( ) );
        BOOST_CHECK( customOcean.calculateDisplacementVelocity( leap, station ).allFinite( ) );
    }
    for( const auto& date : { std::array< int, 3 >{ 2012, 7, 2 }, std::array< int, 3 >{ 1965, 6, 16 } } )
    {
        const double utc = basic_astrodynamics::timeFromDecomposedDateTime< double >( date[ 0 ], date[ 1 ], date[ 2 ], 0, 0, 0.0 );
        const double midnight = sofa_interface::convertUTCtoTT( utc );
        BOOST_CHECK( iers->calculateDisplacementVelocity( midnight, station ).allFinite( ) );
        BOOST_CHECK( ocean.calculateDisplacementVelocity( midnight, station ).allFinite( ) );
    }
}

BOOST_AUTO_TEST_CASE( deformed_station_inertial_velocity_closure )
{
    const Eigen::Vector3d slope( 1.0e-3, -2.0e-3, 3.0e-3 );
    auto model = std::make_shared< FunctionDeformation >(
            [ = ]( double t ) { return ( Eigen::Vector3d( 0.1, 0.2, -0.3 ) + t * slope ).eval( ); } );
    std::vector< std::shared_ptr< BodyDeformationModel > > models = { model };
    auto station = stationWithModels( models );
    const double time = 123.0, rate = 7.292115e-5;
    auto rotation = std::make_shared< ephemerides::SimpleRotationalEphemeris >( Eigen::Quaterniond::Identity( ), rate, 0.0 );
    const Eigen::Vector6d local = station->getCartesianStateInTime( time, "" );
    const Eigen::Vector6d inertial = ephemerides::transformStateToInertialOrientation< double, double >( local, time, rotation );
    const Eigen::Matrix3d r = Eigen::AngleAxisd( rate * time, Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( );
    const Eigen::Vector3d expectedPosition = r * local.head< 3 >( );
    const Eigen::Vector3d expectedVelocity = r * slope + Eigen::Vector3d( 0.0, 0.0, rate ).cross( expectedPosition );
    BOOST_CHECK_SMALL( ( inertial.head< 3 >( ) - expectedPosition ).norm( ), 2.0e-9 );
    BOOST_CHECK_SMALL( ( inertial.tail< 3 >( ) - expectedVelocity ).norm( ), 1.0e-10 );
}

}  // namespace unit_tests
}  // namespace tudat

BOOST_AUTO_TEST_CASE( position_access_preserves_covered_endpoints_and_availability_policy )
{
    using namespace tudat;
    using namespace tudat::unit_tests;
    const Eigen::Vector3d nominal( 6378137.0, 0.0, 0.0 ), slope( 0.001, -0.002, 0.003 );
    auto bounded = std::make_shared< FunctionDeformation >( [ = ]( double t ) {
        if( t < 0.0 || t > 10.0 )
        {
            throw CoverageError( );
        }
        return ( Eigen::Vector3d( 0.1, 0.2, 0.3 ) + t * slope ).eval( );
    } );
    std::vector< std::shared_ptr< BodyDeformationModel > > models{ bounded };
    auto deformation = std::make_shared< ground_stations::BodyDeformationStationMotionModel >(
            [ &models ]( ) -> std::vector< std::shared_ptr< BodyDeformationModel > >& { return models; }, false );
    auto linear = std::make_shared< ground_stations::LinearStationMotionModel >( slope, 0.0 );
    auto combined = std::make_shared< ground_stations::CombinedStationMotionModel >(
            std::vector< std::shared_ptr< ground_stations::StationMotionModel > >{ deformation, linear } );
    auto station = std::make_shared< GroundStationState >( nominal, coordinate_conversions::cartesian_position, nullptr, combined );
    BOOST_CHECK_SMALL( ( station->getCartesianPositionInTime( 0.0, "Earth" ) - nominal - Eigen::Vector3d( 0.1, 0.2, 0.3 ) ).norm( ),
                       5.0e-10 );
    BOOST_CHECK_THROW( station->getCartesianStateInTime( 0.0, "Earth" ), CoverageError );
    BOOST_CHECK_SMALL( ( station->getCartesianStateInTime( 4.0, "Earth" ).tail< 3 >( ) - 2.0 * slope ).norm( ), 1.0e-16 );
    models.push_back( std::make_shared< FunctionDeformation >( []( double ) -> Eigen::Vector3d { throw CoverageError( ); } ) );
    BOOST_CHECK_SMALL( ( station->getCartesianPositionInTime( 0.0, "Earth" ) - nominal - Eigen::Vector3d( 0.1, 0.2, 0.3 ) ).norm( ),
                       5.0e-10 );
    BOOST_CHECK_THROW( station->getCartesianStateInTime( 0.0, "Earth" ), CoverageError );
    auto strict = stationWithModels( models, true );
    BOOST_CHECK_THROW( strict->getCartesianPositionInTime( 4.0, "Earth" ), CoverageError );
    models.back( ) = std::make_shared< FunctionDeformation >(
            []( double ) { return Eigen::Vector3d::Constant( std::numeric_limits< double >::quiet_NaN( ) ).eval( ); } );
    BOOST_CHECK_THROW( station->getCartesianPositionInTime( 4.0, "Earth" ), std::runtime_error );
}
