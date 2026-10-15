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
    explicit FunctionDeformation( const std::function< Eigen::Vector3d( const Time& ) >& function ): function( function ) {}

    Eigen::Vector3d evaluateDisplacement( const Time& time, const std::shared_ptr< GroundStationState > ) override
    {
        lastTime = static_cast< double >( time );
        return function( time );
    }

    std::function< Eigen::Vector3d( const Time& ) > function;
    double lastTime = std::numeric_limits< double >::quiet_NaN( );
};

struct CoverageError : std::out_of_range {
    CoverageError( ): std::out_of_range( "Unavailable derivative sample" ) {}
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

std::shared_ptr< GroundStationState > stationWithModels( std::vector< std::shared_ptr< BodyDeformationModel > >& models,
                                                         const bool throwWhenUnavailable = true )
{
    auto motion = std::make_shared< ground_stations::BodyDeformationStationMotionModel >(
            [ &models ]( ) -> std::vector< std::shared_ptr< BodyDeformationModel > >& { return models; }, throwWhenUnavailable );
    return std::make_shared< GroundStationState >(
            Eigen::Vector3d( 4075578.385, 931852.890, 4801570.154 ), coordinate_conversions::cartesian_position, nullptr, motion );
}
}  // namespace

//! Compare finite-difference velocities with exact constant, linear and sinusoidal derivatives.
BOOST_AUTO_TEST_CASE( analytic_constant_linear_and_sinusoidal_rates )
{
    const Eigen::Vector3d position( 6.4e6, 1.0e6, 2.0e6 );
    FunctionDeformation constant( []( double ) { return Eigen::Vector3d( 0.1, -0.2, 0.3 ); } );

    // Compare finite-difference velocities with exact constant, linear and sinusoidal derivatives.
    BOOST_CHECK_SMALL(
            constant.calculateDisplacementState( 100.0, std::make_shared< GroundStationState >( position ) ).tail< 3 >( ).norm( ),
            1.0e-15 );
    const Eigen::Vector3d slope( 1.0e-4, -3.0e-4, 7.0e-4 );
    FunctionDeformation linear( [ = ]( double t ) { return ( Eigen::Vector3d( 0.1, 0.2, 0.3 ) + t * slope ).eval( ); } );
    BOOST_CHECK_SMALL(
            ( linear.calculateDisplacementState( 100.0, std::make_shared< GroundStationState >( position ) ).tail< 3 >( ) - slope ).norm( ),
            1.0e-16 );
    FunctionDeformation sinusoid( []( double t ) {
        return Eigen::Vector3d( 0.3 * std::sin( 0.01 * t ), 0.1 * std::cos( 0.02 * t ), 0.2 * std::sin( 0.03 * t ) );
    } );
    const double epoch = 43.0;
    const Eigen::Vector3d expected( 0.003 * std::cos( 0.01 * epoch ), -0.002 * std::sin( 0.02 * epoch ), 0.006 * std::cos( 0.03 * epoch ) );
    std::array< double, 3 > errors;
    for( std::size_t i = 0; i < errors.size( ); ++i )
    {
        sinusoid.setDisplacementDerivativeTimeStep( 4.0 / std::pow( 2.0, i ) );
        errors[ i ] =
                ( sinusoid.calculateDisplacementState( epoch, std::make_shared< GroundStationState >( position ) ).tail< 3 >( ) - expected )
                        .norm( );
        BOOST_CHECK_EQUAL( sinusoid.lastTime, epoch );
    }
    BOOST_CHECK_GT( errors[ 0 ] / errors[ 1 ], 15.0 );
    BOOST_CHECK_GT( errors[ 1 ] / errors[ 2 ], 15.0 );
    BOOST_CHECK_SMALL( errors[ 2 ], 2.0e-10 );
    std::cout << "Analytic sinusoid errors h=4,2,1 s: " << errors[ 0 ] << ' ' << errors[ 1 ] << ' ' << errors[ 2 ] << '\n';
}

//! Keep derivative samples distinct when their absolute epochs round to the same double.
BOOST_AUTO_TEST_CASE( precise_time_stencil_retains_offsets_below_double_resolution )
{
    const Time centre( 334368000.125 );
    const Eigen::Vector3d slope( 0.001, -0.002, 0.003 );
    FunctionDeformation model( [ = ]( const Time& time ) { return ( static_cast< double >( time - centre ) * slope ).eval( ); } );
    model.setDisplacementDerivativeTimeStep( std::ldexp( 1.0, -27 ) );
    auto station = std::make_shared< GroundStationState >( Eigen::Vector3d( 6378137, 0, 0 ) );

    // Keep derivative samples distinct when their absolute epochs round to the same double.
    BOOST_CHECK_SMALL( ( model.calculateDisplacementState( centre, station ).tail< 3 >( ) - slope ).norm( ), 1e-17 );
}

//! Recover the original tidal displacement after an interrupted multi-body cache update.
BOOST_AUTO_TEST_CASE( failed_body_refresh_preserves_central_cache )
{
    bool fail = false;
    auto station = std::make_shared< GroundStationState >( Eigen::Vector3d( 6378137, 0, 0 ) );
    basic_astrodynamics::BasicTidalBodyDeformation cached(
            []( const Time& ) { return Eigen::Vector6d::Zero( ).eval( ); },
            { []( const Time& t ) {
                 return ( Eigen::Vector6d( ) << 3.8e8 + static_cast< double >( t ) * 1e6, 1e8, 2e7, 0, 0, 0 ).finished( );
             },
              [ & ]( const Time& t ) {
                  if( fail )
                  {
                      throw CoverageError( );
                  }
                  return ( Eigen::Vector6d( ) << 1.5e11, 5e10, 1e10, 0, 0, 0 ).finished( );
              } },
            []( const Time& ) { return Eigen::Quaterniond::Identity( ); },
            []( ) { return 1.; },
            { []( ) { return 0.0123; }, []( ) { return 332946.; } },
            6378137.,
            { { 2, { 0.6078, 0.0847 } } } );
    const auto central = cached.calculateDisplacement( Time( 0. ), station );
    fail = true;

    // Recover the original tidal displacement after an interrupted multi-body cache update.

    // Interrupting the second body must not leave the first body cached at the wrong epoch.
    BOOST_CHECK_THROW( cached.calculateDisplacement( Time( 1. ), station ), CoverageError );
    fail = false;

    // After the callback recovers, recompute the original central displacement.
    BOOST_CHECK_SMALL( ( cached.calculateDisplacement( Time( 0. ), station ) - central ).norm( ), 1e-18 );
}

//! Preserve the oblate station basis when applying local deformation.
BOOST_AUTO_TEST_CASE( oblate_nominal_station_basis_is_preserved )
{
    struct LocalDeformation : BodyDeformationModel {
        Eigen::Vector3d evaluateDisplacement( const Time& time, const std::shared_ptr< GroundStationState > station ) override
        {
            const double t = static_cast< double >( time );
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

    // Preserve the oblate station basis when applying local deformation.
    BOOST_CHECK_GT( std::fabs( phi - 0.8 ), 0.001 );
    LocalDeformation model;
    BOOST_CHECK_SMALL( ( model.calculateDisplacementState( 123.0, station ).tail< 3 >( ) - 0.001 * east - 0.002 * north ).norm( ),
                       1.0e-16 );
}

//! Check convergence of the complete IERS deformation velocity as the derivative step shrinks.
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
    const Eigen::Vector3d reference = model->calculateDisplacementState( time, station ).tail< 3 >( );
    std::array< double, 3 > errors;
    for( std::size_t i = 0; i < errors.size( ); ++i )
    {
        model->setDisplacementDerivativeTimeStep( 1800.0 / std::pow( 2.0, i ) );
        errors[ i ] = ( model->calculateDisplacementState( time, station ).tail< 3 >( ) - reference ).norm( );
    }

    // Check convergence of the complete IERS deformation velocity as the derivative step shrinks.
    BOOST_CHECK_GT( reference.norm( ), 1.0e-7 );
    BOOST_CHECK_GT( errors[ 0 ] / errors[ 1 ], 10.0 );
    BOOST_CHECK_GT( errors[ 1 ] / errors[ 2 ], 10.0 );
    model->setDisplacementDerivativeTimeStep( 0.5 );
    const double fineDifference = ( model->calculateDisplacementState( time, station ).tail< 3 >( ) - reference ).norm( );
    BOOST_CHECK_SMALL( fineDifference, 2.0e-10 );
    std::cout << "Production IERS velocity " << reference.transpose( ) << " m/s; errors h=1800,900,450 s: " << errors[ 0 ] << ' '
              << errors[ 1 ] << ' ' << errors[ 2 ] << "; h=1 vs0.5 s: " << fineDifference << '\n';
}

//! Include both local deformation velocity and rotation of the displaced station.
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

    // Include both local deformation velocity and rotation of the displaced station.
    BOOST_CHECK_SMALL( ( inertial.head< 3 >( ) - expectedPosition ).norm( ), 2.0e-9 );
    BOOST_CHECK_SMALL( ( inertial.tail< 3 >( ) - expectedVelocity ).norm( ), 1.0e-10 );
}

}  // namespace unit_tests
}  // namespace tudat
