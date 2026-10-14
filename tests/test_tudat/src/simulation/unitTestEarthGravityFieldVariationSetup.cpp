/* Copyright (c) 2026, Delft University of Technology. Modified BSD license. */
#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createEarthGravityFieldVariations.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"
#include "tudat/astro/gravitation/earthSolidBodyTideGravityFieldVariations.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/interface/sofa/fundamentalArguments.h"

using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
namespace
{
constexpr double testTime = ( 57753.5 - 51544.5 ) * 86400.0 + 69.0;
const Eigen::Vector3d probePosition( 7e6, -2e6, 1e6 );
gravitation::OceanPoleGravityConstants oceanConstants( )
{
    return { 7.292115e-5, 6.673e-11, 1025.0, 9.7803278, 3.986004418e14, 6378136.6, { 0.687, 0.0036 } };
}
input_output::OceanPoleGravityData oceanData( )
{
    input_output::OceanPoleGravityData data;
    data.source = "controlled equilibrium integrals for factory wiring";
    data.maximumDegree = 2;
    for( int m = 0; m <= 2; ++m )
    {
        data.coefficients[ { 2, m } ] = { { -0.24 + 0.02 * m, -0.005 }, { 0.005, 0.19 - 0.01 * m }, 1 };
    }
    // The zonal sine integral is identically zero in the published A/B convention.
    data.coefficients.at( { 2, 0 } ).b = 0.0;
    return data;
}
BodyListSettings settings( const std::shared_ptr< ModelInterpolationSettings >& interpolation = nullptr )
{
    input_output::EarthOrientationData data;
    data.source = "controlled selected clock and pole";
    for( int i = 0; i < 4; ++i )
    {
        input_output::EarthOrientationRecord row;
        row.utcModifiedJulianDate = 57752.0 + i;
        row.polarMotion = { 1e-6 + i * 1e-8, -2e-6 };
        row.ut1MinusUtc = i < 2 ? -0.4 : 0.6;
        row.celestialPoleOffsets = std::array< double, 2 >{ 0.0, 0.0 };
        data.records[ row.utcModifiedJulianDate ] = row;
    }
    BodyListSettings result( "SSB", "GCRS" );
    for( const auto& name : { "Earth", "Moon", "Sun", "Probe" } )
    {
        result.addSettings( name );
    }
    Eigen::Vector6d earth = Eigen::Vector6d::Zero( );
    earth( 0 ) = 1e8;
    Eigen::Vector6d moon = Eigen::Vector6d::Zero( );
    moon.head< 3 >( ) = Eigen::Vector3d( 4e8, 2e8, -1e8 );
    Eigen::Vector6d sun = Eigen::Vector6d::Zero( );
    sun.head< 3 >( ) = Eigen::Vector3d( 1e11, -1.3e11, 2e10 );
    result.at( "Earth" )->ephemerisSettings = constantEphemerisSettings( earth, "SSB", "GCRS" );
    result.at( "Moon" )->ephemerisSettings = constantEphemerisSettings( moon, "Earth", "GCRS" );
    result.at( "Sun" )->ephemerisSettings = constantEphemerisSettings( sun, "SSB", "GCRS" );
    earth.head< 3 >( ) += probePosition;
    result.at( "Probe" )->ephemerisSettings = constantEphemerisSettings( earth, "SSB", "GCRS" );
    result.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings( iau_2006, "GCRS", nullptr, nullptr, nullptr, data );
    Eigen::MatrixXd cosine = Eigen::MatrixXd::Zero( 5, 5 );
    cosine( 0, 0 ) = 1.0;
    // Deliberately differ from conventional source GM/R to verify consuming-field scaling.
    result.at( "Earth" )->gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
            4e14, 6.4e6, cosine, Eigen::MatrixXd::Zero( 5, 5 ), "ITRS", TUDAT_NAN, TideSystem::zeroTide );
    result.at( "Moon" )->gravityFieldSettings = centralGravitySettings( 4.902800066e12 );
    result.at( "Sun" )->gravityFieldSettings = centralGravitySettings( 1.32712440018e20 );
    result.at( "Earth" )->gravityFieldVariationSettings = {
        iers2010SolidBodyTideGravityFieldVariationSettings( interpolation ),
        solidPoleTideGravityFieldVariationSettings( interpolation ),
        oceanPoleTideGravityFieldVariationSettings( oceanData( ), { 0, 0, -0.3075 }, oceanConstants( ), 2, 2, interpolation )
    };
    return result;
}
std::shared_ptr< gravitation::TimeDependentSphericalHarmonicsGravityField > field( const SystemOfBodies& bodies )
{
    return std::dynamic_pointer_cast< gravitation::TimeDependentSphericalHarmonicsGravityField >(
            bodies.at( "Earth" )->getGravityFieldModel( ) );
}
void update( const SystemOfBodies& bodies )
{
    for( const auto& name : { "Earth", "Moon", "Sun", "Probe" } )
    {
        bodies.at( name )->setStateFromEphemeris( testTime );
    }
    bodies.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( testTime );
    field( bodies )->update( testTime );
}
}  // namespace

BOOST_AUTO_TEST_CASE( completeAssemblyUsesSelectedClocksPoleScalingAndForceOrientation )
{
    const auto bodies = createSystemOfBodies( settings( ) );
    update( bodies );
    const auto gravity = field( bodies );
    const auto rotation =
            std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >( bodies.at( "Earth" )->getRotationalEphemeris( ) );
    const auto converter = rotation->getAnglesCalculator( )->getTerrestrialTimeScaleConverter( );
    const Time tai = converter->getTaiTime( tdb_scale, Time( testTime ) );
    const Time tt = convertTAItoTT( tai );
    const Eigen::Vector6d arguments = sofa_interface::calculateDoodsonFundamentalArguments(
            static_cast< double >( tt ), static_cast< double >( tt ), static_cast< double >( converter->getUt1FromTai( tai ) ) );
    const double radiansPerArcsecond = mathematical_constants::PI / ( 180 * 3600 );
    Eigen::Vector2d wobble = rotation->getAnglesCalculator( )->getPolarMotionCalculator( )->getPositionOfCipInItrsAtTt( tt ) -
            earth_orientation::getSecularPolePositionInMas( static_cast< double >( tt ) ) * ( radiansPerArcsecond / 1000 );
    wobble.y( ) = -wobble.y( );
    const auto earth = bodies.at( "Earth" );
    gravitation::EarthSolidBodyTideGravityFieldVariations solid(
            std::bind( &Body::getState, earth ),
            std::bind( &Body::getCurrentRotationToLocalFrame, earth ),
            { std::bind( &Body::getState, bodies.at( "Moon" ) ), std::bind( &Body::getState, bodies.at( "Sun" ) ) },
            6378136.3,
            [] { return 3.986004415e14; },
            { [] { return 4.902800066e12; }, [] { return 1.32712440018e20; } },
            { "Moon", "Sun" },
            [ arguments ]( double ) { return arguments; },
            gravity->getReferenceRadius( ),
            gravity->getGravitationalParameter( ),
            TideSystem::zeroTide );
    gravitation::SolidPoleTideGravityFieldVariations pole(
            [ wobble, radiansPerArcsecond ]( double ) -> Eigen::Vector2d { return wobble / radiansPerArcsecond; },
            3.986004415e14,
            6378136.3,
            gravity->getGravitationalParameter( ),
            gravity->getReferenceRadius( ) );
    gravitation::OceanPoleTideGravityFieldVariations ocean(
            oceanData( ),
            [ wobble ]( double ) { return wobble; },
            { 0, 0, -0.3075 },
            oceanConstants( ),
            gravity->getGravitationalParameter( ),
            gravity->getReferenceRadius( ),
            2,
            0,
            2,
            2 );
    Eigen::MatrixXd cosine = Eigen::MatrixXd::Zero( 5, 5 ), sine = cosine;
    for( auto* model : std::vector< gravitation::GravityFieldVariations* >{ &solid, &pole, &ocean } )
    {
        model->addSphericalHarmonicsCorrections( testTime, sine, cosine );
    }
    BOOST_CHECK_SMALL( ( gravity->getCosineCoefficients( ) - gravity->getNominalCosineCoefficients( ) - cosine ).norm( ), 1e-21 );
    BOOST_CHECK_SMALL( ( gravity->getSineCoefficients( ) - sine ).norm( ), 1e-21 );
    BOOST_CHECK_GT( ( arguments - sofa_interface::calculateDoodsonFundamentalArguments( testTime ) ).norm( ), 1e-5 );
    const auto acceleration = createSphericalHarmonicsGravityAcceleration(
            bodies.at( "Probe" ), earth, "Probe", "Earth", sphericalHarmonicAcceleration( 4, 4 ), false, false );
    acceleration->updateMembers( testTime );
    gravitation::SphericalHarmonicsGravityField expected(
            gravity->getGravitationalParameter( ), gravity->getReferenceRadius( ), cosine, sine, "ITRS" );
    const Eigen::Quaterniond toFixed = earth->getCurrentRotationToLocalFrame( );
    const Eigen::Vector3d expectedForce = toFixed.inverse( ) * expected.getGradientOfPotential( toFixed * probePosition );
    BOOST_CHECK_SMALL( ( acceleration->getAcceleration( ) - expectedForce ).norm( ), 1e-20 );
}

BOOST_AUTO_TEST_CASE( interpolatedAssemblyDefersSamplingAndMatchesDirectAtAnInteriorNode )
{
    const auto direct = createSystemOfBodies( settings( ) );
    update( direct );
    const auto interpolation = std::make_shared< ModelInterpolationSettings >(
            testTime - 20,
            testTime + 30,
            10,
            std::make_shared< interpolators::InterpolatorSettings >( interpolators::linear_interpolator ) );
    const auto sampled = createSystemOfBodies( settings( interpolation ) );
    update( sampled );
    BOOST_CHECK_SMALL( ( field( direct )->getCosineCoefficients( ) - field( sampled )->getCosineCoefficients( ) ).norm( ), 1e-21 );
    BOOST_CHECK_SMALL( ( field( direct )->getSineCoefficients( ) - field( sampled )->getSineCoefficients( ) ).norm( ), 1e-21 );
}

BOOST_AUTO_TEST_CASE( groupedOverlapFrameClockTideSystemAndCoverageRejections )
{
    for( int scenario = 0; scenario < 10; ++scenario )
    {
        auto setup = settings( );
        auto& variations = setup.at( "Earth" )->gravityFieldVariationSettings;
        if( scenario == 0 || scenario == 2 || scenario == 9 )
        {
            const int size = scenario == 9 ? 4 : 5;
            Eigen::MatrixXd cosine = Eigen::MatrixXd::Zero( size, size );
            cosine( 0, 0 ) = 1;
            setup.at( "Earth" )->gravityFieldSettings = std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                    4e14,
                    6.4e6,
                    cosine,
                    Eigen::MatrixXd::Zero( size, size ),
                    scenario == 2 ? "IAU_Earth" : "ITRS",
                    TUDAT_NAN,
                    scenario == 0 ? TideSystem::unknown : TideSystem::zeroTide );
        }
        if( scenario == 1 || scenario == 3 )
        {
            const auto selected =
                    std::dynamic_pointer_cast< GcrsToItrsRotationModelSettings >( setup.at( "Earth" )->rotationModelSettings );
            const auto changed = std::make_shared< GcrsToItrsRotationModelSettings >(
                    iau_2006, scenario == 3 ? "J2000" : "GCRS", "", scenario == 1 ? tt_scale : tdb_scale );
            changed->setEopData( *selected->getEopData( ) );
            setup.at( "Earth" )->rotationModelSettings = changed;
        }
        if( scenario >= 4 && scenario <= 6 )
        {
            variations.push_back( variations.at( scenario - 4 ) );
        }
        if( scenario == 7 )
        {
            variations.push_back( fixedSingleDegreeLoveNumberGravityFieldVariationSettings( "Moon", 0.3, 2 ) );
        }
        if( scenario == 8 )
        {
            variations.push_back( std::make_shared< ModeCoupledSolidBodyGravityFieldVariationSettings >(
                    std::vector< std::string >{ "Sun" }, std::map< std::pair< int, int >, std::map< std::pair< int, int >, double > >{} ) );
        }
        BOOST_CHECK_THROW( createSystemOfBodies( setup ), std::runtime_error );
    }
}
