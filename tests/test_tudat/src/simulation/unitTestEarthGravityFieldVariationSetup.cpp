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
#include "tudat/astro/basic_astro/dateTime.h"
#include <sstream>
#include <cmath>
#if TUDAT_BUILD_WITH_SERIALIZATION
#include "tudat/io/serialization/archives.h"
#include "tudat/io/serialization/registrations_environment.h"
#endif

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

//! Read the labelled RL07 provider extract shared by the raw and interpolation tests.
input_output::Aod1bGravityFieldData aodData( )
{
    const std::string source = __FILE__;

    return input_output::readAod1bGravityField( source.substr( 0, source.find_last_of( "/\\" ) ) +
                                                "/../../reference_files/aod1b_gravity/AOD1B_2026-01-01_degree2_extract_07.asc" );
}

//! Supply two adjacent complete C20 pieces for testing the consuming baseline and selected time conversion.
input_output::IcgemGravityFieldData temporalData( )
{
    // The source potential differs from the nominal field; no omitted selected coefficient is synthesized.
    std::istringstream source(
            "begin_of_head\nproduct_type gravity_field\nmodelname controlled_temporal_setup\n"
            "earth_gravity_constant 3e14\nradius 6000000\nmax_degree 2\nerrors no\n"
            "norm fully_normalized\ntide_system zero_tide\nformat icgem2.0\nend_of_head\n"
            "gfct 2 0 -4.8e-4 0 20161231.0000 20170101.0000\n"
            "gfct 2 0 -4.6e-4 0 20170101.0000 20170102.0000\n" );

    return input_output::readIcgemGravityField( source, "controlled temporal setup coefficients" );
}

//! Declare a controlled, fully normalized potential with its source GM and radius.
gravitation::HarmonicGravityMetadata harmonicMetadata( )
{
    return { gravitation::HarmonicGravityComponent::oceanMass, 3e14, 6e6, true, "controlled potential revision",
             gravitation::HarmonicGravityLoading::included };
}

//! Supply one nontrivial sectorial term so both phase quadratures are observable.
gravitation::HarmonicGravityConstituent harmonicWave( const std::string& doodson = "255.555" )
{
    gravitation::HarmonicGravityConstituent wave;

    wave.doodson = doodson;
    wave.multipliers = input_output::parseDoodsonMultipliers( doodson );
    wave.cCos = wave.cSin = wave.sCos = wave.sSin = Eigen::MatrixXd::Zero( 3, 3 );
    wave.cCos( 2, 2 ) = 2e-10;
    wave.cSin( 2, 2 ) = 3e-10;
    wave.sCos( 2, 2 ) = -1e-10;
    wave.sSin( 2, 2 ) = 4e-10;
    return wave;
}

//! Assemble selected Earth rotation and field settings; a UT1 offset tests phase-clock routing.
BodyListSettings settings( const std::shared_ptr< ModelInterpolationSettings >& interpolation = nullptr, double ut1Offset = 0.0 )
{
    input_output::EarthOrientationData data;
    data.source = "controlled selected clock and pole";
    for( int i = 0; i < 4; ++i )
    {
        input_output::EarthOrientationRecord row;
        row.utcModifiedJulianDate = 57752.0 + i;
        row.polarMotion = { 1e-6 + i * 1e-8, -2e-6 };
        row.ut1MinusUtc = ( i < 2 ? -0.4 : 0.6 ) + ut1Offset;
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

BOOST_AUTO_TEST_CASE( harmonicAssemblyUsesSelectedUt1ReferenceScalingAndDeferredNode )
{
    // Compare the assembled field against scalar quadratures and explicit source-to-field scaling.
    const auto metadata = harmonicMetadata( );

    const auto wave = harmonicWave( );

    auto directSettings = settings( );

    directSettings.at( "Earth" )->gravityFieldVariationSettings = { std::make_shared< HarmonicGravityFieldVariationSettings >(
            std::vector< gravitation::HarmonicGravityConstituent >{ wave }, metadata, 2, 2 ) };
    const auto direct = createSystemOfBodies( directSettings );

    update( direct );
    const auto rotation =
            std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >( direct.at( "Earth" )->getRotationalEphemeris( ) );

    const auto converter = rotation->getAnglesCalculator( )->getTerrestrialTimeScaleConverter( );

    const Time tai = converter->getTaiTime( tdb_scale, testTime );

    const double tt = static_cast< double >( convertTAItoTT( tai ) );

    const auto arguments =
            sofa_interface::calculateDoodsonFundamentalArguments( tt, tt, static_cast< double >( converter->getUt1FromTai( tai ) ) );

    const double phase = 2 * arguments( 0 );

    const double scale = ( 3e14 / 4e14 ) * std::pow( 6e6 / 6.4e6, 2 );

    BOOST_CHECK_SMALL(
            field( direct )->getCosineCoefficients( )( 2, 2 ) - scale * ( 2e-10 * std::cos( phase ) + 3e-10 * std::sin( phase ) ), 1e-25 );
    BOOST_CHECK_SMALL( field( direct )->getSineCoefficients( )( 2, 2 ) - scale * ( -1e-10 * std::cos( phase ) + 4e-10 * std::sin( phase ) ),
                       1e-25 );

    // A changed selected UT1 series must change the same coefficients at unchanged TDB.
    auto shiftedSettings = settings( nullptr, 0.2 );

    shiftedSettings.at( "Earth" )->gravityFieldVariationSettings = directSettings.at( "Earth" )->gravityFieldVariationSettings;
    const auto shifted = createSystemOfBodies( shiftedSettings );

    update( shifted );
    BOOST_CHECK_GT( ( field( direct )->getCosineCoefficients( ) - field( shifted )->getCosineCoefficients( ) ).norm( ) +
                            ( field( direct )->getSineCoefficients( ) - field( shifted )->getSineCoefficients( ) ).norm( ),
                    1e-15 );

    // Deferred interpolation samples the same selected clock, including this interior node.
    const auto interpolation = std::make_shared< ModelInterpolationSettings >(
            testTime - 20,
            testTime + 30,
            10,
            std::make_shared< interpolators::InterpolatorSettings >( interpolators::linear_interpolator ) );

    auto sampledSettings = settings( );

    sampledSettings.at( "Earth" )->gravityFieldVariationSettings = { std::make_shared< HarmonicGravityFieldVariationSettings >(
            std::vector< gravitation::HarmonicGravityConstituent >{ wave }, metadata, 2, 2, interpolation ) };
    const auto sampled = createSystemOfBodies( sampledSettings );

    update( sampled );
    BOOST_CHECK_SMALL( ( field( direct )->getCosineCoefficients( ) - field( sampled )->getCosineCoefficients( ) ).norm( ), 1e-25 );
    BOOST_CHECK_SMALL( ( field( direct )->getSineCoefficients( ) - field( sampled )->getSineCoefficients( ) ).norm( ), 1e-25 );
}

BOOST_AUTO_TEST_CASE( rawPotentialSettingsRetainFinalAdmittanceInventoryAndComponentIdentity )
{
    // A controlled major-to-target map tests overlap after admittance, without claiming a provider product.
    std::istringstream source(
            "Coefficients to compute variations in normalized Stokes coefficients (unit = 10^-12)\n"
            "Doodson Darw l m DelC+ DelS+ DelC- DelS-\n"
            "255.555 M2 2 0 0 0 0 0\n255.555 M2 2 1 0 0 0 0\n255.555 M2 2 2 1 2 3 4\n" );

    const auto raw = input_output::readHarmonicPotentialData( source );

    gravitation::HarmonicGravityAdmittance admittance;

    admittance.majorDoodsonIds = { "255.555" };
    admittance.targetDoodsonIds = { "255.555", "245.655" };
    admittance.matrix = Eigen::MatrixXd( 2, 1 );
    admittance.matrix << 1.0, 0.5;
    const auto converted = harmonicPotentialGravityFieldVariationSettings( raw, harmonicMetadata( ), 2, 2, {}, admittance );

    const auto ocean = std::dynamic_pointer_cast< HarmonicGravityFieldVariationSettings >( converted );

    BOOST_REQUIRE( ocean );
    BOOST_REQUIRE_EQUAL( ocean->getConstituents( ).size( ), 2 );
    BOOST_CHECK_EQUAL( ocean->getConstituents( )[ 1 ].doodson, "245.655" );
    BOOST_CHECK_EQUAL( ocean->getConstituents( )[ 1 ].cCos( 2, 2 ), 2e-12 );
    const auto disjointOcean = std::make_shared< HarmonicGravityFieldVariationSettings >(
            std::vector< gravitation::HarmonicGravityConstituent >{ harmonicWave( "273.555" ) }, harmonicMetadata( ), 2, 2 );

    auto airMetadata = harmonicMetadata( );

    airMetadata.component = gravitation::HarmonicGravityComponent::directAtmosphericMass;
    const auto air = std::make_shared< HarmonicGravityFieldVariationSettings >( ocean->getConstituents( ), airMetadata, 2, 2 );

    auto setup = settings( );

    setup.at( "Earth" )->gravityFieldVariationSettings = { ocean, disjointOcean, air };
    const auto bodies = createSystemOfBodies( setup );

    update( bodies );
    const auto models = field( bodies )->getGravityFieldVariationsSet( );

    BOOST_CHECK_NE( ocean->getIdentifier( ), disjointOcean->getIdentifier( ) );
    BOOST_CHECK( models->getGravityFieldVariation( gravitation::ocean_tide, ocean->getIdentifier( ) ).first );
    BOOST_CHECK( models->getGravityFieldVariation( gravitation::ocean_tide, disjointOcean->getIdentifier( ) ).first );
    BOOST_CHECK( models->getGravityFieldVariation( gravitation::atmospheric_tide, air->getIdentifier( ) ).first );

    // A second ocean model containing the inferred N2 target would apply the same mass term twice.
    setup.at( "Earth" )->gravityFieldVariationSettings.push_back( std::make_shared< HarmonicGravityFieldVariationSettings >(
            std::vector< gravitation::HarmonicGravityConstituent >{ harmonicWave( "245.655" ) }, harmonicMetadata( ), 2, 2 ) );
    BOOST_CHECK_THROW( createSystemOfBodies( setup ), std::runtime_error );
}

//! Exercise selected time conversion, additive field scaling and separate component retrieval through Earth setup.
BOOST_AUTO_TEST_CASE( aodAssemblyUsesSelectedConverterScalingAndComponentIdentity )
{
    // A deliberately coarse selected TT/TDB table makes an accidental default converter observable.
    const auto data = aodData( );

    const double gps = data.samples.at( 4 ).gpsSecondsSinceJ2000;

    auto setup = settings( );

    const auto selected = std::dynamic_pointer_cast< GcrsToItrsRotationModelSettings >( setup.at( "Earth" )->rotationModelSettings );

    selected->setTdbToTtInterpolatorSettings( interpolators::interpolatorGenerationSettings< double >(
            std::make_shared< interpolators::InterpolatorSettings >( interpolators::linear_interpolator ),
            gps - 30.0 * 86400.0,
            gps + 150.0 * 86400.0,
            90.0 * 86400.0 ) );
    setup.at( "Earth" )->gravityFieldVariationSettings = {
        aod1bGravityFieldVariationSettings( { data }, input_output::Aod1bComponent::atmosphere, 2, 2 ),
        aod1bGravityFieldVariationSettings( { data }, input_output::Aod1bComponent::ocean, 2, 2 )
    };

    const auto bodies = createSystemOfBodies( setup );

    const auto rotation =
            std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >( bodies.at( "Earth" )->getRotationalEphemeris( ) );

    const auto converter = rotation->getAnglesCalculator( )->getTerrestrialTimeScaleConverter( );

    const double tdb = static_cast< double >( converter->getCurrentTime< Time >( tai_scale, tdb_scale, Time( gps ) + 19.0 ) );

    const auto native = std::make_shared< earth_orientation::TerrestrialTimeScaleConverter >( nullptr, nullptr );

    BOOST_CHECK_GT( std::abs( static_cast< double >( native->getTaiTime( tdb_scale, Time( tdb ) ) -
                                                     converter->getTaiTime( tdb_scale, Time( tdb ) ) ) ),
                    1e-6 );

    // The field receives ATM and OCN as anomalies, with each provider potential scaled exactly once.
    const auto gravity = field( bodies );

    const double scale = data.gravitationalParameter / gravity->getGravitationalParameter( ) *
            std::pow( data.referenceRadius / gravity->getReferenceRadius( ), 2 );

    gravity->update( tdb );
    BOOST_CHECK_SMALL( gravity->getCosineCoefficients( )( 2, 0 ) -
                               scale *
                                       ( data.samples.at( 4 ).coefficients.at( { 2, 0 } ).values[ 0 ] +
                                         data.samples.at( 5 ).coefficients.at( { 2, 0 } ).values[ 0 ] ),
                       1e-25 );
    BOOST_CHECK_SMALL( gravity->getSineCoefficients( )( 2, 1 ) -
                               scale *
                                       ( data.samples.at( 4 ).coefficients.at( { 2, 1 } ).values[ 1 ] +
                                         data.samples.at( 5 ).coefficients.at( { 2, 1 } ).values[ 1 ] ),
                       1e-25 );

    const auto models = gravity->getGravityFieldVariationsSet( );

    const auto atmosphere = models->getGravityFieldVariation( gravitation::aod1b_rl07, "AOD1B_RL07_ATM" );

    const auto ocean = models->getGravityFieldVariation( gravitation::aod1b_rl07, "AOD1B_RL07_OCN" );

    BOOST_REQUIRE( atmosphere.first && ocean.first );
    BOOST_CHECK( atmosphere.second != ocean.second );
    setup.at( "Earth" )->gravityFieldVariationSettings.push_back(
            aod1bGravityFieldVariationSettings( { data }, input_output::Aod1bComponent::global, 2, 2 ) );
    BOOST_CHECK_THROW( createSystemOfBodies( setup ), std::runtime_error );
    setup.at( "Earth" )->gravityFieldVariationSettings = {
        aod1bGravityFieldVariationSettings( { data }, input_output::Aod1bComponent::atmosphere, 2, 2 ),
        aod1bGravityFieldVariationSettings( { data }, input_output::Aod1bComponent::atmosphere, 2, 2 )
    };
    BOOST_CHECK_THROW( createSystemOfBodies( setup ), std::runtime_error );
    setup.at( "Earth" )->gravityFieldVariationSettings = {
        aod1bGravityFieldVariationSettings( { data }, input_output::Aod1bComponent::oceanBottomPressure, 2, 2 ),
        aod1bGravityFieldVariationSettings( { data }, input_output::Aod1bComponent::ocean, 2, 2 )
    };
    BOOST_CHECK_THROW( createSystemOfBodies( setup ), std::runtime_error );
}

//! Assemble temporal anomalies against the actual field, preserving source-clock boundaries and the adopted nominal snapshot.
BOOST_AUTO_TEST_CASE( temporalAssemblyUsesSelectedConverterBaselineAndIntervalOwnership )
{
    const auto data = temporalData( );

    auto setup = settings( );

    const auto selected = std::dynamic_pointer_cast< GcrsToItrsRotationModelSettings >( setup.at( "Earth" )->rotationModelSettings );

    // A coarse selected time conversion makes accidentally constructing a default converter observable.
    selected->setTdbToTtInterpolatorSettings( interpolators::interpolatorGenerationSettings< double >(
            std::make_shared< interpolators::InterpolatorSettings >( interpolators::linear_interpolator ),
            testTime - 30.0 * 86400.0,
            testTime + 150.0 * 86400.0,
            90.0 * 86400.0 ) );

    // Adopt a nonzero nominal C20 so the complete provider baseline must replace, rather than add to, it.
    const auto nominalSettings =
            std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >( setup.at( "Earth" )->gravityFieldSettings );

    Eigen::MatrixXd nominal = nominalSettings->getCosineCoefficients( );

    nominal( 2, 0 ) = -4.7e-4;
    nominalSettings->resetCosineCoefficients( nominal );

    const auto temporal = icgemTemporalGravityFieldVariationSettings( data, tai_scale, gravitation::IcgemYearConvention::julianYear, 2, 0 );

    setup.at( "Earth" )->gravityFieldVariationSettings = { temporal };
    const auto bodies = createSystemOfBodies( setup );

    const auto gravity = field( bodies );

    const auto rotation =
            std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >( bodies.at( "Earth" )->getRotationalEphemeris( ) );

    const auto converter = rotation->getAnglesCalculator( )->getTerrestrialTimeScaleConverter( );

    const Time boundary = timeFromDecomposedDateTime< Time >( 2017, 1, 1, 0, 0, 0.0L );

    const double tdb = static_cast< double >( converter->getCurrentTime< Time >( tai_scale, tdb_scale, boundary ) );

    const double scale = data.gravitationalParameter / gravity->getGravitationalParameter( ) *
            std::pow( data.referenceRadius / gravity->getReferenceRadius( ), 2 );

    const auto native = std::make_shared< earth_orientation::TerrestrialTimeScaleConverter >( nullptr, nullptr );

    BOOST_CHECK_GT( std::abs( static_cast< double >( native->getTaiTime( tdb_scale, Time( tdb ) ) -
                                                     converter->getTaiTime( tdb_scale, Time( tdb ) ) ) ),
                    1e-6 );
    gravity->update( tdb );

    // The exact scalar boundary belongs to the right base. Only binary multiplication/addition rounding remains.
    BOOST_CHECK_SMALL( gravity->getCosineCoefficients( )( 2, 0 ) - scale * -4.6e-4, 2e-19 );
    BOOST_CHECK( gravity->getGravityFieldVariationsSet( )
                         ->getGravityFieldVariation( gravitation::icgem_temporal, "ICGEM_TEMPORAL_controlled_temporal_setup" )
                         .first );

    // An intentional later nominal change adds an offset without redefining the adopted anomaly baseline.
    gravity->setNominalCosineCoefficient( 2, 0, nominal( 2, 0 ) + 1e-5 );
    gravity->update( tdb + 1.0 );
    BOOST_CHECK_SMALL( gravity->getCosineCoefficients( )( 2, 0 ) - ( scale * -4.6e-4 + 1e-5 ), 2e-19 );

    setup.at( "Earth" )->gravityFieldVariationSettings.push_back( temporal );
    BOOST_CHECK_THROW( createSystemOfBodies( setup ), std::runtime_error );

    auto incompatible = data;

    incompatible.tideSystem = TideSystem::tideFree;
    setup.at( "Earth" )->gravityFieldVariationSettings = { icgemTemporalGravityFieldVariationSettings(
            incompatible, tai_scale, gravitation::IcgemYearConvention::julianYear, 2, 0 ) };
    BOOST_CHECK_THROW( createSystemOfBodies( setup ), std::runtime_error );
}

#if TUDAT_BUILD_WITH_SERIALIZATION
BOOST_AUTO_TEST_CASE( earthGravitySettingsRoundTripCompleteInputsAndFinalInventory )
{
    // Retain polymorphic presets, source provenance, physical constants and canonical coefficients exactly.
    const auto interpolation = std::make_shared< ModelInterpolationSettings >(
            testTime - 20,
            testTime + 30,
            10,
            std::make_shared< interpolators::InterpolatorSettings >( interpolators::linear_interpolator ) );

    std::vector< std::shared_ptr< GravityFieldVariationSettings > > original{
        iers2010SolidBodyTideGravityFieldVariationSettings( interpolation ),
        solidPoleTideGravityFieldVariationSettings( ),
        oceanPoleTideGravityFieldVariationSettings( oceanData( ), { 0, 0, -0.3075 }, oceanConstants( ), 2, 2 ),
        std::make_shared< HarmonicGravityFieldVariationSettings >(
                std::vector< gravitation::HarmonicGravityConstituent >{ harmonicWave( ) }, harmonicMetadata( ), 2, 2, interpolation ),
        aod1bGravityFieldVariationSettings( { aodData( ) }, input_output::Aod1bComponent::global, 2, 2 ),
        icgemTemporalGravityFieldVariationSettings( temporalData( ),
                                                    tt_scale,
                                                    gravitation::IcgemYearConvention::fractionalGregorianYear,
                                                    2,
                                                    0,
                                                    { { "20161231.0000", testTime - 100.0 } } )
    };

    std::stringstream stream;

    {
        cereal::BinaryOutputArchive archive( stream );

        archive( original );
    }
    std::vector< std::shared_ptr< GravityFieldVariationSettings > > restored;

    {
        cereal::BinaryInputArchive archive( stream );

        archive( restored );
    }
    BOOST_REQUIRE_EQUAL( original.size( ), restored.size( ) );
    for( std::size_t i = 0; i < original.size( ); ++i )
    {
        BOOST_REQUIRE( restored[ i ] );
        BOOST_CHECK( *original[ i ] == *restored[ i ] );
    }
}
#endif
