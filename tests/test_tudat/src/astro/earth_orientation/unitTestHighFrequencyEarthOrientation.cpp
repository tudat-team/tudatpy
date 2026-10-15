#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/astro/earth_orientation/shortPeriodEarthOrientationCorrectionCalculator.h"
#include "tudat/astro/earth_orientation/earthOrientationRates.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "desai_sibois_reference/expected.h"

// The Earth-specific factory is kept in its own translation unit.
namespace tudat
{
namespace simulation_setup
{
std::shared_ptr< ephemerides::RotationalEphemeris > createGcrsToItrsRotationModel(
        const std::shared_ptr< GcrsToItrsRotationModelSettings >& );
}
}  // namespace tudat
using namespace tudat;
using namespace tudat::earth_orientation;
namespace
{
const double uas = unit_conversions::convertArcSecondsToRadians< double >( 1e-6 );

//! Short-period polar motion and UT1 corrections.
struct ShortPeriodCorrections {
    std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > polarMotion;
    std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > ut1;
};

//! Conventional libration terms, optionally followed by the Desai and Sibois (2016) ocean-tide resource tables.
ShortPeriodCorrections corrections( const bool includeDesaiSibois,
                                    const double minimumPolarMotionAmplitude = 0.0,
                                    const double minimumUt1Amplitude = 0.0 )
{
    const std::string p = paths::getEarthOrientationDataFilesPath( ) + "/";
    std::vector< std::string > polarAmplitudes{ p + "polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt" };
    std::vector< std::string > polarMultipliers{ p + "polarMotionLibrationFundamentalArgumentMultipliersQuasiDiurnalOnly.txt" };
    std::vector< std::string > ut1Amplitudes{ p + "utcLibrationAmplitudes.txt" };
    std::vector< std::string > ut1Multipliers{ p + "utcLibrationFundamentalArgumentMultipliers.txt" };
    if( includeDesaiSibois )
    {
        polarAmplitudes.push_back( p + "desaiSibois2016PolarMotionOceanTidesAmplitudes.txt" );
        polarMultipliers.push_back( p + "desaiSibois2016PolarMotionOceanTidesFundamentalArgumentMultipliers.txt" );
        ut1Amplitudes.push_back( p + "desaiSibois2016Ut1OceanTidesAmplitudes.txt" );
        ut1Multipliers.push_back( p + "desaiSibois2016Ut1OceanTidesFundamentalArgumentMultipliers.txt" );
    }
    return { std::make_shared< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > >(
                     uas, minimumPolarMotionAmplitude, polarAmplitudes, polarMultipliers ),
             std::make_shared< ShortPeriodEarthOrientationCorrectionCalculator< double > >(
                     1e-6, minimumUt1Amplitude, ut1Amplitudes, ut1Multipliers ) };
}
ShortPeriodCorrections libration( )
{
    return corrections( false );
}
input_output::EarthOrientationData dailyData( )
{
    input_output::EarthOrientationData data;
    data.source = "controlled daily EOP";
    const Time origin( 536500837.0 );
    for( int i = -1; i <= 2; ++i )
    {
        input_output::EarthOrientationRecord record;
        record.utcModifiedJulianDate = 57754.0 + i;
        const Time utc( ( record.utcModifiedJulianDate - 51544.5 ) * 86400.0 );
        const Time tai = sofa_interface::convertUTCtoTAI( utc );
        const double dt = static_cast< double >( tai - origin );
        record.polarMotion = { 1e-6 + 2e-12 * dt, 2e-6 - 3e-12 * dt };
        record.ut1MinusUtc = -37.4 - 2e-8 * dt + static_cast< double >( tai - utc );
        record.celestialPoleOffsets = std::array< double, 2 >{ 1e-9, -2e-9 };
        record.lengthOfDay = .1;  // Deliberately inconsistent metadata, never added to UT1 rate.
        data.records.emplace( record.utcModifiedJulianDate, record );
    }
    return data;
}
std::shared_ptr< simulation_setup::GcrsToItrsRotationModelSettings > settings( )
{
    auto result = std::make_shared< simulation_setup::GcrsToItrsRotationModelSettings >(
            basic_astrodynamics::iau_2006, "GCRS", "unused preloaded EOP", basic_astrodynamics::tt_scale );
    result->setEopData( dailyData( ) );
    return result;
}
}  // namespace

BOOST_AUTO_TEST_CASE( explicit_ocean_replacement_matches_independent_phase_values_and_analytic_rates )
{
    const auto ds = corrections( true );
    const auto lib = libration( );
    for( const auto& fixture : desai_sibois_reference::fixtures )
    {
        Eigen::Vector6d a;
        for( int i = 0; i < 6; ++i ) a( i ) = fixture.arguments[ i ];
        const Eigen::Vector2d ocean =
                ds.polarMotion->getCorrectionsFromFundamentalArgument( a ) - lib.polarMotion->getCorrectionsFromFundamentalArgument( a );
        const double oceanUt1 = ds.ut1->getCorrectionsFromFundamentalArgument( a ) - lib.ut1->getCorrectionsFromFundamentalArgument( a );
        BOOST_CHECK_SMALL( ocean.x( ) / uas - fixture.values[ 0 ], 2e-10 );
        BOOST_CHECK_SMALL( ocean.y( ) / uas - fixture.values[ 1 ], 2e-10 );
        BOOST_CHECK_SMALL( oceanUt1 / 1e-6 - fixture.values[ 2 ], 2e-11 );

        // Production phase function uses explicit TT/UT1, tested against the
        // independent high-precision polynomial oracle at 14 epochs (1962--2025).
        const auto physical = sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst( fixture.tt, fixture.tt, fixture.ut1 );
        for( int i = 0; i < 6; ++i ) BOOST_CHECK_SMALL( std::sin( physical( i ) - a( i ) ), 8e-11 );
        const Eigen::Vector2d rate =
                calculateShortPeriodCorrectionRate( ds.polarMotion, Time( fixture.tt ), Time( fixture.ut1 ), .99999998 ) -
                calculateShortPeriodCorrectionRate( lib.polarMotion, Time( fixture.tt ), Time( fixture.ut1 ), .99999998 );
        const double rateUt1 = calculateShortPeriodCorrectionRate( ds.ut1, Time( fixture.tt ), Time( fixture.ut1 ), .99999998 ) -
                calculateShortPeriodCorrectionRate( lib.ut1, Time( fixture.tt ), Time( fixture.ut1 ), .99999998 );
        BOOST_CHECK_SMALL( rate.x( ) / uas - fixture.rates[ 0 ], 2e-10 );
        BOOST_CHECK_SMALL( rate.y( ) / uas - fixture.rates[ 1 ], 2e-10 );
        BOOST_CHECK_SMALL( rateUt1 / 1e-6 - fixture.rates[ 2 ], 2e-11 );
    }
}

BOOST_AUTO_TEST_CASE( amplitude_cutoffs_remove_ocean_and_libration_terms )
{
    const auto none = corrections( true, 1e6, 1e6 );
    for( const auto& f : desai_sibois_reference::fixtures )
    {
        BOOST_CHECK_EQUAL( none.ut1->getCorrections( f.tt ), 0.0 );
        BOOST_CHECK_EQUAL( none.polarMotion->getCorrections( f.tt ).norm( ), 0.0 );
    }
}

BOOST_AUTO_TEST_CASE( factory_opt_in_preserves_daily_eop_cpo_libration_and_native_ut1_rate )
{
    auto configured = settings( );
    BOOST_CHECK( configured->getOceanTideEopModel( ) == OceanTideEopModel::iers_2010 );
    const auto original = simulation_setup::createGcrsToItrsRotationModel( configured );
    configured->setOceanTideEopModel( OceanTideEopModel::desai_sibois_2016 );
    const auto replacement = std::dynamic_pointer_cast< ephemerides::GcrsToItrsRotationModel >(
            simulation_setup::createGcrsToItrsRotationModel( configured ) );
    const auto& f = desai_sibois_reference::fixtures[ 10 ];
    BOOST_CHECK_SMALL( f.tt - 536500869.184, 1e-8 );
    const Time tai = Time( 536500837.0 );
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( dailyData( ) );
    const auto ds = corrections( true );
    TerrestrialTimeScaleConverter converter( nullptr, ds.ut1, nullptr, evaluator );
    const auto result = converter.getUt1FromTaiWithRate( tai );
    const auto lib = libration( );
    const auto args = sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst( f.tt, f.tt, f.ut1 );
    BOOST_CHECK_SMALL( static_cast< double >( result.ut1 - ( Time( 536500800.0 ) - .4 ) ) - f.values[ 2 ] * 1e-6 -
                               lib.ut1->getCorrectionsFromFundamentalArgument( args ),
                       2e-10 );
    const double expectedRate =
            .99999998 + f.rates[ 2 ] * 1e-6 + calculateShortPeriodCorrectionRate( lib.ut1, Time( f.tt ), Time( f.ut1 ), .99999998 );
    BOOST_CHECK_SMALL( result.dUt1DTai - expectedRate, 4e-16 );

    // Independent SOFA CIO assembly with ocean fixture, retained libration,
    // daily PM and CPO; UT1-corrected phase is used for polar motion.
    const auto polarArgs = sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst( f.tt, f.tt, static_cast< double >( result.ut1 ) );

    // Ocean terms, validated against the independent fixture in the first test case.
    const Eigen::Vector2d ocean = ds.polarMotion->getCorrectionsFromFundamentalArgument( polarArgs ) -
            lib.polarMotion->getCorrectionsFromFundamentalArgument( polarArgs );
    const Eigen::Vector2d pm = Eigen::Vector2d( 1e-6, 2e-6 ) + ocean + lib.polarMotion->getCorrectionsFromFundamentalArgument( polarArgs );
    double x, y, s, q[ 3 ][ 3 ], w[ 3 ][ 3 ], reference[ 3 ][ 3 ];
    iauXys06a( 2451545.0, f.tt / 86400.0, &x, &y, &s );
    iauC2ixys( x + 1e-9, y - 2e-9, s, q );
    iauPom00( pm.x( ), pm.y( ), iauSp00( 2451545.0, f.tt / 86400.0 ), w );
    iauC2tcio( q,
               iauEra00( 2451545.0 + result.ut1.fullDaysSinceEpoch( ), static_cast< double >( result.ut1.secondsSinceNoon( ) / 86400.0L ) ),
               w,
               reference );
    const Eigen::Matrix3d actual = replacement->getRotationToBaseFrameFromExtendedTime( tai + 32.184 ).toRotationMatrix( );
    for( int i = 0; i < 3; ++i )
        for( int j = 0; j < 3; ++j ) BOOST_CHECK_SMALL( actual( i, j ) - reference[ j ][ i ], 6e-14 );
    BOOST_CHECK_GT( ( original->getRotationToBaseFrame( f.tt ).toRotationMatrix( ) - actual ).norm( ), 1e-12 );
    const Eigen::Matrix3d numerical = ( 8.0 *
                                                ( replacement->getRotationToBaseFrame( f.tt + 32.0 ).toRotationMatrix( ) -
                                                  replacement->getRotationToBaseFrame( f.tt - 32.0 ).toRotationMatrix( ) ) -
                                        ( replacement->getRotationToBaseFrame( f.tt + 64.0 ).toRotationMatrix( ) -
                                          replacement->getRotationToBaseFrame( f.tt - 64.0 ).toRotationMatrix( ) ) ) /
            384.0;
    BOOST_CHECK_SMALL( ( replacement->getDerivativeOfRotationToBaseFrameFromExtendedTime( Time( f.tt ) ) - numerical ).norm( ), 1e-13 );
}

BOOST_AUTO_TEST_CASE( custom_short_period_files_and_units_cannot_double_count_ocean_terms )
{
    auto configured = settings( );
    configured->setOceanTideEopModel( OceanTideEopModel::desai_sibois_2016 );
    configured->getPolarMotionCorrectionSettings( )->amplitudesFiles_.push_back( "custom-ocean.txt" );
    BOOST_CHECK_THROW( simulation_setup::createGcrsToItrsRotationModel( configured ), std::runtime_error );
    configured = settings( );
    configured->setOceanTideEopModel( OceanTideEopModel::desai_sibois_2016 );
    configured->getUt1CorrectionSettings( )->conversionFactor_ = 1.0;
    BOOST_CHECK_THROW( simulation_setup::createGcrsToItrsRotationModel( configured ), std::runtime_error );
    BOOST_CHECK_THROW( configured->setOceanTideEopModel( static_cast< OceanTideEopModel >( 42 ) ), std::runtime_error );
}
