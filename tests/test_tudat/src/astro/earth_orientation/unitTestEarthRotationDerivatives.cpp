#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include <limits>
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/astro/earth_orientation/shortPeriodEarthOrientationCorrectionCalculator.h"

using namespace tudat;
using namespace tudat::earth_orientation;
using namespace tudat::basic_astrodynamics;

namespace
{
input_output::EarthOrientationData rateData( const double ut1Slope = -2e-8 )
{
    input_output::EarthOrientationData data;
    data.source = "controlled affine TAI EOP rates";
    const Time reference = Time( ( 57754.0 - 51544.5 ) * 86400.0 ) + 37.0;
    for( int i = 0; i < 4; ++i )
    {
        input_output::EarthOrientationRecord record;
        record.utcModifiedJulianDate = 57752.0 + i;
        const Time utc( ( record.utcModifiedJulianDate - 51544.5 ) * 86400.0 );
        const Time tai = sofa_interface::convertUTCtoTAI( utc );
        const double dt = static_cast< double >( tai - reference );
        record.polarMotion = { 1e-6 + 2e-12 * dt, 2e-6 - 3e-12 * dt };
        record.ut1MinusUtc = -36.4 + ut1Slope * dt + static_cast< double >( tai - utc );
        record.celestialPoleOffsets = std::array< double, 2 >{ 1e-9 + 4e-14 * dt, -2e-9 - 5e-14 * dt };
        record.lengthOfDay = 0.004;  // Metadata intentionally differs from the selected UT1 interpolant.
        data.records.emplace( record.utcModifiedJulianDate, record );
    }
    return data;
}
}  // namespace

BOOST_AUTO_TEST_CASE( native_ut1_rate_uses_selected_interpolant_without_independent_lod )
{
    auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( rateData( ) );
    TerrestrialTimeScaleConverter converter( evaluator, nullptr, nullptr );
    const Time midnight = Time( ( 57754.0 - 51544.5 ) * 86400.0 ) + 37.0;
    for( const Time& tai : { evaluator->getStartTai( ), midnight - 0.5, midnight, evaluator->getEndTai( ) } )
    {
        const auto result = converter.getUt1FromTaiWithRate( tai );
        BOOST_CHECK_EQUAL( result.ut1, converter.getUt1FromTai( tai ) );
        BOOST_CHECK_SMALL( result.dUt1DTai - ( 1.0 - 2e-8 ), 2e-16 );
    }
}

BOOST_AUTO_TEST_CASE( native_ut1_rate_matches_complete_position_sequence_through_leap )
{
    auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( rateData( ) );
    TerrestrialTimeScaleConverter converter( evaluator, getDefaultUT1CorrectionCalculator( ), nullptr );
    const Time tai = Time( ( 57754.0 - 51544.5 ) * 86400.0 ) + 36.5;
    const auto result = converter.getUt1FromTaiWithRate( tai );
    BOOST_CHECK_EQUAL( result.ut1, converter.getUt1FromTai( tai ) );
    const Time a = converter.getUt1FromTai( tai + 64.0 ) - converter.getUt1FromTai( tai - 64.0 );
    const Time b = converter.getUt1FromTai( tai + 128.0 ) - converter.getUt1FromTai( tai - 128.0 );
    const double independent = ( 8.0 * static_cast< double >( a ) - static_cast< double >( b ) ) / ( 12.0 * 64.0 );
    BOOST_CHECK_SMALL( result.dUt1DTai - independent, 4e-14 );
    BOOST_CHECK_GT( std::abs( result.dUt1DTai - ( 1.0 - 2e-8 ) ), 1e-11 );
}

BOOST_AUTO_TEST_CASE( ut1_rate_unavailability_retains_covered_and_historical_positions )
{
    auto data = rateData( );
    data.records.erase( std::next( data.records.begin( ) ), std::prev( data.records.end( ) ) );
    auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( data );
    TerrestrialTimeScaleConverter converter( evaluator, nullptr, nullptr );
    BOOST_CHECK_NO_THROW( converter.getUt1FromTai( evaluator->getStartTai( ) ) );
    BOOST_CHECK_THROW( converter.getUt1FromTaiWithRate( evaluator->getStartTai( ) ), std::runtime_error );
    EarthOrientationEvaluationSettings settings;
    settings.historicalUt1 = HistoricalUt1Policy::delta_t;
    auto historicalData = std::make_shared< EarthOrientationDataEvaluator >( rateData( ), settings );
    TerrestrialTimeScaleConverter historical( historicalData, nullptr, nullptr );
    const Time tt( ( 1893.0 - 2000.0 ) * physical_constants::JULIAN_YEAR );
    BOOST_CHECK_NO_THROW( historical.getCurrentTime< Time >( tt_scale, ut1_scale, tt ) );
    BOOST_CHECK_THROW( historical.getUt1FromTaiWithRate( convertTTtoTAI( tt ) ), std::runtime_error );
    const Time transition( ( JULIAN_DAY_OF_EOP_INTRODUCTION - JULIAN_DAY_ON_J2000 ) * physical_constants::JULIAN_DAY );
    const Time before = transition - 1e-10;
    BOOST_CHECK_EQUAL( static_cast< double >( before ), static_cast< double >( transition ) );
    BOOST_CHECK_NO_THROW( historical.getUt1FromTai( before ) );
    BOOST_CHECK_NO_THROW( historical.getCurrentTime< Time >( ut1_scale, tai_scale, before ) );
    BOOST_CHECK_EXCEPTION( historical.getUt1FromTaiWithRate( before ), std::runtime_error, []( const std::runtime_error& error ) {
        return std::string( error.what( ) ).find( "historical Delta T" ) != std::string::npos;
    } );
    BOOST_CHECK_THROW( historical.getUt1FromTai( transition ), std::runtime_error );
    TerrestrialTimeScaleConverter withoutEop( nullptr, nullptr );
    const Time nearby = transition - 1e-6;
    for( const auto scales : { std::make_pair( utc_scale, tt_scale ), std::make_pair( tai_scale, utc_scale ) } )
    {
        const Time reference = withoutEop.getCurrentTime< Time >( scales.first, scales.second, nearby );
        const Time current = withoutEop.getCurrentTime< Time >( scales.first, scales.second, before );
        BOOST_CHECK_SMALL( static_cast< double >( ( current - reference ) - ( before - nearby ) ), 1e-8 );
    }
}

#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"

namespace
{
Time referenceTai( )
{
    return Time( ( 57754.0 - 51544.5 ) * 86400.0 ) + 37.0;
}

std::shared_ptr< EarthOrientationAnglesCalculator > rateCalculator( const input_output::EarthOrientationData& data = rateData( ),
                                                                    const IAUConventions theory = iau_2006,
                                                                    const bool shortPeriods = false,
                                                                    const EarthOrientationEvaluationSettings& settings = {} )
{
    auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( data, settings );
    auto ut1 =
            std::make_shared< TerrestrialTimeScaleConverter >( evaluator, shortPeriods ? getDefaultUT1CorrectionCalculator( ) : nullptr );
    auto polar =
            std::make_shared< PolarMotionCalculator >( evaluator, shortPeriods ? getDefaultPolarMotionCorrectionCalculator( ) : nullptr );
    auto celestial = std::make_shared< PrecessionNutationCalculator >( theory, evaluator, nullptr );
    return std::make_shared< EarthOrientationAnglesCalculator >( polar, celestial, ut1 );
}

Eigen::Matrix3d matrixFromSofa( const double matrix[ 3 ][ 3 ] )
{
    Eigen::Matrix3d result;
    for( int i = 0; i < 3; ++i )
    {
        for( int j = 0; j < 3; ++j )
        {
            result( i, j ) = matrix[ i ][ j ];
        }
    }
    return result;
}

// Independent SOFA CIO recipe, using controlled EOP values linear in TAI, without
// production orientation angles, EOP interpolation or quaternion composition.
// https://www.iausofa.org/s/sofa_pn_c.pdf, sections 5.5/5.6.
Eigen::Matrix3d independentSofaRotation( const Time& tai, const IAUConventions theory = iau_2006, const double ut1Slope = -2e-8 )
{
    const double dt = static_cast< double >( tai - referenceTai( ) );
    const Time tt = tai + 32.184;
    const Time ut1 = tai - 36.4 + ut1Slope * dt;
    double x, y, s, q[ 3 ][ 3 ], w[ 3 ][ 3 ], result[ 3 ][ 3 ];
    const double tt1 = 2451545.0 + tt.fullDaysSinceEpoch( );
    const double tt2 = static_cast< double >( tt.secondsSinceNoon( ) / 86400.0L );
    if( theory == iau_2006 )
    {
        iauXys06a( tt1, tt2, &x, &y, &s );
    }
    else if( theory == iau_2000_a )
    {
        iauXys00a( tt1, tt2, &x, &y, &s );
    }
    else
    {
        iauXys00b( tt1, tt2, &x, &y, &s );
    }
    iauC2ixys( x + 1e-9 + 4e-14 * dt, y - 2e-9 - 5e-14 * dt, s, q );
    iauPom00( 1e-6 + 2e-12 * dt, 2e-6 - 3e-12 * dt, iauSp00( tt1, tt2 ), w );
    iauC2tcio(
            q, iauEra00( 2451545.0 + ut1.fullDaysSinceEpoch( ), static_cast< double >( ut1.secondsSinceNoon( ) / 86400.0L ) ), w, result );
    return matrixFromSofa( result ).transpose( );
}

Eigen::Matrix3d fivePointMatrixRate( const std::function< Eigen::Matrix3d( const Time& ) >& position, const Time& time, const double step )
{
    return ( 8.0 * ( position( time + step ) - position( time - step ) ) -
             ( position( time + 2.0 * step ) - position( time - 2.0 * step ) ) ) /
            ( 12.0 * step );
}

std::map< std::string, Eigen::Matrix3d > independentBaseFrames( const IAUConventions theory )
{
    double bias[ 3 ][ 3 ], precession[ 3 ][ 3 ], combined[ 3 ][ 3 ];
    if( theory == iau_2006 )
    {
        iauBp06( 2451545.0, 0.0, bias, precession, combined );
    }
    else
    {
        iauBp00( 2451545.0, 0.0, bias, precession, combined );
    }

    // NAIF built-in ECLIPJ2000 obliquity 84381.448 arcsec, explicit independent constants.
    Eigen::Matrix3d ecliptic;
    ecliptic << 1.0, 0.0, 0.0, 0.0, 0.9174820620691818, 0.3977771559319137, 0.0, -0.3977771559319137, 0.9174820620691818;

    // J2000 (and its synonym EME2000) includes the frame bias; ECLIPJ2000 is defined with respect to J2000.
    return { { "GCRS", Eigen::Matrix3d::Identity( ) },
             { "ICRF", Eigen::Matrix3d::Identity( ) },
             { "J2000", matrixFromSofa( bias ) },
             { "EME2000", matrixFromSofa( bias ) },
             { "ECLIPJ2000", ecliptic * matrixFromSofa( bias ) } };
}
}  // namespace

BOOST_AUTO_TEST_CASE( complete_rotation_matches_independent_sofa_for_conventions_and_base_frames )
{
    const Time tai = referenceTai( ) + 12345.678901234L;
    for( const auto theory : { iau_2000_a, iau_2000_b, iau_2006 } )
    {
        auto calculator = rateCalculator( rateData( ), theory );
        const auto reference = independentSofaRotation( tai, theory );
        const auto expectedRate = fivePointMatrixRate( [ = ]( const Time& t ) { return independentSofaRotation( t, theory ); }, tai, 32.0 );
        for( const auto& base : independentBaseFrames( theory ) )
        {
            ephemerides::GcrsToItrsRotationModel model( calculator, tai_scale, base.first );
            const auto quaternion = model.getRotationToBaseFrameFromExtendedTime( tai );
            const Eigen::Matrix3d r = quaternion.toRotationMatrix( );
            const Eigen::Matrix3d dr = model.getDerivativeOfRotationToBaseFrameFromExtendedTime( tai );
            BOOST_CHECK_SMALL( ( r - base.second * reference ).norm( ), 3e-14 );
            BOOST_CHECK_SMALL( ( dr - base.second * expectedRate ).norm( ), 8e-15 );
            BOOST_CHECK_SMALL( quaternion.norm( ) - 1.0, 4e-16 );
            BOOST_CHECK_SMALL( ( dr * r.transpose( ) + r * dr.transpose( ) ).norm( ), 2e-19 );
            BOOST_CHECK_SMALL( ( model.getDerivativeOfRotationToTargetFrameFromExtendedTime( tai ) - dr.transpose( ) ).norm( ), 1e-25 );
            const double input = static_cast< double >( tai );
            BOOST_CHECK_SMALL( ( model.getDerivativeOfRotationToBaseFrame( input ) -
                                 model.getDerivativeOfRotationToBaseFrameFromExtendedTime( Time( input ) ) )
                                       .norm( ),
                               1e-25 );
            BOOST_CHECK_SMALL( ( model.getRotationToBaseFrame( input ).toRotationMatrix( ) -
                                 model.getRotationToBaseFrameFromExtendedTime( Time( input ) ).toRotationMatrix( ) )
                                       .norm( ),
                               1e-25 );
        }
    }
}

BOOST_AUTO_TEST_CASE( exact_celestial_matrix_remains_unit_across_supported_era )
{
    for( const auto theory : { iau_2000_a, iau_2000_b, iau_2006 } )
    {
        for( const double years : { -50.0, 0.0, 25.0, 100.0, 500.0 } )
        {
            const auto xys = sofa_interface::getPositionOfCipInGcrs( years * 365.25 * 86400.0, 2451545.0, theory );
            double q[ 3 ][ 3 ];
            iauC2ixys( xys.x( ), xys.y( ), xys.z( ), q );
            const auto result = calculateRotationFromCirsToGcrs( xys.x( ), xys.y( ), xys.z( ) );
            BOOST_CHECK_SMALL( ( result.toRotationMatrix( ) - matrixFromSofa( q ).transpose( ) ).norm( ), 6e-16 );
            BOOST_CHECK_SMALL( result.norm( ) - 1.0, 4e-16 );
        }
    }
}

BOOST_AUTO_TEST_CASE( full_sofa_derivative_convergence_and_station_velocity_closure )
{
    auto calculator = rateCalculator( );
    ephemerides::GcrsToItrsRotationModel model( calculator, tai_scale );
    const Time tai = referenceTai( ) + 1234.5678;
    const Eigen::Vector3d station( 4075578.385, 931852.890, 4801570.154 );
    const Eigen::Matrix3d derivative = model.getDerivativeOfRotationToBaseFrameFromExtendedTime( tai );
    std::vector< double > errors;
    for( const double step : { 512.0, 128.0, 32.0, 8.0 } )
    {
        const auto numerical = fivePointMatrixRate( []( const Time& t ) { return independentSofaRotation( t ); }, tai, step );
        errors.push_back( ( ( numerical - derivative ) * station ).norm( ) );
        BOOST_TEST_MESSAGE( "Independent complete SOFA station derivative h=" << step << " s error=" << errors.back( ) << " m/s" );
    }
    BOOST_CHECK_LT( errors[ 1 ], errors[ 0 ] / 100.0 );
    BOOST_CHECK_LT( errors[ 2 ], 3e-8 );
    BOOST_CHECK_LT( errors[ 3 ], 1e-7 );
    const auto closure = fivePointMatrixRate(
            [ &model ]( const Time& t ) { return model.getRotationToBaseFrameFromExtendedTime( t ).toRotationMatrix( ); }, tai, 32.0 );
    BOOST_CHECK_SMALL( ( ( closure - derivative ) * station ).norm( ), 1e-7 );
}

BOOST_AUTO_TEST_CASE( affine_ut1_polar_and_celestial_slopes_contribute_without_phase_offset )
{
    const Time tai = referenceTai( );
    auto moving = rateCalculator( rateData( 4e-5 ) );
    auto staticData = rateData( 0.0 );
    for( auto& record : staticData.records )
    {
        record.second.polarMotion = { 1e-6, 2e-6 };
        record.second.celestialPoleOffsets = std::array< double, 2 >{ 1e-9, -2e-9 };
    }
    auto stationary = rateCalculator( staticData );
    const auto values = moving->getRotationAnglesAndRatesFromItrsToGcrs( tai, tai_scale );
    const auto fixed = stationary->getRotationAnglesAndRatesFromItrsToGcrs( tai, tai_scale );
    BOOST_CHECK_SMALL( ( values.angles - fixed.angles ).norm( ), 1e-20 );
    BOOST_CHECK_EQUAL( values.ut1, fixed.ut1 );
    BOOST_CHECK_SMALL( values.dUt1DInput - 1.00004, 2e-16 );
    BOOST_CHECK_SMALL( values.angleRates[ 3 ] - 2e-12, 1e-26 );
    BOOST_CHECK_SMALL( values.angleRates[ 4 ] + 3e-12, 1e-26 );
    BOOST_CHECK_SMALL( values.angleRates[ 0 ] - fixed.angleRates[ 0 ] - 4e-14, 1e-26 );
    BOOST_CHECK_SMALL( values.angleRates[ 1 ] - fixed.angleRates[ 1 ] + 5e-14, 1e-26 );
    auto onlySpin = fixed;
    onlySpin.dUt1DInput = values.dUt1DInput;
    const Eigen::Matrix3d change = calculateRotationRateFromItrsToGcrs( onlySpin ) - calculateRotationRateFromItrsToGcrs( fixed );
    BOOST_CHECK_SMALL( change.norm( ) - std::sqrt( 2.0 ) * 2.0 * mathematical_constants::PI * 1.00273781191135448 / 86400.0 * 4e-5, 5e-20 );
    BOOST_CHECK_GT( ( calculateRotationRateFromItrsToGcrs( values ) - calculateRotationRateFromItrsToGcrs( onlySpin ) ).norm( ), 1e-12 );
}

BOOST_AUTO_TEST_CASE( continuous_rates_cross_leap_seconds_and_share_clock_scaling )
{
    auto calculator = rateCalculator( rateData( ), iau_2006, true );
    auto converter = calculator->getTerrestrialTimeScaleConverter( );
    ephemerides::GcrsToItrsRotationModel taiModel( calculator, tai_scale );
    for( const double delta : { -1.5, -0.5, 0.0, 0.5, 12345.0 } )
    {
        const Time tai = referenceTai( ) + delta;
        const Eigen::Matrix3d r = taiModel.getRotationToBaseFrameFromExtendedTime( tai ).toRotationMatrix( );
        const Eigen::Matrix3d dr = taiModel.getDerivativeOfRotationToBaseFrameFromExtendedTime( tai );
        const auto fd = fivePointMatrixRate(
                [ &taiModel ]( const Time& t ) { return taiModel.getRotationToBaseFrameFromExtendedTime( t ).toRotationMatrix( ); },
                tai,
                32.0 );
        BOOST_CHECK_SMALL( ( dr - fd ).norm( ), 1e-14 );
        for( const auto scale : { tt_scale, tdb_scale, ut1_scale } )
        {
            const Time input = converter->getCurrentTime< Time >( tai_scale, scale, tai );
            ephemerides::GcrsToItrsRotationModel model( calculator, scale );
            const auto state = calculator->getRotationAnglesAndRatesFromItrsToGcrs( input, scale );
            BOOST_CHECK_SMALL( ( model.getRotationToBaseFrameFromExtendedTime( input ).toRotationMatrix( ) - r ).norm( ), 1e-14 );
            BOOST_CHECK_SMALL( ( model.getDerivativeOfRotationToBaseFrameFromExtendedTime( input ) - dr * state.dTaiDInput ).norm( ),
                               3e-19 );
            const auto inputDerivative = fivePointMatrixRate(
                    [ &model ]( const Time& time ) { return model.getRotationToBaseFrameFromExtendedTime( time ).toRotationMatrix( ); },
                    input,
                    32.0 );
            BOOST_CHECK_SMALL( ( model.getDerivativeOfRotationToBaseFrameFromExtendedTime( input ) - inputDerivative ).norm( ), 1e-14 );
            if( scale == ut1_scale )
            {
                BOOST_CHECK_SMALL( state.dUt1DInput - 1.0, 2e-16 );
            }

            BOOST_CHECK_SMALL( ( model.getRotationBetweenIntermediateFrames( ephemerides::EarthOrientationIntermediateFrame::itrs,
                                                                             ephemerides::EarthOrientationIntermediateFrame::icrs,
                                                                             static_cast< double >( input ) )
                                         .toRotationMatrix( ) -
                                 model.getRotationToBaseFrame( static_cast< double >( input ) ).toRotationMatrix( ) )
                                       .norm( ),
                               1e-15 );
        }
    }
    const Time utc = referenceTai( ) - 37.0;
    ephemerides::GcrsToItrsRotationModel utcModel( calculator, utc_scale );
    BOOST_CHECK_NO_THROW( utcModel.getRotationToBaseFrameFromExtendedTime( utc ) );
    BOOST_CHECK_THROW( utcModel.getDerivativeOfRotationToBaseFrameFromExtendedTime( utc ), std::runtime_error );
    BOOST_CHECK_NO_THROW( utcModel.getDerivativeOfRotationToBaseFrameFromExtendedTime( utc - 0.1 ) );
    BOOST_CHECK_NO_THROW( utcModel.getDerivativeOfRotationToBaseFrameFromExtendedTime( utc + 0.1 ) );
}

BOOST_AUTO_TEST_CASE( pre_1972_utc_drift_scales_the_full_rotation_derivative )
{
    auto data = rateData( );
    data.records.clear( );
    for( int i = 0; i < 4; ++i )
    {
        input_output::EarthOrientationRecord record;
        record.utcModifiedJulianDate = 39000.0 + i;  // 1965, SOFA drift 0.001296 seconds per UTC day.
        record.polarMotion = { 1e-6, 2e-6 };
        record.celestialPoleOffsets = std::array< double, 2 >{ 1e-9, -2e-9 };
        const Time utc( ( record.utcModifiedJulianDate - 51544.5 ) * 86400.0 );
        record.ut1MinusUtc = -4.0 + static_cast< double >( sofa_interface::convertUTCtoTAI( utc ) - utc );
        data.records.emplace( record.utcModifiedJulianDate, record );
    }
    auto calculator = rateCalculator( data );
    const Time utc( ( 39001.25 - 51544.5 ) * 86400.0 );
    const Time tai = sofa_interface::convertUTCtoTAI( utc );
    const auto utcRate = calculator->getRotationAnglesAndRatesFromItrsToGcrs( utc, utc_scale );
    const auto taiRate = calculator->getRotationAnglesAndRatesFromItrsToGcrs( tai, tai_scale );
    BOOST_CHECK_SMALL( utcRate.dTaiDInput - ( 1.0 + 0.001296 / 86400.0 ), 2e-16 );
    BOOST_CHECK_SMALL( ( calculateRotationRateFromItrsToGcrs( utcRate ) -
                         ( 1.0 + 0.001296 / 86400.0 ) * calculateRotationRateFromItrsToGcrs( taiRate ) )
                               .norm( ),
                       3e-20 );
}

BOOST_AUTO_TEST_CASE( strict_rates_retain_endpoints_knot_policy_and_missing_rate_failures )
{
    auto calculator = rateCalculator( );
    const auto converter = calculator->getTerrestrialTimeScaleConverter( );
    const auto evaluator = converter->getEarthOrientationDataEvaluator( );
    for( const Time& tai : { evaluator->getStartTai( ), evaluator->getEndTai( ) } )
    {
        for( const auto scale : { tai_scale, tt_scale, tdb_scale, ut1_scale } )
        {
            const Time input = converter->getCurrentTime< Time >( tai_scale, scale, tai );
            const auto rate = calculator->getRotationAnglesAndRatesFromItrsToGcrs( input, scale );
            const auto doubleRate = calculator->getRotationAnglesAndRatesFromItrsToGcrs( static_cast< double >( input ), scale );
            if( std::numeric_limits< long double >::digits > std::numeric_limits< double >::digits )
            {
                BOOST_CHECK_EQUAL( rate.tai, tai );
            }
            else
            {
                // Four fractional-hour roundings bound the native clock round trip on binary64 long-double platforms.
                const long double epochTolerance = 2.0L * std::numeric_limits< long double >::epsilon( ) * TIME_NORMALIZATION_INTEGER_TERM;

                BOOST_CHECK_SMALL( static_cast< long double >( rate.tai - tai ), epochTolerance );
            }

            // A rounded double end point lies within the EOP coverage tolerance of the exact end point.
            BOOST_CHECK_SMALL( static_cast< double >( doubleRate.tai - tai ), 1e-6 );
            BOOST_CHECK_SMALL( ( rate.angleRates - doubleRate.angleRates ).norm( ), 1e-22 );
        }
    }
    auto gap = rateData( );
    gap.records.erase( std::next( gap.records.begin( ) ), std::prev( gap.records.end( ) ) );
    auto gapCalculator = rateCalculator( gap );
    const Time start = gapCalculator->getTerrestrialTimeScaleConverter( )->getEarthOrientationDataEvaluator( )->getStartTai( );
    BOOST_CHECK_NO_THROW( gapCalculator->getRotationAnglesAndTimesFromItrsToGcrs( start, tai_scale ) );
    BOOST_CHECK_THROW( gapCalculator->getRotationAnglesAndRatesFromItrsToGcrs( start, tai_scale ), std::runtime_error );
    auto missing = rateData( );
    std::next( missing.records.begin( ) )->second.celestialPoleOffsets.reset( );
    auto missingCalculator = rateCalculator( missing );
    BOOST_CHECK_THROW( missingCalculator->getRotationAnglesAndRatesFromItrsToGcrs( start, tai_scale ), std::runtime_error );
    auto knotData = rateData( );
    std::prev( knotData.records.end( ) )->second.polarMotion[ 0 ] += 1e-6;
    auto knotCalculator = rateCalculator( knotData );
    const auto at = knotCalculator->getRotationAnglesAndRatesFromItrsToGcrs( referenceTai( ), tai_scale );
    const auto before = knotCalculator->getRotationAnglesAndRatesFromItrsToGcrs( referenceTai( ) - 0.01, tai_scale );
    BOOST_CHECK_SMALL( at.angleRates[ 3 ] - ( 2e-12 + 1e-6 / 86400.0 ), 1e-25 );

    // Fractional-hour clock roundoff on either side of the knot must retain the right-hand slope.
    const long double roundoff = std::numeric_limits< long double >::epsilon( ) * TIME_NORMALIZATION_INTEGER_TERM;

    for( const long double shift : { -roundoff, roundoff } )
    {
        const auto nearby = knotCalculator->getRotationAnglesAndRatesFromItrsToGcrs( referenceTai( ) + shift, tai_scale );

        BOOST_CHECK_EQUAL( nearby.angleRates[ 3 ], at.angleRates[ 3 ] );
    }
    BOOST_CHECK_SMALL( before.angleRates[ 3 ] - 2e-12, 1e-25 );
}

BOOST_AUTO_TEST_CASE( nominal_cip_rate_matches_independent_sofa_and_ignores_position_interpolation )
{
    const Time tt = referenceTai( ) + 32.184 + 12345.0;
    const auto eopData = std::make_shared< EarthOrientationDataEvaluator >( rateData( ) );
    PrecessionNutationCalculator direct( iau_2006, eopData );

    // Unit-scale rotation-matrix arithmetic is allowed eight epsilon of absolute
    // error per returned coordinate. The five-point weights amplify that by
    // 3/(2h); add the two stencils' budgets, in three dimensions. Truncation
    // allowance: amplitude 1e-3 rad and the shortest nominal CIP period 2 days,
    // giving amplitude*(2*pi/period)^5*h^4/30 for a five-point derivative.
    const auto budget = []( const double h1, const double h2 ) {
        const double rounding = 1.5 * 8.0 * std::numeric_limits< double >::epsilon( ) * ( 1.0 / h1 + 1.0 / h2 );
        const double frequency = 2.0 * mathematical_constants::PI / ( 2.0 * 86400.0 );
        const double truncation = 1e-3 * std::pow( frequency, 5 ) * ( std::pow( h1, 4 ) + std::pow( h2, 4 ) ) / 30.0;
        return std::sqrt( 3.0 ) * ( rounding + truncation );
    };
    for( const double year : { 1962.0, 2000.0, 2017.0, 2026.0 } )
    {
        const Time sample( ( year - 2000.0 ) * 365.25 * 86400.0 + 12345.6789 );
        const auto nominal = direct.getNominalPositionAndRate( sample );
        const auto sofaValue = []( const Time& epoch ) -> Eigen::Vector3d {
            double x, y, s;
            iauXys06a( 2451545.0 + epoch.fullDaysSinceEpoch( ), static_cast< double >( epoch.secondsSinceNoon( ) / 86400.0L ), &x, &y, &s );
            return Eigen::Vector3d( x, y, s );
        };
        const double h = 64.0;
        const Eigen::Vector3d independent = ( 8.0 * ( sofaValue( sample + h ) - sofaValue( sample - h ) ) -
                                              ( sofaValue( sample + 2 * h ) - sofaValue( sample - 2 * h ) ) ) /
                ( 12.0 * h );
        BOOST_CHECK_SMALL( ( nominal.second - independent ).norm( ), budget( 32.0, 64.0 ) );
    }

    // Interpolated nominal positions retain the rate of the direct model.
    const double initial = static_cast< double >( tt ) - 20000.0;
    const double final = initial + 40000.0;
    for( const auto& settings :
         { interpolators::linearInterpolation( ), interpolators::cubicSplineInterpolation( ), interpolators::lagrangeInterpolation( 6 ) } )
    {
        PrecessionNutationCalculator interpolated(
                iau_2006, eopData, interpolators::interpolatorGenerationSettings( settings, initial, final, 2000.0 ) );
        for( const double epoch : { initial + 1000.0, static_cast< double >( tt ) } )
        {
            const auto result = interpolated.getNominalPositionAndRate( Time( epoch ) );
            BOOST_CHECK_EQUAL( ( result.second - direct.getNominalPositionAndRate( Time( epoch ) ).second ).norm( ), 0.0 );
        }
    }
}

BOOST_AUTO_TEST_CASE( short_period_tangent_rates_converge_without_probing_daily_eop )
{
    const Time tt = referenceTai( ) + 32.184;
    const Time ut1 = referenceTai( ) - 36.4;
    auto polar = getDefaultPolarMotionCorrectionCalculator( );
    auto clock = getDefaultUT1CorrectionCalculator( );
    const auto polar32 = polar->getCorrectionRate( tt, ut1, 1.0002, 32.0 );
    const double clock32 = clock->getCorrectionRate( tt, ut1, 1.0002, 32.0 );
    for( const double h : { 16.0, 64.0, 128.0 } )
    {
        const auto polarOther = polar->getCorrectionRate( tt, ut1, 1.0002, h );
        const double clockOther = clock->getCorrectionRate( tt, ut1, 1.0002, h );
        BOOST_TEST_MESSAGE( "Short-period h=" << h << " PM tangent discrepancy=" << ( polarOther - polar32 ).norm( )
                                              << " rad/s, UT1 tangent discrepancy=" << std::abs( clockOther - clock32 ) << " s/s" );
        BOOST_CHECK_SMALL( ( polarOther - polar32 ).norm( ), 2e-20 );
        BOOST_CHECK_SMALL( clockOther - clock32, 5e-16 );
    }
}
