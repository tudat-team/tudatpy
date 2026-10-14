#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

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
        record.celestialPoleOffsets = std::array< double, 2 >{ 1e-9, -2e-9 };
        record.lengthOfDay = 0.004;  // Metadata intentionally differs from the selected UT1 interpolant.
        data.records.emplace( record.utcModifiedJulianDate, record );
    }
    return data;
}
}  // namespace

BOOST_AUTO_TEST_CASE( native_ut1_rate_uses_selected_interpolant_without_independent_lod )
{
    auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( rateData( ) );
    TerrestrialTimeScaleConverter converter( nullptr, nullptr, nullptr, evaluator );
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
    TerrestrialTimeScaleConverter converter( nullptr, getDefaultUT1CorrectionCalculator( ), nullptr, evaluator );
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
    TerrestrialTimeScaleConverter converter( nullptr, nullptr, nullptr, evaluator );
    BOOST_CHECK_NO_THROW( converter.getUt1FromTai( evaluator->getStartTai( ) ) );
    BOOST_CHECK_THROW( converter.getUt1FromTaiWithRate( evaluator->getStartTai( ) ), std::runtime_error );
    EarthOrientationEvaluationSettings settings;
    settings.historicalUt1 = HistoricalUt1Policy::delta_t;
    auto historicalData = std::make_shared< EarthOrientationDataEvaluator >( rateData( ), settings );
    TerrestrialTimeScaleConverter historical( nullptr, nullptr, nullptr, historicalData );
    const Time tt( ( 1893.0 - 2000.0 ) * physical_constants::JULIAN_YEAR );
    BOOST_CHECK_NO_THROW( historical.getCurrentTime< Time >( tt_scale, ut1_scale, tt ) );
    BOOST_CHECK_THROW( historical.getUt1FromTaiWithRate( convertTTtoTAI( tt ) ), std::runtime_error );
}
