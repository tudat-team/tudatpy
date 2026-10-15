/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "tudat/math/basic/numericalDerivative.h"

#include <limits>

namespace tudat
{

namespace earth_orientation
{

// Load the existing historical Delta-T values only when the fallback is needed.
void TerrestrialTimeScaleConverter::loadHistoricalDeltaTInterpolator( )
{
    if( historicalDeltaTInterpolator_ == nullptr )
    {
        const auto values = input_output::readFloatingPointMapFromFile< double, double >(
                paths::getEarthOrientationDataFilesPath( ) + "/historicalDeltaT.txt", ",", "#" );
        historicalDeltaTInterpolator_ =
                interpolators::createOneDimensionalInterpolator( values, interpolators::cubicSplineInterpolation( ) );
    }
}

// Throw an exception for an unknown time scale.
void TerrestrialTimeScaleConverter::checkTimeScale( const basic_astrodynamics::TimeScales timeScale )
{
    switch( timeScale )
    {
        case basic_astrodynamics::tai_scale:
        case basic_astrodynamics::tt_scale:
        case basic_astrodynamics::tdb_scale:
        case basic_astrodynamics::utc_scale:
        case basic_astrodynamics::ut1_scale:
            return;
        default:
            throw std::runtime_error( "Unknown time scale in time scale conversion" );
    }
}

// UT1 from historical Delta T (TT - UT1), evaluated at the year computed from TT.
Time TerrestrialTimeScaleConverter::getHistoricalUt1FromTai( const Time& tai )
{
    loadHistoricalDeltaTInterpolator( );
    const Time tt = basic_astrodynamics::convertTAItoTT( tai );
    return tt - historicalDeltaTInterpolator_->interpolate( static_cast< double >( tt ) / physical_constants::JULIAN_YEAR + 2000.0 );
}

// Inverse of getHistoricalUt1FromTai: four fixed-point updates of TT = UT1 + Delta T( TT ).
Time TerrestrialTimeScaleConverter::getTaiFromHistoricalUt1( const Time& ut1 )
{
    loadHistoricalDeltaTInterpolator( );
    Time tt = ut1;
    for( int i = 0; i < 4; ++i )
    {
        tt = ut1 + historicalDeltaTInterpolator_->interpolate( static_cast< double >( tt ) / physical_constants::JULIAN_YEAR + 2000.0 );
    }
    return basic_astrodynamics::convertTTtoTAI( tt );
}

// TDB - TT, interpolated if requested; at a point away from the geocentre, this requires UT1 when Earth orientation data is
// provided.
double TerrestrialTimeScaleConverter::getTdbMinusTt( const Time& tt, const Eigen::Vector3d& position )
{
    const auto key = std::make_tuple( position.x( ), position.y( ), position.z( ) );
    if( tdbToTtInterpolators_.count( key ) )
    {
        return tdbToTtInterpolators_.at( key )->interpolate( static_cast< double >( tt ) );
    }
    if( !eopData_ && !position.isZero( 0.0 ) )
    {
        return sofa_interface::getTDBminusTT( static_cast< double >( tt ), position );
    }

    // At the geocentre, the terms that depend on UT1 vanish, so TT is used instead of UT1. This avoids requiring Earth
    // orientation data for geocentric TDB.
    const Time ut1 = position.isZero( 0.0 ) ? tt : getUt1FromTai( basic_astrodynamics::convertTTtoTAI( tt ) );
    const double fraction = static_cast< double >( ut1.secondsSinceMidnight( ) / 86400.0L );
    return sofa_interface::getTDBminusTT( static_cast< double >( tt ), fraction, position );
}

// TT - TDB (seconds) at the geocentre, at a given TDB.
double TerrestrialTimeScaleConverter::getGeocentricTtMinusTdb( const Time& tdb )
{
    return static_cast< double >( basic_astrodynamics::convertTAItoTT( getTaiTime( basic_astrodynamics::tdb_scale, tdb ) ) - tdb );
}

// dTT/dTDB at the geocentre, from TT - TDB (interpolated, if selected) with a fourth-order central difference over 32 s.
double TerrestrialTimeScaleConverter::getGeocentricTtRateWithRespectToTdb( const Time& tdb )
{
    if( !std::isfinite( static_cast< double >( tdb ) ) )
    {
        throw std::invalid_argument( "Non-finite TDB clock-rate epoch" );
    }

    const std::function< double( const Time& ) > ttMinusTdb =
            std::bind( &TerrestrialTimeScaleConverter::getGeocentricTtMinusTdb, this, std::placeholders::_1 );
    const double rate = 1.0 + numerical_derivatives::computeFourthOrderCentralDifference( ttMinusTdb, tdb, 32.0 );

    if( !std::isfinite( rate ) || rate <= 0.0 )
    {
        throw std::runtime_error( "Invalid geocentric TT/TDB clock rate" );
    }
    return rate;
}

// UT1 without requiring its rate.
Time TerrestrialTimeScaleConverter::getUt1FromTai( const Time& tai )
{
    return evaluateUt1FromTai( tai, false ).ut1;
}

// UT1 and its derivative with respect to TAI.
Ut1ValueAndRate TerrestrialTimeScaleConverter::getUt1FromTaiWithRate( const Time& tai )
{
    return evaluateUt1FromTai( tai, true );
}

// UT1 from the daily data and short-period variations or, before 1962, from historical Delta T.
Ut1ValueAndRate TerrestrialTimeScaleConverter::evaluateUt1FromTai( const Time& tai, const bool requireRate )
{
    // Historical Delta T is used before 1962: always without Earth orientation data, and with these data if they start
    // later than the requested epoch and the settings select Delta T.
    const bool useHistoricalDeltaT = tai < Time( utcIntroductionEpochInTai_ ) &&
            ( !eopData_ || ( tai < eopData_->getStartTai( ) && eopData_->getSettings( ).historicalUt1 == HistoricalUt1Policy::delta_t ) );
    if( useHistoricalDeltaT )
    {
        if( requireRate )
        {
            throw std::runtime_error( "The rate of UT1 is not available from historical Delta T" );
        }
        return { getHistoricalUt1FromTai( tai ) };
    }
    if( !eopData_ )
    {
        throw std::runtime_error( "UT1 after 1962 requires Earth orientation data" );
    }

    // Daily UT1-TAI and, if available, its rate.
    const auto daily = eopData_->evaluate( tai );

    if( requireRate && !daily.ut1MinusTaiRate )
    {
        throw std::runtime_error( "Rate of the interpolated UT1-TAI is unavailable" );
    }

    // Fundamental arguments use TT and daily UT1 before adding any short-period correction.
    const Time tt = basic_astrodynamics::convertTAItoTT( tai );

    const Time baselineUt1 = tai + daily.ut1MinusTai;

    Ut1ValueAndRate result{ baselineUt1 };

    if( requireRate )
    {
        result.dUt1DTai = 1.0 + *daily.ut1MinusTaiRate;
    }
    if( shortPeriodUt1CorrectionCalculator_ )
    {
        // Short-period values and rates are evaluated at the daily UT1, without short-period variations.
        const auto arguments = sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst(
                static_cast< double >( tt ), static_cast< double >( tt ), static_cast< double >( baselineUt1 ) );

        result.ut1 += shortPeriodUt1CorrectionCalculator_->getCorrectionsFromFundamentalArgument( arguments );
        if( requireRate )
        {
            result.dUt1DTai += shortPeriodUt1CorrectionCalculator_->getCorrectionRate( tt, baselineUt1, result.dUt1DTai );
        }
    }
    if( requireRate && ( !std::isfinite( result.dUt1DTai ) || result.dUt1DTai <= 0.0 ) )
    {
        throw std::runtime_error( "Invalid UT1 rate" );
    }
    return result;
}

// Convert an epoch in the given time scale to TAI, without reducing Time precision.
Time TerrestrialTimeScaleConverter::getTaiTime( const basic_astrodynamics::TimeScales inputScale,
                                                const Time& input,
                                                const Eigen::Vector3d& position )
{
    switch( inputScale )
    {
        case basic_astrodynamics::tai_scale:
            return input;
        case basic_astrodynamics::tt_scale:
            return basic_astrodynamics::convertTTtoTAI( input );
        case basic_astrodynamics::tdb_scale: {
            // TDB away from the geocentre depends on UT1, and is inverted within the Earth orientation data.
            if( eopData_ && !position.isZero( 0.0 ) )
            {
                return invertTimeFromTai( inputScale, input, position );
            }

            // TDB - TT changes slowly; three fixed-point updates give TT to Time precision.
            Time tt = input;

            for( int i = 0; i < 3; ++i )
            {
                tt = input - getTdbMinusTt( tt, position );
            }
            return basic_astrodynamics::convertTTtoTAI( tt );
        }
        case basic_astrodynamics::utc_scale:
            // Before 1962, UTC is taken to be equal to UT1.
            if( input < Time( utcIntroductionEpochInTai_ ) )
            {
                return getTaiFromHistoricalUt1( input );
            }
            return sofa_interface::convertUTCtoTAI( input );
        case basic_astrodynamics::ut1_scale:
            if( !eopData_ )
            {
                if( input < Time( utcIntroductionEpochInTai_ ) )
                {
                    return getTaiFromHistoricalUt1( input );
                }
                throw std::runtime_error( "UT1 after 1962 requires Earth orientation data" );
            }
            return invertTimeFromTai( inputScale, input, position );
        default:
            throw std::runtime_error( "Unknown input time scale" );
    }
}

// Convert TAI to the given time scale.
Time TerrestrialTimeScaleConverter::getTimeFromTai( const basic_astrodynamics::TimeScales outputScale,
                                                    const Time& tai,
                                                    const Eigen::Vector3d& position )
{
    const Time tt = basic_astrodynamics::convertTAItoTT( tai );
    switch( outputScale )
    {
        case basic_astrodynamics::tai_scale:
            return tai;
        case basic_astrodynamics::tt_scale:
            return tt;
        case basic_astrodynamics::tdb_scale:
            return tt + getTdbMinusTt( tt, position );
        case basic_astrodynamics::utc_scale:
            // Before 1962, UTC is taken to be equal to UT1.
            if( tai < Time( utcIntroductionEpochInTai_ ) )
            {
                return getHistoricalUt1FromTai( tai );
            }
            return EarthOrientationDataEvaluator::convertTaiToScalarUtc( tai );
        case basic_astrodynamics::ut1_scale:
            return getUt1FromTai( tai );
        default:
            throw std::runtime_error( "Unknown output time scale" );
    }
}

// Find the TAI at which UT1 (or TDB away from the geocentre) has the requested value, by bisection between data points.
Time TerrestrialTimeScaleConverter::invertTimeFromTai( const basic_astrodynamics::TimeScales inputScale,
                                                       const Time& requested,
                                                       const Eigen::Vector3d& position )
{
    // Forward clock values at the raw endpoints define this inverse's coverage.
    const auto& epochs = eopData_->getEpochsTai( );

    const Time first = getTimeFromTai( inputScale, epochs.front( ), position );

    const Time last = getTimeFromTai( inputScale, epochs.back( ), position );

    // A request at, or within the EOP coverage tolerance outside, an end point maps to that observation.
    const double tolerance = EarthOrientationDataEvaluator::coverageToleranceSeconds;

    if( requested <= first && first - requested <= tolerance )
    {
        return epochs.front( );
    }
    if( requested >= last && requested - last <= tolerance )
    {
        return epochs.back( );
    }
    if( requested < first || requested > last )
    {
        // Before 1962, UT1 can be computed from historical Delta T (if selected in the settings).
        if( inputScale == basic_astrodynamics::ut1_scale && requested < first &&
            eopData_->getSettings( ).historicalUt1 == HistoricalUt1Policy::delta_t && requested < Time( utcIntroductionEpochInTai_ ) )
        {
            return getTaiFromHistoricalUt1( requested );
        }
        throw std::runtime_error( inputScale == basic_astrodynamics::ut1_scale ? "Requested UT1 outside corrected raw EOP coverage"
                                                                               : "Requested topocentric TDB outside raw EOP coverage" );
    }

    // Locate the raw TAI interval whose forward clock values bracket the request.
    const auto right = std::lower_bound(
            epochs.begin( ), epochs.end( ), requested, [ this, inputScale, &position ]( const Time& epoch, const Time& input ) {
                return getTimeFromTai( inputScale, epoch, position ) < input;
            } );

    if( right == epochs.end( ) || right == epochs.begin( ) )
    {
        throw std::runtime_error( inputScale == basic_astrodynamics::ut1_scale ? "Cannot bracket UT1 in EOP coverage"
                                                                               : "Cannot bracket topocentric TDB in EOP coverage" );
    }

    if( getTimeFromTai( inputScale, *right, position ) == requested )
    {
        return *right;
    }

    // Bisect in Time; only the residual comparison uses double seconds.
    Time lo = *( right - 1 ), hi = *right;

    for( int i = 0; i < 100; ++i )
    {
        const Time middle = lo + ( hi - lo ) * 0.5;

        const Time residual = getTimeFromTai( inputScale, middle, position ) - requested;

        if( std::abs( static_cast< double >( residual ) ) < 1e-12 )
        {
            return middle;
        }
        if( residual < Time( 0.0 ) )
        {
            lo = middle;
        }
        else
        {
            hi = middle;
        }
        if( inputScale == basic_astrodynamics::ut1_scale && middle == lo && middle == hi )
        {
            return middle;
        }
    }
    throw std::runtime_error( inputScale == basic_astrodynamics::ut1_scale
                                      ? "UT1 inversion did not converge within its EOP interval"
                                      : "Topocentric TDB inversion did not converge within its EOP interval" );
}

//! Function to create the default Earth time scales conversion object
std::shared_ptr< TerrestrialTimeScaleConverter > createDefaultTimeConverter(
        const std::shared_ptr< EOPReader > eopReader,
        const bool reuseExistingConverter,
        const std::optional< EarthOrientationEvaluationSettings >& settings )
{
    // Only the default data and default evaluation settings are shared between callers.
    if( reuseExistingConverter && eopReader == nullptr && !settings.has_value( ) )
    {
        static std::shared_ptr< TerrestrialTimeScaleConverter > defaultTimeConverter = createDefaultTimeConverter( nullptr, false );
        return defaultTimeConverter;
    }

    std::shared_ptr< EOPReader > eopReaderToUse = ( eopReader == nullptr ) ? std::make_shared< EOPReader >( ) : eopReader;
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( eopReaderToUse->getRawEarthOrientationData( ),
                                                                              settings.value_or( EarthOrientationEvaluationSettings( ) ) );
    return std::make_shared< TerrestrialTimeScaleConverter >( evaluator, getDefaultUT1CorrectionCalculator( ) );
}

}  // namespace earth_orientation

}  // namespace tudat
