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
#include "tudat/astro/earth_orientation/earthOrientationRates.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "tudat/math/interpolators/jumpDataLinearInterpolator.h"

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

double TerrestrialTimeScaleConverter::getContinuousTdbMinusTt( const Time& tt, const Eigen::Vector3d& position )
{
    const auto key = std::make_tuple( position.x( ), position.y( ), position.z( ) );
    if( tdbToTtInterpolators_.count( key ) ) return tdbToTtInterpolators_.at( key )->interpolate( static_cast< double >( tt ) );
    if( !continuousEopData_ && !position.isZero( 0.0 ) ) return sofa_interface::getTDBminusTT( static_cast< double >( tt ), position );
    // Geocentric TT/TDB needs no EOP. The topocentric SOFA terms require
    // actual UT1 and therefore share the strict TAI-axis coverage policy.
    // The UT1 fraction is irrelevant when all observer terms vanish. Using TT
    // in that case avoids UTC/DAT lookups, including SOFA's future-DAT warning.
    const Time ut1 = position.isZero( 0.0 ) ? tt : getUt1FromTai( basic_astrodynamics::convertTTtoTAI( tt ) );
    const double fraction = static_cast< double >( ut1.secondsSinceMidnight( ) / 86400.0L );
    return sofa_interface::getTDBminusTT( static_cast< double >( tt ), fraction, position );
}

double TerrestrialTimeScaleConverter::getGeocentricTtRateWithRespectToTdb( const Time& tdb )
{
    if( tdbToTtInterpolators_.count( std::make_tuple( 0.0, 0.0, 0.0 ) ) )
        throw std::runtime_error( "TT/TDB interpolator clock rates are unavailable; position conversion remains supported" );
    if( !std::isfinite( static_cast< double >( tdb ) ) ) throw std::invalid_argument( "Non-finite TDB clock-rate epoch" );
    constexpr double h = 32.0;
    const auto offset = [ this, &tdb ]( const double delta ) {
        const Time shifted = tdb + delta;
        const long double actualSeparation = static_cast< long double >( shifted - tdb );
        const long double tolerance = 16.0L * std::numeric_limits< double >::epsilon( ) * std::abs( delta );
        if( !std::isfinite( actualSeparation ) || std::abs( actualSeparation - delta ) > tolerance )
            throw std::invalid_argument( "Unrepresentable geocentric clock-rate stencil" );
        return static_cast< double >( basic_astrodynamics::convertTAItoTT( getTaiTime( basic_astrodynamics::tdb_scale, shifted ) ) -
                                      shifted );
    };
    const double rate = 1.0 + ( 8.0 * ( offset( h ) - offset( -h ) ) - ( offset( 2.0 * h ) - offset( -2.0 * h ) ) ) / ( 12.0 * h );
    if( !std::isfinite( rate ) || rate <= 0.0 ) throw std::runtime_error( "Invalid native geocentric TT/TDB clock rate" );
    return rate;
}

// Return the modeled UT1 position without requiring a daily interpolation slope.
Time TerrestrialTimeScaleConverter::getUt1FromTai( const Time& tai )
{
    return evaluateUt1FromTai( tai, false ).ut1;
}

// Return the same UT1 value together with the derivative of its actual daily and short-period models.
Ut1ValueAndRate TerrestrialTimeScaleConverter::getUt1FromTaiWithRate( const Time& tai )
{
    return evaluateUt1FromTai( tai, true );
}

// Share the value calculation while keeping unsupported rates independent of position availability.
Ut1ValueAndRate TerrestrialTimeScaleConverter::evaluateUt1FromTai( const Time& tai, const bool requireRate )
{
    if( !continuousEopData_ )
    {
        throw std::runtime_error( requireRate ? "UT1 rates require a supported continuous EOP model"
                                             : "Continuous UT1 evaluation requires EOP data" );
    }
    if( tai < continuousEopData_->getStartTai( ) && continuousEopData_->getSettings( ).historicalUt1 == HistoricalUt1Policy::delta_t &&
        tai < Time( utcIntroductionEpochInTai_ ) )
    {
        if( requireRate )
        {
            throw std::runtime_error( "Historical Delta-T clock rates are unavailable; position conversion remains supported" );
        }
        loadHistoricalDeltaTInterpolator( );

        // Historical Delta-T is TT minus UT1 in seconds; it supplies no supported derivative.
        const Time tt = basic_astrodynamics::convertTAItoTT( tai );

        return { tt - historicalDeltaTInterpolator_->interpolate(
                              static_cast< double >( tt ) / physical_constants::JULIAN_YEAR + 2000.0 ) };
    }

    // Daily UT1-TAI and its optional slope come from the selected interpolation policy.
    const auto daily = continuousEopData_->evaluate( tai );

    if( requireRate && !daily.ut1MinusTaiRate )
    {
        throw std::runtime_error( "Selected UT1-TAI interpolation slope is unavailable" );
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
        // Short-period values and tangents are evaluated at the same baseline clocks.
        const auto arguments = sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst(
                static_cast< double >( tt ), static_cast< double >( tt ), static_cast< double >( baselineUt1 ) );

        result.ut1 += shortPeriodUt1CorrectionCalculator_->getCorrectionsFromFundamentalArgument( arguments );
        if( requireRate )
        {
            // The correction's tangent uses the baseline daily rate, not the final corrected UT1 rate.
            result.dUt1DTai +=
                    calculateShortPeriodCorrectionRate< double >( shortPeriodUt1CorrectionCalculator_, tt, baselineUt1, result.dUt1DTai );
        }
    }
    if( requireRate && ( !std::isfinite( result.dUt1DTai ) || result.dUt1DTai <= 0.0 ) )
    {
        throw std::runtime_error( "Invalid modeled UT1 rate" );
    }
    return result;
}

// Bracket and bisect the requested continuous clock using the original raw TAI nodes.
Time TerrestrialTimeScaleConverter::invertContinuousTime( const basic_astrodynamics::TimeScales inputScale,
                                                          const Time& requested,
                                                          const Eigen::Vector3d& position,
                                                          const bool roundedDoubleInput )
{
    // Forward clock values at the raw endpoints define this inverse's coverage.
    const auto& epochs = continuousEopData_->getEpochsTai( );

    const Time first = getContinuousOutputTime( inputScale, epochs.front( ), position );

    const Time last = getContinuousOutputTime( inputScale, epochs.back( ), position );

    if( requested == first || ( roundedDoubleInput && static_cast< double >( requested ) == static_cast< double >( first ) ) )
    {
        return epochs.front( );
    }
    if( requested == last || ( roundedDoubleInput && static_cast< double >( requested ) == static_cast< double >( last ) ) )
    {
        return epochs.back( );
    }
    if( requested < first || requested > last )
    {
        // Only UT1 has the existing pre-UTC historical Delta-T inverse.
        if( inputScale == basic_astrodynamics::ut1_scale && requested < first &&
            continuousEopData_->getSettings( ).historicalUt1 == HistoricalUt1Policy::delta_t &&
            requested < Time( utcIntroductionEpochInTai_ ) )
        {
            loadHistoricalDeltaTInterpolator( );

            // Retain the four fixed-point updates used by the historical UT1 route.
            Time tt = requested;

            for( int i = 0; i < 4; ++i )
            {
                tt = requested +
                        historicalDeltaTInterpolator_->interpolate( static_cast< double >( tt ) / physical_constants::JULIAN_YEAR +
                                                                    2000.0 );
            }
            return basic_astrodynamics::convertTTtoTAI( tt );
        }
        throw std::runtime_error( inputScale == basic_astrodynamics::ut1_scale ? "Requested UT1 outside corrected raw EOP coverage"
                                                                             : "Requested topocentric TDB outside raw EOP coverage" );
    }

    // Locate the raw TAI interval whose forward clock values bracket the request.
    const auto right = std::lower_bound( epochs.begin( ), epochs.end( ), requested,
                                         [ this, inputScale, &position ]( const Time& epoch, const Time& input ) {
                                             return getContinuousOutputTime( inputScale, epoch, position ) < input;
                                         } );

    if( right == epochs.end( ) || right == epochs.begin( ) )
    {
        throw std::runtime_error( inputScale == basic_astrodynamics::ut1_scale ? "Cannot bracket UT1 in EOP coverage"
                                                                             : "Cannot bracket topocentric TDB in EOP coverage" );
    }

    if( getContinuousOutputTime( inputScale, *right, position ) == requested )
    {
        return *right;
    }

    // Bisect in Time; only the residual comparison uses double seconds.
    Time lo = *( right - 1 ), hi = *right;

    for( int i = 0; i < 100; ++i )
    {
        const Time middle = lo + ( hi - lo ) * 0.5;

        const Time residual = getContinuousOutputTime( inputScale, middle, position ) - requested;

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
    throw std::runtime_error( inputScale == basic_astrodynamics::ut1_scale ? "UT1 inversion did not converge within its EOP interval"
                                                                         : "Topocentric TDB inversion did not converge within its EOP interval" );
}

//! Function to get current time list at double precision
template<>
CurrentTimes< double >& TerrestrialTimeScaleConverter::getCurrentTimeList< double >( )
{
    return currentTimes_;
}

//! Function to get current time list at Time precision
template<>
CurrentTimes< Time >& TerrestrialTimeScaleConverter::getCurrentTimeList< Time >( )
{
    return currentTimesSplit_;
}

//! Function to get ground station position used on last call to updateTimes function at double precision
template<>
Eigen::Vector3d& TerrestrialTimeScaleConverter::getPreviousGroundStationPosition< double >( )
{
    return previousEarthFixedPosition_;
}

//! Function to get ground station position used on last call to updateTimes function at Time precision
template<>
Eigen::Vector3d& TerrestrialTimeScaleConverter::getPreviousGroundStationPosition< Time >( )
{
    return previousEarthFixedPositionSplit_;
}

//! Function to reset ground station position used on last call to updateTimes function at double precision
template<>
void TerrestrialTimeScaleConverter::setCurrentGroundStation< double >( const Eigen::Vector3d& currentGroundStation )
{
    previousEarthFixedPosition_ = currentGroundStation;
}

//! Function to reset ground station position used on last call to updateTimes function at Time precision
template<>
void TerrestrialTimeScaleConverter::setCurrentGroundStation< Time >( const Eigen::Vector3d& currentGroundStation )
{
    previousEarthFixedPositionSplit_ = currentGroundStation;
}

//! Function to create the default Earth time scales conversion object
std::shared_ptr< TerrestrialTimeScaleConverter > createDefaultTimeConverter( const std::shared_ptr< EOPReader > eopReader,
                                                                             const bool reuseExistingConverter,
                                                                             const EarthOrientationEvaluationSettings& settings )
{
    if( reuseExistingConverter && eopReader == nullptr && settings.maximumGapSeconds == 36.0 * 3600.0 &&
        settings.missingCelestialPole == MissingCelestialPolePolicy::reject && settings.historicalUt1 == HistoricalUt1Policy::reject &&
        settings.interpolation == EarthOrientationInterpolation::linear )
    {
        static std::shared_ptr< TerrestrialTimeScaleConverter > defaultTimeConverter = createDefaultTimeConverter( nullptr, false );
        return defaultTimeConverter;
    }

    std::shared_ptr< EOPReader > eopReaderToUse = ( eopReader == nullptr ) ? std::make_shared< EOPReader >( ) : eopReader;
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( eopReaderToUse->getRawEarthOrientationData( ), settings );
    return std::make_shared< TerrestrialTimeScaleConverter >( nullptr, getDefaultUT1CorrectionCalculator( ), nullptr, evaluator );
}

}  // namespace earth_orientation

}  // namespace tudat
