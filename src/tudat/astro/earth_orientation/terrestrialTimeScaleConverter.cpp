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

Time TerrestrialTimeScaleConverter::getUt1FromTai( const Time& tai )
{
    if( !continuousEopData_ ) throw std::runtime_error( "Continuous UT1 evaluation requires EOP data" );
    if( tai < continuousEopData_->getStartTai( ) && continuousEopData_->getSettings( ).historicalUt1 == HistoricalUt1Policy::delta_t &&
        static_cast< double >( tai ) < utcIntroductionEpochInTai_ )
    {
        loadHistoricalDeltaTInterpolator( );
        const Time tt = basic_astrodynamics::convertTAItoTT( tai );
        return tt - historicalDeltaTInterpolator_->interpolate( static_cast< double >( tt ) / physical_constants::JULIAN_YEAR + 2000.0 );
    }
    const double baseline = continuousEopData_->evaluate( tai ).ut1MinusTai;
    const Time tt = basic_astrodynamics::convertTAItoTT( tai );
    Time ut1 = tai + baseline;
    if( shortPeriodUt1CorrectionCalculator_ )
    {
        const Eigen::Vector6d arguments = sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst(
                static_cast< double >( tt ), static_cast< double >( tt ), static_cast< double >( ut1 ) );
        ut1 += shortPeriodUt1CorrectionCalculator_->getCorrectionsFromFundamentalArgument( arguments );
    }
    return ut1;
}

Ut1ValueAndRate TerrestrialTimeScaleConverter::getUt1FromTaiWithRate( const Time& tai )
{
    if( !continuousEopData_ ) throw std::runtime_error( "UT1 rates require a supported continuous EOP model" );
    if( tai < continuousEopData_->getStartTai( ) && continuousEopData_->getSettings( ).historicalUt1 == HistoricalUt1Policy::delta_t &&
        static_cast< double >( tai ) < utcIntroductionEpochInTai_ )
        throw std::runtime_error( "Historical Delta-T clock rates are unavailable; position conversion remains supported" );
    const auto daily = continuousEopData_->evaluate( tai );
    if( !daily.ut1MinusTaiRate ) throw std::runtime_error( "Selected UT1-TAI interpolation slope is unavailable" );
    const Time tt = basic_astrodynamics::convertTAItoTT( tai );
    const Time baselineUt1 = tai + daily.ut1MinusTai;
    Ut1ValueAndRate result{ baselineUt1, 1.0 + *daily.ut1MinusTaiRate };
    if( shortPeriodUt1CorrectionCalculator_ )
    {
        const auto arguments = sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst(
                static_cast< double >( tt ), static_cast< double >( tt ), static_cast< double >( baselineUt1 ) );
        result.ut1 += shortPeriodUt1CorrectionCalculator_->getCorrectionsFromFundamentalArgument( arguments );
        // The position path evaluates this correction once at baseline UT1.
        // Its tangent therefore uses the baseline daily rate, not final corrected UT1.
        result.dUt1DTai +=
                calculateShortPeriodCorrectionRate< double >( shortPeriodUt1CorrectionCalculator_, tt, baselineUt1, result.dUt1DTai );
    }
    if( !std::isfinite( result.dUt1DTai ) || result.dUt1DTai <= 0.0 ) throw std::runtime_error( "Invalid modeled UT1 rate" );
    return result;
}

Time TerrestrialTimeScaleConverter::invertContinuousUt1( const Time& requested, const bool roundedDoubleInput )
{
    const auto& epochs = continuousEopData_->getEpochsTai( );
    const Time firstUt1 = getUt1FromTai( epochs.front( ) );
    const Time lastUt1 = getUt1FromTai( epochs.back( ) );
    if( requested == firstUt1 || ( roundedDoubleInput && static_cast< double >( requested ) == static_cast< double >( firstUt1 ) ) )
        return epochs.front( );
    if( requested == lastUt1 || ( roundedDoubleInput && static_cast< double >( requested ) == static_cast< double >( lastUt1 ) ) )
        return epochs.back( );
    if( requested < firstUt1 || requested > lastUt1 )
    {
        if( requested < firstUt1 && continuousEopData_->getSettings( ).historicalUt1 == HistoricalUt1Policy::delta_t &&
            static_cast< double >( requested ) < utcIntroductionEpochInTai_ )
        {
            loadHistoricalDeltaTInterpolator( );
            Time tt = requested;
            for( int i = 0; i < 4; ++i )
                tt = requested +
                        historicalDeltaTInterpolator_->interpolate( static_cast< double >( tt ) / physical_constants::JULIAN_YEAR +
                                                                    2000.0 );
            return basic_astrodynamics::convertTTtoTAI( tt );
        }
        throw std::runtime_error( "Requested UT1 outside corrected raw EOP coverage" );
    }
    const auto right = std::lower_bound( epochs.begin( ), epochs.end( ), requested, [ this ]( const Time& epoch, const Time& ut1 ) {
        return getUt1FromTai( epoch ) < ut1;
    } );
    if( right == epochs.end( ) || right == epochs.begin( ) ) throw std::runtime_error( "Cannot bracket UT1 in EOP coverage" );
    if( getUt1FromTai( *right ) == requested ) return *right;
    Time lo = *( right - 1 ), hi = *right;
    for( int i = 0; i < 100; ++i )
    {
        const Time middle = lo + ( hi - lo ) * 0.5;
        const Time residual = getUt1FromTai( middle ) - requested;
        if( std::abs( static_cast< double >( residual ) ) < 1e-12 ) return middle;
        if( residual < Time( 0.0 ) )
            lo = middle;
        else
            hi = middle;
        if( middle == lo && middle == hi ) return middle;
    }
    throw std::runtime_error( "UT1 inversion did not converge within its EOP interval" );
}

Time TerrestrialTimeScaleConverter::invertContinuousTdb( const Time& requested,
                                                         const Eigen::Vector3d& position,
                                                         const bool roundedDoubleInput )
{
    const auto tdbFromTai = [ this, &position ]( const Time& tai ) {
        const Time tt = basic_astrodynamics::convertTAItoTT( tai );
        return tt + getContinuousTdbMinusTt( tt, position );
    };
    const auto& epochs = continuousEopData_->getEpochsTai( );
    const Time firstTdb = tdbFromTai( epochs.front( ) );
    const Time lastTdb = tdbFromTai( epochs.back( ) );
    if( requested == firstTdb || ( roundedDoubleInput && static_cast< double >( requested ) == static_cast< double >( firstTdb ) ) )
        return epochs.front( );
    if( requested == lastTdb || ( roundedDoubleInput && static_cast< double >( requested ) == static_cast< double >( lastTdb ) ) )
        return epochs.back( );
    if( requested < firstTdb || requested > lastTdb ) throw std::runtime_error( "Requested topocentric TDB outside raw EOP coverage" );
    const auto right = std::lower_bound( epochs.begin( ), epochs.end( ), requested, [ &tdbFromTai ]( const Time& tai, const Time& tdb ) {
        return tdbFromTai( tai ) < tdb;
    } );
    if( right == epochs.end( ) || right == epochs.begin( ) ) throw std::runtime_error( "Cannot bracket topocentric TDB in EOP coverage" );
    if( tdbFromTai( *right ) == requested ) return *right;
    Time lo = *( right - 1 ), hi = *right;
    for( int i = 0; i < 100; ++i )
    {
        const Time middle = lo + ( hi - lo ) * 0.5;
        const Time residual = tdbFromTai( middle ) - requested;
        if( std::abs( static_cast< double >( residual ) ) < 1e-12 ) return middle;
        if( residual < Time( 0.0 ) )
            lo = middle;
        else
            hi = middle;
    }
    throw std::runtime_error( "Topocentric TDB inversion did not converge within its EOP interval" );
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
    const auto view = std::make_shared< EarthOrientationDataInterpolator< double > >( evaluator, EarthOrientationDataField::ut1_minus_utc );
    return std::make_shared< TerrestrialTimeScaleConverter >( view, getDefaultUT1CorrectionCalculator( ), nullptr, evaluator );
}

}  // namespace earth_orientation

}  // namespace tudat
