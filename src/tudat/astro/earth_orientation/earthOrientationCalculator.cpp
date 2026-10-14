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

#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"
#include "tudat/astro/earth_orientation/earthOrientationRates.h"

namespace tudat
{

namespace earth_orientation
{

//! Function to get the value of the TIO locator, approximated as described in IERS 2010 conventions
double getApproximateTioLocator( const double secondsSinceJ2000 )
{
    return iauSp00( basic_astrodynamics::JULIAN_DAY_ON_J2000, secondsSinceJ2000 / physical_constants::JULIAN_DAY );
}

//! Function to calculate rotation from CIRS to GCRS, i.e. applying CIO-based rotations due to nutation and precession
Eigen::Quaterniond calculateRotationFromCirsToGcrs( const double celestialPoleXPosition,
                                                    const double celestialPoleYPosition,
                                                    const double cioLocator )
{
    // Pre-compute reused parameters
    double xParameterSquared = celestialPoleXPosition * celestialPoleXPosition;
    double yParameterSquared = celestialPoleYPosition * celestialPoleYPosition;
    double xyCrossTerm = celestialPoleXPosition * celestialPoleYPosition;
    const double zSquared = 1.0 - xParameterSquared - yParameterSquared;
    if( !std::isfinite( zSquared ) || zSquared <= 0.0 ) throw std::invalid_argument( "Invalid celestial pole coordinates" );
    double parameterA = 1.0 / ( 1.0 + std::sqrt( zSquared ) );

    // Set up rotation.
    Eigen::Matrix3d rotationMatrix;
    rotationMatrix << 1.0 - parameterA * xParameterSquared, -parameterA * xyCrossTerm, celestialPoleXPosition, -parameterA * xyCrossTerm,
            1.0 - parameterA * yParameterSquared, celestialPoleYPosition, -celestialPoleXPosition, -celestialPoleYPosition,
            1.0 - parameterA * ( xParameterSquared + yParameterSquared );
    return Eigen::Quaterniond( rotationMatrix ) * Eigen::Quaterniond( Eigen::AngleAxisd( -cioLocator, Eigen::Vector3d::UnitZ( ) ) );
}

//! Calculate rotation from TIRS to CIRS.
Eigen::Quaterniond calculateRotationFromTirsToCirs( const double earthRotationAngle )
{
    return Eigen::Quaterniond( Eigen::AngleAxisd( earthRotationAngle, Eigen::Vector3d::UnitZ( ) ) );
}

//! Calculate rotation from ITRS to TIRS.
Eigen::Quaterniond calculateRotationFromItrsToTirs( const double xPolePosition, const double yPolePosition, const double tioLocator )
{
    return Eigen::Quaterniond( Eigen::AngleAxisd( tioLocator, Eigen::Vector3d::UnitZ( ) ) *
                               Eigen::AngleAxisd( -xPolePosition, Eigen::Vector3d::UnitY( ) ) *
                               Eigen::AngleAxisd( -yPolePosition, Eigen::Vector3d::UnitX( ) ) );
}

namespace
{
Eigen::Matrix3d rotationGenerator( const int axis )
{
    Eigen::Matrix3d result = Eigen::Matrix3d::Zero( );
    result( ( axis + 1 ) % 3, ( axis + 2 ) % 3 ) = -1.0;
    result( ( axis + 2 ) % 3, ( axis + 1 ) % 3 ) = 1.0;
    return result;
}

// SOFA DAT is affine in UTC within each pre-1972 segment and constant afterwards.
// Scalar UTC derivatives at an offset/rate change are undefined, including the
// midnight following a leap second. Continuous TAI/TT derivatives remain defined.
double taiRateWithRespectToUtc( const Time& utc )
{
    int year, month, day;
    double fraction, dat0, datHalf;
    if( iauJd2cal( 2451545.0 + utc.fullDaysSinceEpoch( ),
                   static_cast< double >( utc.secondsSinceNoon( ) / 86400.0L ),
                   &year,
                   &month,
                   &day,
                   &fraction ) ||
        iauDat( year, month, day, 0.0, &dat0 ) != 0 || iauDat( year, month, day, 0.5, &datHalf ) != 0 )
        throw std::runtime_error( "UTC rotation rate outside supported SOFA DAT dates" );
    const double driftPerDay = 2.0 * ( datHalf - dat0 );
    if( utc.secondsSinceMidnight( ) == 0.0L )
    {
        double origin, mjd, previous0, previousHalf;
        iauCal2jd( year, month, day, &origin, &mjd );
        iauJd2cal( origin, mjd - 1.0, &year, &month, &day, &fraction );
        if( iauDat( year, month, day, 0.0, &previous0 ) != 0 || iauDat( year, month, day, 0.5, &previousHalf ) != 0 )
            throw std::runtime_error( "UTC rotation rate has no supported left-hand clock" );
        const double previousDrift = 2.0 * ( previousHalf - previous0 );
        if( std::abs( dat0 - previous0 - previousDrift ) > 1e-12 || std::abs( driftPerDay - previousDrift ) > 1e-12 )
            throw std::runtime_error( "UTC rotation rate is undefined at a UTC offset or drift boundary; use TAI/TT" );
    }
    return 1.0 + driftPerDay / 86400.0;
}

template< typename Value >
void requireLegacyStencilInSingleInterval(
        const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Value > >& interpolator,
        const double first,
        const double last )
{
    if( !interpolator ) throw std::runtime_error( "Legacy Earth orientation rate requires its daily interpolator" );
    const auto type = interpolator->getInterpolatorType( );
    if( type != interpolators::linear_interpolator && type != interpolators::discrete_jump_linear_interpolator &&
        type != interpolators::cubic_spline_interpolator && type != interpolators::hermite_spline_interpolator &&
        type != interpolators::lagrange_interpolator && type != interpolators::piecewise_constant_interpolator )
        throw std::runtime_error( "Unsupported legacy Earth orientation interpolator rate" );
    const auto nodes = interpolator->getIndependentValues( );
    if( nodes.size( ) < 2 || first <= nodes.front( ) || last >= nodes.back( ) ||
        std::lower_bound( nodes.begin( ), nodes.end( ), first ) != std::upper_bound( nodes.begin( ), nodes.end( ), last ) )
        throw std::runtime_error( "Legacy Earth orientation rate stencil touches an interpolation knot or coverage boundary" );
}
}  // namespace

Eigen::Matrix3d calculateRotationRateFromItrsToGcrs( const EarthOrientationAnglesAndRates& orientation )
{
    const double x = orientation.angles[ 0 ], y = orientation.angles[ 1 ], s = orientation.angles[ 2 ];
    const double dx = orientation.angleRates[ 0 ], dy = orientation.angleRates[ 1 ], ds = orientation.angleRates[ 2 ];
    const double z = std::sqrt( 1.0 - x * x - y * y );
    if( !orientation.angles.allFinite( ) || !orientation.angleRates.allFinite( ) || !std::isfinite( z ) || z <= 0.0 ||
        !std::isfinite( orientation.dTtDInput ) || !std::isfinite( orientation.dUt1DInput ) )
        throw std::invalid_argument( "Invalid complete Earth rotation derivative inputs" );
    const double a = 1.0 / ( 1.0 + z );
    const double da = a * a * ( x * dx + y * dy ) / z;
    Eigen::Matrix3d b, db;
    b << 1.0 - a * x * x, -a * x * y, x, -a * x * y, 1.0 - a * y * y, y, -x, -y, z;
    db << -da * x * x - 2.0 * a * x * dx, -da * x * y - a * ( dx * y + x * dy ), dx, -da * x * y - a * ( dx * y + x * dy ),
            -da * y * y - 2.0 * a * y * dy, dy, -dx, -dy, -( x * dx + y * dy ) / z;
    const Eigen::Matrix3d generatorZ = rotationGenerator( 2 );
    const Eigen::Matrix3d spinS = Eigen::AngleAxisd( -s, Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( );
    const Eigen::Matrix3d q = b * spinS;
    const Eigen::Matrix3d dq = db * spinS - ds * b * generatorZ * spinS;
    const Eigen::Matrix3d spin =
            calculateRotationFromTirsToCirs( sofa_interface::calculateEarthRotationAngleTemplated< Time >( orientation.ut1 ) )
                    .toRotationMatrix( );
    const double spinRate = 2.0 * mathematical_constants::PI * 1.00273781191135448 / 86400.0 * orientation.dUt1DInput;
    const Eigen::Matrix3d tio =
            Eigen::AngleAxisd( getApproximateTioLocator( static_cast< double >( orientation.tt ) ), Eigen::Vector3d::UnitZ( ) )
                    .toRotationMatrix( );
    const Eigen::Matrix3d polarX = Eigen::AngleAxisd( -orientation.angles[ 3 ], Eigen::Vector3d::UnitY( ) ).toRotationMatrix( );
    const Eigen::Matrix3d polarY = Eigen::AngleAxisd( -orientation.angles[ 4 ], Eigen::Vector3d::UnitX( ) ).toRotationMatrix( );
    // SOFA iauSp00: -47 microarcseconds per Julian century of TT (IERS 2010, 5.13).
    const double tioRate = -47e-6 * mathematical_constants::PI / ( 180.0 * 3600.0 ) / ( 36525.0 * 86400.0 ) * orientation.dTtDInput;
    const Eigen::Matrix3d w = tio * polarX * polarY;
    const Eigen::Matrix3d dw = tioRate * generatorZ * w - orientation.angleRates[ 3 ] * tio * rotationGenerator( 1 ) * polarX * polarY -
            orientation.angleRates[ 4 ] * tio * polarX * rotationGenerator( 0 ) * polarY;
    return dq * spin * w + spinRate * q * generatorZ * spin * w + q * spin * dw;
}

EarthOrientationAnglesAndRates EarthOrientationAnglesCalculator::getRotationAnglesAndRatesFromItrsToGcrs(
        const double input,
        const basic_astrodynamics::TimeScales scale )
{
    return calculateRates( getRotationAnglesAndTimesFromItrsToGcrs( input, scale ), scale, Time( input ) );
}

EarthOrientationAnglesAndRates EarthOrientationAnglesCalculator::getRotationAnglesAndRatesFromItrsToGcrs(
        const Time& input,
        const basic_astrodynamics::TimeScales scale )
{
    return calculateRates( getRotationAnglesAndTimesFromItrsToGcrs( input, scale ), scale, input );
}

EarthOrientationAnglesAndRates EarthOrientationAnglesCalculator::calculateRates( const EarthOrientationAngles& orientation,
                                                                                 const basic_astrodynamics::TimeScales inputScale,
                                                                                 const Time& input )
{
    using namespace basic_astrodynamics;
    EarthOrientationAnglesAndRates result;
    static_cast< EarthOrientationAngles& >( result ) = orientation;
    if( terrestrialTimeScaleConverter_->usesContinuousEop( ) )
    {
        const auto ut1 = terrestrialTimeScaleConverter_->getUt1FromTaiWithRate( result.tai );
        const auto polarEvaluator = polarMotionCalculator_->getEarthOrientationDataEvaluator( );
        const auto cpoEvaluator = precessionNutationCalculator_->getEarthOrientationDataEvaluator( );
        if( !polarEvaluator || !cpoEvaluator ) throw std::runtime_error( "Inconsistent continuous Earth orientation models" );
        const auto polar = polarEvaluator->evaluate( result.tai );
        const auto cpo = cpoEvaluator->evaluate( result.tai, true );
        if( !polar.polarMotionRate || !cpo.celestialPoleOffsetRate )
            throw std::runtime_error( "Selected polar motion or celestial pole interpolation rate is unavailable" );
        const auto nominal = precessionNutationCalculator_->getNominalPositionAndRate( result.tt );
        result.angleRates.head< 3 >( ) = nominal.second;
        result.angleRates.head< 2 >( ) += *cpo.celestialPoleOffsetRate;
        result.angleRates.tail< 2 >( ) = *polar.polarMotionRate;
        const auto shortPeriodPolar = polarMotionCalculator_->getShortPeriodPolarMotionCalculator( );
        if( shortPeriodPolar )
            result.angleRates.tail< 2 >( ) +=
                    calculateShortPeriodCorrectionRate< Eigen::Vector2d >( shortPeriodPolar, result.tt, result.ut1, ut1.dUt1DTai );
        result.dTaiDInput = 1.0;
        switch( inputScale )
        {
            case tai_scale:
            case tt_scale:
                break;
            case ut1_scale:
                result.dTaiDInput = 1.0 / ut1.dUt1DTai;
                break;
            case tdb_scale:
                result.dTaiDInput = terrestrialTimeScaleConverter_->getGeocentricTtRateWithRespectToTdb( input );
                break;
            case utc_scale:
                result.dTaiDInput = taiRateWithRespectToUtc( input );
                break;
            default:
                throw std::runtime_error( "Unknown Earth rotation input time scale" );
        }
        result.dTtDInput = result.dTaiDInput;
        result.dUt1DInput = ut1.dUt1DTai * result.dTaiDInput;
        result.angleRates *= result.dTaiDInput;
    }
    else
    {
        // The legacy constructor may use arbitrary interpolators and approximate
        // clock inverses. Differentiate exactly that supplied position path, only
        // on a stencil inside a single daily interpolation interval and UTC era.
        constexpr double h = 32.0;
        if( static_cast< double >( result.tai ) < ( 41317.0 - 51544.5 ) * 86400.0 + 10.0 )
            throw std::runtime_error( "Legacy Earth rotation rates before 1972 are unsupported; use continuous EOP" );
        const auto first = getRotationAnglesAndTimesFromItrsToGcrs( input - 2.0 * h, inputScale );
        const auto last = getRotationAnglesAndTimesFromItrsToGcrs( input + 2.0 * h, inputScale );
        const Time firstUtc = EarthOrientationDataEvaluator::convertTaiToScalarUtc( first.tai );
        const Time lastUtc = EarthOrientationDataEvaluator::convertTaiToScalarUtc( last.tai );
        if( std::abs( static_cast< double >( ( last.tai - lastUtc ) - ( first.tai - firstUtc ) ) ) > 1e-12 )
            throw std::runtime_error( "Legacy Earth rotation rate stencil crosses a UTC discontinuity" );
        requireLegacyStencilInSingleInterval( polarMotionCalculator_->getDailyIersValueInterpolator( ),
                                              static_cast< double >( firstUtc ),
                                              static_cast< double >( lastUtc ) );
        requireLegacyStencilInSingleInterval( precessionNutationCalculator_->getDailyCorrectionInterpolator( ),
                                              static_cast< double >( firstUtc ),
                                              static_cast< double >( lastUtc ) );
        requireLegacyStencilInSingleInterval( terrestrialTimeScaleConverter_->getDailyUtcUt1CorrectionInterpolator( ),
                                              static_cast< double >( firstUtc ),
                                              static_cast< double >( lastUtc ) );
        if( precessionNutationCalculator_->getNominalCipInterpolator( ) )
            requireLegacyStencilInSingleInterval( precessionNutationCalculator_->getNominalCipInterpolator( ),
                                                  static_cast< double >( first.tt ),
                                                  static_cast< double >( last.tt ) );
        const auto before = getRotationAnglesAndTimesFromItrsToGcrs( input - h, inputScale );
        const auto after = getRotationAnglesAndTimesFromItrsToGcrs( input + h, inputScale );
        result.angleRates = ( 8.0 * ( after.angles - before.angles ) - ( last.angles - first.angles ) ) / ( 12.0 * h );
        result.dTtDInput =
                ( 8.0 * static_cast< double >( after.tt - before.tt ) - static_cast< double >( last.tt - first.tt ) ) / ( 12.0 * h );
        result.dTaiDInput = result.dTtDInput;
        result.dUt1DInput =
                ( 8.0 * static_cast< double >( after.ut1 - before.ut1 ) - static_cast< double >( last.ut1 - first.ut1 ) ) / ( 12.0 * h );
        getRotationAnglesAndTimesFromItrsToGcrs( input, inputScale );
    }
    if( !result.angleRates.allFinite( ) || !std::isfinite( result.dTaiDInput ) || !std::isfinite( result.dUt1DInput ) ||
        result.dTaiDInput <= 0.0 || result.dUt1DInput <= 0.0 )
        throw std::runtime_error( "Invalid Earth orientation rate" );
    return result;
}

//! Function to create an EarthOrientationAnglesCalculator object, with default settings
std::shared_ptr< EarthOrientationAnglesCalculator > createStandardEarthOrientationCalculator(
        const std::shared_ptr< EOPReader > eopReader,
        const EarthOrientationEvaluationSettings& settings )
{
    const auto evaluator = std::make_shared< EarthOrientationDataEvaluator >( eopReader->getRawEarthOrientationData( ), settings );
    const auto ut1 = std::make_shared< EarthOrientationDataInterpolator< double > >( evaluator, EarthOrientationDataField::ut1_minus_utc );
    const auto timeConverter =
            std::make_shared< TerrestrialTimeScaleConverter >( ut1, getDefaultUT1CorrectionCalculator( ), nullptr, evaluator );
    const auto polar =
            std::make_shared< EarthOrientationDataInterpolator< Eigen::Vector2d > >( evaluator, EarthOrientationDataField::polar_motion );
    const auto cpo =
            std::make_shared< EarthOrientationDataInterpolator< Eigen::Vector2d > >( evaluator, EarthOrientationDataField::celestial_pole );
    const auto polarModel = std::make_shared< PolarMotionCalculator >(
            polar, getDefaultPolarMotionCorrectionCalculator( ), evaluator, [ timeConverter ]( const Time& tai ) {
                return timeConverter->getUt1FromTai( tai );
            } );
    const auto precession = std::make_shared< PrecessionNutationCalculator >( basic_astrodynamics::iau_2006, cpo, nullptr, evaluator );
    return std::make_shared< EarthOrientationAnglesCalculator >( polarModel, precession, timeConverter );
}

////! Function to create an interpolator for the Earth orientation angles
// std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 6,1 > > >
// createInterpolatorForItrsToGcrsAngles(
//         const double intervalStart, const double intervalEnd, const double timeStep,
//         const basic_astrodynamics::TimeScales timeScale,
//         const std::shared_ptr< EarthOrientationAnglesCalculator > earthOrientationCalculator,
//         const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings )
//{
//     // Interpolate Earth orientation angles
//     std::map< double, Eigen::Matrix< double, 6,1 > > orientationMap;
//     double currentTime = intervalStart;
//     while( currentTime < intervalEnd )
//     {
//         orientationMap[ currentTime ] = earthOrientationCalculator->getRotationAnglesFromItrsToGcrs(
//                     currentTime, timeScale );
//         currentTime += timeStep;
//     }

//    return interpolators::createOneDimensionalInterpolator(
//                orientationMap, interpolatorSettings );
//}

}  // namespace earth_orientation

}  // namespace tudat
