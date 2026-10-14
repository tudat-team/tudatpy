#include "tudat/astro/ground_stations/earthDeformationTime.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

namespace tudat
{
namespace basic_astrodynamics
{

std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > getGeocentricTimeConverter( )
{
    static const auto converter = std::make_shared< earth_orientation::TerrestrialTimeScaleConverter >( nullptr, nullptr );
    return converter;
}

Time normalizeScalarTdb( const double tdb, const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter )
{
    if( !std::isfinite( tdb ) ) throw std::invalid_argument( "Station TDB epoch must be finite." );
    if( converter && converter->getEarthOrientationDataEvaluator( ) )
    {
        // Reuse the converter's original-double candidate filtering: an interior
        // clock interpolator must not be probed at unrelated EOP endpoints.
        const Time tai = converter->getTaiTime( tdb_scale, tdb );
        const auto& evaluator = converter->getEarthOrientationDataEvaluator( );
        const Time start = evaluator->getStartTai( ), end = evaluator->getEndTai( );
        if( tai == start || tai == end )
        {
            const Time image = converter->getCurrentTime< Time >( tai_scale, tdb_scale, tai );
            if( static_cast< double >( image ) == tdb )
            {
                const Time other = tai == start ? end : start;
                const double candidateWindow =
                        std::max( 1e-6, 8 * std::numeric_limits< double >::epsilon( ) * std::abs( static_cast< double >( tai ) ) );
                if( other != tai && std::abs( static_cast< double >( other - tai ) ) <= candidateWindow &&
                    static_cast< double >( converter->getCurrentTime< Time >( tai_scale, tdb_scale, other ) ) == tdb )
                    throw std::invalid_argument( "Double TDB endpoint image is ambiguous; supply precise Time." );
                return image;
            }
        }
    }
    return Time( tdb );
}

void validateEarthDeformationUtcInterval( const double lowerTt, const double upperTt )
{
    if( !std::isfinite( lowerTt ) || !std::isfinite( upperTt ) || lowerTt > upperTt )
    {
        throw std::invalid_argument( "Earth deformation derivative interval must have finite ordered TT epochs." );
    }
    const double intervalWidth = upperTt - lowerTt;
    if( !std::isfinite( intervalWidth ) )
    {
        throw std::invalid_argument( "Earth deformation derivative interval width must be finite." );
    }
    double minimumUtc = std::numeric_limits< double >::infinity( );
    double maximumUtc = -minimumUtc;
    double epochTolerance = 1.0e-6;
    for( int i = 0; i < 5; ++i )
    {
        const double tt = lowerTt + intervalWidth * ( static_cast< double >( i ) / 4.0 );
        const double tai = convertTTtoTAI( tt );
        const double utc = sofa_interface::convertTTtoUTC( tt );
        const double reconstructedTai = sofa_interface::convertUTCtoTAI( utc );
        if( !std::isfinite( tai ) || !std::isfinite( utc ) || !std::isfinite( reconstructedTai ) )
        {
            throw std::runtime_error( "Non-finite Earth deformation TT/UTC conversion." );
        }
        // Bound the rounding tolerance independently of epoch magnitude, so it cannot
        // hide the historical 0.05-0.1 s steps or the modern 1 s leap-second offsets.
        const double tolerance =
                std::min( 1.0e-5, std::max( 1.0e-6, 16.0 * std::numeric_limits< double >::epsilon( ) * std::fabs( tai ) ) );
        epochTolerance = std::max( epochTolerance, tolerance );
        if( std::fabs( reconstructedTai - tai ) > tolerance )
        {
            throw std::domain_error(
                    "Earth deformation derivative interval includes ambiguous scalar UTC; "
                    "this is a current UTC-phase representation limitation." );
        }
        minimumUtc = std::min( minimumUtc, utc );
        maximumUtc = std::max( maximumUtc, utc );
    }

    const double day = physical_constants::JULIAN_DAY;
    double midnight = std::floor( ( minimumUtc + day / 2.0 ) / day ) * day - day / 2.0;
    if( !std::isfinite( midnight ) )
    {
        throw std::domain_error( "Earth deformation UTC midnight cannot be represented." );
    }
    while( midnight <= maximumUtc + epochTolerance )
    {
        if( midnight >= minimumUtc - epochTolerance )
        {
            const double beforeOne = sofa_interface::getDeltaAtFromUtc( ( midnight - 1.0 ) / day );
            const double beforeTwo = sofa_interface::getDeltaAtFromUtc( ( midnight - 2.0 ) / day );
            const double atMidnight = sofa_interface::getDeltaAtFromUtc( midnight / day );
            if( !std::isfinite( beforeOne ) || !std::isfinite( beforeTwo ) || !std::isfinite( atMidnight ) )
            {
                throw std::runtime_error( "Non-finite Earth deformation UTC offset." );
            }
            // Extrapolate the local smooth drift to midnight. Constancy of DAT is
            // deliberately not required: UTC drifted smoothly before 1972.
            if( std::fabs( atMidnight - ( 2.0 * beforeOne - beforeTwo ) ) > 1.0e-7 )
            {
                throw std::domain_error(
                        "Earth deformation derivative interval touches a UTC offset step; "
                        "this is a current UTC-phase representation limitation." );
            }
        }
        const double nextMidnight = midnight + day;
        if( !std::isfinite( nextMidnight ) || nextMidnight <= midnight )
        {
            throw std::domain_error( "Earth deformation UTC midnight iteration cannot progress." );
        }
        midnight = nextMidnight;
    }
}

}  // namespace basic_astrodynamics
}  // namespace tudat
