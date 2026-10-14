#ifndef TUDAT_EARTHORIENTATIONRATES_H
#define TUDAT_EARTHORIENTATIONRATES_H

#include <array>
#include <cmath>
#include <limits>
#include <type_traits>
#include "tudat/astro/earth_orientation/shortPeriodEarthOrientationCorrectionCalculator.h"

namespace tudat
{
namespace earth_orientation
{

//! Time derivative of short-period corrections with respect to TT seconds.
//! Sample the harmonic model at four nearby TT times. Shift UT1 by dUT1/dTT
//! times the TT shift, using the daily interpolation slope at the requested time.
//! This does not query daily EOP observations at the nearby sample times.
template< typename Output >
Output calculateShortPeriodCorrectionRate( const std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Output > >& calculator,
                                           const Time& tt,
                                           const Time& ut1,
                                           const double dUt1DTt,
                                           const double stepSeconds = 32.0 )
{
    if( !calculator || !std::isfinite( stepSeconds ) || stepSeconds <= 0.0 || !std::isfinite( dUt1DTt ) )
    {
        throw std::invalid_argument( "Invalid short-period orientation rate inputs" );
    }
    // Four TT offsets for the fourth-order derivative, in seconds.
    const std::array< double, 4 > offsets = { -2.0 * stepSeconds, -stepSeconds, stepSeconds, 2.0 * stepSeconds };

    // Harmonic corrections at each shifted time, in radians or UT1 seconds.
    std::array< Output, 4 > corrections;

    for( std::size_t i = 0; i < offsets.size( ); ++i )
    {
        // Apply this sample's TT shift and check that Time can represent it.
        const double offset = offsets[ i ];

        const Time shiftedTt = tt + offset;

        // The accepted spacing error is proportional to the requested shift.
        const long double actualOffset = static_cast< long double >( shiftedTt - tt );

        const long double tolerance = 16.0L * std::numeric_limits< double >::epsilon( ) * std::abs( offset );

        if( !std::isfinite( actualOffset ) || std::abs( actualOffset - offset ) > tolerance )
        {
            throw std::invalid_argument( "Unrepresentable short-period orientation rate stencil" );
        }
        // UT1 advances at the local daily-interpolation rate relative to TT.
        const Time shiftedUt1 = ut1 + offset * dUt1DTt;

        // Evaluate the same astronomical argument convention at the shifted clocks.
        const Eigen::Vector6d arguments = sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst(
                static_cast< double >( shiftedTt ), static_cast< double >( shiftedTt ), static_cast< double >( shiftedUt1 ) );

        corrections[ i ] = calculator->getCorrectionsFromFundamentalArgument( arguments );
        if constexpr( std::is_arithmetic< Output >::value )
        {
            if( !std::isfinite( corrections[ i ] ) )
            {
                throw std::runtime_error( "Non-finite short-period orientation correction" );
            }
        }
        else if( !corrections[ i ].allFinite( ) )
        {
            throw std::runtime_error( "Non-finite short-period orientation correction" );
        }
    }

    // Pair opposite samples before weighting, to avoid subtracting large sums.
    return ( 8.0 * ( corrections[ 2 ] - corrections[ 1 ] ) - ( corrections[ 3 ] - corrections[ 0 ] ) ) / ( 12.0 * stepSeconds );
}

}  // namespace earth_orientation
}  // namespace tudat
#endif
