#ifndef TUDAT_EARTHORIENTATIONRATES_H
#define TUDAT_EARTHORIENTATIONRATES_H

#include <cmath>
#include <limits>
#include <type_traits>
#include "tudat/astro/earth_orientation/shortPeriodEarthOrientationCorrectionCalculator.h"

namespace tudat
{
namespace earth_orientation
{

//! Differentiate small harmonic sums along the local TT/UT1 argument tangent.
//! No EOP samples are queried outside the central selected interpolation stencil.
template< typename Output >
Output calculateShortPeriodCorrectionRate( const std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Output > >& calculator,
                                           const Time& tt,
                                           const Time& ut1,
                                           const double dUt1DTt,
                                           const double h = 32.0 )
{
    if( !calculator || !std::isfinite( h ) || h <= 0.0 || !std::isfinite( dUt1DTt ) )
        throw std::invalid_argument( "Invalid short-period orientation rate inputs" );
    const auto sample = [ & ]( const double delta ) -> Output {
        const Time shiftedTt = tt + delta;
        const long double separation = static_cast< long double >( shiftedTt - tt );
        if( !std::isfinite( separation ) ||
            std::abs( separation - delta ) > 16.0L * std::numeric_limits< double >::epsilon( ) * std::abs( delta ) )
            throw std::invalid_argument( "Unrepresentable short-period orientation rate stencil" );
        const Time shiftedUt1 = ut1 + delta * dUt1DTt;
        const auto arguments = sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst(
                static_cast< double >( shiftedTt ), static_cast< double >( shiftedTt ), static_cast< double >( shiftedUt1 ) );
        const Output result = calculator->getCorrectionsFromFundamentalArgument( arguments );
        if constexpr( std::is_arithmetic< Output >::value )
        {
            if( !std::isfinite( result ) ) throw std::runtime_error( "Non-finite short-period orientation correction" );
        }
        else if( !result.allFinite( ) )
            throw std::runtime_error( "Non-finite short-period orientation correction" );
        return result;
    };
    return ( 8.0 * ( sample( h ) - sample( -h ) ) - ( sample( 2.0 * h ) - sample( -2.0 * h ) ) ) / ( 12.0 * h );
}

}  // namespace earth_orientation
}  // namespace tudat
#endif
