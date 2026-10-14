#ifndef TUDAT_EARTHORIENTATIONRATES_H
#define TUDAT_EARTHORIENTATIONRATES_H

#include <cmath>
#include <functional>
#include "tudat/astro/earth_orientation/shortPeriodEarthOrientationCorrectionCalculator.h"
#include "tudat/math/basic/numericalDerivative.h"

namespace tudat
{
namespace earth_orientation
{

//! Time derivative of short-period corrections with respect to TT seconds.
//! The harmonic model is differentiated with a fourth-order central difference in TT. UT1 is shifted by dUT1/dTT
//! times the TT shift, using the daily interpolation slope at the requested time, so that no daily EOP
//! observations are queried at the nearby sample times.
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

    // Corrections (radians or UT1 seconds) at a shifted TT, with UT1 shifted at the local daily rate.
    const std::function< Output( const Time& ) > corrections = [ & ]( const Time& shiftedTt ) {
        const Time shiftedUt1 = ut1 + static_cast< double >( shiftedTt - tt ) * dUt1DTt;
        const Eigen::Vector6d arguments = sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst(
                static_cast< double >( shiftedTt ), static_cast< double >( shiftedTt ), static_cast< double >( shiftedUt1 ) );
        return calculator->getCorrectionsFromFundamentalArgument( arguments );
    };
    return numerical_derivatives::computeFourthOrderCentralDifference( corrections, tt, stepSeconds );
}

}  // namespace earth_orientation
}  // namespace tudat
#endif
