#include "tudat/math/interpolators/chebyshevApproximation.h"

namespace tudat
{

namespace interpolators
{

//! Convert point to Chebyshev-independent variable in [-1, 1] domain
double determineChebyshevPolynomialIndependentVariable(
        const double intervalStart,
        const double intervalEnd,
        const double evaluationPoint )
{
    return 2.0 * ( evaluationPoint - intervalStart ) / ( intervalEnd - intervalStart ) - 1.0;
}

//! Evaluate Chebyshev polynomials T_0 through T_{n-1} at physical coordinate x
std::vector< double > evaluateChebyshevPolynomials(
        const double intervalStart,
        const double intervalEnd,
        const double evaluationPoint,
        const int maximumOrder )
{
    const double chebyshevX = determineChebyshevPolynomialIndependentVariable(
                intervalStart, intervalEnd, evaluationPoint );
    return evaluateChebyshevPolynomials( chebyshevX, maximumOrder );
}

//! Evaluate Chebyshev polynomials T_0 through T_{n-1} at normalized [-1,1] coordinate
std::vector< double > evaluateChebyshevPolynomials(
        const double independentChebyshevVariable,
        const int maximumOrder )
{
    std::vector< double > chebyshevPolynomials;
    chebyshevPolynomials.reserve( maximumOrder );

    if ( maximumOrder >= 1 )
    {
        chebyshevPolynomials.push_back( 1.0 );
    }

    if ( maximumOrder >= 2 )
    {
        chebyshevPolynomials.push_back( independentChebyshevVariable );
    }

    for ( int i = 2; i < maximumOrder; ++i )
    {
        chebyshevPolynomials.push_back(
            2.0 * independentChebyshevVariable * chebyshevPolynomials[ i - 1 ] - chebyshevPolynomials[ i - 2 ] );
    }

    return chebyshevPolynomials;
}

} // namespace interpolators

} // namespace tudat
