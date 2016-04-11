/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120926    E. Dekens         File created.
 *      121218    S. Billemont      Added output fuctions to display Legendre polynomial data,
 *                                  for debugging.
 *
 *    References
 *
 *    Notes
 *
 */

#include <sstream>
#include <stdexcept>

#include <boost/exception/all.hpp>

#include "Tudat/Mathematics/BasicMathematics/legendrePolynomials.h"

namespace tudat
{
namespace basic_mathematics
{

//! Get Legendre polynomial from cache when possible, and from direct computation otherwise.
void LegendreCache::update( const double polynomialParameter,
                            const LegendrePolynomialFunction legendrePolynomialFunction )
{
    for( int i = 0; i <= maximumDegree_; i++ )
    {
        for( int j = 0; ( ( j <= i ) && ( j <= maximumOrder_ ) ) ; j++ )
        {
            legendreValues_[ i * ( maximumOrder_ + 1 ) + j ] = legendrePolynomialFunction( i, j, polynomialParameter, this );
        }
    }
}

void LegendreCache::resetMaximumDegreeAndOrder( const int degree, const int order )
{
    maximumDegree_ = degree;
    maximumOrder_ = order;
    legendreValues_.resize( ( maximumDegree_ + 1 ) * ( maximumOrder_ + 1 ) );
    sinesOfLongitude_.resize( maximumOrder_ + 3 );
    cosinesOfLongitude_.resize( maximumOrder_ + 3 );
    referenceRadiusRatioPowers_.resize( maximumDegree_ + 3 );

    currentPolynomialParameter_ = -1.0E100;
    currentPolynomialParameterComplement_ = -1.0E100;
}

//! Get Legendre polynomial from cache when possible, and from direct computation otherwise.
double LegendreCache::getOrElseUpdate(
        const int degree, const int order, const double polynomialParameter,
        const LegendrePolynomialFunction& legendrePolynomialFunction )
{
    if( polynomialParameter != currentPolynomialParameter_ )
    {
        currentPolynomialParameter_ = polynomialParameter;
        currentPolynomialParameterComplement_ = std::sqrt( 1.0 - currentPolynomialParameter_ * currentPolynomialParameter_ );

        update( polynomialParameter, legendrePolynomialFunction );
    }

    if( degree > maximumDegree_ || order > maximumOrder_ )
    {
        std::cerr<<"Error when requesting legendre cache, maximum degree or order exceeded "<<
                   degree<<" "<<maximumDegree_<<" "<<order<<" "<<maximumOrder_<<std::endl;
    }
    else if( order > degree )
    {
        returnValue_ = 0.0;
    }
    else
    {
        returnValue_ = legendreValues_[ degree * ( maximumOrder_ + 1  ) + order ];
    }

    return returnValue_;

}

//! Compute unnormalized associated Legendre polynomial.
double computeLegendrePolynomial( const int degree,
                                  const int order,
                                  const double polynomialParameter,
                                  basic_mathematics::LegendreCache* legendreCache )
{
    // If degree or order is negative...
    if ( degree < 0 || order < 0 )
    {
        // Set error message.
        std::stringstream errorMessage;
        errorMessage << "Error: the Legendre polynomial of = " << degree << " and order = "
                     << order << " is undefined." << std::endl;

        // Throw a run-time error.
        boost::throw_exception( boost::enable_error_info( std::runtime_error(
               errorMessage.str( ) ) ) );
    }

    // Else if order is greater than degree...
    else if ( order > degree && degree >= 0 )
    {
        // Return zero.
        return 0.0;
    }

    // Else if order and degree are lower than 2...
    else if ( degree <= 1 && order <= 1 )
    {
        // Compute polynomial explicitly.
        return computeLegendrePolynomialExplicit( degree, order, polynomialParameter );
    }

    // Else if degree and order are sectoral...
    else if ( degree == order )
    {
        // Obtain polynomial of degree one and order one.
        const double degreeOneOrderOnePolynomial = legendreCache->getOrElseUpdate(
                    1, 1, polynomialParameter, legendrePolynomialFunction);

        // Obtain prior sectoral polynomial.
        const double priorSectoralPolynomial = legendreCache->getOrElseUpdate(
                    degree - 1, order - 1, polynomialParameter, legendrePolynomialFunction);

        // Compute polynomial.
        return computeLegendrePolynomialDiagonal(
                    degree, degreeOneOrderOnePolynomial, priorSectoralPolynomial );
    }

    // Else degree and order are zonal/tessoral...
    else
    {
        // Obtain prior degree polynomial.
        const double oneDegreePriorPolynomial = legendreCache->getOrElseUpdate(
                    degree - 1, order, polynomialParameter, legendrePolynomialFunction);

        // Obtain two degrees prior polynomial.
        const double twoDegreesPriorPolynomial = legendreCache->getOrElseUpdate(
                    degree - 2, order, polynomialParameter, legendrePolynomialFunction);

        // Compute polynomial.
        return computeLegendrePolynomialVertical( degree,
                                                  order,
                                                  polynomialParameter,
                                                  oneDegreePriorPolynomial,
                                                  twoDegreesPriorPolynomial );
    }
}

//! Compute geodesy-normalized associated Legendre polynomial.
double computeGeodesyLegendrePolynomial( const int degree,
                                         const int order,
                                         const double polynomialParameter,
                                         basic_mathematics::LegendreCache* geodesyLegendreCache )
{
    // If degree or order is negative...
    if ( degree < 0 || order < 0 )
    {
        // Set error message.
        std::stringstream errorMessage;
        errorMessage << "Error: the Legendre polynomial of = " << degree << " and order = "
                     << order << " is undefined." << std::endl;

        // Throw a run-time error.
        boost::throw_exception( boost::enable_error_info( std::runtime_error(
               errorMessage.str( ) ) ) );
    }

    // Else if order is greater than degree...
    else if ( order > degree && degree >= 0 )
    {
        // Return zero.
        return 0.0;
    }

    // Else if order and degree are lower than 2...
    else if ( degree <= 1 && order <= 1 )
    {
        // Compute polynomial explicitly.
        return computeGeodesyLegendrePolynomialExplicit( degree, order, polynomialParameter );
    }

    // Else if degree and order are sectoral...
    else if ( degree == order )
    {
        // Obtain polynomial of degree one and order one.
        double degreeOneOrderOnePolynomial = geodesyLegendreCache->getOrElseUpdate(
                    1, 1, polynomialParameter, geodesyNormalizedLegendrePolynomialFunction );

        // Obtain prior sectoral polynomial.
        double priorSectoralPolynomial = geodesyLegendreCache->getOrElseUpdate(
                    degree - 1, order - 1, polynomialParameter, geodesyNormalizedLegendrePolynomialFunction);

        // Compute polynomial.
        return computeGeodesyLegendrePolynomialDiagonal(
                    degree, degreeOneOrderOnePolynomial, priorSectoralPolynomial );
    }

    // Else degree and order are zonal/tessoral...
    else
    {
        // Obtain prior degree polynomial.
        double oneDegreePriorPolynomial = geodesyLegendreCache->getOrElseUpdate(
                    degree - 1, order, polynomialParameter, geodesyNormalizedLegendrePolynomialFunction);

        // Obtain two degrees prior polynomial.
        double twoDegreesPriorPolynomial = geodesyLegendreCache->getOrElseUpdate(
                    degree - 2, order, polynomialParameter, geodesyNormalizedLegendrePolynomialFunction);

        // Compute polynomial.
        return computeGeodesyLegendrePolynomialVertical( degree,
                                                         order,
                                                         polynomialParameter,
                                                         oneDegreePriorPolynomial,
                                                         twoDegreesPriorPolynomial );
    }
}

//! Compute derivative of unnormalized Legendre polynomial.
double computeLegendrePolynomialDerivative( const int order,
                                            const double polynomialParameter,
                                            const double currentLegendrePolynomial,
                                            const double incrementedLegendrePolynomial )
{
    // Return polynomial derivative.
    return incrementedLegendrePolynomial
            / std::sqrt( 1.0 - polynomialParameter * polynomialParameter )
            - static_cast< double >( order ) * polynomialParameter
            / ( 1.0 - polynomialParameter * polynomialParameter )
            * currentLegendrePolynomial;
}

//! Compute derivative of geodesy-normalized Legendre polynomial.
double computeGeodesyLegendrePolynomialDerivative( const int degree,
                                                   const int order,
                                                   const double polynomialParameter,
                                                   const double currentLegendrePolynomial,
                                                   const double incrementedLegendrePolynomial )
{
    // Compute normalization correction factor.
    double normalizationCorrection = std::sqrt( ( static_cast< double >( degree )
                                                  + static_cast< double >( order ) + 1.0 )
                                                * ( static_cast< double >( degree - order ) ) );

    // If order is zero apply multiplication factor.
    if ( order == 0 )
    {
        normalizationCorrection *= std::sqrt( 0.5 );
    }

    // Return polynomial derivative.
    return normalizationCorrection * incrementedLegendrePolynomial
            / std::sqrt( 1.0 - polynomialParameter * polynomialParameter )
            - static_cast< double >( order ) * polynomialParameter
            / ( 1.0 - polynomialParameter * polynomialParameter )
            * currentLegendrePolynomial;
}

//! Compute low degree/order unnormalized Legendre polynomial explicitly.
double computeLegendrePolynomialExplicit( const int degree,
                                          const int order,
                                          const double polynomialParameter )
{
    // If 0,0 term is requested return polynomial value.
    if ( degree == 0 && order == 0 )
    {
        return 1.0;
    }

    // Else if 1,0 term is requested return polynomial value.
    else if ( degree == 1 && order == 0 )
    {
        return polynomialParameter;
    }

    // Else if 1,1 term is requested return Legendre polynomial value.
    else if ( degree == 1 && order == 1 )
    {
        return std::sqrt( 1 - polynomialParameter * polynomialParameter );
    }

    // Else the requested term cannot be computed; throw a run-time error.
    else
    {
        // Set error message.
        std::stringstream errorMessage;
        errorMessage << "Error: computation of Legendre polynomial of = " << degree
                     << " and order = " << order << " is not supported." << std::endl;

        // Throw a run-time error.
        boost::throw_exception( boost::enable_error_info( std::runtime_error(
               errorMessage.str( ) ) ) );
    }
}

//! Compute low degree/order geodesy-normalized Legendre polynomials explicitly.
double computeGeodesyLegendrePolynomialExplicit( const int degree,
                                                 const int order,
                                                 const double polynomialParameter )
{
    // If 0,0 term is requested return Legendre polynomial value.
    if ( degree == 0 && order == 0 )
    {
        return 1.0;
    }

    // Else if 1,0 term is requested return polynomial value.
    else if ( degree == 1 && order == 0 )
    {
        return std::sqrt( 3.0 ) * polynomialParameter;
    }

    // Else if 1,1 term is requested return polynomial value.
    else if ( degree == 1 && order == 1 )
    {
        return std::sqrt( 3.0 - 3.0 * polynomialParameter * polynomialParameter );
    }

    // Else the requested term cannot be computed; throw a run-time error.
    else
    {
        // Set error message.
        std::stringstream errorMessage;
        errorMessage << "Error: computation of Legendre polynomial of = " << degree
                     << " and order = " << order << " is not supported." << std::endl;

        // Throw a run-time error.
        boost::throw_exception( boost::enable_error_info( std::runtime_error(
               errorMessage.str( ) ) ) );
    }
}

//! Compute unnormalized Legendre polynomial through sectoral recursion.
double computeLegendrePolynomialDiagonal( const int degree,
                                          const double degreeOneOrderOnePolynomial,
                                          const double priorSectoralPolynomial )
{
    // Return polynomial.
    return ( 2.0 * static_cast< double >( degree ) - 1.0 )
            * degreeOneOrderOnePolynomial * priorSectoralPolynomial;

}

//! Compute geodesy-normalized Legendre polynomial through sectoral recursion.
double computeGeodesyLegendrePolynomialDiagonal( const int degree,
                                                 const double degreeOneOrderOnePolynomial,
                                                 const double priorSectoralPolynomial )
{
    // Return polynomial.
    return std::sqrt( ( 2.0 * static_cast< double >( degree ) + 1.0 )
                      / ( 6.0 * static_cast< double >( degree ) ) )
            * degreeOneOrderOnePolynomial * priorSectoralPolynomial;
}

//! Compute unnormalized Legendre polynomial through degree recursion.
double computeLegendrePolynomialVertical( const int degree,
                                          const int order,
                                          const double polynomialParameter,
                                          const double oneDegreePriorPolynomial,
                                          const double twoDegreesPriorPolynomial )
{
    // Return polynomial.
    return ( ( 2.0 * static_cast< double >( degree ) - 1.0 ) * polynomialParameter
             * oneDegreePriorPolynomial - ( static_cast< double >( degree + order ) - 1.0 )
             * twoDegreesPriorPolynomial ) / ( static_cast< double >( degree - order ) );
}

//! Compute geodesy-normalized Legendre polynomial through degree recursion.
double computeGeodesyLegendrePolynomialVertical( const int degree,
                                                 const int order,
                                                 const double polynomialParameter,
                                                 const double oneDegreePriorPolynomial,
                                                 const double twoDegreesPriorPolynomial )
{
    // Return polynomial.
    return std::sqrt( ( 2.0 * static_cast< double >( degree ) + 1.0 )
                      / ( ( static_cast< double >( degree + order ) ) * ( static_cast< double >( degree - order ) ) ) )
            * ( std::sqrt( 2.0 * static_cast< double >( degree ) - 1.0 ) * polynomialParameter * oneDegreePriorPolynomial
                - std::sqrt( ( static_cast< double >( degree + order ) - 1.0 )
                             * ( static_cast< double >( degree - order ) - 1.0 )
                             / ( 2.0 * static_cast< double >( degree ) - 3.0 ) )
                * twoDegreesPriorPolynomial );
}

} // namespace basic_mathematics
} // namespace tudat
