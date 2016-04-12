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
 *      130121    K. Kumar          Added shared-ptr typedefs.
 *
 *    References
 *      Eberly, D. Spherical Harmonics. Help documentation of Geometric Tools, 2008. Available at
 *        URL http://www.geometrictools.com/Documentation/SphericalHarmonics.pdf. Last access:
 *        09-09-2012.
 *      Heiskanen, W.A., Moritz, H. Physical geodesy. Freeman, 1967.
 *      Holmes, S.A., Featherstone, W.E. A unified approach to the Clenshaw summation and the
 *        recursive computation of very high degree and order normalised associated Legendre
 *        functions. Journal of Geodesy, 76(5):279-299, 2002.
 *      Vallado, D. and McClain, W. Fundamentals of astrodynammics and applications. Microcosm
 *        Press, 2001.
 *      Weisstein, E.W. Associated Legendre Polynomial, 2012.
 *        URL http://mathworld.wolfram.com/AssociatedLegendrePolynomial.html. Last access:
 *        12-09-2012.
 *
 *    Notes
 *      For information on how the caching mechanism works, please contact S. Billemont
 *      (S.Billemont@studelft.tudelft.nl).
 *
 */

#ifndef TUDAT_LEGENDRE_POLYNOMIALS_H
#define TUDAT_LEGENDRE_POLYNOMIALS_H

#include <cstddef>
#include <iostream>

#include <boost/bind.hpp>

#include <boost/circular_buffer.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <boost/enable_shared_from_this.hpp>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{
namespace basic_mathematics
{

//! Class for creating and accessing a back-end cache of Legendre polynomials.
class LegendreCache
{
private:


public:

    //! Define Legendre polynomial function pointer.
    typedef boost::function< double ( int, int, LegendreCache& ) > LegendrePolynomialFunction;

    //! Initialize LegendreCache instance.
    LegendreCache( const bool useGeodesyNormalization = 1 );

    LegendreCache( const int maximumDegree, const int maximumOrder, const bool useGeodesyNormalization = 1 );

    void resetMaximumDegreeAndOrder( const int degree, const int order );


    void update( const double polynomialParameter );

    double getCurrentPolynomialParameter( )
    {
        return currentPolynomialParameter_;
    }

    double getCurrentPolynomialParameterComplement( )
    {
        return currentPolynomialParameterComplement_;
    }


    //! Get Legendre polynomial value from either cache or from computation.
    /*!
    * \param degree Degree of requested Legendre polynomial.
    * \param order Order of requested Legendre polynomial.
    * \param polynomialParameter Free variable  of requested Legendre polynomial.
    * \param legendrePolynomialFunction Function which takes degree, order and
    *          polynomialParameter as arguments. The function must return the corresponding
    *          Legendre polynomial value.
    * \return Legendre polynomial value.
    */
    double getLegendrePolynomial( const int degree, const int order );

    int getMaximumDegree( ){ return maximumDegree_; }

    int getMaximumOrder( ){ return maximumOrder_; }

    bool getUseGeodesyNormalization( )
    {
        return useGeodesyNormalization_;
    }


private:
    int maximumDegree_;

    int maximumOrder_;

    double currentPolynomialParameter_;

    double currentPolynomialParameterComplement_;

    std::vector< double > legendreValues_;

    double returnValue_ ;

    LegendrePolynomialFunction legendrePolynomialFunction_;

    bool useGeodesyNormalization_;

};



//! Compute unnormalized associated Legendre polynomial.
/*!
 * This function returns an unnormalized associated Legendre polynomial \f$ P _{ n, m }( u ) \f$
 * with degree \f$ n \f$, order \f$ m \f$ and polynomial parameter \f$ u \f$.
 * \f$ P _{ n, m }( u ) \f$ obeys the definition given by Vallado [2001]:
 * \f[
 *     P _{ n, m } ( u ) = ( 1 - u ^ 2 ) ^ { m / 2 } \frac{ d ^ m }{ du ^ m } \left[ P_n ( u )
 *     \right] \textrm{ for }n \geq 0, 0 \leq m \leq n
 * \f]
 * in which \f$ P_{ n, m }( u ) \f$ is the ordinary Legendre polynomial with degree \f$ n \f$,
 * order \f$ m \f$ and polynomial parameter \f$ u \f$.
 *
 * For \f$ n \geq 0, m \geq n \f$ the \f$ P _{ n, m }( u ) \f$ has been defined here as
 * follows:
 * \f[
 *     P _{ n, m } ( u ) = 0 \textrm{ for }n \geq 0, m \geq n
 * \f]
 *
 * This function has been optimized for repeated calls with varying 'degree' and 'order' arguments
 * (but with identical 'polynomialParameter' argument). To this end the function maintains a
 * back-end cache with intermediate results which is automatically carried over between calls.
 * \param degree Degree of requested Legendre polynomial.
 * \param order Order of requested Legendre polynomial.
 * \param polynomialParameter Free variable  of requested Legendre polynomial.
 * \return Unnormalized Legendre polynomial.
*/
double computeLegendrePolynomial( const int degree,
                                  const int order,
                                  LegendreCache& legendreCache );

double computeLegendrePolynomial( const int degree,
                                  const int order,
                                  const double legendreParameter );

//! Compute geodesy-normalized associated Legendre polynomial.
/*!
 * This function returns a normalized associated Legendre polynomial
 * \f$ \bar{ P }_{ n, m }( u ) \f$ with degree \f$ n \f$, order \f$ m \f$ and polynomial
 * parameter \f$ u \f$. The normalization obeys the definition:
 * \f[
 *     \bar{ P }_{ n, m } ( u ) = \Pi_{ n, m } P_{ n, m } ( u )
 * \f]
 * in which \f$ \Pi_{ n, m } \f$ is the normalization factor which is commonly used in geodesy and
 * is given by Heiskanen & Moritz [1967] as:
 * \f[
 *     \Pi_{ n, m } = \sqrt{ \frac{ ( 2 - \delta_{ 0, m } ) ( 2 n + 1 ) ( n - m )! }
 *     { ( n + m )! } }
 * \f]
 * in which \f$ n \f$ is the degree, \f$ m \f$ is the order and \f$ \delta_{ 0, m } \f$ is the
 * Kronecker delta.

 * \f$ P _{ n, m }( u ) \f$ obeys the definition given by Vallado [2001]:
 * \f[
 *     P _{ n, m } ( u ) = ( 1 - u ^ 2 ) ^ { m / 2 } \frac{ d ^ m }{ du ^ m } \left[ P_n ( u )
 *     \right]
 * \f]
 * in which \f$ P_{ n, m }( u ) \f$ is the ordinary Legendre polynomial with degree \f$ n \f$,
 *  order \f$ m \f$ and polynomial parameter \f$ u \f$.
 *
 * For \f$ n \geq 0, m \geq n \f$ the \f$ P _{ n, m }( u ) \f$ has been defined here as
 * follows:
 * \f[
 *     P _{ n, m } ( u ) = 0 \textrm{ for }n \geq 0, m \geq n
 * \f]
 *
 * This function has been optimized for repeated calls with varying 'degree' and 'order' arguments
 * (but with identical 'polynomialParameter' argument). To this end the function maintains a
 * back-end cache with intermediate results which is automatically carried over between calls.
 * \param degree Degree of requested Legendre polynomial.
 * \param order Order of requested Legendre polynomial.
 * \param polynomialParameter Free variable of requested Legendre polynomial.
 * \return Geodesy-normalized Legendre polynomial.
*/
double computeGeodesyLegendrePolynomial( const int degree,
                                         const int order,
                                         LegendreCache& geodesyLegendreCache );

double computeGeodesyLegendrePolynomial( const int degree,
                                         const int order,
                                         const double legendreParameter );

//! Compute derivative of unnormalized Legendre polynomial.
/*!
 * The derivative is computed as:
 * \f[
 *     \frac{ \mathrm{ d } P_{ n, m } ( u ) }{ \mathrm{ d } u } = \frac{ \sqrt{ 1 - u ^ 2 }
 *     P_{ n, m + 1 }( u ) - m u P_{ n, m }( u ) }{ 1 - u ^ 2 }
 * \f]
 * in which \f$ n \f$ is the degree, \f$ m \f$ is the order, \f$ u \f$ is the polynomial parameter,
 * \f$ P_{ n, m } ( u ) \f$ is an associated Legendre polynomial, and \f$ P_{ n, m + 1 } ( u ) \f$
 * is an associated Legendre polynomial.
 * \param order Order of requested Legendre polynomial derivative.
 * \param polynomialParameter Free variable  of requested Legendre polynomial derivative.
 * \param currentLegendrePolynomial Unnormalized Legendre polynomial with the same degree, order
 *          and polynomial parameter as the requested Legendre polynomial derivative.
 * \param incrementedLegendrePolynomial Unnormalized Legendre polynomial with the same degree and
 *          polynomial parameter as the requested Legendre polynomial derivative, but with an order
 *          of one more.
 * \return Unnormalized Legendre polynomial derivative with respect to the polynomial parameter.
*/
double computeLegendrePolynomialDerivative( const int order,
                                            const double polynomialParameter,
                                            const double currentLegendrePolynomial,
                                            const double incrementedLegendrePolynomial );

//! Compute derivative of geodesy-normalized Legendre polynomial.
/*!
 * The derivative is computed as:
 * \f[
 *     \frac{ \mathrm{ d } \bar{ P }_{ n, m } ( u ) }{ \mathrm{ d } u } = \sqrt{ \frac{
 *     ( n + m + 1 )( n - m ) }{ ( 1 + \delta_{ 0, m } ) ( 1 - u^2 ) } } \bar{ P }_{ n, m + 1 }
 *     - \frac{ m u }{ 1 - u^2 } \bar{ P }_{ n, m }
 * \f]
 * in which \f$ n \f$ is the degree, \f$ m \f$ is the order, \f$ u \f$ is the polynomial parameter,
 * \f$ \delta_{ 0, m } \f$ is the Kronecker delta, \f$ P_{ n, m } ( u ) \f$ is an associated
 * Legendre polynomial, and \f$ P_{ n, m + 1 } ( u ) \f$ is an associated Legendre polynomial.
 * \param degree Degree of requested Legendre polynomial derivative.
 * \param order Order of requested Legendre polynomial derivative.
 * \param polynomialParameter Free variable  of requested Legendre polynomial derivative.
 * \param currentLegendrePolynomial Geodesy-normalized Legendre polynomial with the same degree, order
 *          and polynomial parameter as the requested Legendre polynomial derivative.
 * \param incrementedLegendrePolynomial Geodesy-normalized Legendre polynomial with the same degree and
 *          polynomial parameter as the requested Legendre polynomial derivative, but with an order
 *          of one more.
 * \return Geodesy-normalized Legendre polynomial derivative with respect to the polynomial parameter.
*/
double computeGeodesyLegendrePolynomialDerivative( const int degree,
                                                   const int order,
                                                   const double polynomialParameter,
                                                   const double currentLegendrePolynomial,
                                                   const double incrementedLegendrePolynomial );

//! Compute low degree/order unnormalized Legendre polynomial explicitly.
/*!
 * The associated Legendre polynomial \f$ P_{ n, m }(u) \f$ with degree \f$ n \f$, order
 * \f$ m \f$ and polynomial parameter \f$ u \f$ is calculated using explicit formulas given by
 * Eberly [2008]:
 * \f{eqnarray*}{
 *     P_{ 0, 0 }( u ) = 1 \\
 *     P_{ 1, 0 }( u ) = u \\
 *     P_{ 1, 1 }( u ) = \sqrt{ 1 - u^2 }
 * \f}
 * Calculation up to \f$ n = 1 \f$ and \f$ m = 1 \f$ is supported by this function.
 * \param degree Degree of requested Legendre polynomial.
 * \param order Order of requested Legendre polynomial.
 * \param polynomialParameter Free variable of requested Legendre polynomial.
 * \return Unnormalized Legendre polynomial.
*/
double computeLegendrePolynomialExplicit( const int degree,
                                          const int order,
                                          const double polynomialParameter );

//! Compute low degree/order geodesy-normalized Legendre polynomial explicitly.
/*!
 * The normalized associated Legendre polynomial \f$ \bar{ P }_{ n, m }(u) \f$ with degree
 * \f$ n \f$, order \f$ m \f$ and polynomial parameter \f$ u \f$ is calculated using explicit
 * formulas:
 * \f{eqnarray*}{
 *     \bar{ P }_{ 0, 0 }( u ) = 1 \\
 *     \bar{ P }_{ 1, 0 }( u ) = \sqrt{ 3 } u \\
 *     \bar{ P }_{ 1, 1 }( u ) = \sqrt{ 3 - 3 u^2 }
 * \f}
 * Calculation up to \f$ n = 1 \f$ and \f$ m = 1 \f$ is supported by this function.
 * \param degree Degree of requested Legendre polynomial.
 * \param order Order of requested Legendre polynomial.
 * \param polynomialParameter Free variable of requested Legendre polynomial.
 * \return Geodesy-normalized Legendre polynomial.
*/
double computeGeodesyLegendrePolynomialExplicit( const int degree,
                                                 const int order,
                                                 const double polynomialParameter );

//! Compute unnormalized Legendre polynomial through sectoral recursion.
/*!
 * The associated Legendre polynomial \f$ P_{ n, m }( u ) \f$ with degree \f$ n \f$, order
 * \f$ m \f$ and polynomial parameter \f$ u \f$ is calculated through degree recursion
 * as given by Vallado [2001]:
 * \f[
 *     P_{ n, m }( u ) = ( 2 n - 1 ) P_{ 1, 1 }( u ) P_{ n - 1, m - 1 }( u )
 * \f]
 * \param degree Degree of requested Legendre polynomial.
 * \param degreeOneOrderOnePolynomial Unnormalized Legendre polynomial with the same polynomial
 *          parameter as the requested Legendre polynomial, but with a degree and order of one.
 * \param priorSectoralPolynomial Unnormalized Legendre polynomial with the same polynomial
 *          parameter as the requested Legendre polynomial, but with a degree and order of one
 *          less.
 * \return Unnormalized Legendre polynomial.
*/
double computeLegendrePolynomialDiagonal( const int degree,
                                          const double degreeOneOrderOnePolynomial,
                                          const double priorSectoralPolynomial );

//! Compute geodesy-normalized Legendre polynomial through sectoral recursion.
/*!
 * The associated Legendre polynomial \f$ P_{ n, m }( u ) \f$ with degree \f$ n \f$, order
 * \f$ m \f$ and polynomial parameter \f$ u \f$ is calculated through degree recursion
 * as given by Holmes & Featherstone [2002]:
 * \f[
 *     \bar{ P }_{ n, m }( u ) = \sqrt{ \frac{ 2 n - 1 }{ 6 n } } \bar{ P }_{ 1, 1 }( u )
 *     \bar{ P }_{ n - 1, m - 1 }( u )
 * \f]
 * \param degree Degree of requested Legendre polynomial.
 * \param degreeOneOrderOnePolynomial Unnormalized Legendre polynomial with the same polynomial
 *          parameter as the requested Legendre polynomial, but with a degree and order of one.
 * \param priorSectoralPolynomial Legendre polynomial with the same polynomial parameter as the
 *          requested Legendre polynomial, but with a degree and order of one less.
 * \return Geodesy-normalized Legendre polynomial.
*/
double computeGeodesyLegendrePolynomialDiagonal( const int degree,
                                                 const double degreeOneOrderOnePolynomial,
                                                 const double priorSectoralPolynomial );

//! Compute unnormalized Legendre polynomial through degree recursion.
/*!
 * The associated Legendre polynomial \f$ P_{ n, m }( u ) \f$ with degree \f$ n \f$, order
 * \f$ m \f$ and polynomial parameter \f$ u \f$ is calculated through degree recursion
 * as given by Weisstein [2012]:
 * \f[
 *     P_{ n, m }( u ) = \frac{ ( 2 n - 1 ) u P_{ n - 1, m }( u )
 *     - ( n + m - 1 ) P_{ n - 2, m }( u ) }{ n - m }
 * \f]
 * \param degree Degree of requested Legendre polynomial.
 * \param order Order of requested Legendre polynomial.
 * \param polynomialParameter Free variable of requested Legendre polynomial.
 * \param oneDegreePriorPolynomial Unnormalized Legendre polynomial with the same order and
 *          polynomial parameter as the requested Legendre polynomial, but with a degree of one
 *          less.
 * \param twoDegreesPriorPolynomial Geodesy-normalized Legendre polynomial with the same order
 *          and polynomial parameter as the requested Legendre polynomial, but with a degree of
 *          two less.
 * \return Unnormalized Legendre polynomial.
*/
double computeLegendrePolynomialVertical( const int degree,
                                          const int order,
                                          const double polynomialParameter,
                                          const double oneDegreePriorPolynomial,
                                          const double twoDegreesPriorPolynomial );

//! Compute geodesy-normalized Legendre polynomial through degree recursion.
/*!
 * The normalized associated Legendre polynomial \f$ \bar{ P }_{ n, m }( u ) \f$ with degree
 * \f$ n \f$, order \f$ m \f$ and polynomial parameter \f$ u \f$ is calculated through degree
 * recursion as given by Holmes & Featherstone [2002]:
 * \f[
 *     \bar{ P }_{ n, m }( u ) = \sqrt{ \frac{ ( 2 n - 1 )( 2 n + 1 ) }{ ( n - m ) ( n + m ) } }
 *     u \bar{ P }_{ n - 1, m }( u ) - \sqrt{ \frac{ ( 2 n + 1 )( n + m - 1 ) ( n - m - 1 ) }
 *     { ( n - m)( n + m )( 2 n - 3 ) } } \bar{ P }_{ n - 2, m }( u )
 * \f]
 * \param degree Degree of requested Legendre polynomial.
 * \param order Order of requested Legendre polynomial.
 * \param polynomialParameter Free variable of requested Legendre polynomial.
 * \param oneDegreePriorPolynomial Geodesy-normalized Legendre polynomial with the same polynomial
 *          parameter as the requested Legendre polynomial, but with a degree of one less.
 * \param twoDegreesPriorPolynomial Geodesy-normalized Legendre polynomial with the same polynomial
 *          parameter as the requested Legendre polynomial, but with a degree of two less.
 * \return Geodesy-normalized Legendre polynomial.
*/
double computeGeodesyLegendrePolynomialVertical( const int degree,
                                                 const int order,
                                                 const double polynomialParameter,
                                                 const double oneDegreePriorPolynomial,
                                                 const double twoDegreesPriorPolynomial );



static const LegendreCache::LegendrePolynomialFunction geodesyNormalizedLegendrePolynomialFunction =
        boost::bind( static_cast< double(&)( const int,  const int,  LegendreCache&  )>( &computeGeodesyLegendrePolynomial ), _1, _2, _3 );

static const LegendreCache::LegendrePolynomialFunction regularLegendrePolynomialFunction =
        boost::bind( static_cast< double(&)( const int,  const int,  LegendreCache& )>( &computeLegendrePolynomial ), _1, _2, _3 );

} // namespace basic_mathematics
} // namespace tudat

#endif // TUDAT_LEGENDRE_POLYNOMIALS_H
