/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
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

#include <Eigen/Core>




#include <boost/circular_buffer.hpp>
#include <functional>
#include <memory>
#include <vector>



#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace basic_mathematics
{

//! Class for creating and accessing a back-end cache of Legendre polynomials.
class LegendreCache
{

public:

    //! Define Legendre polynomial function pointer.
    typedef std::function< double ( int, int, LegendreCache& ) > LegendrePolynomialFunction;

    //! Default constructor, initializes cache object with 0 maximum degree and order.
    /*!
     * Default constructor, initializes cache object with 0 mazimum degree and order.
     * \param useGeodesyNormalization Parameter defining whether the cache is used for a normalized or unnormalized
     * gravity field.
     */
    LegendreCache( const bool useGeodesyNormalization = 1 );

    //! Constructor
    /*!
     * Constructor
     * \param maximumDegree Maximum degree to which to update cache
     * \param maximumOrder Maximum order to which to update cache
     * \param useGeodesyNormalization Parameter defining whether the cache is used for a normalized or unnormalized
     * gravity field.
     */
    LegendreCache( const int maximumDegree, const int maximumOrder, const bool useGeodesyNormalization = 1 );

    ~LegendreCache( ){ }
    //! Update maximum degree and order of cache
    /*!
     * Update maximum degree and order of cache
     * \param maximumDegree Maximum degree to which to update cache
     * \param maximumOrder Maximum order to which to update cache
     */
    void resetMaximumDegreeAndOrder( const int maximumDegree, const int maximumOrder );

    //! Update cache with new polynomial parameter (sine of latitude)
    /*!
     * Update cache with new polynomial parameter (sine of latitude).
     * \param polynomialParameter Parameter used as input argument for Legendre polynomials, in astro
     * applications, this is typically the sine of the body-fixed latitude.
     */
   void update( const double polynomialParameter );

    //! Function to return the current polynomial parameter (typically sine of latitude)
    /*!
     * Function to return the current polynomial parameter (typically sine of latitude)
     * \return Current polynomial parameter
     */
    double getCurrentPolynomialParameter( )
    {
        return currentPolynomialParameter_;
    }

    //! Function to return the complement to the current polynomial parameter (typically cosine of latitude)
    /*!
     *  Function to return the complement to the current polynomial parameter (typically cosine of latitude)
     * \return Complement to the current polynomial parameter
     */
    double getCurrentPolynomialParameterComplement( )
    {
        return currentPolynomialParameterComplement_;
    }

    //! Get Legendre polynomial value from the cache.
    /*!
    * Get Legendre polynomial value from the cache, as computed by last call to update function.
    * \param degree Degree of requested Legendre polynomial.
    * \param order Order of requested Legendre polynomial.
    * \return Legendre polynomial value.
    */
    double getLegendrePolynomial( const int degree, const int order );

    //! Get first derivative of Legendre polynomial value from the cache.
    /*!
    * Get first derivative of Legendre polynomial value from the cache, as computed by last call to update function.
    * \param degree Degree of requested Legendre polynomial.
    * \param order Order of requested Legendre polynomial.
    * \return First derivative of Legendre polynomial value.
    */
    double getLegendrePolynomialDerivative( const int degree, const int order );

    //! Get second derivative of Legendre polynomial value from the cache.
    /*!
    * Get second derivative of Legendre polynomial value from the cache, as computed by last call to update function.
    * \param degree Degree of requested Legendre polynomial.
    * \param order Order of requested Legendre polynomial.
    * \return Second derivative of Legendre polynomial value.
    */
    double getLegendrePolynomialSecondDerivative( const int degree, const int order );

    //! Function to get the maximum degree of cache.
    /*!
     * Function to get the maximum degree of cache
     * \return Maximum degree of cache.
     */
    int getMaximumDegree( )
    {
        return maximumDegree_;
    }

    //! Function to get the maximum order of cache.
    /*!
     * Function to get the maximum order of cache.
     * \return Maximum order of cache.
     */
    int getMaximumOrder( )
    {
        return maximumOrder_;
    }

    //! Function to get whether the Legendre polynomials are geodesy-normalized or unnormalized
    /*!
     * Function to get whether the Legendre polynomials are geodesy-normalized or unnormalized
     * \return Boolean denoting to get whether the Legendre polynomials are geodesy-normalized or unnormalized
     */
    bool getUseGeodesyNormalization( )
    {
        return useGeodesyNormalization_;
    }

    //! Function to reset whether the first derivatives are to be computed when calling update function
    /*!
     * Function to reset whether the first derivatives are to be computed when calling update function
     * \param computeFirstDerivatives Boolean denoting whether the first derivatives of the Legendre polynomials are
     * to be computed when calling update function.
     */
    void setComputeFirstDerivatives( const bool computeFirstDerivatives )
    {
        computeFirstDerivatives_ = computeFirstDerivatives;
        currentPolynomialParameter_ = TUDAT_NAN;
    }

    //! Function to reset whether the second derivatives are to be computed when calling update function
    /*!
     * Function to reset whether the second derivatives are to be computed when calling update function
     * \param computeSecondDerivatives Boolean denoting whether the second derivatives of the Legendre polynomials are
     * to be computed when calling update function.
     */
    void setComputeSecondDerivatives( const bool computeSecondDerivatives )
    {
        computeSecondDerivatives_ = computeSecondDerivatives;
        currentPolynomialParameter_ = TUDAT_NAN;
    }

    double getVerticalLegendreValuesComputationMultipliersOne( const int degree, const int order );

    double getVerticalLegendreValuesComputationMultipliersTwo( const int degree, const int order );

private:

    //! Maximum degree of cache.
    int maximumDegree_;

    //! Maximum order of cache.
    int maximumOrder_;

    //! Current polynomial parameter (sine of latitude).
    double currentPolynomialParameter_;

    //! Current 'complement' to polynomial parameter (cosine of latitude).
    double currentPolynomialParameterComplement_;

    double currentOneOverPolynomialParameterComplement_;

    //! List of current values of Legendre polynomials at degree and order (n,m)
    /*!
     * List of current values of Legendre polynomials at degree and order (n,m). The corresponding polynomial is at entry
     * n * ( maximumOrder_ + 1 ) + m.
     */
    std::vector< double > legendreValues_;

    std::vector< double > verticalLegendreValuesComputationMultipliersOne_;

    std::vector< double > verticalLegendreValuesComputationMultipliersTwo_;

    //! List of current values of first derivatives of Legendre polynomials at degree and order (n,m)
    /*!
     * List of current values of first derivatives of Legendre polynomials at degree and order (n,m).
     * The corresponding polynomial is at entry n * ( maximumOrder_ + 1 ) + m.
     */
    std::vector< double > legendreDerivatives_;

    //! List of current values of second derivatives of Legendre polynomials at degree and order (n,m)
    /*!
     * List of current values of second derivatives of Legendre polynomials at degree and order (n,m).
     * The corresponding polynomial is at entry n * ( maximumOrder_ + 1 ) + m.
     */
    std::vector< double > legendreSecondDerivatives_;

    //! Function from which to compute the Legendre polynomials.
    LegendrePolynomialFunction legendrePolynomialFunction_;

    //! Boolean denoting whether the Legendre polynomials are geodesy-normalized or unnormalized
    bool useGeodesyNormalization_;

    //! Vector of ratio of reference radius over current radius to power i, with i the entry in the vector.
    std::vector< double > referenceRadiusRatioPowers_;

    //! Prec-computed normalization factors that are to be used for computation fo Legendre polynomial derivative
    std::vector< double > derivativeNormalizations_;

    //! Boolean denoting whether the first derivatives of the Legendre polynomials are to be computed when calling
    //! update function.
    bool computeFirstDerivatives_{true};

    //! Boolean denoting whether the second derivatives of the Legendre polynomials are to be computed when calling
    //! update function.
    bool computeSecondDerivatives_;


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
 * \param legendreCache Legendre cache from which to retrieve Legendre polynomial.
 * \return Unnormalized Legendre polynomial.
*/
double computeLegendrePolynomialFromCache( const int degree,
                                           const int order,
                                           LegendreCache& legendreCache );


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
 * \param legendreParameter Free variable  of requested Legendre polynomial.
 * \return Unnormalized Legendre polynomial.
*/
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
 * \param geodesyLegendreCache Legendre cache from which to retrieve Legendre polynomial.
 * \return Geodesy-normalized Legendre polynomial.
*/
double computeGeodesyLegendrePolynomialFromCache( const int degree,
                                                  const int order,
                                                  LegendreCache& geodesyLegendreCache );

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
 * \param legendreParameter Free variable of requested Legendre polynomial.
 * \return Geodesy-normalized Legendre polynomial.
*/
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

//! Compute derivative of geodesy-normalized associated Legendre polynomial.
/*!
 * Compute derivative of geodesy-normalized associated Legendre polynomial.
 * \param degree Degree of requested Legendre polynomial derivative.
 * \param order Order of requested Legendre polynomial derivative.
 * \param polynomialParameter Free variable  of requested Legendre polynomial derivative.
 * \param currentLegendrePolynomial Geodesy-normalized Legendre polynomial with the same degree, order
 *          and polynomial parameter as the requested Legendre polynomial derivative.
 * \param incrementedLegendrePolynomial Geodesy-normalized Legendre polynomial with the same degree and
 *          polynomial parameter as the requested Legendre polynomial derivative, but with an order
 *          of one more.
 * \param normalizationCorrection Pre-computed scaling term used for part of computations.
 * \return Geodesy-normalized Legendre polynomial derivative with respect to the polynomial parameter.
*/
double computeGeodesyLegendrePolynomialDerivative( const int order,
                                                   const double polynomialParameter,
                                                   const double oneOverPolynomialParameterComplement,
                                                   const double currentLegendrePolynomial,
                                                   const double incrementedLegendrePolynomial,
                                                   const double normalizationCorrection );


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

//! Compute second derivative of geodesy-normalized associated Legendre polynomial.
/*!
 * Compute second derivative of geodesy-normalized associated Legendre polynomial.
 * \param degree Degree of requested Legendre polynomial derivative.
 * \param order Order of requested Legendre polynomial derivative.
 * \param polynomialParameter Free variable  of requested Legendre polynomial derivative.
 * \param currentLegendrePolynomial Geodesy-normalized Legendre polynomial with the same degree, order
 *          and polynomial parameter as the requested Legendre polynomial derivative.
 * \param incrementedLegendrePolynomial Geodesy-normalized Legendre polynomial with the same degree and
 *          polynomial parameter as the requested Legendre polynomial derivative, but with an order
 *          of one more.
 * \param currentLegendrePolynomialDerivative Geodesy-normalized derivative of Legendre polynomial with the same degree,
 *          order and polynomial parameter as the requested Legendre polynomial derivative.
 * \param incrementedLegendrePolynomialDerivative Geodesy-normalized derivative of Legendre polynomial with the same degree
 *          and polynomial parameter as the requested Legendre polynomial derivative, but with an order
 *          of one more.
 * \param normalizationCorrection Pre-computed scaling term used for part of computations.
 * \return Geodesy-normalized Legendre polynomial derivative with respect to the polynomial parameter.
*/
double computeGeodesyLegendrePolynomialSecondDerivative( const int order,
                                                         const double polynomialParameter,
                                                         const double oneOverPolynomialParameterComplement,
                                                         const double currentLegendrePolynomial,
                                                         const double incrementedLegendrePolynomial,
                                                         const double currentLegendrePolynomialDerivative,
                                                         const double incrementedLegendrePolynomialDerivative,
                                                         const double normalizationCorrection );

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

double computeGeodesyLegendrePolynomialVertical( const int degree,
                                                 const int order,
                                                 const double polynomialParameter,
                                                 const double firstMultiplier,
                                                 const double secondMultiplier,
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


//! Predefine boost function for geodesy-normalized Legendre polynomial.
static const LegendreCache::LegendrePolynomialFunction geodesyNormalizedLegendrePolynomialFunction =
        std::bind( &computeGeodesyLegendrePolynomialFromCache, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 );

//! Function to calculate the normalization factor for Legendre polynomials to geodesy-normalized.
/*!
 *  Function to calculate the normalization factor K_nm to convert regular Legendre polynomials P_nm
 *  to geodesy-normalized Legendre polynomials bar{P}_nm, so that bar{P}_nm=K_nm*P_nm.
 *  Note that the same conversion for spherical harmonic coefficients from regular to geodesy-
 *  normalized requires the inverse factor: 1/K_nm: C_nm=K_nm* bar{C}_nm.
 * \param degree Degree of Legendre polynomial
 * \param order Order of Legendre polynomial
 * \return Normalization factor K_nm
 */
double calculateLegendreGeodesyNormalizationFactor( const int degree, const int order );

//! Predefine boost function for unnormalized Legendre polynomial.
const LegendreCache::LegendrePolynomialFunction regularLegendrePolynomialFunction =
        std::bind( &computeLegendrePolynomialFromCache, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3 );

//! Function to convert unnormalized to geodesy-normalized (4-pi normalized) spherical harmonic coefficients
/*!
 * Function to convert unnormalized to geodesy-normalized (4-pi normalized) spherical harmonic coefficients
 * \param unnormalizedCosineCoefficients Original, unnormalized, cosine coefficients
 * \param unnormalizedSineCoefficients Original, unnormalized, sine coefficients
 * \param normalizedCosineCoefficients Transformed, normalized, cosine coefficients (returned by reference)
 * \param normalizedSineCoefficients Transformed, normalized, sine coefficients (returned by reference)
 */
void convertUnnormalizedToGeodesyNormalizedCoefficients(
        const Eigen::MatrixXd& unnormalizedCosineCoefficients,
        const Eigen::MatrixXd& unnormalizedSineCoefficients,
        Eigen::MatrixXd& normalizedCosineCoefficients,
        Eigen::MatrixXd& normalizedSineCoefficients );

std::pair< Eigen::MatrixXd, Eigen::MatrixXd > convertUnnormalizedToGeodesyNormalizedCoefficients(
        const Eigen::MatrixXd& unnormalizedCosineCoefficients,
        const Eigen::MatrixXd& unnormalizedSineCoefficients );

//! Function to convert geodesy-normalized (4-pi normalized) to unnormalized spherical harmonic coefficients
/*!
 * Function to convert geodesy-normalized (4-pi normalized) to unnormalized spherical harmonic coefficients
 * \param normalizedCosineCoefficients Original, normalized, cosine coefficients
 * \param normalizedSineCoefficients Original, normalized, sine coefficients
 * \param unnormalizedCosineCoefficients Transformed, unnormalized, cosine coefficients (returned by reference)
 * \param unnormalizedSineCoefficients Transformed, unnormalized, sine coefficients (returned by reference)
 */
void convertGeodesyNormalizedToUnnormalizedCoefficients(
        const Eigen::MatrixXd& normalizedCosineCoefficients,
        const Eigen::MatrixXd& normalizedSineCoefficients,
        Eigen::MatrixXd& unnormalizedCosineCoefficients,
        Eigen::MatrixXd& unnormalizedSineCoefficients );

std::pair< Eigen::MatrixXd, Eigen::MatrixXd > convertGeodesyNormalizedToUnnormalizedCoefficients(
        const Eigen::MatrixXd& normalizedCosineCoefficients,
        const Eigen::MatrixXd& normalizedSineCoefficients );

//! Function to convert unnormalized to geodesy-normalized (4-pi normalized) spherical harmonic coefficients
/*!
 * Function to convert unnormalized to geodesy-normalized (4-pi normalized) spherical harmonic coefficients, by modifying input
 * coefficients
 * \param cosineCoefficients Original, unnormalized, cosine coefficients, to be transformed to normalized coefficients (returned
 * by reference)
 * \param sineCoefficients Original, unnormalized, cosine coefficients, to be transformed to normalized coefficients (returned
 * by reference)
 */
void geodesyNormalizeUnnormalizedCoefficients(
        Eigen::MatrixXd& cosineCoefficients,
        Eigen::MatrixXd& sineCoefficients );

} // namespace basic_mathematics
} // namespace tudat

#endif // TUDAT_LEGENDRE_POLYNOMIALS_H
