/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_SPHERICAL_HARMONICS_H
#define TUDAT_SPHERICAL_HARMONICS_H

#include <Eigen/Core>

#include <boost/make_shared.hpp>

#include "Tudat/Mathematics/BasicMathematics/legendrePolynomials.h"

namespace tudat
{
namespace basic_mathematics
{

//! Cache object in which variables that are required for the computation of spherical harmonic potential are stored.
/*!
 *  Cache object in which variables that are required for the computation of spherical harmonic potential are stored.
 *  The variables are the Legendre polynomials at the required degree and order, the cosine of the latitude, the
 *  sine and cosine of the order times the longitude, and the ratio of the distance and the reference radius to the
 *  power degree + 1.
 */
class SphericalHarmonicsCache
{
public:

    //! Default constructor, initializes cache object with 0 maximum degree and order.
    /*!
     * Default constructor, initializes cache object with 0 mazimum degree and order.
     * \param useGeodesyNormalization Parameter defining whether the cache is used for a normalized or unnormalized
     * gravity field.
     */
    SphericalHarmonicsCache( const bool useGeodesyNormalization = 1 )
    {
        legendreCache_ = std::make_shared< LegendreCache >( useGeodesyNormalization );
        currentLongitude_ = TUDAT_NAN;
        referenceRadiusRatio_ = TUDAT_NAN;

        resetMaximumDegreeAndOrder( 0, 0 );
    }

    //! Constructor
    /*!
     * Constructor
     * \param maximumDegree Maximum degree to which to update cache
     * \param maximumOrder Maximum order to which to update cache
     * \param useGeodesyNormalization Parameter defining whether the cache is used for a normalized or unnormalized
     * gravity field.
     */
    SphericalHarmonicsCache( const int maximumDegree, const int maximumOrder, const bool useGeodesyNormalization = 1 )
    {
        legendreCache_ = std::make_shared< LegendreCache >( maximumDegree, maximumOrder, useGeodesyNormalization );

        currentLongitude_ = TUDAT_NAN;
        referenceRadiusRatio_ = TUDAT_NAN;

        resetMaximumDegreeAndOrder( maximumDegree, maximumOrder );
    }

    //! Update maximum degree and order of cache
    /*!
     * Update maximum degree and order of cache
     * \param maximumDegree Maximum degree to which to update cache
     * \param maximumOrder Maximum order to which to update cache
     */
    void resetMaximumDegreeAndOrder( const int maximumDegree, const int maximumOrder );


    //! Update cached variables to current state.
    /*!
     * Update cached variables to current state.
     * \param radius Distance from origin
     * \param polynomialParameter Input parameter to Legendre polynomials (sine of latitude)
     * \param longitude Current latitude
     * \param referenceRadius Reference (typically equatorial) radius of gravity field.
     */
    void update( const double radius, const double polynomialParameter,
                 const double longitude, const double referenceRadius )
    {
        legendreCache_->update( polynomialParameter );
        updateSines( longitude );
        updateRadiusPowers( referenceRadius / radius );
    }

    //! Function to retrieve the current sine of m times the longitude.
    /*!
     * Function to retrieve the current sine of order times the longitude.
     * \param order Order as input to sine( order * longitude )
     * \return Sine( order * longitude )
     */
    double getSineOfMultipleLongitude( const int order )
    {
        return sinesOfLongitude_[ order ];
    }

    //! Function to retrieve the current cosine of m times the longitude.
    /*!
     * Function to retrieve the current cosine of order times the longitude.
     * \param order Order as input to cosine( order * longitude )
     * \return Cosine( order * longitude )
     */
    double getCosineOfMultipleLongitude( const int order )
    {
        return cosinesOfLongitude_[ order ];
    }

    //! Function to get an integer power of the distance divided by the reference radius.
    /*!
     * Function to get an integer power of the distance divided by the reference radius.
     * \param degreePlusOne Power to which the ratio of distance and reference radius is to be computed (typically
     * degree + 1).
     * \return Ratio of distance and reference radius to power of input argument.
     */
    double getReferenceRadiusRatioPowers( const int degreePlusOne )
    {
        return referenceRadiusRatioPowers_[ degreePlusOne ];
    }

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

    //! Function to get current longitude
    /*!
     * Function to get current longitude
     * \return Current longitude
     */
    double getCurrentLongitude( )
    {
        return currentLongitude_;
    }

    //! Function to get object for caching and computing Legendre polynomials.
    /*!
     * Function to get object for caching and computing Legendre polynomials.
     * \return Object for caching and computing Legendre polynomials.
     */
    std::shared_ptr< LegendreCache > getLegendreCache( )
    {
        return legendreCache_;
    }

private:

    //! Update cached values of sines and cosines of longitude/
    /*!
     * Update cached values of sines and cosines of longitude/
     * \param longitude Current longitude.
     */
    void updateSines( const double longitude )
    {
        //! Check if update is needed.
        if( !( currentLongitude_ == longitude ) )
        {
            currentLongitude_ = longitude;
            for( unsigned int i = 0; i < sinesOfLongitude_.size( ); i++ )
            {
                sinesOfLongitude_[ i ] = std::sin( static_cast< double >( i ) * longitude );
                cosinesOfLongitude_[ i ] = std::cos( static_cast< double >( i ) * longitude );
            }
        }
    }

    //! Update cached values of powers of distance over reference radius.
    /*!
     * Update cached values of powers of distance over reference radius.
     * \param referenceRadiusRatio Distance divided by reference radius.
     */
    void updateRadiusPowers( const double referenceRadiusRatio )
    {
        //! Check if update is needed.
        if( !( referenceRadiusRatio_ == referenceRadiusRatio ) )
        {
            referenceRadiusRatio_ = referenceRadiusRatio;
            double currentRatioPower = 1.0;
            for( int i = 0; i <= maximumDegree_ + 1; i++ )
            {
                referenceRadiusRatioPowers_[ i ] = currentRatioPower;
                currentRatioPower *= referenceRadiusRatio_;
            }
        }
    }

    //! Maximum degree of cache.
    int maximumDegree_;

    //! Maximum order of cache.
    int maximumOrder_;

    //! Current longitude.
    double currentLongitude_;

    //! Current ratio of distance to reference radius
    double referenceRadiusRatio_;

    //! List of sines of order times longitude.
    /*!
     *  List of sines of order times longitude. Entry i denotes sin(i times longitude).
     */
    std::vector< double > sinesOfLongitude_;

    //! List of cosines of order times longitude.
    /*!
     *  List of cosines of order times longitude. Entry i denotes cos(i times longitude).
     */
    std::vector< double > cosinesOfLongitude_;

    //! List of powers of distance divided by reference radius.
    /*!
     * List of powers of distance divided by reference radius. Entry i denoted (distance/reference radius) to the power i.
     */
    std::vector< double > referenceRadiusRatioPowers_;

    //! Object for caching and computing Legendre polynomials.
    std::shared_ptr< LegendreCache > legendreCache_;

};

//! Spherical coordinate indices.
enum SphericalCoordinatesIndices{ radiusIndex, latitudeIndex, longitudeIndex };

//! Compute the gradient of a single term of a spherical harmonics potential field.
/*!
 * This function returns a vector with the derivatives of a generic potential field (defined by
 * spherical harmonics) from pre-computed quatities.
 * \param distance Distance to center of body with gravity field at which the potential gradient is to be calculated
 * \param radiusPowerTerm Distance divided by the reference radius of the gravity field, to the power (degree + 1)
 * \param cosineOfOrderLongitude Cosine of order times the longitude at which the potential is to be calculated
 * \param sineOfOrderLongitude Sine of order times the longitude at which the potential is to be calculated
 * \param cosineOfLatitude Cosine of the latitude at which the potential is to be calculated
 * \param preMultiplier Generic multiplication factor.
 * \param degree Degree of the harmonic for which the gradient is to be computed.
 * \param order Order of the harmonic for which the gradient is to be computed.
 * \param cosineHarmonicCoefficient Coefficient which characterizes relative strengh of a harmonic
 *          term.
 * \param sineHarmonicCoefficient Coefficient which characterizes relative strengh of a harmonic
 *          term.
 * \param legendrePolynomial Value of associated Legendre polynomial with the same degree and order
 *          as the to be computed harmonic, and with the sine of the latitude coordinate as
 *          polynomial parameter. Make sure that the Legendre polynomial has the same
 *          normalization as the harmonic coefficients.
 * \param legendrePolynomialDerivative Value of the derivative of parameter 'legendrePolynomial'
 *          with respect to the sine of the latitude angle.
 * \return Vector with derivatives of potential field.
 *          The order is important!
 *          gradient( 0 ) = derivative with respect to radial distance,
 *          gradient( 1 ) = derivative with respect to latitude angle,
 *          gradient( 2 ) = derivative with respect to longitude angle.
 */
Eigen::Vector3d computePotentialGradient(
        const double distance,
        const double radiusPowerTerm,
        const double cosineOfOrderLongitude,
        const double sineOfOrderLongitude,
        const double cosineOfLatitude,
        const double preMultiplier,
        const int degree,
        const int order,
        const double cosineHarmonicCoefficient,
        const double sineHarmonicCoefficient,
        const double legendrePolynomial,
        const double legendrePolynomialDerivative );



//! Compute the gradient of a single term of a spherical harmonics potential field.
/*!
 * This function returns a vector with the derivatives of a generic potential field (defined by
 * spherical harmonics). *
 * \param sphericalPosition Vector with spherical coordinates.
 *          The order is important!
 *          sphericalPosition( 0 ) = radial coordinate,
 *          sphericalPosition( 1 ) = latitude coordinate,
 *          sphericalPosition( 2 ) = longitude coordinate.
 * \param referenceRadius Radius of harmonics reference sphere.
 * \param preMultiplier Generic multiplication factor.
 * \param degree Degree of the harmonic for which the gradient is to be computed.
 * \param order Order of the harmonic for which the gradient is to be computed.
 * \param cosineHarmonicCoefficient Coefficient which characterizes relative strengh of a harmonic
 *          term.
 * \param sineHarmonicCoefficient Coefficient which characterizes relative strengh of a harmonic
 *          term.
 * \param legendrePolynomial Value of associated Legendre polynomial with the same degree and order
 *          as the to be computed harmonic, and with the sine of the latitude coordinate as
 *          polynomial parameter. Make sure that the Legendre polynomial has the same
 *          normalization as the harmonic coefficients.
 * \param legendrePolynomialDerivative Value of the derivative of parameter 'legendrePolynomial'
 *          with respect to the sine of the latitude angle.
 * \return Vector with derivatives of potential field.
 *          The order is important!
 *          gradient( 0 ) = derivative with respect to radial distance,
 *          gradient( 1 ) = derivative with respect to latitude angle,
 *          gradient( 2 ) = derivative with respect to longitude angle.
 */
Eigen::Vector3d computePotentialGradient( const Eigen::Vector3d& sphericalPosition,
                                          const double referenceRadius,
                                          const double preMultiplier,
                                          const int degree,
                                          const int order,
                                          const double cosineHarmonicCoefficient,
                                          const double sineHarmonicCoefficient,
                                          const double legendrePolynomial,
                                          const double legendrePolynomialDerivative );

//! Compute the gradient of a single term of a spherical harmonics potential field.
/*!
 * This function returns a vector with the derivatives of a generic potential field (defined by
 * spherical harmonics). It is assumed that the potential field of a single harmonic is
 * characterized by:
 * \f[
 *     U = A \left( \frac{ R }{ r} \right) ^{ n + 1}
 *     P _{ n, m } ( \sin \phi ) \left[ C _{ n, m } \cos( m \lambda ) + S _{ n, m } \sin( \lambda )
 *     \right]
 * \f]
 * in which \f$ A \f$ is the generic multiplication factor, \f$ R \f$ is the radius of the harmonics
 * reference sphere, \f$ P _{ n, m }( \sin \phi ) \f$ is the associated Legendre polynomial with
 * \f$ \sin \phi \f$  as polynomial parameter, \f$ r \f$ is the radial coordinate, \f$ \phi \f$ is
 * the latitude coordinate, \f$ \lambda \f$ is the longitude coordinate, \f$ n \f$ is the harmonics
 * degree, \f$ m \f$ is the harmonics order, \f$ C _{ n, m } \f$ is the cosine harmonics
 * coefficient, and \f$ S _{ n, m } \f$ is the sine harmonics coefficient.
 *
 * The potential derivatives are calculated through:
 * \f{eqnarray*}{
 *     \frac{ \mathrm{ d } U }{ \mathrm{ d } r } = -\frac{ A }{ r } \left( \frac{ R }{ r } \right)
 *     ^{ n + 1 } ( n + 1 ) P_{ n, m }( \sin \phi )[ C_{ n, m } \cos( m \lambda ) + S_{ n,m }
 *     \sin( m \lambda ) ] \\
 *     \frac{ \mathrm{ d } U }{ \mathrm{ d } \phi } = A \left( \frac{ R }{ r } \right)^{ n + 1 }
 *     \frac{ \mathrm{ d } [ P( \sin \phi ) ] }{ \mathrm{ d } [ \sin \phi ] } \cos \phi
 *     [ C_{ n, m } \cos( m \lambda ) + S_{ n, m } \sin( m \lambda ) ] \\
 *     \frac{ \mathrm{ d } U }{ \mathrm{ d } \lambda } = A \left( \frac{ R }{ r } \right)^{ n + 1 }
 *     m P_{ n, m }( \sin \phi ) [ S_{ n, m } \cos( m \lambda ) - C_{ n, m }
 *     \sin( m \lambda ) ]
 * \f}
 *
 * \param sphericalPosition Vector with spherical coordinates.
 *          The order is important!
 *          sphericalPosition( 0 ) = radial coordinate,
 *          sphericalPosition( 1 ) = latitude coordinate,
 *          sphericalPosition( 2 ) = longitude coordinate.
 * \param preMultiplier Generic multiplication factor.
 * \param degree Degree of the harmonic for which the gradient is to be computed.
 * \param order Order of the harmonic for which the gradient is to be computed.
 * \param cosineHarmonicCoefficient Coefficient which characterizes relative strengh of a harmonic
 *          term.
 * \param sineHarmonicCoefficient Coefficient which characterizes relative strengh of a harmonic
 *          term.
 * \param legendrePolynomial Value of associated Legendre polynomial with the same degree and order
 *          as the to be computed harmonic, and with the sine of the latitude coordinate as
 *          polynomial parameter. Make sure that the Legendre polynomial has the same
 *          normalization as the harmonic coefficients.
 * \param legendrePolynomialDerivative Value of the derivative of parameter 'legendrePolynomial'
 *          with respect to the sine of the latitude angle.
 * \param sphericalHarmonicsCache Cache object containing current values of trigonometric funtions of latitude anf longitude,
 *          as well as legendre polynomials at current state.
 * \return Vector with derivatives of potential field.
 *          The order is important!
 *          gradient( 0 ) = derivative with respect to radial distance,
 *          gradient( 1 ) = derivative with respect to latitude angle,
 *          gradient( 2 ) = derivative with respect to longitude angle.
 */
Eigen::Vector3d computePotentialGradient( const Eigen::Vector3d& sphericalPosition,
                                          const double preMultiplier,
                                          const int degree,
                                          const int order,
                                          const double cosineHarmonicCoefficient,
                                          const double sineHarmonicCoefficient,
                                          const double legendrePolynomial,
                                          const double legendrePolynomialDerivative,
                                          const std::shared_ptr< SphericalHarmonicsCache > sphericalHarmonicsCache );

} // namespace basic_mathematics
} // namespace tudat

#endif // TUDAT_SPHERICAL_HARMONICS_H
