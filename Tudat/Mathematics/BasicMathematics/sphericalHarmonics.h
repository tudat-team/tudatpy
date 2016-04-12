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
 *
 *    References
 *
 *    Notes
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

class SphericalHarmonicsCache
{
public:

    SphericalHarmonicsCache( const bool useGeodesyNormalization = 1 )
    {
        legendreCache_ = boost::make_shared< LegendreCache >( useGeodesyNormalization );

        currentLongitude_ = TUDAT_NAN;
        referenceRadiusRatio_ = TUDAT_NAN;

        resetMaximumDegreeAndOrder( 0, 0 );
    }

    SphericalHarmonicsCache( const int maximumDegree, const int maximumOrder, const bool useGeodesyNormalization = 1 )
    {
        legendreCache_ = boost::make_shared< LegendreCache >( maximumDegree, maximumOrder, useGeodesyNormalization );

        currentLongitude_ = TUDAT_NAN;
        referenceRadiusRatio_ = TUDAT_NAN;

        resetMaximumDegreeAndOrder( maximumDegree, maximumOrder );
    }

    void resetMaximumDegreeAndOrder( const int degree, const int order );


    double getSineOfMultipleLongitude( const int i )
    {
        return sinesOfLongitude_[ i ];
    }

    double getCosineOfMultipleLongitude( const int i )
    {
        return cosinesOfLongitude_[ i ];
    }

    double getReferenceRadiusRatioPowers( const int i )
    {
        return referenceRadiusRatioPowers_[ i ];
    }

    int getMaximumDegree( ){ return maximumDegree_; }

    int getMaximumOrder( ){ return maximumOrder_; }

    double getCurrentLongitude( ){ return currentLongitude_; }

    boost::shared_ptr< LegendreCache > getLegendreCache( )
    {
        return legendreCache_;
    }

    void update( const double radius, const double polynomialParameter, const double longitude, const double referenceRadius )
    {
        legendreCache_->update( polynomialParameter );
        updateSines( longitude );
        updateRadiusPowers( referenceRadius / radius );
    }

private:

    void updateSines( const double longitude )
    {
        if( !(currentLongitude_ == longitude ) )
        {
            currentLongitude_ = longitude;
            for( unsigned int i = 0; i < sinesOfLongitude_.size( ); i++ )
            {
                sinesOfLongitude_[ i ] = std::sin( static_cast< double >( i ) * longitude );
                cosinesOfLongitude_[ i ] = std::cos( static_cast< double >( i ) * longitude );
            }
        }
    }
    void updateRadiusPowers( const double referenceRadiusRatio )
    {
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

    int maximumDegree_;

    int maximumOrder_;

    double currentLongitude_;

    double referenceRadiusRatio_;

    std::vector< double > sinesOfLongitude_;

    std::vector< double > cosinesOfLongitude_;

    std::vector< double > referenceRadiusRatioPowers_;

    double returnValue_ ;

    boost::shared_ptr< LegendreCache > legendreCache_;



};

//! Spherical coordinate indices.
enum SphericalCoordinatesIndices{ radiusIndex, latitudeIndex, longitudeIndex };

Eigen::Vector3d computePotentialGradient(
        const double distance,
        const double radiusPowerTerm,
        const double cosineOfOrderLongitude, const double sineOfOrderLongitude,
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

Eigen::Vector3d computePotentialGradient( const Eigen::Vector3d& sphericalPosition,
                                          const double preMultiplier,
                                          const int degree,
                                          const int order,
                                          const double cosineHarmonicCoefficient,
                                          const double sineHarmonicCoefficient,
                                          const double legendrePolynomial,
                                          const double legendrePolynomialDerivative,
                                          boost::shared_ptr< SphericalHarmonicsCache > sphericalHarmonicsCache );

} // namespace basic_mathematics
} // namespace tudat

#endif // TUDAT_SPHERICAL_HARMONICS_H
