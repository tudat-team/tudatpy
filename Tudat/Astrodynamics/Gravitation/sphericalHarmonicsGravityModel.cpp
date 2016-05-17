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
 *      121017    E. Dekens         Code created.
 *
 *    References
 *
 *    Notes
 *
 */

#include <cmath>
#include <limits>
#include <stdexcept>

#include <boost/exception/all.hpp>
#include <boost/math/constants/constants.hpp>

#include "Eigen/Core"

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2GravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2J3GravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModel.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/legendrePolynomials.h"
#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonics.h"

namespace tudat
{
namespace gravitation
{

//! Compute gravitational acceleration due to multiple spherical harmonics terms, defined using
//! geodesy-normalization.
Eigen::Vector3d computeGeodesyNormalizedGravitationalAccelerationSum(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameter,
        const double equatorialRadius,
        const Eigen::MatrixXd& cosineHarmonicCoefficients,
        const Eigen::MatrixXd& sineHarmonicCoefficients,
        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache )
{
    // Set highest degree and order.
    const int highestDegree = cosineHarmonicCoefficients.rows( );
    const int highestOrder = cosineHarmonicCoefficients.cols( );

    // Declare spherical position vector.
    Eigen::Vector3d sphericalpositionOfBodySubjectToAcceleration = Eigen::Vector3d::Zero( );

    // Convert Cartesian coordinates to cylindrical.
    const Eigen::Vector3d cylindricalCoordinates = coordinate_conversions::
            convertCartesianToCylindrical( positionOfBodySubjectToAcceleration );

    // Compute radius coordinate.
    sphericalpositionOfBodySubjectToAcceleration( 0 )
            = std::sqrt( cylindricalCoordinates( 0 ) * cylindricalCoordinates( 0 )
                       + cylindricalCoordinates( 2 ) * cylindricalCoordinates( 2 ) );

    // If radius coordinate is smaller than planetary radius...
    if ( sphericalpositionOfBodySubjectToAcceleration( 0 ) < equatorialRadius )
    {
        // ...throw runtime error.
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Distance to origin is smaller than the size of the main body." ) ) );
    }

    // If radius coordinate is zero...
    if ( std::fabs( cylindricalCoordinates( 0 ) ) < std::numeric_limits< double >::epsilon( ) )
    {
        // ...set latitude coordinate to 90 degrees.
        sphericalpositionOfBodySubjectToAcceleration( 1 )
                = mathematical_constants::PI / 2.0;
    }

    // Else...
    else
    {
        // ...compute latitude coordinate.
        sphericalpositionOfBodySubjectToAcceleration( 1 )
                = std::atan( cylindricalCoordinates( 2 ) / cylindricalCoordinates( 0 ) );
    }

    // Compute longitude coordinate.
    sphericalpositionOfBodySubjectToAcceleration( 2 ) = cylindricalCoordinates( 1 );

    double sineOfAngle = std::sin( sphericalpositionOfBodySubjectToAcceleration( 1 ) );
    sphericalHarmonicsCache->update( sphericalpositionOfBodySubjectToAcceleration( 0 ),
                                     sineOfAngle,
                                     sphericalpositionOfBodySubjectToAcceleration( 2 ),
                                     equatorialRadius );
    basic_mathematics::LegendreCache& legendreCacheReference = *( sphericalHarmonicsCache->getLegendreCache( ) );

    // Compute gradient premultiplier.
    const double preMultiplier = gravitationalParameter / equatorialRadius;

    // Initialize gradient vector.
    Eigen::Vector3d sphericalGradient = Eigen::Vector3d::Zero( );

    // Loop through all degrees.
    for ( int degree = 0; degree < highestDegree; degree++ )
    {
        // Loop through all orders.
        for ( int order = 0; ( order <= degree ) && ( order < highestOrder ); order++ )
        {
            // Compute geodesy-normalized Legendre polynomials.
            const double legendrePolynomial = basic_mathematics::computeGeodesyLegendrePolynomial(
                        degree, order, legendreCacheReference );
            const double incrementedLegendrePolynomial =
                    basic_mathematics::computeGeodesyLegendrePolynomial(
                        degree, order + 1, legendreCacheReference );

            // Compute geodesy-normalized Legendre polynomial derivative.
            const double legendrePolynomialDerivative =
                    basic_mathematics::computeGeodesyLegendrePolynomialDerivative(
                        degree, order,
                        sineOfAngle,
                        legendrePolynomial, incrementedLegendrePolynomial );

            // Compute the potential gradient of a single spherical harmonic term.
            sphericalGradient += basic_mathematics::computePotentialGradient(
                        sphericalpositionOfBodySubjectToAcceleration,
                        preMultiplier,
                        degree,
                        order,
                        cosineHarmonicCoefficients( degree, order ),
                        sineHarmonicCoefficients( degree, order ),
                        legendrePolynomial,
                        legendrePolynomialDerivative, sphericalHarmonicsCache );
        }
    }

    // Convert from spherical gradient to Cartesian gradient (which equals acceleration vector) and
    // return the resulting acceleration vector.
    return coordinate_conversions::convertSphericalToCartesianGradient(
                sphericalGradient, positionOfBodySubjectToAcceleration );
}

//! Compute gravitational acceleration due to single spherical harmonics term.
Eigen::Vector3d computeSingleGeodesyNormalizedGravitationalAcceleration(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameter,
        const double equatorialRadius,
        const int degree,
        const int order,
        const double cosineHarmonicCoefficient,
        const double sineHarmonicCoefficient,
        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache )
{
    // Declare spherical position vector.
    Eigen::Vector3d sphericalpositionOfBodySubjectToAcceleration;

    // Convert Cartesian coordinates to cylindrical.
    const Eigen::Vector3d cylindricalCoordinates = coordinate_conversions::
            convertCartesianToCylindrical( positionOfBodySubjectToAcceleration );

    // Compute radius coordinate.
    sphericalpositionOfBodySubjectToAcceleration( 0 )
            = std::sqrt(cylindricalCoordinates( 0 ) * cylindricalCoordinates( 0 )
                       + cylindricalCoordinates( 2 ) * cylindricalCoordinates( 2 ) );

    // If radius coordinate is smaller than planetary radius...
    if (sphericalpositionOfBodySubjectToAcceleration( 0 ) < equatorialRadius)
    {
        // ...trow runtime error.
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Distance to origin is smaller than the size of the main body." ) ) );
    }

    // If radius coordinate is zero...
    if ( std::fabs( cylindricalCoordinates( 0 ) ) < std::numeric_limits< double >::epsilon( ) )
    {
        // ...set latitude coordinate to 90 degrees.
        sphericalpositionOfBodySubjectToAcceleration( 1 )
                = mathematical_constants::PI / 2.0;
    }

    // Else...
    else
    {
        // ...compute latitude coordinate.
        sphericalpositionOfBodySubjectToAcceleration( 1 )
                = std::atan( cylindricalCoordinates( 2 ) / cylindricalCoordinates( 0 ) );
    }

    // Compute longitude coordinate.
    sphericalpositionOfBodySubjectToAcceleration( 2 ) = cylindricalCoordinates( 1 );


    double sineOfAngle = std::sin( sphericalpositionOfBodySubjectToAcceleration( 1 ) );
    sphericalHarmonicsCache->update( sphericalpositionOfBodySubjectToAcceleration( 0 ),
                                     sineOfAngle,
                                     sphericalpositionOfBodySubjectToAcceleration( 2 ),
                                     equatorialRadius );
    basic_mathematics::LegendreCache& legendreCacheReference = *( sphericalHarmonicsCache->getLegendreCache( ) );

    // Compute gradient premultiplier.
    const double preMultiplier = gravitationalParameter / equatorialRadius;

    // Compute geodesy-normalized Legendre polynomials.
    const double legendrePolynomial = basic_mathematics::computeGeodesyLegendrePolynomial(
                degree, order, legendreCacheReference );
    const double incrementedLegendrePolynomial =
            basic_mathematics::computeGeodesyLegendrePolynomial(
                degree, order + 1,  legendreCacheReference );

    // Compute geodesy-normalized Legendre polynomial derivative.x
    const double legendrePolynomialDerivative =
            basic_mathematics::computeGeodesyLegendrePolynomialDerivative(
                degree, order,
                sineOfAngle, legendrePolynomial,
                incrementedLegendrePolynomial );

    // Compute the potential gradient resulting from the spherical harmonic term.
    const Eigen::Vector3d sphericalGradient = basic_mathematics::computePotentialGradient(
                sphericalpositionOfBodySubjectToAcceleration,
                preMultiplier,
                degree,
                order,
                cosineHarmonicCoefficient,
                sineHarmonicCoefficient,
                legendrePolynomial,
                legendrePolynomialDerivative, sphericalHarmonicsCache );

    // Convert from spherical gradient to Cartesian gradient (which equals acceleration vector),
    // and return resulting acceleration vector.
    return coordinate_conversions::convertSphericalToCartesianGradient(
                sphericalGradient, positionOfBodySubjectToAcceleration );
}

} // namespace gravitation
} // namespace tudat
