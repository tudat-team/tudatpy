/*    Copyright (c) 2010-2012, Delft University of Technology
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
#include <stdexcept>

#include <boost/exception/all.hpp>
#include <boost/math/constants/constants.hpp>

#include "Eigen/Core"

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravitationalAccelerationModel.h"
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
        const Eigen::Vector3d &position,
        const int highestDegree,
        const int highestOrder,
        const Eigen::MatrixXd& cosineHarmonicCoefficients,
        const Eigen::MatrixXd& sineHarmonicCoefficients,
        const double gravitationalParameter,
        const double planetaryRadius )
{
    // Declare spherical position vector.
    Eigen::Vector3d sphericalPosition;

    // Convert Cartesian coordinates to cylindrical.
    const Eigen::Vector3d cylindricalCoordinates = basic_mathematics::coordinate_conversions::
            convertCartesianToCylindrical( position );

    // Compute radius coordinate.
    sphericalPosition( 0 ) = std::sqrt( cylindricalCoordinates( 0 ) * cylindricalCoordinates( 0 )
                       + cylindricalCoordinates( 2 ) * cylindricalCoordinates( 2 ) );

    // If radius coordinate is smaller than planetary radius...
    if ( sphericalPosition( 0 ) < planetaryRadius )
    {
        // ...throw runtime error.
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Distance to origin is smaller than the size of the main body." ) ) );
    }

    // If radius coordinate is zero...
    if ( cylindricalCoordinates( 0 ) == 0 )
    {
        // ...set latitude coordinate to 90 degrees.
        sphericalPosition( 1 ) = mathematics::PI / 2.0;
    }

    // Else...
    else
    {
        // ...compute latitude coordinate.
        sphericalPosition( 1 ) = std::atan( cylindricalCoordinates( 2 )
                                            / cylindricalCoordinates( 0 ) );
    }

    // Compute longitude coordinate.
    sphericalPosition( 2 ) = cylindricalCoordinates( 1 );

    // Compute gradient premultiplier.
    const double preMultiplier = gravitationalParameter / planetaryRadius;

    // Initialize gradient vector.
    Eigen::Vector3d sphericalGradient = Eigen::VectorXd::Zero( 3 );

    // Loop through all degrees.
    for ( int degree = 0; degree <= highestDegree; degree++ )
    {
        // Loop through all orders.
        for ( int order = 0; order <= degree && order <= highestOrder; order++ )
        {
            // Compute geodesy-normalized Legendre polynomials.
            const double legendrePolynomial = basic_mathematics::computeGeodesyLegendrePolynomial(
                        degree, order, std::sin( sphericalPosition( 1 ) ) );
            const double incrementedLegendrePolynomial =
                    basic_mathematics::computeGeodesyLegendrePolynomial(
                        degree, order + 1, std::sin( sphericalPosition( 1 ) ) );

            // Compute geodesy-normalized Legendre polynomial derivative.
            const double legendrePolynomialDerivative =
                    basic_mathematics::computeGeodesyLegendrePolynomialDerivative(
                        degree, order, std::sin( sphericalPosition( 1 ) ), legendrePolynomial,
                        incrementedLegendrePolynomial );

            // Compute the potential gradient of a single spherical harmonic term.
            sphericalGradient += basic_mathematics::computePotentialGradient(
                        sphericalPosition,
                        planetaryRadius,
                        preMultiplier,
                        degree,
                        order,
                        cosineHarmonicCoefficients( degree, order),
                        sineHarmonicCoefficients( degree, order),
                        legendrePolynomial,
                        legendrePolynomialDerivative );
        }
    }

    // Convert from spherical gradient to Cartesian gradient (which equals acceleration vector) and
    // return the resulting acceleration vector.
    return basic_mathematics::coordinate_conversions::convertSphericalToCartesianGradient(
                sphericalGradient, position );
}

//! Compute gravitational acceleration due to single spherical harmonics term.
Eigen::Vector3d computeSingleGeodesyNormalizedGravitationalAcceleration(
        const Eigen::Vector3d &position,
        const int degree,
        const int order,
        const double cosineHarmonicCoefficient,
        const double sineHarmonicCoefficient,
        const double gravitationalParameter,
        const double planetaryRadius )
{
    // Declare spherical position vector.
    Eigen::Vector3d sphericalPosition;

    // Convert Cartesian coordinates to cylindrical.
    const Eigen::Vector3d cylindricalCoordinates = basic_mathematics::coordinate_conversions::
            convertCartesianToCylindrical( position );

    // Compute radius coordinate.
    sphericalPosition( 0 ) = std::sqrt(cylindricalCoordinates( 0 ) * cylindricalCoordinates( 0 )
                       + cylindricalCoordinates( 2 ) * cylindricalCoordinates( 2 ) );

    // If radius coordinate is smaller than planetary radius...
    if (sphericalPosition( 0 ) < planetaryRadius)
    {
        // ...trow runtime error.
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Distance to origin is smaller than the size of the main body." ) ) );
    }

    // If radius coordinate is zero...
    if ( cylindricalCoordinates( 0 ) == 0 )
    {
        // ...set latitude coordinate to 90 degrees.
        sphericalPosition( 1 ) = mathematics::PI / 2.0;
    }

    // Else...
    else
    {
        // ...compute latitude coordinate.
        sphericalPosition( 1 ) = std::atan( cylindricalCoordinates( 2 ) /
                                          cylindricalCoordinates( 0 ) );
    }

    // Compute longitude coordinate.
    sphericalPosition( 2 ) = cylindricalCoordinates( 1 );

    // Compute gradient premultiplier.
    const double preMultiplier = gravitationalParameter / planetaryRadius;

    // Compute geodesy-normalized Legendre polynomials.
    const double legendrePolynomial = basic_mathematics::computeGeodesyLegendrePolynomial(
                degree, order, std::sin( sphericalPosition( 1 ) ) );
    const double incrementedLegendrePolynomial =
            basic_mathematics::computeGeodesyLegendrePolynomial(
                degree, order + 1, std::sin( sphericalPosition( 1 ) ) );

    // Compute geodesy-normalized Legendre polynomial derivative.
    const double legendrePolynomialDerivative =
            basic_mathematics::computeGeodesyLegendrePolynomialDerivative(
                degree, order, std::sin( sphericalPosition( 1 ) ), legendrePolynomial,
                incrementedLegendrePolynomial );

    // Compute the potential gradient resulting from the spherical harmonic term.
    const Eigen::Vector3d sphericalGradient = basic_mathematics::computePotentialGradient(
                sphericalPosition,
                planetaryRadius,
                preMultiplier,
                degree,
                order,
                cosineHarmonicCoefficient,
                sineHarmonicCoefficient,
                legendrePolynomial,
                legendrePolynomialDerivative );

    // Convert from spherical gradient to Cartesian gradient (which equals acceleration vector),
    // and return resulting acceleration vector.
    return basic_mathematics::coordinate_conversions::convertSphericalToCartesianGradient(
                sphericalGradient, position );
}

} // namespace gravitation
} // namespace tudat
