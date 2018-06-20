/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPHERICALHARMONICPARTIALFUNCTIONS_H
#define TUDAT_SPHERICALHARMONICPARTIALFUNCTIONS_H

#include <memory>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonics.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function to compute the spherical Hessian of a single term of a spherical harmonic potential
/*!
 *  Function to compute the spherical Hessian (i.e. matrix of second derivatives w.r.t. spherical components radius, latitude
 *  and longitude) of a single term of a spherical harmonic potential.
 *  \param distance Distance to center of body with gravity field at which the partials are to be calculated
 *  \param radiusPowerTerm Distance divided by the reference radius of the gravity field, to the power (degree + 1)
 *  \param cosineOfOrderLongitude Cosine of order times the longitude at which the potential is to be calculated
 *  \param sineOfOrderLongitude Sine of order times the longitude at which the potential is to be calculated
 *  \param cosineOfLatitude Cosine of the latitude at which the potential is to be calculated
 *  \param sineOfLatitude Sine of the latitude at which the potential is to be calculated
 *  \param preMultiplier Pre-multiplier of potential (gravitational parametere divided by reference radius in normal
 *  representation)
 *  \param degree Degree of the harmonic for which the gradient is to be computed.
 *  \param order Order of the harmonic for which the gradient is to be computed.
 *  \param cosineHarmonicCoefficient Coefficient which characterizes relative strengh of a harmonic term.
 *  \param sineHarmonicCoefficient Coefficient which characterizes relative strengh of a harmonic term.
 *  \param legendrePolynomial Value of associated Legendre polynomial with the same degree and order as the to be computed
 *  harmonic, and with the sine of the latitude coordinate as polynomial parameter.
 *  \param legendrePolynomialDerivative Value of the derivative of parameter 'legendrePolynomial' with respect to the sine
 *  of the latitude angle.
 *  \param legendrePolynomialSecondDerivative Value of second derivative of parameter 'legendrePolynomial' with respect to
 *  the sine of the latitude angle.
 *  \param sphericalHessian Hessian of potential term in spherical coordinates (returned by reference).
 */
void computePotentialSphericalHessian(
        const double distance,
        const double radiusPowerTerm,
        const double cosineOfOrderLongitude,
        const double sineOfOrderLongitude,
        const double cosineOfLatitude,
        const double sineOfLatitude,
        const double preMultiplier,
        const int degree,
        const int order,
        const double cosineHarmonicCoefficient,
        const double sineHarmonicCoefficient,
        const double legendrePolynomial,
        const double legendrePolynomialDerivative,
        const double legendrePolynomialSecondDerivative,
        Eigen::Matrix3d& sphericalHessian );

//! Function to compute the spherical Hessian of a single term of a spherical harmonic potential
/*!
 *  Function to compute the spherical Hessian (i.e. matrix of second derivatives w.r.t. spherical components radius, latitude
 *  and longitude) of a single term of a spherical harmonic potential.
 *  \param sphericalPosition Spherical position (radius, ,latitude, longitude) at which potential partials are to be
 *  evaluated
 *  \param referenceRadius Reference radius of spherical harmonic potential.
 *  \param preMultiplier Pre-multiplier of potential (gravitational parametere divided by reference radius in normal
 *  representation)
 *  \param degree Degree of the harmonic for which the gradient is to be computed.
 *  \param order Order of the harmonic for which the gradient is to be computed.
 *  \param cosineHarmonicCoefficient Coefficient which characterizes relative strengh of a harmonic term.
 *  \param sineHarmonicCoefficient Coefficient which characterizes relative strengh of a harmonic term.
 *  \param legendrePolynomial Value of associated Legendre polynomial with the same degree and order as the to be computed
 *  harmonic, and with the sine of the latitude coordinate as polynomial parameter.
 *  \param legendrePolynomialDerivative Value of the derivative of parameter 'legendrePolynomial' with respect to the sine
 *  of the latitude angle.
 *  \param legendrePolynomialSecondDerivative Value of second derivative of parameter 'legendrePolynomial' with respect to
 *  the sine of the latitude angle.
 *  \param sphericalHessian Hessian of potential term in spherical coordinates (returned by reference).
 */
void computePotentialSphericalHessian(
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double preMultiplier,
        const int degree,
        const int order,
        const double cosineHarmonicCoefficient,
        const double sineHarmonicCoefficient,
        const double legendrePolynomial,
        const double legendrePolynomialDerivative,
        const double legendrePolynomialSecondDerivative,
        Eigen::Matrix3d& sphericalHessian );

//! Function to compute the spherical Hessian of a single term of a spherical harmonic potential
/*!
 *  Function to compute the spherical Hessian (i.e. matrix of second derivatives w.r.t. spherical components radius, latitude
 *  and longitude) of a single term of a spherical harmonic potential.
 *  \param sphericalPosition Spherical position (radius, ,latitude, longitude) at which potential partials are to be
 *  evaluated
 *  \param preMultiplier Pre-multiplier of potential (gravitational parametere divided by reference radius in normal
 *  representation)
 *  \param degree Degree of the harmonic for which the gradient is to be computed.
 *  \param order Order of the harmonic for which the gradient is to be computed.
 *  \param cosineHarmonicCoefficient Coefficient which characterizes relative strengh of a harmonic term.
 *  \param sineHarmonicCoefficient Coefficient which characterizes relative strengh of a harmonic term.
 *  \param sphericalHarmonicsCache Cache object containing precomputed spherical harmonics terms.
 *  \param sphericalHessian Hessian of potential term in spherical coordinates (returned by reference).
 */
void computePotentialSphericalHessian(
        const Eigen::Vector3d& sphericalPosition,
        const double preMultiplier,
        const int degree,
        const int order,
        const double cosineHarmonicCoefficient,
        const double sineHarmonicCoefficient,
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        Eigen::Matrix3d& sphericalHessian );

//! Function to compute the spherical Hessian of a full spherical harmonic potential
/*!
 *  Function to compute the spherical Hessian (i.e. matrix of second derivatives w.r.t. spherical components radius, latitude
 *  and longitude) of a full spherical harmonic potential.
 *  \param sphericalPosition Spherical position (radius, ,latitude, longitude) at which potential partials are to be
 *  evaluated
 *  \param referenceRadius Reference radius of spherical harmonic potential.
 *  \param gravitionalParameter Gravitational parameter used for spherical harmonic expansion
 *  \param cosineHarmonicCoefficients Matrix of coefficient which characterize the relative strengh of cosine harmonic
 *  terms.
 *  \param sineHarmonicCoefficients Matrix of coefficient which characterize the relative strengh of sine harmonic terms.
 *  \param sphericalHarmonicsCache Cache object containing precomputed spherical harmonics terms.
 *  \return Hessian of potential in spherical coordinates (returned by reference).
 */
Eigen::Matrix3d computeCumulativeSphericalHessian(
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const Eigen::MatrixXd cosineHarmonicCoefficients,
        const Eigen::MatrixXd sineHarmonicCoefficients,
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache );

//! Calculate partial of spherical harmonic acceleration w.r.t. position of body undergoing acceleration
//! (in the body-fixed frame)
/*!
 *  Calculate partial of spherical harmonic acceleration w.r.t. position of body undergoing acceleration
 * (in the body-fixed frame)
 * \param cartesianPosition Cartesian position  at which potential partials are to be evaluated
 * \param sphericalPosition Spherical position (radius, ,latitude, longitude) at which potential partials are to be
 * evaluated
 * \param referenceRadius Reference radius of spherical harmonic potential.
 * \param gravitionalParameter Gravitational parameter used for spherical harmonic expansion
 * \param cosineHarmonicCoefficients Cosine spherical harmonic coefficients.
 * \param sineHarmonicCoefficients Sine spherical harmonic coefficients
 * \param sphericalHarmonicsCache Cache object containing precomputed spherical harmonics terms.
 * \param sphericalPotentialGradient Potential gradient in spherical coordinates
 * \param sphericalToCartesianGradientMatrix Matrix to convert (by premultiplication) a spherical gradient to a Cartesian
 * gradient
 * \return Partial of spherical harmonic acceleration w.r.t. position of body undergoinng acceleration (equals minus
 * partial of spherical harmonic acceleration w.r.t. position of body exerting acceleration) with both acceleration and
 * position in body-fixed frame.
 */
Eigen::Matrix3d computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
        const Eigen::Vector3d& cartesianPosition,
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const Eigen::MatrixXd cosineHarmonicCoefficients,
        const Eigen::MatrixXd sineHarmonicCoefficients,
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const Eigen::Vector3d& sphericalPotentialGradient,
        const Eigen::Matrix3d& sphericalToCartesianGradientMatrix );

//! Calculate partial of spherical harmonic acceleration w.r.t. position of body undergoing acceleration
//! (in the body-fixed frame)
/*!
 *  Calculate partial of spherical harmonic acceleration w.r.t. position of body undergoing acceleration
 * (in the body-fixed frame)
 * \param cartesianPosition Cartesian position  at which potential partials are to be evaluated
 * \param referenceRadius Reference radius of spherical harmonic potential.
 * \param gravitionalParameter Gravitational parameter used for spherical harmonic expansion
 * \param cosineHarmonicCoefficients Cosine spherical harmonic coefficients.
 * \param sineHarmonicCoefficients Sine spherical harmonic coefficients
 * \param sphericalHarmonicsCache Cache object containing precomputed spherical harmonics terms.
 * \return Partial of spherical harmonic acceleration w.r.t. position of body undergoinng acceleration (equals minus
 * partial of spherical harmonic acceleration w.r.t. position of body exerting acceleration) with both acceleration and
 * position in body-fixed frame.
 */
Eigen::Matrix3d computePartialDerivativeOfBodyFixedSphericalHarmonicAcceleration(
        const Eigen::Vector3d& cartesianPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const Eigen::MatrixXd cosineHarmonicCoefficients,
        const Eigen::MatrixXd sineHarmonicCoefficients,
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache );

//! Calculate partial of spherical harmonic acceleration w.r.t. a set of cosine coefficients
/*!
 *  Calculate partial of spherical harmonic acceleration w.r.t. a set of cosine coefficients
 *  \param sphericalPosition Spherical coordinate of body undergoing acceleration in frame fixed to body exerting
 *  acceleration, as radius, latitude, longitude.
 *  \param referenceRadius Reference radius of spherical harmonic potential.
 *  \param gravitionalParameter Gravitational parameter used for spherical harmonic expansion
 *  \param sphericalHarmonicsCache Cache object containing precomputed spherical harmonics terms.
 *  \param blockIndices List of cosine coefficient indices wrt which the partials are to be taken (first and second
 *  are degree and order for each vector entry).
 *  \param sphericalToCartesianGradientMatrix Matrix to convert (by premultiplication) a spherical gradient to a Cartesian
 *  gradient
 *  \param bodyFixedToIntegrationFrame Matrix to rotate from body-fixed to integration frame.
 *  \param partialsMatrix Partials of spherical harmonic acceleration w.r.t. to requested set of cosine coefficients
 *  (returned by reference).
 */
void calculateSphericalHarmonicGravityWrtCCoefficients(
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const std::vector< std::pair< int, int > >& blockIndices,
        const Eigen::Matrix3d& sphericalToCartesianGradientMatrix,
        const Eigen::Matrix3d& bodyFixedToIntegrationFrame,
        Eigen::MatrixXd& partialsMatrix );

//! Calculate partial of spherical harmonic acceleration w.r.t. a set of sine coefficients
/*!
 *  Calculate partial of spherical harmonic acceleration w.r.t. a set of sine coefficients
 *  \param sphericalPosition Spherical coordinate of body undergoing acceleration in frame fixed to body exerting
 *  acceleration, as radius, latitude, longitude.
 *  \param referenceRadius Reference radius of spherical harmonic potential.
 *  \param gravitionalParameter Gravitational parameter used for spherical harmonic expansion
 *  \param sphericalHarmonicsCache Cache object containing precomputed spherical harmonics terms.
 *  \param blockIndices List of cosine coefficient indices wrt which the partials are to be taken (first and second
 *  are degree and order for each vector entry).
 *  \param sphericalToCartesianGradientMatrix Matrix to convert (by premultiplication) a spherical gradient to a Cartesian
 *  gradient
 *  \param bodyFixedToIntegrationFrame Matrix to rotate from body-fixed to integration frame.
 *  \param partialsMatrix Partials of spherical harmonic acceleration w.r.t. to requested set of sine coefficients
 *  (returned by reference).
 */
void calculateSphericalHarmonicGravityWrtSCoefficients(
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const std::vector< std::pair< int, int > >& blockIndices,
        const Eigen::Matrix3d& sphericalToCartesianGradientMatrix,
        const Eigen::Matrix3d& bodyFixedToIntegrationFrame,
        Eigen::MatrixXd& partialsMatrix  );

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_SPHERICALHARMONICPARTIALFUNCTIONS_H
