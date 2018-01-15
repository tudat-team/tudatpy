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

#ifndef TUDAT_TRIAXIALELLIPSOIDGRAVITY_H
#define TUDAT_TRIAXIALELLIPSOIDGRAVITY_H

#include <Eigen/Core>

namespace tudat
{

namespace gravitation
{

//! Function to calculate (non-normalized) cosine spherical harmonic coefficient for a homogeneous
//! triaxial ellipsoid
/*!
 *  Function to calculate (non-normalized) cosine spherical harmonic coefficient for a homogeneous
 *  triaxial ellipsoid, from method of Boyce (1997) at a given degree and order. X-axis is alligned
 *  with largest axis, y-axis with middle axis and z-axis with smallest axis
 *  \param aSquaredMinusCSquared Square of largest axis minus square of shortest axis
 *  \param bSquaredMinusCSquared Square of middle axis minus square of shortest axis
 *  \param referenceRadius Reference radius of spherical harmonic gravity field.
 *  \param degree Degree of coefficient
 *  \param order Order of coefficient
 *  \return Spherical harmonic cosine coefficient (non-normalized) of triaxial ellipsoid at
 *  requested degree and order.
 */
double calculateCosineTermForTriaxialEllipsoidSphericalHarmonicGravity(
        const double aSquaredMinusCSquared, const double bSquaredMinusCSquared,
        const double referenceRadius, const int degree, const int order );

//! Function to calculate triaxial ellipsoid reference radius
/*!
 *  Function to calculate triaxial ellipsoid reference radius of spherical harmonic expansion,
 *  per Balmino (1994)
 *  \param axisA Largest axis of triaxial ellipsoid
 *  \param axisB Middle axis of triaxial ellipsoid
 *  \param axisC Smallest axis of triaxial ellipsoid
 *  \return Triaxial ellipsoid reference radius of spherical harmonic expansion
 */
double calculateTriAxialEllipsoidReferenceRadius(
        const double axisA, const double axisB, const double axisC );

//! Function to calculate triaxial ellipsoid volume
/*!
 *  Function to calculate triaxial ellipsoid volume
 *  \param axisA Largest axis of triaxial ellipsoid
 *  \param axisB Middle axis of triaxial ellipsoid
 *  \param axisC Smallest axis of triaxial ellipsoid
 *  \return Volume of triaxial ellipsoid.
 */
double calculateTriAxialEllipsoidVolume(
        const double axisA, const double axisB, const double axisC );

//! Function to calculate (non-normalized) cosine spherical harmonic coefficients for a homogeneous
//! triaxial ellipsoid
/*!
 *  Function to calculate (non-normalized) cosine spherical harmonic coefficients for a homogeneous
 *  triaxial ellipsoid, from method of Boyce (1997) up  to given degree and order. X-axis is
 *  aligned with largest axis, y-axis with middle axis and z-axis with smallest axis
 *  \param axisA Largest axis of triaxial ellipsoid
 *  \param axisB Middle axis of triaxial ellipsoid
 *  \param axisC Smallest axis of triaxial ellipsoid
 *  \param maximumDegree Maximum degree of expansion
 *  \param maximumOrder Maximum oredr of expansion
 *  \return Spherical harmonic cosine coefficient matrix (non-normalized) of triaxial ellipsoid up
 *  to requested degree and order.
 */
Eigen::MatrixXd createTriAxialEllipsoidSphericalHarmonicCosineCoefficients(
        const double axisA, const double axisB, const double axisC,
        const int maximumDegree, const int maximumOrder );

//! Function to calculate (non-normalized) cosine and sine spherical harmonic coefficients for a
//! homogeneous triaxial ellipsoid
/*!
 *  Function to calculate (non-normalized) cosine and sinespherical harmonic coefficients for a
 *  homogeneous triaxial ellipsoid, from method of Boyce (1997) up to given degree and order.
 *  X-axis is alligned with largest axis, y-axis with middle axis and z-axis with smallest axis.
 *  Note that all sine coefficients are zero for all homogeneous triaxila ellipsoids
 *  (in given frame).
 *  \param axisA Largest axis of triaxial ellipsoid
 *  \param axisB Middle axis of triaxial ellipsoid
 *  \param axisC Smallest axis of triaxial ellipsoid
 *  \param maximumDegree Maximum degree of expansion
 *  \param maximumOrder Maximum oredr of expansion
 *  \return Spherical harmonic cosine and coefficient matrix (non-normalized) pair of triaxial
 *  ellipsoid up to requested degree and order.
 */
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > createTriAxialEllipsoidSphericalHarmonicCoefficients(
        const double axisA, const double axisB, const double axisC,
        const int maximumDegree, const int maximumOrder );

//! Function to calculate (normalized) cosine and sine spherical harmonic coefficients for a
//! homogeneous triaxial ellipsoid
/*!
 *  Function to calculate (normalized) cosine and sinespherical harmonic coefficients for a
 *  homogeneous triaxial ellipsoid, from method of Boyce (1997) up to given degree and order.
 *  X-axis is alligned with largest axis, y-axis with middle axis and z-axis with smallest axis.
 *  Note that all sine coefficients are zero for all homogeneous triaxila ellipsoids
 *  (in given frame).
 *  \param axisA Largest axis of triaxial ellipsoid
 *  \param axisB Middle axis of triaxial ellipsoid
 *  \param axisC Smallest axis of triaxial ellipsoid
 *  \param maximumDegree Maximum degree of expansion
 *  \param maximumOrder Maximum oredr of expansion
 *  \return Spherical harmonic cosine and coefficient matrix (normalized) pair of triaxial
 *  ellipsoid up to requested degree and order.
 */
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > createTriAxialEllipsoidNormalizedSphericalHarmonicCoefficients(
        const double axisA, const double axisB, const double axisC,
        const int maximumDegree, const int maximumOrder );
}

}

#endif // TUDAT_TRIAXIALELLIPSOIDGRAVITY_H
