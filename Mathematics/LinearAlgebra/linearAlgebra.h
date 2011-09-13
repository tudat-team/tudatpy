/*! \file linearAlgebra.h
 *    This file contains include statements for common Eigen components,
 *    typedef for vectors and a number of useful vector operation definitions.
 *
 *    Path              : /Mathematics/LinearAlgebra/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Date created      : 7 August, 2009
 *    Last modified     : 30 September, 2010
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      090807    J. Melman         First creation of code.
 *      100930    D. Dirkx          Modified to comply with Tudat standards.
 *      100930    J. Melman         Implemented namespace, minor comment
 *                                  changes.
 */

#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H

// Notice that coefficient access methods in Eigen have assertions
// checking the ranges. So if you do a lot of coefficient access, 
// these assertions can have an important cost. If you want to 
// save cost, define EIGEN_NO_DEBUG, and it won't check assertions.
//#ifndef EIGEN_NO_DEBUG
//#define EIGEN_NO_DEBUG
//#endif

// Include statements.
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/QR>
#include <Eigen/Cholesky>

// Import most common Eigen types.
USING_PART_OF_NAMESPACE_EIGEN

//! Linear algebra namespace.
/*!
 *  Linear algebra namespace.
 */
namespace linear_algebra
{

//! Determine the cosine of the angle between two vectors.
/*!
 * Function to determine the cosine of the angle between two vectors;
 * both vectors must have non-zero norm.
 */
double determineCosineOfAngleBetweenVectors( const Vector3d& vector0,
                                             const Vector3d& vector1 );

//! Determine the angle between two vectors.
/*!
 * Function to determine the angle between two vectors;
 * both vectors must have non-zero norm.
 */
double determineAngleBetweenVectors( const Vector3d& vector0,
                                     const Vector3d& vector1 );

//! Determine the average of the components of a vector.
/*!
 * Function to determine the average (arithmetic mean) of the components of a
 * vector.
 */
double determineAverageOfVectorComponents( const VectorXd& vector0 );

//! Determine the standard deviation of the components of a vector.
/*!
 * Function to determine the standard deviation of the components of a vector.
 */
double determineStandardDeviationOfVectorComponents( const VectorXd& vector0 );
    
}

#endif // LINEAR_ALGEBRA_H

// End of file.

