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

#ifndef TUDAT_LINEAR_ALGEBRA_H
#define TUDAT_LINEAR_ALGEBRA_H

#include <map>

#include <boost/function.hpp>

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Geometry>

#include "Tudat/Basics/basicTypedefs.h"
namespace tudat
{

namespace linear_algebra
{

//! Function to put a quaternion in 'vector format', e.g. a Vector4d with entries (w,x,y,z) of the quaternion
/*!
 * Function to put a quaternion in 'vector format', e.g. a Vector4d with entries (w,x,y,z) of the quaternion
 * \param quaternion Quaternion that is to be put into vector format.
 * \return Vector format of input quaternion
 */
Eigen::Vector4d convertQuaternionToVectorFormat( const Eigen::Quaterniond& quaternion );

//! Function that returns that 'cross-product matrix'
/*!
 *  Function that returns that 'cross-product matrix', i.e. for vectors a,b and c, with c = a x b, the matrix A such that
 *  c = Ab.
 *  \param leftHandVector The left-multiplying vector (a in above example)
 *  \return The matrix by which to premultiply the right-multiplying vector to obtain the cross product of the two matrices.
 */
Eigen::Matrix3d getCrossProductMatrix( const Eigen::Vector3d& leftHandVector );

//! Compute cosine of the angle between two vectors.
/*!
 * Computes the cosine of the angle between two vectors; both vectors must have non-zero norm.
 * \param vector0 First vector.
 * \param vector1 Second vector.
 * \return Cosine of angle between vectors.
 */
double computeCosineOfAngleBetweenVectors( const Eigen::VectorXd& vector0,
                                           const Eigen::VectorXd& vector1 );

//! Compute angle between two vectors.
/*!
 * Computes the angle between two vectors; both vectors must have non-zero norm.
 * \param vector0 First vector.
 * \param vector1 Second vector.
 * \return Angle between vectors.
 */
double computeAngleBetweenVectors( const Eigen::VectorXd& vector0,
                                   const Eigen::VectorXd& vector1 );

//! Computes the difference between two 3d vectors.
/*!
 * Computes the difference between two 3d vectors (first input minus second input, i.e vector from second input to
 * first input).
 * \param vector0 First vector.
 * \param vector1 Second vector.
 * \return Difference between vectors
 */
template< int VectorSize >
Eigen::Matrix< double, VectorSize, 1 > computeVectorDifference(
        const Eigen::Matrix< double, VectorSize, 1 >& vector0,
        const Eigen::Matrix< double, VectorSize, 1 >& vector1 )
{
    return ( vector0 - vector1 );
}

//! Computes norm of the the difference between two 3d vectors.
/*!
 * Computes the norm of the difference between two 3d vectors (i.e. distance between vectors)
 * \param vector0 First vector.
 * \param vector1 Second vector.
 * \return Norm of difference between vectors
 */
double computeNormOfVectorDifference( const Eigen::Vector3d& vector0,
                                      const Eigen::Vector3d& vector1 );

//! Computes the norm of a 3d vector
/*!
 * Computes the norm of a 3d vector
 * \param vector Vector for which the norm is to be computed
 * \return Vector norm
 */
double getVectorNorm( const Eigen::Vector3d& vector );

//! Computes the norm of a 3d vector from a vector-returning function.
/*!
 * Computes the norm of a 3d vector from a vector-returning function.
 * \param vectorFunction Function returning the vector for which the norm is to be computed
 * \return Vector norm
 */
double getVectorNormFromFunction( const boost::function< Eigen::Vector3d( ) > vectorFunction );

//! Flip matrix rows.
/*!
 * Flips all rows of an Eigen-matrix, i.e., order of rows is reversed.
 * \param matrixToFlip Matrix that should be flipped. The flipping is done in place, i.e. the input
 *          matrix is modified directly (no copies are made).
 */
static inline void flipMatrixRows( Eigen::MatrixXd& matrixToFlip )
{
    // Loop through each row.
    for ( int i = 0; i < static_cast< int >( matrixToFlip.rows( ) / 2 ); i++ )
    {
        // Save top row in a temporary matrix.
        const Eigen::MatrixXd temporaryRow = matrixToFlip.row( i );

        // Set top row to bottom row.
        matrixToFlip.row( i ) = matrixToFlip.row( matrixToFlip.rows( ) - 1 - i );

        // Set bottom row to temporarily stored matrix.
        matrixToFlip.row( matrixToFlip.rows( ) - 1 - i ) = temporaryRow;
    }
}

Eigen::Vector3d evaluateSecondBlockInStateVector(
        const boost::function< Eigen::Vector6d( const double ) > stateFunction,
        const double time );

double computeNormOfVectorDifference( const Eigen::Vector3d& vector0,
                                      const Eigen::Vector3d& vector1 );

//! Function to calculate the jacobian of a normalized vector, from the jacobian of the unnormalized vector.
/*!
 *  Function to calculate the jacobian (partial matrix) of a normalized vector, from the jacobian
 *  (partial matrix) of the unnormalized vector and the unnormalized vector itself, i.e. d/dp(x/|x|) from d/dp(x) and x, with p and x
 *  3d vectors.
 *  \param partialOfUnnormalizedVector The jacobian of the unnormalized vector
 *  \param unnormalizedVector Unnormalized vector wrt which partialOfUnnormalizedVector is taken
 *  \return The jacobian of the normalized vector.
 */
Eigen::Matrix3d calculatePartialOfNormalizedVector(
        const Eigen::Matrix3d& partialOfUnnormalizedVector,
        const Eigen::Vector3d& unnormalizedVector );

//! Function to check whether an Eigen Matrix has any NaN entries
/*!
 *  Function to check whether an Eigen Matrix has any NaN entries
 *  \param matrixToCheck Eigen Matrix to check for any NaN entries
 *  \return True if matrix has NaN entries, false otherwise
 */
template< typename StateScalarType, int NumberOfRows, int NumberOfColumns >
bool doesMatrixHaveNanEntries( const Eigen::Matrix< StateScalarType, NumberOfRows, NumberOfColumns > matrixToCheck )
{
    bool areNanEntriesPresent = false;
    for( int i = 0; i < matrixToCheck.rows( ); i++ )
    {
        for( int j = 0; j < matrixToCheck.cols( ); j++ )
        {
            if( !( matrixToCheck( i, j ) == matrixToCheck( i, j ) ) )
            {
                areNanEntriesPresent = true;
            }
        }
    }
    return areNanEntriesPresent;
}

//! Function to compute the root mean square value of the entries in an Eigen vector
/*!
 *  Function to compute the root mean square (RMS) value of the entries in an Eigen vector
 * \param inputVector Vector for which the RMS is to be computed
 * \return RMS of input vector
 */
double getVectorEntryRootMeanSquare( const Eigen::VectorXd& inputVector );


} // namespace linear_algebra

} // namespace tudat

#endif // TUDAT_LINEAR_ALGEBRA_H
