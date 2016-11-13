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
 *      130116    K. Kumar          File created.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_LINEAR_ALGEBRA_H
#define TUDAT_LINEAR_ALGEBRA_H

#include <boost/function.hpp>

#include <Eigen/Core>

namespace tudat
{

namespace linear_algebra
{

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
Eigen::Vector3d computeVectorDifference( const Eigen::Vector3d& vector0,
                                         const Eigen::Vector3d& vector1 );

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

double computeNormOfVectorDifference( const Eigen::Vector3d& vector0,
                                      const Eigen::Vector3d& vector1 );
} // namespace linear_algebra

} // namespace tudat

#endif // TUDAT_LINEAR_ALGEBRA_H
