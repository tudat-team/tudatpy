/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <Tudat/Basics/basicTypedefs.h>

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

Eigen::Vector3d evaluateSecondBlockInStateVector(
        const boost::function< Eigen::Vector6d( const double ) > stateFunction,
        const double time );

double computeNormOfVectorDifference( const Eigen::Vector3d& vector0,
                                      const Eigen::Vector3d& vector1 );

double getConditionNumberOfInformationMatrix( const Eigen::MatrixXd informationMatrix );

double getConditionNumberOfDecomposedMatrix( const Eigen::JacobiSVD< Eigen::MatrixXd >& singularValueDecomposition );

Eigen::JacobiSVD< Eigen::MatrixXd > getSVDDecompositionOfInformationMatrix( const Eigen::MatrixXd& informationMatrix );

Eigen::VectorXd solveSystemOfEquationsWithSvd( const Eigen::MatrixXd matrixToInvert,
                                               const Eigen::VectorXd rightHandSideVector,
                                               const bool checkConditionNumber = 1,
                                               const double maximumAllowedConditionNumber = 1.0E-8 );

Eigen::MatrixXd multiplyInformationMatrixByDiagonalWeightMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix );

Eigen::MatrixXd calculateCovarianceMatrixWithConsiderParameters(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix,
        const Eigen::MatrixXd& considerInformationMatrix,
        const Eigen::MatrixXd& considerCovarianceMatrix );

Eigen::MatrixXd calculateInverseOfUpdatedCovarianceMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix );

Eigen::MatrixXd calculateInverseOfUpdatedCovarianceMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix );

std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromInformationMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& observationResiduals,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix,
        const Eigen::VectorXd& aPrioriAdjustmentEstimate,
        const bool checkConditionNumber = 1,
        const double maximumAllowedConditionNumber = 1.0E8  );

std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromInformationMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& observationResiduals,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const bool checkConditionNumber = 1,
        const double maximumAllowedConditionNumber = 1.0E8  );

std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromInformationMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& observationResiduals,
        const bool checkConditionNumber = 1,
        const double maximumAllowedConditionNumber = 1.0E8 );

Eigen::VectorXd getLeastSquaresPolynomialFit(
        const Eigen::VectorXd& independentValues,
        const Eigen::VectorXd& dependentValues,
        const std::vector< double >& polynomialPowers );

std::vector< double > getLeastSquaresPolynomialFit(
        const std::map< double, double >& independentDependentValueMap,
        const std::vector< double >& polynomialPowers );

} // namespace linear_algebra

} // namespace tudat

#endif // TUDAT_LINEAR_ALGEBRA_H
