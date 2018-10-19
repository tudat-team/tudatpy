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

#include <cmath>

#include <Eigen/LU>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Basics/utilities.h"

namespace tudat
{

namespace linear_algebra
{

//! Function to put a quaternion in 'vector format', e.g. a Vector4d with entries (w,x,y,z) of the quaternion
Eigen::Vector4d convertQuaternionToVectorFormat( const Eigen::Quaterniond& quaternion )
{
    Eigen::Vector4d vector;
    vector( 0 ) = quaternion.w( );
    vector( 1 ) = quaternion.x( );
    vector( 2 ) = quaternion.y( );
    vector( 3 ) = quaternion.z( );
    return vector;
}

//! Function to put a vector in 'quaternion format', i.e. a Quaterniond.
Eigen::Quaterniond convertVectorToQuaternionFormat( const Eigen::Vector4d& vector )
{
    Eigen::Quaterniond quaternion;
    quaternion.w( ) = vector( 0 );
    quaternion.x( ) = vector( 1 );
    quaternion.y( ) = vector( 2 );
    quaternion.z( ) = vector( 3 );

    return quaternion;
}


//! Function to take the product of two quaternions.
Eigen::Vector4d quaternionProduct( const Eigen::Vector4d& firstQuaternion, const Eigen::Vector4d& secondQuaternion )
{
    Eigen::Vector4d resultantQuaternion;
    resultantQuaternion[ 0 ] = firstQuaternion[ 0 ] * secondQuaternion[ 0 ] -
            firstQuaternion.segment( 1, 3 ).dot( secondQuaternion.segment( 1, 3 ) );
    resultantQuaternion.segment( 1, 3 ) = firstQuaternion[ 0 ] * secondQuaternion.segment( 1, 3 ) +
            secondQuaternion[ 0 ] * firstQuaternion.segment( 1, 3 ) +
            getCrossProductMatrix( firstQuaternion.segment( 1, 3 ) ) * secondQuaternion.segment( 1, 3 );
    return resultantQuaternion;
}

//! Function to invert a quaternion.
void invertQuaternion( Eigen::Vector4d& quaternionVector )
{
    quaternionVector.segment( 1, 3 ) *= -1.0;
}

//! Function that returns that 'cross-product matrix'
Eigen::Matrix3d getCrossProductMatrix( const Eigen::Vector3d& vector )
{
    Eigen::Matrix3d crossProductMatrix = Eigen::Matrix3d::Zero( );
    crossProductMatrix( 1, 0 ) = vector.z( );
    crossProductMatrix( 0, 1 ) = -vector.z( );
    crossProductMatrix( 2, 0 ) = -vector.y( );
    crossProductMatrix( 0, 2 ) = vector.y( );
    crossProductMatrix( 2, 1 ) = vector.x( );
    crossProductMatrix( 1, 2 ) = -vector.x( );
    return crossProductMatrix;
}

//! Compute cosine of the angle between two vectors.
double computeCosineOfAngleBetweenVectors( const Eigen::VectorXd& vector0,
                                           const Eigen::VectorXd& vector1 )
{
    if( !( vector0.size( ) == vector1.size( ) ) )
    {
        throw std::runtime_error( "Error when computing angle between vectors; size is incompatible" );
    }

    // Get the cosine of the angle by dotting the normalized vectors.
    double dotProductOfNormalizedVectors = vector0.normalized( ).dot( vector1.normalized( ) );

    // Explicitly define the extreme cases, which can give problems with the acos function.
    if ( dotProductOfNormalizedVectors >= 1.0 )
    {
        return 1.0;
    }

    else if ( dotProductOfNormalizedVectors <= -1.0 )
    {
        return -1.0;
    }
    // Determine the actual angle.
    else
    {
        return dotProductOfNormalizedVectors;
    }
}

//! Compute angle between two vectors.
double computeAngleBetweenVectors( const Eigen::VectorXd& vector0, const Eigen::VectorXd& vector1 )
{
    // Determine the cosine of the angle by using another routine.
    double dotProductOfNormalizedVectors = computeCosineOfAngleBetweenVectors( vector0, vector1 );

    // Return arccosine of the above, which is effectively the angle.
    return std::acos( dotProductOfNormalizedVectors );
}

//! Computes norm of the the difference between two 3d vectors.
double computeNormOfVectorDifference( const Eigen::Vector3d& vector0,
                                      const Eigen::Vector3d& vector1 )
{
    return ( vector0 - vector1 ).norm( );
}

//! Computes the norm of a 3d vector
double getVectorNorm( const Eigen::Vector3d& vector )
{
    return vector.norm( );
}

Eigen::Vector3d evaluateSecondBlockInStateVector(
        const std::function< Eigen::Vector6d( const double ) > stateFunction,
        const double time )
{
    return stateFunction( time ).segment( 3, 3 );
}

//! Computes the norm of a 3d vector from a vector-returning function.
double getVectorNormFromFunction( const std::function< Eigen::Vector3d( ) > vectorFunction )
{
    return getVectorNorm( vectorFunction( ) );
}

//! Function to calculate the jacobian of a normalized vector, from the partial of the unnormalized vector.
Eigen::Matrix3d calculatePartialOfNormalizedVector( const Eigen::Matrix3d& partialOfUnnormalizedVector,
                                                    const Eigen::Vector3d& unnormalizedVector )
{
    double normOfVector = unnormalizedVector.norm( );

    return ( Eigen::Matrix3d::Identity( ) / normOfVector - unnormalizedVector * unnormalizedVector.transpose( ) /
             ( normOfVector * normOfVector * normOfVector ) ) * partialOfUnnormalizedVector;
}

//! Function to compute the root mean square value of the entries in an Eigen vector
double getVectorEntryRootMeanSquare( const Eigen::VectorXd& inputVector )
{
    // Calculate RMS for vector
    double vectorRms = 0.0;
    for( int i = 0; i < inputVector.rows( ); i++ )
    {
        vectorRms += inputVector( i ) * inputVector( i );
    }
    vectorRms = std::sqrt( vectorRms / inputVector.rows( ) );

    return vectorRms;
}

//! Function to compute the partial derivative of a rotation matrix w.r.t. its associated quaterion elements
void computePartialDerivativeOfRotationMatrixWrtQuaternion(
        const Eigen::Vector4d quaternionVector,
        std::vector< Eigen::Matrix3d >& partialDerivatives )
{
    partialDerivatives[ 0 ]<< quaternionVector( 0 ), -quaternionVector( 3 ), quaternionVector( 2 ),
            quaternionVector( 3 ), quaternionVector( 0 ), -quaternionVector( 1 ),
            -quaternionVector( 2 ), quaternionVector( 1 ), quaternionVector( 0 );
     partialDerivatives[ 0 ] *= 2.0;

    partialDerivatives[ 1 ]<< quaternionVector( 1 ), quaternionVector( 2 ), quaternionVector( 3 ),
            quaternionVector( 2 ), -quaternionVector( 1 ), -quaternionVector( 0 ),
            quaternionVector( 3 ), quaternionVector( 0 ), -quaternionVector( 1 );
    partialDerivatives[ 1 ] *= 2.0;

    partialDerivatives[ 2 ]<< -quaternionVector( 2 ), quaternionVector( 1 ), quaternionVector( 0 ),
            quaternionVector( 1 ), quaternionVector( 2 ), quaternionVector( 3 ),
            -quaternionVector( 0 ), quaternionVector( 3 ), -quaternionVector( 2 );
    partialDerivatives[ 2 ] *= 2.0;

    partialDerivatives[ 3 ]<< -quaternionVector( 3 ), -quaternionVector( 0 ), quaternionVector( 1 ),
            quaternionVector( 0 ), -quaternionVector( 3 ), quaternionVector( 2 ),
            quaternionVector( 1 ), quaternionVector( 2 ), quaternionVector( 3 );
    partialDerivatives[ 3 ] *= 2.0;

}


} // namespace linear_algebra

} // namespace tudat
