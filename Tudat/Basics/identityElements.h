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

#ifndef TUDAT_IDENTITY_ELEMENTS_H
#define TUDAT_IDENTITY_ELEMENTS_H

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{

//! Interface class to output the identity elements for any type.
/*!
 *  Interface class to output the identity elements for any type. This class is then specialized for
 *  each type needed in Tudat.
 */
template< typename VariableType >
class IdentityElement
{
public:

    //! Function to output the zero value (i.e., the addition identity) for the specified type.
    /*!
     *  Function to output the zero value (i.e., the addition identity) for the specified type.
     *  \return Addition identity of the specified type.
     */
    static VariableType getAdditionIdentity( );

    //! Function to output the unit value (i.e., the multiplication identity) for the specified type.
    /*!
     *  Function to output the unit value (i.e., the multiplication identity) for the specified type.
     *  \return Multiplication identity of the specified type.
     */
    static VariableType getMultiplicationIdentity( );

    //! Function to output the NaN value (i.e., the null identity) for the specified type.
    /*!
     *  Function to output the NaN value (i.e., the null identity) for the specified type.
     *  \return Null identity of the specified type.
     */
    static VariableType getNullIdentity( );

};

//! Classes for signed integer types.
/*!
 *  Classes for signed integer types.
 */
// Class for scalar of integer type.
template< >
class IdentityElement< int >
{
public:

    static int getAdditionIdentity( )
    {
        return 0;
    }

    static int getMultiplicationIdentity( )
    {
        return 1;
    }

    static int getNullIdentity( )
    {
        return static_cast< int >( TUDAT_NAN );
    }

};

// Class for scalar of long integer type.
template< >
class IdentityElement< long >
{
public:

    static long getAdditionIdentity( )
    {
        return 0L;
    }

    static long getMultiplicationIdentity( )
    {
        return 1L;
    }

    static long getNullIdentity( )
    {
        return static_cast< long >( TUDAT_NAN );
    }

};

// Class for scalar of long long integer type.
template< >
class IdentityElement< long long >
{
public:

    static long long getAdditionIdentity( )
    {
        return 0LL;
    }

    static long long getMultiplicationIdentity( )
    {
        return 1LL;
    }

    static long long getNullIdentity( )
    {
        return static_cast< long long >( TUDAT_NAN );
    }

};

//! Classes for unsigned integer types.
/*!
 *  Classes for unsigned integer types.
 */
// Class for scalar of unsigned integer type.
template< >
class IdentityElement< unsigned int >
{
public:

    static unsigned int getAdditionIdentity( )
    {
        return 0;
    }

    static unsigned int getMultiplicationIdentity( )
    {
        return 1;
    }

    static unsigned int getNullIdentity( )
    {
        return static_cast< unsigned int >( TUDAT_NAN );
    }

};

// Class for scalar of unsigned long integer type.
template< >
class IdentityElement< unsigned long >
{
public:

    static unsigned long getAdditionIdentity( )
    {
        return 0L;
    }

    static unsigned long getMultiplicationIdentity( )
    {
        return 1L;
    }

    static unsigned long getNullIdentity( )
    {
        return static_cast< unsigned long >( TUDAT_NAN );
    }

};

// Class for scalar of unsigned long long integer type.
template< >
class IdentityElement< unsigned long long >
{
public:

    static unsigned long long getAdditionIdentity( )
    {
        return 0LL;
    }

    static unsigned long long getMultiplicationIdentity( )
    {
        return 1LL;
    }

    static unsigned long long getNullIdentity( )
    {
        return static_cast< unsigned long long >( TUDAT_NAN );
    }

};

//! Classes for floating-point types.
/*!
 *  Classes for floating-point types.
 */
// Class for scalar of float type.
template< >
class IdentityElement< float >
{
public:

    static float getAdditionIdentity( )
    {
        return 0.0;
    }

    static float getMultiplicationIdentity( )
    {
        return 1.0;
    }

    static float getNullIdentity( )
    {
        return static_cast< float >( TUDAT_NAN );
    }

};

// Class for scalar of double type.
template< >
class IdentityElement< double >
{
public:

    static double getAdditionIdentity( )
    {
        return 0.0;
    }

    static double getMultiplicationIdentity( )
    {
        return 1.0;
    }

    static double getNullIdentity( )
    {
        return static_cast< double >( TUDAT_NAN );
    }

};

// Class for scalar of long double type.
template< >
class IdentityElement< long double >
{
public:

    static long double getAdditionIdentity( )
    {
        return 0.0L;
    }

    static long double getMultiplicationIdentity( )
    {
        return 1.0L;
    }

    static long double getNullIdentity( )
    {
        return static_cast< long double >( TUDAT_NAN );
    }

};

//! Classes for vector types.
/*!
 *  Classes for vector types.
 */
// Class for column vector of length 1 of double type.
template< >
class IdentityElement< Eigen::Vector1d >
{
public:

    static Eigen::Vector1d getAdditionIdentity( )
    {
        return Eigen::Vector1d::Zero( );
    }

    static Eigen::Vector1d getNullIdentity( )
    {
        return Eigen::Vector1d::Constant( TUDAT_NAN );
    }

};

// Class for column vector of length 2 of double type.
template< >
class IdentityElement< Eigen::Vector2d >
{
public:

    static Eigen::Vector2d getAdditionIdentity( )
    {
        return Eigen::Vector2d::Zero( );
    }

    static Eigen::Vector2d getNullIdentity( )
    {
        return Eigen::Vector2d::Constant( TUDAT_NAN );
    }

};

// Class for column vector of length 3 of double type.
template< >
class IdentityElement< Eigen::Vector3d >
{
public:

    static Eigen::Vector3d getAdditionIdentity( )
    {
        return Eigen::Vector3d::Zero( );
    }

    static Eigen::Vector3d getNullIdentity( )
    {
        return Eigen::Vector3d::Constant( TUDAT_NAN );
    }

};

// Class for column vector of length 4 of double type.
template< >
class IdentityElement< Eigen::Vector4d >
{
public:

    static Eigen::Vector4d getAdditionIdentity( )
    {
        return Eigen::Vector4d::Zero( );
    }

    static Eigen::Vector4d getNullIdentity( )
    {
        return Eigen::Vector4d::Constant( TUDAT_NAN );
    }

};

// Class for column vector of length 5 of double type.
template< >
class IdentityElement< Eigen::Vector5d >
{
public:

    static Eigen::Vector5d getAdditionIdentity( )
    {
        return Eigen::Vector5d::Zero( );
    }

    static Eigen::Vector5d getNullIdentity( )
    {
        return Eigen::Vector5d::Constant( TUDAT_NAN );
    }

};

// Class for column vector of length 6 of integer type.
template< >
class IdentityElement< Eigen::Vector6i >
{
public:

    static Eigen::Vector6i getAdditionIdentity( )
    {
        return Eigen::Vector6i::Zero( );
    }

    static Eigen::Vector6i getNullIdentity( )
    {
        return Eigen::Vector6i::Constant( TUDAT_NAN );
    }

};

// Class for column vector of length 6 of float type.
template< >
class IdentityElement< Eigen::Vector6f >
{
public:

    static Eigen::Vector6f getAdditionIdentity( )
    {
        return Eigen::Vector6f::Zero( );
    }

    static Eigen::Vector6f getNullIdentity( )
    {
        return Eigen::Vector6f::Constant( TUDAT_NAN );
    }

};

// Class for column vector of length 6 of double type.
template< >
class IdentityElement< Eigen::Vector6d >
{
public:

    static Eigen::Vector6d getAdditionIdentity( )
    {
        return Eigen::Vector6d::Zero( );
    }

    static Eigen::Vector6d getNullIdentity( )
    {
        return Eigen::Vector6d::Constant( TUDAT_NAN );
    }

};

// Class for column vector of length 6 of long double type.
template< >
class IdentityElement< Eigen::Vector6ld >
{
public:

    static Eigen::Vector6ld getAdditionIdentity( )
    {
        return Eigen::Vector6ld::Zero( );
    }

    static Eigen::Vector6ld getNullIdentity( )
    {
        return Eigen::Vector6ld::Constant( TUDAT_NAN );
    }

};

// Class for column vector of length 7 of double type.
template< >
class IdentityElement< Eigen::Vector7d >
{
public:

    static Eigen::Vector7d getAdditionIdentity( )
    {
        return Eigen::Vector7d::Zero( );
    }

    static Eigen::Vector7d getNullIdentity( )
    {
        return Eigen::Vector7d::Constant( TUDAT_NAN );
    }

};

// Class for column vector of length 7 of long double type.
template< >
class IdentityElement< Eigen::Vector7ld >
{
public:

    static Eigen::Vector7ld getAdditionIdentity( )
    {
        return Eigen::Vector7ld::Zero( );
    }

    static Eigen::Vector7ld getNullIdentity( )
    {
        return Eigen::Vector7ld::Constant( TUDAT_NAN );
    }

};

//! Classes for matrix types.
/*!
 *  Classes for matrix types.
 */
// Class for square 2-by-2 matrix of double type.
template< >
class IdentityElement< Eigen::Matrix2d >
{
public:

    static Eigen::Matrix2d getAdditionIdentity( )
    {
        return Eigen::Matrix2d::Zero( );
    }

    static Eigen::Matrix2d getMultiplicationIdentity( )
    {
        return Eigen::Matrix2d::Identity( );
    }

    static Eigen::Matrix2d getNullIdentity( )
    {
        return Eigen::Matrix2d::Constant( TUDAT_NAN );
    }

};

// Class for square 3-by-3 matrix of double type.
template< >
class IdentityElement< Eigen::Matrix3d >
{
public:

    static Eigen::Matrix3d getAdditionIdentity( )
    {
        return Eigen::Matrix3d::Zero( );
    }

    static Eigen::Matrix3d getMultiplicationIdentity( )
    {
        return Eigen::Matrix3d::Identity( );
    }

    static Eigen::Matrix3d getNullIdentity( )
    {
        return Eigen::Matrix3d::Constant( TUDAT_NAN );
    }

};

// Class for square 6-by-6 matrix of double type.
template< >
class IdentityElement< Eigen::Matrix6d >
{
public:

    static Eigen::Matrix6d getAdditionIdentity( )
    {
        return Eigen::Matrix6d::Zero( );
    }

    static Eigen::Matrix6d getMultiplicationIdentity( )
    {
        return Eigen::Matrix6d::Identity( );
    }

    static Eigen::Matrix6d getNullIdentity( )
    {
        return Eigen::Matrix6d::Constant( TUDAT_NAN );
    }

};

// Class for square 6-by-6 matrix of integer type.
template< >
class IdentityElement< Eigen::Matrix6i >
{
public:

    static Eigen::Matrix6i getAdditionIdentity( )
    {
        return Eigen::Matrix6i::Zero( );
    }

    static Eigen::Matrix6i getMultiplicationIdentity( )
    {
        return Eigen::Matrix6i::Identity( );
    }

    static Eigen::Matrix6i getNullIdentity( )
    {
        return Eigen::Matrix6i::Constant( TUDAT_NAN );
    }

};

// Class for square 6-by-6 matrix of float type.
template< >
class IdentityElement< Eigen::Matrix6f >
{
public:

    static Eigen::Matrix6f getAdditionIdentity( )
    {
        return Eigen::Matrix6f::Zero( );
    }

    static Eigen::Matrix6f getMultiplicationIdentity( )
    {
        return Eigen::Matrix6f::Identity( );
    }

    static Eigen::Matrix6f getNullIdentity( )
    {
        return Eigen::Matrix6f::Constant( TUDAT_NAN );
    }

};

//! Classes for dynamic matrix types.
/*!
 *  Classes for dynamic matrix types.
 */
// Class for dynamic square matrix of integer type.
template< >
class IdentityElement< Eigen::MatrixXi >
{
public:

    static Eigen::MatrixXi getAdditionIdentity( )
    {
        Eigen::MatrixXi zeroMatrix;
        return zeroMatrix.setZero( );
    }

    static Eigen::MatrixXi getMultiplicationIdentity( )
    {
        Eigen::MatrixXi identityMatrix;
        return identityMatrix.setIdentity( );
    }

    static Eigen::MatrixXi getNullIdentity( )
    {
        Eigen::MatrixXi nanMatrix;
        return nanMatrix.setConstant( TUDAT_NAN );
    }

};

// Class for dynamic square matrix of long integer type.
template< >
class IdentityElement< Eigen::MatrixXl >
{
public:

    static Eigen::MatrixXl getAdditionIdentity( )
    {
        Eigen::MatrixXl zeroMatrix;
        return zeroMatrix.setZero( );
    }

    static Eigen::MatrixXl getMultiplicationIdentity( )
    {
        Eigen::MatrixXl identityMatrix;
        return identityMatrix.setIdentity( );
    }

    static Eigen::MatrixXl getNullIdentity( )
    {
        Eigen::MatrixXl nanMatrix;
        return nanMatrix.setConstant( TUDAT_NAN );
    }

};

// Class for dynamic square matrix of long long integer type.
template< >
class IdentityElement< Eigen::MatrixXll >
{
public:

    static Eigen::MatrixXll getAdditionIdentity( )
    {
        Eigen::MatrixXll zeroMatrix;
        return zeroMatrix.setZero( );
    }

    static Eigen::MatrixXll getMultiplicationIdentity( )
    {
        Eigen::MatrixXll identityMatrix;
        return identityMatrix.setIdentity( );
    }

    static Eigen::MatrixXll getNullIdentity( )
    {
        Eigen::MatrixXll nanMatrix;
        return nanMatrix.setConstant( TUDAT_NAN );
    }

};

// Class for dynamic square matrix of float type.
template< >
class IdentityElement< Eigen::MatrixXf >
{
public:

    static Eigen::MatrixXf getAdditionIdentity( )
    {
        Eigen::MatrixXf zeroMatrix;
        return zeroMatrix.setZero( );
    }

    static Eigen::MatrixXf getMultiplicationIdentity( )
    {
        Eigen::MatrixXf identityMatrix;
        return identityMatrix.setIdentity( );
    }

    static Eigen::MatrixXf getNullIdentity( )
    {
        Eigen::MatrixXf nanMatrix;
        return nanMatrix.setConstant( TUDAT_NAN );
    }

};

// Class for dynamic square matrix of double type.
template< >
class IdentityElement< Eigen::MatrixXd >
{
public:

    static Eigen::MatrixXd getAdditionIdentity( )
    {
        Eigen::MatrixXd zeroMatrix;
        return zeroMatrix.setZero( );
    }

    static Eigen::MatrixXd getMultiplicationIdentity( )
    {
        Eigen::MatrixXd identityMatrix;
        return identityMatrix.setIdentity( );
    }

    static Eigen::MatrixXd getNullIdentity( )
    {
        Eigen::MatrixXd nanMatrix;
        return nanMatrix.setConstant( TUDAT_NAN );
    }

};

// Class for dynamic square matrix of long double type.
template< >
class IdentityElement< Eigen::MatrixXld >
{
public:

    static Eigen::MatrixXld getAdditionIdentity( )
    {
        Eigen::MatrixXld zeroMatrix;
        return zeroMatrix.setZero( );
    }

    static Eigen::MatrixXld getMultiplicationIdentity( )
    {
        Eigen::MatrixXld identityMatrix;
        return identityMatrix.setIdentity( );
    }

    static Eigen::MatrixXld getNullIdentity( )
    {
        Eigen::MatrixXld nanMatrix;
        return nanMatrix.setConstant( TUDAT_NAN );
    }

};

// Class for vector of double type.
template< >
class IdentityElement< Eigen::VectorXd >
{
public:

    static Eigen::VectorXd getAdditionIdentity( )
    {
        Eigen::VectorXd zeroMatrix;
        return zeroMatrix.setZero( );
    }

    static Eigen::VectorXd getNullIdentity( )
    {
        Eigen::VectorXd nanMatrix;
        return nanMatrix.setConstant( TUDAT_NAN );
    }

};

// Class for vector of long double type.
template< >
class IdentityElement< Eigen::VectorXld >
{
public:

    static Eigen::VectorXld getAdditionIdentity( )
    {
        Eigen::VectorXld zeroMatrix;
        return zeroMatrix.setZero( );
    }

    static Eigen::VectorXld getNullIdentity( )
    {
        Eigen::VectorXld nanMatrix;
        return nanMatrix.setConstant( TUDAT_NAN );
    }

};

// Class for row vector of double type.
template< >
class IdentityElement< Eigen::RowVectorXd >
{
public:

    static Eigen::RowVectorXd getAdditionIdentity( )
    {
        Eigen::RowVectorXd zeroMatrix;
        return zeroMatrix.setZero( );
    }

    static Eigen::RowVectorXd getNullIdentity( )
    {
        Eigen::RowVectorXd nanMatrix;
        return nanMatrix.setConstant( TUDAT_NAN );
    }

};

} // namespace tudat

#endif // TUDAT_IDENTITY_ELEMENTS_H
