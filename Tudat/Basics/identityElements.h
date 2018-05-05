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

};

// Class for semi-dynamic (1 column) matrix of double type.
template< >
class IdentityElement< Eigen::MatrixX1d >
{
public:

    static Eigen::MatrixX1d getAdditionIdentity( )
    {
        Eigen::MatrixX1d zeroMatrix;
        return zeroMatrix.setZero( );
    }

};

// Class for semi-dynamic (1 column) matrix of long double type.
template< >
class IdentityElement< Eigen::MatrixX1ld >
{
public:

    static Eigen::MatrixX1ld getAdditionIdentity( )
    {
        Eigen::MatrixX1ld zeroMatrix;
        return zeroMatrix.setZero( );
    }

};

// Class for semi-dynamic (2 columns) matrix of double type.
template< >
class IdentityElement< Eigen::MatrixX2d >
{
public:

    static Eigen::MatrixX2d getAdditionIdentity( )
    {
        Eigen::MatrixX2d zeroMatrix;
        return zeroMatrix.setZero( );
    }

};

// Class for semi-dynamic (3 columns) matrix of double type.
template< >
class IdentityElement< Eigen::MatrixX3d >
{
public:

    static Eigen::MatrixX3d getAdditionIdentity( )
    {
        Eigen::MatrixX3d zeroMatrix;
        return zeroMatrix.setZero( );
    }

};

// Class for semi-dynamic (1 row) matrix of double type.
template< >
class IdentityElement< Eigen::Matrix1Xd >
{
public:

    static Eigen::Matrix1Xd getAdditionIdentity( )
    {
        Eigen::Matrix1Xd zeroMatrix;
        return zeroMatrix.setZero( );
    }

};

// Class for semi-dynamic (2 rows) matrix of double type.
template< >
class IdentityElement< Eigen::Matrix2Xd >
{
public:

    static Eigen::Matrix2Xd getAdditionIdentity( )
    {
        Eigen::Matrix2Xd zeroMatrix;
        return zeroMatrix.setZero( );
    }

};

// Class for semi-dynamic (3 rows) matrix of double type.
template< >
class IdentityElement< Eigen::Matrix3Xd >
{
public:

    static Eigen::Matrix3Xd getAdditionIdentity( )
    {
        Eigen::Matrix3Xd zeroMatrix;
        return zeroMatrix.setZero( );
    }

};

} // namespace tudat

#endif // TUDAT_IDENTITY_ELEMENTS_H
