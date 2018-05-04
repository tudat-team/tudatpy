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

#ifndef TUDAT_ADDITION_IDENTITIES_H
#define TUDAT_ADDITION_IDENTITIES_H

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{

//! Interface class to output the zero value, or addition identity, for any type.
/*!
 *  Interface class to output the zero value, or addition identity, for any type. This class is then specialized for
 *  each type needed in Tudat.
 */
template< typename VariableType >
class AdditionIdentity
{
public:

    //! Function to output the zero value (i.e., the addition identity) for the specified type.
    /*!
     *  Function to output the zero value (i.e., the addition identity) for the specified type.
     *  \return Addition identity of the specified type.
     */
    static VariableType getZeroValue( );

};

//! Classes for signed integer types.
/*!
 *  Classes for signed integer types.
 */
template< >
class AdditionIdentity< int >
{
public:

    static int getZeroValue( )
    {
        return 0;
    }

};

template< >
class AdditionIdentity< long >
{
public:

    static long getZeroValue( )
    {
        return 0L;
    }

};

template< >
class AdditionIdentity< long long >
{
public:

    static long long getZeroValue( )
    {
        return 0LL;
    }

};

//! Classes for unsigned integer types.
/*!
 *  Classes for unsigned integer types.
 */
template< >
class AdditionIdentity< unsigned int >
{
public:

    static unsigned int getZeroValue( )
    {
        return 0;
    }

};

template< >
class AdditionIdentity< unsigned long >
{
public:

    static unsigned long getZeroValue( )
    {
        return 0L;
    }

};

template< >
class AdditionIdentity< unsigned long long >
{
public:

    static unsigned long long getZeroValue( )
    {
        return 0LL;
    }

};

//! Classes for floating-point types.
/*!
 *  Classes for floating-point types.
 */
template< >
class AdditionIdentity< float >
{
public:

    static float getZeroValue( )
    {
        return 0.0;
    }

};

template< >
class AdditionIdentity< double >
{
public:

    static double getZeroValue( )
    {
        return 0.0;
    }

};

template< >
class AdditionIdentity< long double >
{
public:

    static long double getZeroValue( )
    {
        return 0.0L;
    }

};

//! Classes for vector types.
/*!
 *  Classes for vector types.
 */
template< >
class AdditionIdentity< Eigen::Vector1d >
{
public:

    static Eigen::Vector1d getZeroValue( )
    {
        return Eigen::Vector1d::Zero( );
    }

};

template< >
class AdditionIdentity< Eigen::Vector2d >
{
public:

    static Eigen::Vector2d getZeroValue( )
    {
        return Eigen::Vector2d::Zero( );
    }

};

template< >
class AdditionIdentity< Eigen::Vector3d >
{
public:

    static Eigen::Vector3d getZeroValue( )
    {
        return Eigen::Vector3d::Zero( );
    }

};

template< >
class AdditionIdentity< Eigen::Vector4d >
{
public:

    static Eigen::Vector4d getZeroValue( )
    {
        return Eigen::Vector4d::Zero( );
    }

};

template< >
class AdditionIdentity< Eigen::Vector5d >
{
public:

    static Eigen::Vector5d getZeroValue( )
    {
        return Eigen::Vector5d::Zero( );
    }

};

template< >
class AdditionIdentity< Eigen::Vector6i >
{
public:

    static Eigen::Vector6i getZeroValue( )
    {
        return Eigen::Vector6i::Zero( );
    }

};

template< >
class AdditionIdentity< Eigen::Vector6f >
{
public:

    static Eigen::Vector6f getZeroValue( )
    {
        return Eigen::Vector6f::Zero( );
    }

};

template< >
class AdditionIdentity< Eigen::Vector6d >
{
public:

    static Eigen::Vector6d getZeroValue( )
    {
        return Eigen::Vector6d::Zero( );
    }

};

template< >
class AdditionIdentity< Eigen::Matrix< long double, 6, 1 > >
{
public:

    static Eigen::Matrix< long double, 6, 1 > getZeroValue( )
    {
        return Eigen::Matrix< long double, 6, 1 >::Zero( );
    }

};

template< >
class AdditionIdentity< Eigen::Vector7d >
{
public:

    static Eigen::Vector7d getZeroValue( )
    {
        return Eigen::Vector7d::Zero( );
    }

};

template< >
class AdditionIdentity< Eigen::Matrix< long double, 7, 1 > >
{
public:

    static Eigen::Matrix< long double, 7, 1 > getZeroValue( )
    {
        return Eigen::Matrix< long double, 7, 1 >::Zero( );
    }

};

//! Classes for matrix types.
/*!
 *  Classes for matrix types.
 */
template< >
class AdditionIdentity< Eigen::Matrix2d >
{
public:

    static Eigen::Matrix2d getZeroValue( )
    {
        return Eigen::Matrix2d::Zero( );
    }

};

template< >
class AdditionIdentity< Eigen::Matrix3d >
{
public:

    static Eigen::Matrix3d getZeroValue( )
    {
        return Eigen::Matrix3d::Zero( );
    }

};

template< >
class AdditionIdentity< Eigen::Matrix6d >
{
public:

    static Eigen::Matrix6d getZeroValue( )
    {
        return Eigen::Matrix6d::Zero( );
    }

};

template< >
class AdditionIdentity< Eigen::Matrix6i >
{
public:

    static Eigen::Matrix6i getZeroValue( )
    {
        return Eigen::Matrix6i::Zero( );
    }

};

template< >
class AdditionIdentity< Eigen::Matrix6f >
{
public:

    static Eigen::Matrix6f getZeroValue( )
    {
        return Eigen::Matrix6f::Zero( );
    }

};

//! Classes for dynamic matrix types.
/*!
 *  Classes for dynamic matrix types.
 */
template< >
class AdditionIdentity< Eigen::MatrixXd >
{
public:

    static Eigen::MatrixXd getZeroValue( )
    {
        Eigen::MatrixXd zeroMatrix;
        return zeroMatrix.setZero( );
    }

};

template< >
class AdditionIdentity< Eigen::Matrix< double, Eigen::Dynamic, 1 > >
{
public:

    static Eigen::Matrix< double, Eigen::Dynamic, 1 > getZeroValue( )
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > zeroMatrix;
        return zeroMatrix.setZero( );
    }

};

template< >
class AdditionIdentity< Eigen::Matrix< double, Eigen::Dynamic, 2 > >
{
public:

    static Eigen::Matrix< double, Eigen::Dynamic, 2 > getZeroValue( )
    {
        Eigen::Matrix< double, Eigen::Dynamic, 2 > zeroMatrix;
        return zeroMatrix.setZero( );
    }

};

template< >
class AdditionIdentity< Eigen::Matrix< double, Eigen::Dynamic, 3 > >
{
public:

    static Eigen::Matrix< double, Eigen::Dynamic, 3 > getZeroValue( )
    {
        Eigen::Matrix< double, Eigen::Dynamic, 3 > zeroMatrix;
        return zeroMatrix.setZero( );
    }

};

template< >
class AdditionIdentity< Eigen::Matrix< double, 1, Eigen::Dynamic > >
{
public:

    static Eigen::Matrix< double, 1, Eigen::Dynamic > getZeroValue( )
    {
        Eigen::Matrix< double, 1, Eigen::Dynamic > zeroMatrix;
        return zeroMatrix.setZero( );
    }

};

template< >
class AdditionIdentity< Eigen::Matrix< double, 2, Eigen::Dynamic > >
{
public:

    static Eigen::Matrix< double, 2, Eigen::Dynamic > getZeroValue( )
    {
        Eigen::Matrix< double, 2, Eigen::Dynamic > zeroMatrix;
        return zeroMatrix.setZero( );
    }

};

template< >
class AdditionIdentity< Eigen::Matrix< double, 3, Eigen::Dynamic > >
{
public:

    static Eigen::Matrix< double, 3, Eigen::Dynamic > getZeroValue( )
    {
        Eigen::Matrix< double, 3, Eigen::Dynamic > zeroMatrix;
        return zeroMatrix.setZero( );
    }

};

} // namespace tudat

#endif // TUDAT_ADDITION_IDENTITIES_H
