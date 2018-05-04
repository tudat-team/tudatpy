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

// Template function
template< typename VariableType >
class AdditionIdentity
{
public:

    static VariableType getZeroValue( );

};

// Signed integers types
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

// Unsigned integers types
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

// Floating-point types
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

// Vector types
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
class AdditionIdentity< Eigen::Vector5d >
{
public:

    static Eigen::Vector5d getZeroValue( )
    {
        return Eigen::Vector5d::Zero( );
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
class AdditionIdentity< Eigen::Vector7d >
{
public:

    static Eigen::Vector7d getZeroValue( )
    {
        return Eigen::Vector7d::Zero( );
    }

};

// Matrix types
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

} // namespace tudat

#endif // TUDAT_ADDITION_IDENTITIES_H
