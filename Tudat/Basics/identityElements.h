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
#include "Tudat/Basics/tudatTypeTraits.h"
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
    template< typename std::enable_if< is_eigen_matrix< VariableType >::value, int >::type = 0 >
    static VariableType getAdditionIdentity( )
    {        
        return VariableType::Zero( VariableType::RowsAtCompileTime, VariableType::ColsAtCompileTime );
    }

    //! Function to output the unit value (i.e., the multiplication identity) for the specified type.
    /*!
     *  Function to output the unit value (i.e., the multiplication identity) for the specified type.
     *  \return Multiplication identity of the specified type.
     */
    template< typename std::enable_if< is_eigen_matrix< VariableType >::value, int >::type = 0 >
    static VariableType getMultiplicationIdentity( )
    {
        if( VariableType::RowsAtCompileTime == VariableType::ColsAtCompileTime )
        {
        return VariableType::Identity( VariableType::RowsAtCompileTime, VariableType::ColsAtCompileTime );
        }
        else
        {
            throw std::runtime_error( "Error, multiplication identity not defined for non-square matrix" );
        }
    }

    //! Function to output the NaN value (i.e., the null identity) for the specified type.
    /*!
     *  Function to output the NaN value (i.e., the null identity) for the specified type.
     *  \return Null identity of the specified type.
     */
    template< typename std::enable_if< is_eigen_matrix< VariableType >::value, int >::type = 0 >
    static VariableType getNullIdentity( )
    {
        return VariableType::Constant( VariableType::RowsAtCompileTime, VariableType::ColsAtCompileTime , TUDAT_NAN );
    }

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

} // namespace tudat

#endif // TUDAT_IDENTITY_ELEMENTS_H
