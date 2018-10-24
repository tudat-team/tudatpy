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

#include <type_traits>

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
class IdentityElement
{
public:

    //! Function to output the zero value (i.e., the addition identity) for Eigen types.
    /*!
     *  Function to output the zero value (i.e., the addition identity) for Eigen types.
     *  \return Addition identity of Eigen types.
     */
    template< typename VariableType, typename std::enable_if< is_eigen_matrix< VariableType >::value, int >::type = 0 >
    static VariableType getAdditionIdentity( )
    {
        return VariableType::Zero(
                    ( VariableType::RowsAtCompileTime > 0 ) ? VariableType::RowsAtCompileTime : 0,
                    ( VariableType::ColsAtCompileTime > 0 ) ? VariableType::ColsAtCompileTime : 0 );
    }

    //! Function to output the zero value (i.e., the addition identity) for integer and floating point types.
    /*!
     *  Function to output the zero value (i.e., the addition identity) for integer and floating point types.
     *  \return Addition identity of integer and floating point types.
     */
    template< typename VariableType, typename std::enable_if< ( std::is_integral< VariableType >::value ||
                                                              std::is_floating_point< VariableType >::value ), int >::type = 0 >
    static VariableType getAdditionIdentity( )
    {
        return tudat::mathematical_constants::getFloatingInteger< VariableType >( 0 );
    }

    //! Function to output the unit value (i.e., the multiplication identity) for Eigen types.
    /*!
     *  Function to output the unit value (i.e., the multiplication identity) for Eigen types.
     *  \return Multiplication identity of Eigen types.
     */
    template< typename VariableType, typename std::enable_if< is_eigen_matrix< VariableType >::value, int >::type = 0 >
    static VariableType getMultiplicationIdentity( )
    {
        if( VariableType::RowsAtCompileTime == VariableType::ColsAtCompileTime )
        {
            return VariableType::Identity(
                        VariableType::RowsAtCompileTime < 0 ? 0 : VariableType::RowsAtCompileTime,
                        VariableType::ColsAtCompileTime < 0 ? 0 : VariableType::ColsAtCompileTime );
        }
        else
        {
            throw std::runtime_error( "Error, multiplication identity not defined for non-square matrix" );
        }
    }

    //! Function to output the unit value (i.e., the multiplication identity) for integer and floating point types.
    /*!
     *  Function to output the unit value (i.e., the multiplication identity) for integer and floating point types.
     *  \return Multiplication identity of integer and floating point types.
     */
    template< typename VariableType, typename std::enable_if< ( std::is_integral< VariableType >::value ||
                                                              std::is_floating_point< VariableType >::value ), int >::type = 0 >
    static VariableType getMultiplicationIdentity( )
    {
        return tudat::mathematical_constants::getFloatingInteger< VariableType >( 1 );
    }

    //! Function to output the NaN value (i.e., the null identity) for Eigen types.
    /*!
     *  Function to output the NaN value (i.e., the null identity) for Eigen types.
     *  \return Null identity of Eigen types.
     */
    template< typename VariableType, typename std::enable_if< is_eigen_matrix< VariableType >::value, int >::type = 0 >
    static VariableType getNullIdentity( )
    {
        return VariableType::Constant( ( VariableType::RowsAtCompileTime > 0 ) ? VariableType::RowsAtCompileTime : 0,
                                       ( VariableType::ColsAtCompileTime > 0 ) ? VariableType::ColsAtCompileTime : 0, TUDAT_NAN );
    }

    //! Function to output the NaN value (i.e., the null identity) for floating point types.
    /*!
     *  Function to output the NaN value (i.e., the null identity) for floating point types.
     *  \return Null identity of floating point types.
     */
    template< typename VariableType, typename std::enable_if< std::is_floating_point< VariableType >::value, int >::type = 0 >
    static VariableType getNullIdentity( )
    {
        return TUDAT_NAN;
    }

};

} // namespace tudat

#endif // TUDAT_IDENTITY_ELEMENTS_H
