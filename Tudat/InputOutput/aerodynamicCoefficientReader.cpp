/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/InputOutput/aerodynamicCoefficientReader.h"

namespace tudat
{

namespace input_output
{

//! Function to merge three double multi-arrays of 1 dimension into a single Vector3d multi-array
boost::multi_array< Eigen::Vector3d, 1 > mergeOneDimensionalCoefficients(
        const boost::multi_array< double, 1 > xComponents,
        const boost::multi_array< double, 1 > yComponents,
        const boost::multi_array< double, 1 > zComponents )
{
    boost::multi_array< Eigen::Vector3d, 1 > vectorArray;

    // Check input consistency
    if( !( xComponents.shape( )[ 0 ] == yComponents.shape( )[ 0 ] ) ||
            !( xComponents.shape( )[ 0 ] == zComponents.shape( )[ 0 ] ) )
    {
        throw std::runtime_error( "Error when creating 1-D merged multi-array, input sizes are inconsistent" );
    }

    // Create Vector3d multi array
    vectorArray.resize( boost::extents[ xComponents.shape( )[ 0 ] ] );
    for( unsigned int i = 0; i < xComponents.shape( )[ 0 ]; i++ )
    {
        vectorArray[ i ] = ( Eigen::Vector3d( ) << xComponents[ i ], yComponents[ i ], zComponents[ i ] ).finished( );
    }

    return vectorArray;
}

//! Function to merge three double multi-arrays of 2 dimensions into a single Vector3d multi-array
boost::multi_array< Eigen::Vector3d, 2 > mergeTwoDimensionalCoefficients(
        const boost::multi_array< double, 2 > xComponents,
        const boost::multi_array< double, 2 > yComponents,
        const boost::multi_array< double, 2 > zComponents )
{
    boost::multi_array< Eigen::Vector3d, 2 > vectorArray;

    // Check input consistency
    for( unsigned int i = 0; i < 2; i++ )
    {
        if( !( xComponents.shape( )[ i ] == yComponents.shape( )[ i ] ) ||
                !( xComponents.shape( )[ i ] == zComponents.shape( )[ i ] ) )
        {
            throw std::runtime_error( "Error when creating 1-D merged multi-array, input sizes are inconsistent" );
        }
    }

    // Create Vector3d multi array
    vectorArray.resize( boost::extents[ xComponents.shape( )[ 0 ] ][ xComponents.shape( )[ 1 ] ] );
    for( unsigned int i = 0; i < xComponents.shape( )[ 0 ]; i++ )
    {
        for( unsigned int j = 0; j < xComponents.shape( )[ 1 ]; j++ )
        {
            vectorArray[ i ][ j ] = ( Eigen::Vector3d( ) << xComponents[ i ][ j ], yComponents[ i ][ j ], zComponents[ i ][ j ] ).finished( );
        }
    }

    return vectorArray;
}

//! Function to merge three double multi-arrays of 3 dimensions into a single Vector3d multi-array
boost::multi_array< Eigen::Vector3d, 3 > mergeThreeDimensionalCoefficients(
        const boost::multi_array< double, 3 > xComponents,
        const boost::multi_array< double, 3 > yComponents,
        const boost::multi_array< double, 3 > zComponents )
{
    boost::multi_array< Eigen::Vector3d, 3 > vectorArray;

    // Check input consistency
    for( unsigned int i = 0; i < 3; i++ )
    {
        if( !( xComponents.shape( )[ i ] == yComponents.shape( )[ i ] ) ||
                !( xComponents.shape( )[ i ] == zComponents.shape( )[ i ] ) )
        {
            throw std::runtime_error( "Error when creating 1-D merged multi-array, input sizes are inconsistent" );
        }
    }

    // Create Vector3d multi array
    vectorArray.resize( boost::extents[ xComponents.shape( )[ 0 ] ][ xComponents.shape( )[ 1 ] ][ xComponents.shape( )[ 2 ] ] );
    for( unsigned int i = 0; i < xComponents.shape( )[ 0 ]; i++ )
    {
        for( unsigned int j = 0; j < xComponents.shape( )[ 1 ]; j++ )
        {
            for( unsigned int k = 0; k < xComponents.shape( )[ 2 ]; k++ )
            {
                vectorArray[ i ][ j ][ k ] = ( Eigen::Vector3d( ) << xComponents[ i ][ j ][ k ], yComponents[ i ][ j ][ k ], zComponents[ i ][ j ][ k ] ).finished( );
            }
        }
    }
    return vectorArray;
}

//! Function to compare if two lists of aerodynamic coefficient independent variables are equal
bool compareIndependentVariables( const std::vector< std::vector< double > >& list1,
                                  const std::vector< std::vector< double > >& list2 )
{
    // Check if sizes of outer vectors are equal
    bool isEqual = 1;
    if( list1.size( ) != list2.size( ) )
    {
        isEqual = 0;
    }

    // Check if sizes of inner vectors are equal
    if( isEqual )
    {
        for( unsigned int i = 0; i < list1.size( ); i++ )
        {
            if( list1.at( i ).size( ) != list2.at( i ).size( ) )
            {
                isEqual = 0;
            }
            else
            {
                for( unsigned int j = 0; j < list1.at( i ).size( ); j++ )
                {
                    if( list1.at( i ).at( j ) != list2.at( i ).at( j ) )
                    {
                        isEqual = 0;
                    }
                }
            }
        }
    }

    return isEqual;
}

}

}
