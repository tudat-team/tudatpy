/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_UTILITIES_H
#define TUDAT_UTILITIES_H

#include <vector>

#include <iostream>
#include <map>

#include <boost/function.hpp>
#include <boost/multi_array.hpp>

#include <Eigen/Core>

namespace tudat
{

namespace utilities
{

//! Function to create a vector from the values of a map
/*!
 *  Function to create a vector from the values of a map. The output vector is in the order of the map entries, i.e. as provided by a forward iterator.
 *  The map keys are not used for the return vector.
 *  \param inputMap Original map from which the vector is to be created
 *  \return Vector created from the map values
 */
template< typename VectorArgument, typename KeyType >
std::vector< VectorArgument > createVectorFromMapValues( const std::map< KeyType, VectorArgument >& inputMap )
{
    // Create and size return vector.
    std::vector< VectorArgument > outputVector;
    outputVector.resize( inputMap.size( ) );

    // Iterate over all map entries and create vector
    int currentIndex = 0;
    for( typename std::map< KeyType, VectorArgument >::const_iterator mapIterator = inputMap.begin( );
         mapIterator != inputMap.end( ); mapIterator++ )
    {
        outputVector[ currentIndex ] = mapIterator->second;
        currentIndex++;
    }

    return outputVector;

}

//! Function to create a vector from the keys of a map
/*!
 *  Function to create a vector from the keys of a map. The output vector is in the order of the map entries, i.e. as provided by a forward iterator.
 *  The map values are not used for the return vector.
 *  \param inputMap Original map from which the vector is to be created
 *  \return Vector created from the map keys
 */
template< typename VectorArgument, typename KeyType >
std::vector< KeyType > createVectorFromMapKeys( const std::map< KeyType, VectorArgument >& inputMap )
{
    // Create and size return vector.
    std::vector< KeyType > outputVector;
    outputVector.resize( inputMap.size( ) );

    // Iterate over all map entries and create vector
    int currentIndex = 0;
    for( typename std::map< KeyType, VectorArgument >::const_iterator mapIterator = inputMap.begin( );
         mapIterator != inputMap.end( ); mapIterator++ )
    {
        outputVector[ currentIndex ] = mapIterator->first;
        currentIndex++;
    }

    return outputVector;
}

//! Function to sum the return values of two boost function with empty input argument list.
/*!
 * Function to sum the return values of two boost function with empty input argument list.
 * \param function1 First function to be added.
 * \param function2 Second function to be added.
 * \return Sum of return values of function1 and function2
 */
template< typename S >
S sumFunctionReturn( const boost::function< S( ) > function1, const boost::function< S( ) > function2 )
{
    return function1( ) + function2( );
}

//! Function to create a vector block history from full matrix history.
/*!
 *  Function to create a vector matrix block history from full matrix history.
 *  \param matrixHistory Full matrix history
 *  \param blockMatrixHistory Block vector history (return by reference).
 *  \param startIndices Starting point (row,column) in matrix of return vector blocks.
 *  \param segmentSize Number of rows in vector.
 */
template< typename S, typename T >
void createVectorBlockMatrixHistory(
        const std::map< S, Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > >& matrixHistory,
        std::map< S, Eigen::Matrix< T, Eigen::Dynamic, 1 > >& blockMatrixHistory,
        const std::pair< int, int > startIndices, const int segmentSize )
{
    blockMatrixHistory.clear( );

    // Loop over integration output and put results in corresponding data structures.
    for( typename std::map< S, Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > >::const_iterator matrixIterator = matrixHistory.begin( );
         matrixIterator != matrixHistory.end( ); matrixIterator++ )
    {
        // Set numerically integrated states of bodies.
        blockMatrixHistory[ matrixIterator->first ] =
                matrixIterator->second.block( startIndices.first, startIndices.second, segmentSize, 1 );
    }
}

//! Function to print the contents of a map, line by line
/*!
 *  Function to print the contents of a map, line by line. Both the key and value types must have the << operator defined
 *  \param mapToPrint Map that is to be printed.
 */
template< typename S, typename T >
void printMapContents( const std::map< S, T >& mapToPrint)
{
    for( typename std::map< S, T >::const_iterator mapIterator = mapToPrint.begin( );
         mapIterator != mapToPrint.end( ); mapIterator++ )
    {
        std::cout<<mapIterator->first<<", "<<mapIterator->second<<std::endl;
    }
}

template< typename S, int NumberOfDimensions >
void copyMultiArray( const boost::multi_array< S, NumberOfDimensions >& arrayToCopy,
                     boost::multi_array< S, NumberOfDimensions >& targetArray )
{
    std::vector< size_t > ex;
    const size_t* shape = arrayToCopy.shape( );
    ex.assign( shape, shape + arrayToCopy.num_dimensions( ) );
    targetArray.resize( ex );
    targetArray = arrayToCopy;
}


template< unsigned int NumberOfDimensions >
typename boost::multi_array< double ,NumberOfDimensions >::index getMultiArrayIndex(
        const typename boost::multi_array< double, NumberOfDimensions >& m, const double* requestedElement,
        const unsigned short int direction)
{
    int offset = requestedElement - m.origin( );
    return( offset / m.strides( )[ direction] % m.shape( )[ direction ] +  m.index_bases( )[direction] );
}

boost::array< boost::multi_array< double, 1 >::index, 1 > getMultiArrayIndexArray(
        const boost::multi_array< double, 1 >& m, const double* requestedElement );

boost::array< boost::multi_array< double, 2 >::index, 2 > getMultiArrayIndexArray(
        const boost::multi_array< double, 2 >& m, const double* requestedElement );

boost::array< boost::multi_array< double, 3 >::index, 3 > getMultiArrayIndexArray(
        const boost::multi_array< double, 3 >& m, const double* requestedElement );

} // namespace utilities

} // namespace tudat

#endif // TUDAT_UTILITIES_H
