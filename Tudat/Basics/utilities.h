#ifndef UTILITIES_H
#define UTILITIES_H

#include <map>
#include <boost/function.hpp>
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

}

}

#endif // UTILITIES_H
