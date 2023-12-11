/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <unordered_map>
#include <map>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <iomanip>

#include <functional>
#include <boost/multi_array.hpp>
#include <memory>

#include <Eigen/Core>

namespace tudat
{

namespace utilities
{

//! Function to recalculate map keys as a linear function of original map keys.
/*!
 *  Function to recalculate map keys as a linear function of original map keys, i.e. new map key is constant * old key - offset or
 *  ( old key - offset ) * constant, where the choise between these two is provided by an input boolean.
 *  \param originalMap Orignal, unscaled map
 *  \param offset Offset that is to be applied to (subtracted from) map keys
 *  \param scale Value by which existing map keys are to be scaled (either before or agter application of offset variabled, depending on value
 *  of isOffsetAppliedFirst input variable)
 *  \param isOffsetAppliedFirst Boolean denoting order in which offset and scale are to be applied to existing map keys.
 */
template< typename S, typename T >
std::map< S, T > linearlyScaleKeyOfMap(
        const std::map< S, T >& originalMap, S offset, S scale, bool isOffsetAppliedFirst = true )
{
    // Declare new map
    std::map< S, T > scaledMap;
    double newKey;

    // Iterate over old map and calculate new key values.
    for( typename std::map< S, T >::const_iterator mapIterator = originalMap.begin( ); mapIterator != originalMap.end( ); mapIterator++ )
    {
        // Determine order in which modifications are to be applied
        if( isOffsetAppliedFirst )
        {
            newKey = ( mapIterator->first - offset ) * scale;
        }
        else
        {
            newKey = mapIterator->first * scale - offset;
        }

        // Set scalted map key with corresponding value in new map
        scaledMap[ newKey ] = mapIterator->second;
    }
    return scaledMap;
}

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
         mapIterator != inputMap.end( ); mapIterator++, currentIndex++ )
    {
        outputVector[ currentIndex ] = mapIterator->second;
    }

    return outputVector;
}

template< typename VectorArgument, typename KeyType >
std::vector< KeyType > createVectorFromUnorderedMapKeys( const std::unordered_map< KeyType, VectorArgument >& inputMap )
{
    // Create and size return vector.
    std::vector< KeyType > outputVector;
    outputVector.resize( inputMap.size( ) );

    // Iterate over all map entries and create vector
    int currentIndex = 0;
    for( typename std::unordered_map< KeyType, VectorArgument >::const_iterator mapIterator = inputMap.begin( );
         mapIterator != inputMap.end( ); mapIterator++, currentIndex++ )
    {
        outputVector[ currentIndex ] = mapIterator->first;
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
         mapIterator != inputMap.end( ); mapIterator++, currentIndex++ )
    {
        outputVector[ currentIndex ] = mapIterator->first;
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
S sumFunctionReturn( const std::function< S( ) > function1, const std::function< S( ) > function2 )
{
    return function1( ) + function2( );
}

//! Function to subtract the return values of two boost function with empty input argument list.
/*!
 * Function to subtract the return values of two boost function with empty input argument list.
 * \param function1 First function to be subtracted from.
 * \param function2 Second function to be subtracted.
 * \return Return values of function1 - return value of function2
 */
template< typename S >
S subtractFunctionReturn( const std::function< S( ) > function1, const std::function< S( ) > function2 )
{
    return function1( ) - function2( );
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
void printMapContents( const std::map< S, T >& mapToPrint )
{
    for( typename std::map< S, T >::const_iterator mapIterator = mapToPrint.begin( );
         mapIterator != mapToPrint.end( ); mapIterator++ )
    {
        std::cout << mapIterator->first << ", " << mapIterator->second << std::endl;
    }
}

//! Function to cast a map of Eigen matrices from one key/matrix scalar type set to another set
/*!
 *  Function to produce a map of Eigen matrices, cast from one set of key/matrix scalar type set to another set.
 *  \param originalMap Map in original types
 *  \param newTypesMap Map that is to be created (returned by reference).
 */
template< typename S, typename T, typename U, typename V, int Rows, int Columns >
void castMatrixMap( const std::map< S, Eigen::Matrix< T, Rows, Columns > >& originalMap,
                          std::map< U, Eigen::Matrix< V, Rows, Columns > >& newTypesMap )
{
    newTypesMap.clear( );
    for( typename std::map< S, Eigen::Matrix< T, Rows, Columns > >::const_iterator mapIterator = originalMap.begin( );
         mapIterator != originalMap.end( ); mapIterator++ )
    {
        newTypesMap[ static_cast< U >( mapIterator->first ) ] = mapIterator->second.template cast< V >( );
    }
}

//! Function to dynamic cast vector of shared pointers of one type to shared pointers of another type.
/*!
 *  Function to dynamic cast vector of shared pointers of one type (S) to shared pointers of another type (T). The dynamic cast must be permissible,
 *  i.e. an S pointer must succesfully dynamic cast to a T pointer (T shoudl typically derive from S).
 *  \param originalVector Vector of S shared pointers.
 *  \return Dynamic casted vector of T shared pointers.
 */
template< typename S, typename T >
std::vector< std::shared_ptr< T > >dynamicCastSVectorToTVector( const std::vector< std::shared_ptr< S > >& originalVector )
{
    std::vector< std::shared_ptr< T > > castVector;

    // Iterate over all entries and perform dynamic cast for each entry.
    for( unsigned int i = 0; i < originalVector.size( ); i++ )
    {
        castVector.push_back( std::dynamic_pointer_cast< T >( originalVector.at( i ) ) );
    }

    return castVector;
}

//! Function to concatenate matrix values of map.
template< typename KeyType, typename ScalarType, int NumberOfRows, int NumberOfColumns = 1 >
Eigen::Matrix< ScalarType, Eigen::Dynamic, NumberOfColumns > createConcatenatedEigenMatrixFromMapValues(
        const std::map< KeyType, Eigen::Matrix< ScalarType, NumberOfRows, NumberOfColumns > >& inputMap )
{
    // Create and size return vector.
    Eigen::Matrix< ScalarType, Eigen::Dynamic, NumberOfColumns > outputVector;

    int columns;
    if( NumberOfColumns != Eigen::Dynamic )
    {
        columns = NumberOfColumns;
    }
    else
    {
        columns = inputMap.begin( )->second.cols( );
    }

    if( NumberOfRows != Eigen::Dynamic )
    {
        int counter = 0;

        outputVector.resize( inputMap.size( ) * NumberOfRows, columns );
        for( typename std::map< KeyType, Eigen::Matrix< ScalarType, NumberOfRows, NumberOfColumns > >::const_iterator mapIterator =
             inputMap.begin( ); mapIterator != inputMap.end( ); mapIterator++ )
        {
            outputVector.block( counter, 0, NumberOfRows, columns ) = mapIterator->second;
            counter += NumberOfRows;
        }
    }
    else
    {
        int concatenatedSize = 0;

        for( typename std::map< KeyType, Eigen::Matrix< ScalarType, NumberOfRows, NumberOfColumns > >::const_iterator mapIterator =
             inputMap.begin( ); mapIterator != inputMap.end( ); mapIterator++ )
        {
            concatenatedSize += mapIterator->second.rows( );
        }

        outputVector.resize( concatenatedSize, columns );
        int currentSize;
        int counter = 0;

        for( typename std::map< KeyType, Eigen::Matrix< ScalarType, NumberOfRows, NumberOfColumns > >::const_iterator mapIterator =
             inputMap.begin( ); mapIterator != inputMap.end( ); mapIterator++ )
        {
            currentSize = mapIterator->second.rows( );
            outputVector.block( counter, 0, currentSize, columns ) = mapIterator->second;
            counter += currentSize;
        }
    }

    return outputVector;
}

//! Function to extract both keys and values from map, and output them as a pair.
template< typename KeyType, typename ScalarType, int NumberOfRows >
std::pair< Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 >, Eigen::Matrix< ScalarType, Eigen::Dynamic, Eigen::Dynamic > >
extractKeyAndValuesFromMap( const std::map< KeyType, Eigen::Matrix< ScalarType, NumberOfRows, 1 > >& inputMap )
{
    // Declare eventual output variables
    Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > keyValuesVector;
    Eigen::Matrix< ScalarType, Eigen::Dynamic, Eigen::Dynamic > mappedValuesMatrix;

    // Make compatible with Eigen::Dynamic vectors
    bool dynamicInput = ( NumberOfRows == Eigen::Dynamic );
    unsigned int resizingDimension = dynamicInput ? inputMap.begin( )->second.rows( ) : NumberOfRows;

    // Assign size to matrices
    unsigned int numberOfKeys = inputMap.size( );
    keyValuesVector.resize( numberOfKeys, 1 );
    mappedValuesMatrix.resize( resizingDimension, numberOfKeys );

    // Loop over map and save elements
    int i = 0;
    for ( typename std::map< KeyType, Eigen::Matrix< ScalarType, NumberOfRows, 1 > >::const_iterator
          mapIterator = inputMap.begin( ); mapIterator != inputMap.end( ); mapIterator++, i++ )
    {
        keyValuesVector[ i ] = mapIterator->first;
        mappedValuesMatrix.col( i ) = mapIterator->second;
    }

    // Give output
    return std::make_pair( keyValuesVector, mappedValuesMatrix );
}

//! Function to convert Eigen::Vector to std::vector.
template< typename T >
std::vector< T > convertEigenVectorToStlVector( const Eigen::Matrix< T, Eigen::Dynamic, 1 >& eigenVector )
{
    std::vector< T > stlVector;
    stlVector.resize( eigenVector.rows( ) );
    for( int i = 0; i < eigenVector.rows( ); i++ )
    {
        stlVector[ i ] = eigenVector( i );
    }
    return stlVector;
}

//! Function to convert std::vector to Eigen::Vector.
template< typename T >
Eigen::Matrix< T, Eigen::Dynamic, 1 > convertStlVectorToEigenVector( const std::vector< T >& stlVector )
{
    Eigen::Matrix< T, Eigen::Dynamic, 1 > eigenVector = Eigen::Matrix< T, Eigen::Dynamic, 1 >::Zero( stlVector.size( ) );
    for( unsigned int i = 0; i < stlVector.size( ); i++ )
    {
        eigenVector[ i ] = stlVector[ i ];
    }
    return eigenVector;
}

//! Function to convert std::vector to Eigen::Matrix.
template< typename T, int Rows = Eigen::Dynamic >
Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > convertStlVectorToEigenMatrix(
        const std::vector< Eigen::Matrix< T, Rows, 1 > >& stlVector )
{
    // Create empty vector
    int numberOfRows = Rows;
    if ( numberOfRows == Eigen::Dynamic )
    {
        numberOfRows = stlVector.front( ).rows( );
    }
    Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > eigenMatrix =
            Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic >::Zero( numberOfRows, stlVector.size( ) );

    // Assign values
    for( unsigned int i = 0; i < stlVector.size( ); i++ )
    {
        eigenMatrix.col( i ) = stlVector.at( i );
    }
    return eigenMatrix;
}

//! Function to add a double to all entries in an STL vector.
/*!
 *  Function to add a double to all entries in an STL vector (addition of a double must be defined for T type).
 *  \param vector Vector to which a double is to be added.
 *  \param scalar Value that is to be added to vector
 *  \return New vector with scalar added to all entries of input vector.
 */
template< typename T >
std::vector< T > addScalarToVector( const std::vector< T >& vector, const double scalar )
{
    // Declare and resize return vector.
    std::vector< T > addedVector;
    addedVector.resize( vector.size( ) );

    // Add scalar to all entries of input vector
    for( unsigned int i = 0; i < vector.size( ); i++ )
    {
        addedVector[ i ] = vector[ i ] + scalar;
    }

    return addedVector;
}

//! Function to copy a multi-array into another multi-array
/*!
 *  Function to copy a multi-array into another multi-array, resizing the new multi-array accordingly
 *  \param arrayToCopy Multi-array that is to be copied
 *  \param targetArray New multi-array into which arrayToCopy is to be copied (returned by reference).
 */
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

//! Get index in nth direction of pointer to single entry in multi-array of doubles
/*!
 *  Get index in nth direction of pointer to single entry in multi-array of doubles
 *  \param multiArray Multi-array for which the index is to be retrieved
 *  \param requestedElement Pointer to element for which index is to be retrieved
 *  \param direction Dimension of multi-array for which index is to be retrieved
 *  \return Index in nth direction of pointer to single entry in multi-array of doubles
 */
template< unsigned int NumberOfDimensions >
typename boost::multi_array< double, NumberOfDimensions >::index getMultiArrayIndex(
        const typename boost::multi_array< double, NumberOfDimensions >& multiArray, const double* requestedElement,
        const unsigned short int direction )
{
    int offset = requestedElement - multiArray.origin( );
    return( offset / multiArray.strides( )[ direction ] % multiArray.shape( )[ direction ] +
            multiArray.index_bases( )[ direction ] );
}

//! Get indices of pointer to single entry in multi-array (size 1) of doubles
/*!
 *  Get indices of pointer to single entry in multi-array (size 1) of doubles
 *  \param multiArray Multi-array for which the index is to be retrieved
 *  \param requestedElement Pointer to element for which index is to be retrieved
 *  \return Indices of pointer to single entry in multi-array of doubles
 */
boost::array< boost::multi_array< double, 1 >::index, 1 > getMultiArrayIndexArray(
        const boost::multi_array< double, 1 >& multiArray, const double* requestedElement );

//! Get indices of pointer to single entry in multi-array (size 2) of doubles
/*!
 *  Get indices of pointer to single entry in multi-array (size 2) of doubles
 *  \param multiArray Multi-array for which the index is to be retrieved
 *  \param requestedElement Pointer to element for which index is to be retrieved
 *  \return Indices of pointer to single entry in multi-array of doubles
 */
boost::array< boost::multi_array< double, 2 >::index, 2 > getMultiArrayIndexArray(
        const boost::multi_array< double, 2 >& multiArray, const double* requestedElement );

//! Get indices of pointer to single entry in multi-array (size 3) of doubles
/*!
 *  Get indices of pointer to single entry in multi-array (size 3) of doubles
 *  \param multiArray Multi-array for which the index is to be retrieved
 *  \param requestedElement Pointer to element for which index is to be retrieved
 *  \return Indices of pointer to single entry in multi-array of doubles
 */
boost::array< boost::multi_array< double, 3 >::index, 3 > getMultiArrayIndexArray(
        const boost::multi_array< double, 3 >& multiArray, const double* requestedElement );

template< typename S, typename T >
std::vector< S > createVectorFromVectorOfPairFirsts( const std::vector< std::pair< S, T > > inputVector )
{
    std::vector< S > outputVector;

    for( unsigned int i = 0; i < inputVector.size( ); i++ )
    {
        outputVector.push_back( inputVector.at( i ).first );
    }
    return outputVector;
}

template< typename T, typename S >
T evaluateFunctionWithoutInputArgumentDependency( std::function< T( ) > inputFreeFunction, const S dummyInput )
{
    return inputFreeFunction( );
}

//! Function to get the order in which the input vector would be sorted (in ascending order)
/*!
 *  Function to get the order in which the input vector would be sorted (in ascending order). Example: for inout vector
 *  (5,2,6,7,4,0), output would be (5,1,4,0,2,3).
 *  \param unsortedVector Vector of which the sort order is to be determined
 *  \return Order in which the input vector would be sorted (in ascending order)
 */
template< typename T >
std::vector< int > getSortOrderOfVector( const std::vector< T > unsortedVector )
{
    return getSortOrderOfVectorAndSortedVector( unsortedVector ).first;
}

//! Function to create a vector from the values of a multimap
/*!
 *  Function to create a vector from the values of a multimap. The output vector is in the order of the multimap entries, i.e. as provided by a
 *  forward iterator. The multimap keys are not used for the return vector.
 *  \param inputMap Original multimap from which the vector is to be created
 *  \return Vector created from the multimap values
 */
template< typename VectorArgument, typename KeyType >
std::vector< VectorArgument > createVectorFromMultiMapValues( const std::multimap< KeyType, VectorArgument >& inputMap )
{
    // Create and size return vector.
    std::vector< VectorArgument > outputVector;
    outputVector.resize( inputMap.size( ) );

    // Iterate over all map entries and create vector
    int currentIndex = 0;
    for( typename std::multimap< KeyType, VectorArgument >::const_iterator mapIterator = inputMap.begin( );
         mapIterator != inputMap.end( ); mapIterator++ )
    {
        outputVector[ currentIndex ] = mapIterator->second;
        currentIndex++;
    }

    return outputVector;
}

//! Function to create a vector from the keys of a multimap
/*!
 *  Function to create a vector from the keys of a multimap. The output vector is in the order of the multimap entries, i.e. as provided by a
 *  forward iterator. The multimap values are not used for the return vector.
 *  \param inputMap Original multimap from which the vector is to be created
 *  \return Vector created from the multimap keys
 */
template< typename VectorArgument, typename KeyType >
std::vector< KeyType > createVectorFromMultiMapKeys( const std::multimap< KeyType, VectorArgument >& inputMap )
{
    // Create and size return vector.
    std::vector< KeyType > outputVector;
    outputVector.resize( inputMap.size( ) );

    // Iterate over all map entries and create vector
    int currentIndex = 0;
    for( typename std::multimap< KeyType, VectorArgument >::const_iterator mapIterator = inputMap.begin( );
         mapIterator != inputMap.end( ); mapIterator++ )
    {
        outputVector[ currentIndex ] = mapIterator->first;
        currentIndex++;
    }

    return outputVector;
}

//! Function to get sorted vector of an input vector, as well as the order in which this input has been be sorted (ascending)
/*!
 *  Function to get sorted vector of an input vector, as well as the order in which this input has been be sorted (ascending)).
 *  Example: for inout vector (5,2,6,7,4,0), output would be [(5,1,4,0,2,3), (0,2,4,5,6,7)].
 *  \param unsortedVector Vector that is to be sorted
 *  \return Parit, with first: order in which the input vector is sorted (in ascending order), second: sorted input vector
 */
template< typename T >
std::pair< std::vector< int >, std::vector< T > > getSortOrderOfVectorAndSortedVector( const std::vector< T > unsortedVector )
{
    std::multimap< T, int > sortMap;
    for( unsigned int i = 0; i < unsortedVector.size( ); i++ )
    {
        sortMap.insert( std::pair< T, int >( unsortedVector[ i ], i ) );
    }

    return std::make_pair( createVectorFromMultiMapValues( sortMap ), createVectorFromMultiMapKeys( sortMap ) );
}

template< typename T >
bool doStlVectorContentsMatch(
        const std::vector< T >& vectorA, const std::vector< T >& vectorB )
{
    bool doVectorsMatch = true;
    if( vectorA.size( ) != vectorB.size( ) )
    {
        doVectorsMatch = false;
    }
    else
    {
        for( unsigned int i = 0; i < vectorA.size( ); i++ )
        {
            if( std::count( vectorB.begin( ), vectorB.end( ), vectorA.at( i ) ) !=
                    std::count( vectorA.begin( ), vectorA.end( ), vectorA.at( i ) ))
            {
                doVectorsMatch = false;
            }
        }
    }

    return doVectorsMatch;
}

//! Transform from map of std::vector (output of text file reader) to map of Eigen::Array
template< typename MapKey, typename ScalarType >
std::map< MapKey, Eigen::Array< ScalarType, Eigen::Dynamic, 1 > > convertStlVectorMapToEigenVectorMap(
        std::map< MapKey, std::vector< ScalarType > > stlVectorMap )
{
    std::map< MapKey, Eigen::Array< ScalarType, Eigen::Dynamic, 1 > > eigenMap;
    for ( auto ent: stlVectorMap )
    {
        Eigen::Array< ScalarType, Eigen::Dynamic, 1 > array( ent.second.size( ) );
        for ( int i = 0; i < array.rows( ); i++ )
        {
            array.row( i ) = ent.second.at( i );
        }
        eigenMap[ ent.first ] = array;
    }
    return eigenMap;
}

//! Function to slice standard library vector, given an optional initial and final slicing values.
template< typename T >
std::vector< T > sliceStlVector( const std::vector< T >& vectorToBeSliced, unsigned int startIndex = 0,
                                 unsigned int endIndex = std::numeric_limits< unsigned int >::signaling_NaN( ) )
{
    // Declare output vector
    std::vector< T > slicedVector;

    // Give value to end index
    if ( endIndex == std::numeric_limits< unsigned int >::signaling_NaN( ) )
    {
        endIndex = vectorToBeSliced.size( ) - 1;
    }

    // Check that boundaries make sense
    if ( startIndex > endIndex )
    {
        // Warn user of inconsistency
        std::cerr << "Warning in slicing of std::vector. The starting index is greater than the end index. "
                     "The indices will be swapped." << std::endl;

        // Swap indices
        unsigned int temporaryIndex = startIndex;
        startIndex = endIndex;
        endIndex = temporaryIndex;
    }

    // Transfer values to sliced array
    for ( unsigned int i = startIndex; i < ( endIndex + 1 ); i++ )
    {
        slicedVector.push_back( vectorToBeSliced.at( i ) );
    }
    return slicedVector;
}

template< typename KeyType, typename ScalarType, int FixedSize >
void castDynamicToFixedSizeEigenVectorMap(
        std::map< KeyType, Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > >& originalMap,
        std::map< KeyType, Eigen::Matrix< ScalarType, FixedSize, 1 > >& fixedSizeMap )
{
    for( auto mapIterator : originalMap )
    {
        fixedSizeMap[ mapIterator.first ] = mapIterator.second;
    }
}

//! Function to return the sign (+1 or -1) of a variable of type T
/*!
 *  \param val Variable for which sign is to be determined
 *  \return Sign of variable
 */
template < typename T > int sgn( const T& val )
{
    return ( T( 0 ) < val ) - ( val < T( 0 ) );
}

//! From https://stackoverflow.com/questions/27028226/python-linspace-in-c
template<typename T>
std::vector< T > linspace(T start_in, T end_in, int num_in)
{

    std::vector<double> linspaced;

    T start = static_cast< T >(start_in);
    T end = static_cast< T >(end_in);
    T num = static_cast< T >(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end);
    return linspaced;
}

template< typename KeyType, typename ScalarType >
std::map< KeyType, Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > > sliceMatrixHistory(
        const std::map< KeyType, Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > >& fullHistory,
        const std::pair< int, int > sliceStartIndexAndSize )
{
    std::map< KeyType, Eigen::Matrix< ScalarType, Eigen::Dynamic, 1 > > slicedHistory;

    for( auto mapIterator : fullHistory )
    {
        slicedHistory[ mapIterator.first ] = mapIterator.second.segment(
                    sliceStartIndexAndSize.first, sliceStartIndexAndSize.second );
    }
    return slicedHistory;
}

template< typename KeyType, typename ValueType >
std::map< KeyType, ValueType > concatenateMaps(
        const std::vector< std::map< KeyType, ValueType > >& mapList )
{
    std::map< KeyType, ValueType > concatenatedMaps;

    for( unsigned int i = 0; i < mapList.size( ); i++ )
    {
        std::map< KeyType, ValueType > mapToInsert = mapList.at( i );
        concatenatedMaps.insert( mapToInsert.begin( ), mapToInsert.end( ) );
    }
    return concatenatedMaps;
}

template< typename KeyType, typename ValueType >
std::map< KeyType, ValueType > createMapFromVectors(
        const std::vector< KeyType >& keyList,
        const std::vector< ValueType >& valueList )
{
    if( keyList.size( ) != valueList.size( ) )
    {
        throw std::runtime_error( "Error when creating map for keys and values, sizes are incompatible" );
    }
    std::map< KeyType, ValueType > createdMap;
    for( unsigned int i = 0; i < keyList.size( ); i++ )
    {
        createdMap[ keyList.at( i ) ] = valueList.at( i );
    }
    return createdMap;

}

template< typename KeyType, typename ScalarType, int SizeType >
std::map< KeyType, ScalarType > getSingleVectorEntryHistory(
        const std::map< KeyType, Eigen::Matrix< ScalarType, SizeType, 1 > > originalMap, int vectorEntry )
{
    std::map< KeyType, ScalarType > extractedMap;
    for( auto mapIterator : originalMap )
    {
        extractedMap[ mapIterator.first ] = mapIterator.second( vectorEntry );
    }
    return extractedMap;
}

template< typename ValueType >
std::vector< ValueType > getVectorEntries( const std::vector< ValueType >& fullVector, const std::vector< int > indices )
{
    std::vector< ValueType > subVector;
    subVector.resize( indices.size( ) );
    for( unsigned int i = 0; i < indices.size( ); i++ )
    {
        subVector[ i ] = fullVector.at( indices.at( i ) );
    }
    return subVector;
}


template< typename A >
std::shared_ptr< A > deepcopyPointer(
        const std::shared_ptr< A > originalPointer )
{
    return std::make_shared< A >( *originalPointer );
}

template< typename A >
std::vector< std::shared_ptr< A > > deepcopyDuplicatePointers(
        const std::vector< std::shared_ptr< A > > originalPointers )
{
    std::vector< std::shared_ptr< A > > noDuplicatePointers;
    for( unsigned int i = 0; i < originalPointers.size( ); i++ )
    {
        if( std::find( noDuplicatePointers.begin( ), noDuplicatePointers.end( ), originalPointers.at( i ) ) == noDuplicatePointers.end( ) )
        {
            noDuplicatePointers.push_back( originalPointers.at( i ) );
        }
        else
        {
            noDuplicatePointers.push_back( std::make_shared< A >( *originalPointers.at( i ) ) );
        }
    }
    return noDuplicatePointers;
}

template< typename A >
std::vector< std::shared_ptr< A > > cloneDuplicatePointers(
        const std::vector< std::shared_ptr< A > > originalPointers )
{
    std::vector< std::shared_ptr< A > > noDuplicatePointers;
    for( unsigned int i = 0; i < originalPointers.size( ); i++ )
    {
        if( std::find( noDuplicatePointers.begin( ), noDuplicatePointers.end( ), originalPointers.at( i ) ) == noDuplicatePointers.end( ) )
        {
            noDuplicatePointers.push_back( originalPointers.at( i ) );
        }
        else
        {
            noDuplicatePointers.push_back( originalPointers.at( i )->clone( ) );
        }
    }
    return noDuplicatePointers;
}

template< typename A, typename B >
std::vector< std::pair< A, B > > mergeVectorsIntoVectorOfPairs(
        const std::vector< A >& firstVector,
        const std::vector< B >& secondVector )
{
    std::vector< std::pair< A, B > > mergedVector;
    if( firstVector.size( ) != secondVector.size( ) )
    {
        throw std::runtime_error( "Error when merging vectors in vector of pairs; size is inconsistent" );
    }
    for( unsigned  int i = 0; i < firstVector.size( ); i++ )
    {
        mergedVector.push_back( std::make_pair( firstVector.at( i ), secondVector.at( i ) ) );
    }
    return mergedVector;
}

template< typename A, typename... ArgTypes >
std::vector< A > getFirstTupleEntryVector( const std::vector< std::tuple< A, ArgTypes... > >& tupleVector )
{
    std::vector< A > vectorOfFirstEntries;
    for( unsigned int i = 0; i < tupleVector.size( ); i++ )
    {
        vectorOfFirstEntries.push_back( std::get< 0 >( tupleVector.at( i ) ) );
    }
    return vectorOfFirstEntries;
}

template< typename T >
constexpr int constexpr_int_floor(const T& f)
{
    const int i = static_cast<int>(f);
    return f < i ? i - 1 : i;
}

template< typename IntType >
std::string paddedZeroIntString(
    const IntType valueToConvert,
    const int numberOfDigits )
{
    std::stringstream intStream;
    intStream << std::setw( numberOfDigits ) << std::setfill( '0' ) << valueToConvert;
    std::string intString = intStream.str( );
    return intString;
}

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return std::move(out).str();
}
template <typename T>
std::vector< T > getStlVectorSegment( const std::vector< T > originalVector, const int startIndex, const int size )
{
    typename std::vector< T >::const_iterator first = originalVector.begin() + startIndex;
    typename std::vector< T >::const_iterator last = originalVector.begin() + startIndex + size;
    return std::vector< T >(first, last);
}

template< typename T, typename S >
std::vector< T > staticCastVector( const std::vector< S > originalVector )
{
    std::vector< T > castVector;
    for( unsigned int i = 0; i < originalVector.size( ); i++ )
    {
        castVector.push_back( static_cast< T >( originalVector.at( i ) ) );
    }
    return castVector;
}

template <typename T>
Eigen::Matrix< T, Eigen::Dynamic, 1 > getSuccesivelyConcatenatedVector(
    const Eigen::Matrix< T, Eigen::Dynamic, 1 > baseVector, const unsigned int numberOfConcatenations )
{
    int singleSize = baseVector.rows( );
    Eigen::Matrix< T, Eigen::Dynamic, 1 > concatenatedVector = Eigen::Matrix< T, Eigen::Dynamic, 1 >::Zero(
        numberOfConcatenations * singleSize );
    for( unsigned int i = 0; i < numberOfConcatenations; i++ )
    {
        concatenatedVector.segment( i * singleSize, singleSize ) = baseVector;
    }
    return concatenatedVector;
}

template <typename T>
int countNumberOfOccurencesInVector( const std::vector< T >& vector, const T& value )
{
    int counter = 0;
    for( unsigned int i = 0; i < vector.size( ); i++ )
    {
        if( vector.at( i ) == value )
        {
            counter++;
        }
    }
    return counter;
}

} // namespace utilities

} // namespace tudat

#endif // TUDAT_UTILITIES_H
