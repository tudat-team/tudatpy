/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <map>
#include "Tudat/Basics/utilities.h"

#include "Tudat/InputOutput/multiDimensionalArrayReader.h"

namespace tudat
{

namespace input_output
{

//! Function to merge three double multi-arrays of N dimension into a single Vector3d multi-array
/*!
 *  Function to merge three double multi-arrays of N dimension into a single Vector3d multi-array, where the three
 *  double multi-arrays represent the x-, y- and z-components of the Vector3ds.
 *  \param xComponents Multi-array containing the x-components of the Vector3d
 *  \param yComponents Multi-array containing the y-components of the Vector3d
 *  \param zComponents Multi-array containing the z-components of the Vector3d
 *  \return Single multi-array containing Vector3ds according to double multi-arrays.
 */
template< unsigned int NumberOfDimensions >
boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > mergeNDimensionalCoefficients(
        boost::multi_array< double, static_cast< size_t >( NumberOfDimensions ) > xComponents,
        boost::multi_array< double, static_cast< size_t >( NumberOfDimensions ) > yComponents,
        boost::multi_array< double, static_cast< size_t >( NumberOfDimensions ) > zComponents )
{
    boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > vectorArray;

    // Check input consistency
    for( unsigned int i = 0; i < NumberOfDimensions; i++ )
    {
        if( !( xComponents.shape( )[ i ] == yComponents.shape( )[ i ] ) ||
                !( xComponents.shape( )[ i ] == zComponents.shape( )[ i ] ) )
        {
            throw std::runtime_error( "Error when creating N-D merged multi-array, input sizes are inconsistent" );
        }
    }

    // Retrieve multi-array shape/size
    std::vector< size_t > sizeVector;
    const size_t* arrayShape = xComponents.shape( );
    sizeVector.assign( arrayShape, arrayShape + xComponents.num_dimensions( ) );

    // Resize coefficient multi-array
    vectorArray.resize( sizeVector );

    // Iterate over all elements and combine x,y and z-components into vector3d of associated entry in vectorVector
    int numberOfEntries = xComponents.num_elements( );
    Eigen::Vector3d* vectorVector = new Eigen::Vector3d[ numberOfEntries ] ;

    typedef typename boost::multi_array< double, NumberOfDimensions >::index tIndex;
    typedef boost::array< tIndex, NumberOfDimensions > tIndexArray;

    double* p = xComponents.data( );
    tIndexArray index;
    for( int i = 0; i < numberOfEntries; i++ )
    {
        index = utilities::getMultiArrayIndexArray( xComponents, p );

        vectorVector[ i ] = ( Eigen::Vector3d( ) << xComponents( index ), yComponents( index ), zComponents( index ) ).finished( );
        ++p;
    }

    // Assign array of entries to output multi-array
    vectorArray.assign( vectorVector, vectorVector + numberOfEntries );

    delete[ ] vectorVector;

    return vectorArray;
}

//! Function to compare if two lists of aerodynamic coefficient independent variables are equal
/*!
 * Function to compare if two lists of aerodynamic coefficient independent variables (vector of vector of doubles) are equal
 * \param list1 First list that is to be compared.
 * \param list2 Second list that is to be compared.
 * \return True of the two lists are completely equal in size and contents, false otherwise.
 */
bool compareIndependentVariables( const std::vector< std::vector< double > >& list1,
                                  const std::vector< std::vector< double > >& list2 );

//! Function to read a list of aerodynamic coefficients and associated independent variables from a set of files
/*!
 *  Function to read a list of aerodynamic coefficients of NumberOfDimensions independent variables and associated
 *  independent variables from a set of files.
 *  \param fileNames Vector of size 3, containing the file names for the x-, y- and z- components of the aerodynamic
 *  coefficients. Note that the independent variables for each components must be identical.
 *  \return  Pair: first entry containing multi-array of aerodynamic coefficients, second containing list of independent
 *  variables at which coefficients are defined.
 */
template< unsigned int NumberOfDimensions >
std::pair< boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) >,
std::vector< std::vector< double > > >
readAerodynamicCoefficients( const std::vector< std::string >& fileNames )
{
    if( fileNames.size( ) != 3 )
    {
        throw std::runtime_error( "Error when reading 1-Dimensional aeroynamic coefficients, wrong number of files" );
    }

    std::map< int, std::string > fileNameMap;
    for( unsigned int i = 0; i < 3; i++ )
    {
        fileNameMap[ i ] = fileNames.at( i );
    }

    return readAerodynamicCoefficients< NumberOfDimensions >( fileNameMap );
}

//! Function to read a list of aerodynamic coefficients and associated independent variables from a set of files
/*!
 *  Function to read a list of aerodynamic coefficients of 2 independent variables and associated independent variables
 *  from a set of files.
 *  \param fileNames Map of file names, with the key  required to be 0, 1 and/or 2. These indices denote the  x-, y- and z-
 *  components of the aerodynamic coefficients. All indices that are not provided in this map are assumed to have associated
 *  coefficients equal to zero for all values of the independent variables
 *  Note that the independent variables for each components must be identical.
 *  \return  Pair: first entry containing multi-array of aerodynamic coefficients, second containing list of independent
 *  variables at which coefficients are defined.
 */
template< unsigned int NumberOfDimensions >
std::pair< boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) >,
std::vector< std::vector< double > > >
readAerodynamicCoefficients( const std::map< int, std::string >& fileNames )
{
    std::map< int, boost::multi_array< double, static_cast< size_t >( NumberOfDimensions ) > > rawCoefficientArrays;

    std::vector< boost::multi_array< double, static_cast< size_t >( NumberOfDimensions ) > > coefficientArrays;
    std::vector< std::vector< double > > independentVariables;

    // Iterate over files and read the contents into rawCoefficientArrays/independentVariables.
    for( std::map< int, std::string >::const_iterator fileIterator = fileNames.begin( ); fileIterator != fileNames.end( );
         fileIterator++ )
    {
        // Read current coefficients/independent variables
        std::pair< boost::multi_array< double, static_cast< size_t >( NumberOfDimensions ) >,
                std::vector< std::vector< double > > > currentCoefficients =
                MultiArrayFileReader< NumberOfDimensions >::readMultiArrayAndIndependentVariables( fileIterator->second );

        // Save/check consistency of independent variables
        if( rawCoefficientArrays.size( ) == 0 )
        {
            independentVariables = currentCoefficients.second;
        }
        else
        {
            bool areIndependentVariablesEqual = compareIndependentVariables(
                        independentVariables, currentCoefficients.second );

            if( !areIndependentVariablesEqual )
            {
                throw std::runtime_error( "Error when reading 1-Dimensional aeroynamic coefficients, inconsistent independent variables." );
            }
        }

        // Save file contents into rawCoefficientArrays
        utilities::copyMultiArray< double, NumberOfDimensions >(
                    currentCoefficients.first, rawCoefficientArrays[ fileIterator->first ] );

    }

    // Check if anything has been read from files.
    if( rawCoefficientArrays.size( ) == 0 )
    {
        throw std::runtime_error( "Error when reading aerodynamic coefficients, no files read" );
    }
    else
    {
        coefficientArrays.resize( 3 );
        boost::multi_array< double, static_cast< size_t >( NumberOfDimensions ) > firstMultiArray =
                rawCoefficientArrays.begin( )->second;

        // Iterate over all 3 coefficient entries
        for( unsigned int i = 0; i < 3; i++ )
        {
            // Copy contents for current index into coefficientArrays read from file.
            if( rawCoefficientArrays.count( i ) != 0 )
            {
                utilities::copyMultiArray< double, NumberOfDimensions >(
                            rawCoefficientArrays.at( i ), coefficientArrays[ i ] );
            }
            // Set zero multi-array for current index.
            else
            {
                std::vector< size_t > sizeVector;
                const size_t* arrayShape = firstMultiArray.shape( );
                sizeVector.assign( arrayShape, arrayShape+ firstMultiArray.num_dimensions( ) );

                coefficientArrays[ i ].resize( sizeVector );

                std::fill( coefficientArrays[ i ].data( ),
                           coefficientArrays[ i ].data() + coefficientArrays[ i ].num_elements( ), 0.0 );
            }
        }
    }

    // Merge coefficient entries
    return std::make_pair(
                mergeNDimensionalCoefficients< NumberOfDimensions >(
                    coefficientArrays.at( 0 ), coefficientArrays.at( 1 ), coefficientArrays.at( 2 ) ), independentVariables );
}

} // namespace input_output

} // namespace tudat
