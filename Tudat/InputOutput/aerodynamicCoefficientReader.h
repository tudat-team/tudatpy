/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <map>
#include <iostream>

#include "Tudat/Basics/utilities.h"

#include "Tudat/InputOutput/multiDimensionalArrayReader.h"

namespace tudat
{

namespace input_output
{

//! Function to merge three double multi-arrays of 1 dimension into a single Vector3d multi-array
/*!
 *  Function to merge three double multi-arrays of 1 dimension into a single Vector3d multi-array, where the three
 *  double multi-arrays represent the x-, y- and z-components of the Vector3ds.
 *  \param xComponents Multi-array containing the x-components of the Vector3d
 *  \param yComponents Multi-array containing the y-components of the Vector3d
 *  \param zComponents Multi-array containing the z-components of the Vector3d
 *  \return Single multi-array containing Vector3ds according to double multi-arrays.
 */
boost::multi_array< Eigen::Vector3d, 1 > mergeOneDimensionalCoefficients(
        const boost::multi_array< double, 1 > xComponents,
        const boost::multi_array< double, 1 > yComponents,
        const boost::multi_array< double, 1 > zComponents );

//! Function to merge three double multi-arrays of 2 dimension into a single Vector3d multi-array
/*!
 *  Function to merge three double multi-arrays of 2 dimension into a single Vector3d multi-array, where the three
 *  double multi-arrays represent the x-, y- and z-components of the Vector3ds.
 *  \param xComponents Multi-array containing the x-components of the Vector3d
 *  \param yComponents Multi-array containing the y-components of the Vector3d
 *  \param zComponents Multi-array containing the z-components of the Vector3d
 *  \return Single multi-array containing Vector3ds according to double multi-arrays.
 */
boost::multi_array< Eigen::Vector3d, 2 > mergeTwoDimensionalCoefficients(
        const boost::multi_array< double, 2 > xComponents,
        const boost::multi_array< double, 2 > yComponents,
        const boost::multi_array< double, 2 > zComponents );

//! Function to merge three double multi-arrays of 3 dimension into a single Vector3d multi-array
/*!
 *  Function to merge three double multi-arrays of 3 dimension into a single Vector3d multi-array, where the three
 *  double multi-arrays represent the x-, y- and z-components of the Vector3ds.
 *  \param xComponents Multi-array containing the x-components of the Vector3d
 *  \param yComponents Multi-array containing the y-components of the Vector3d
 *  \param zComponents Multi-array containing the z-components of the Vector3d
 *  \return Single multi-array containing Vector3ds according to double multi-arrays.
 */
boost::multi_array< Eigen::Vector3d, 3 > mergeThreeDimensionalCoefficients(
        const boost::multi_array< double, 3 > xComponents,
        const boost::multi_array< double, 3 > yComponents,
        const boost::multi_array< double, 3 > zComponents );

template< int NumberOfDimensions >
typename boost::multi_array< double ,NumberOfDimensions >::index getIndex(
        const typename boost::multi_array< double, NumberOfDimensions >& m, const double* requestedElement,
        const unsigned short int direction)
{
    int offset = requestedElement - m.origin( );
    return( offset / m.strides( )[direction] % m.shape( )[ direction ] +  m.index_bases( )[direction] );
}

template< int NumberOfDimensions >
boost::array< boost::multi_array< double, 2 >::index, 2 > getMultiArrayIndexArray(
        const typename boost::multi_array< double ,2 >& m, const double* requestedElement )
{
    typename boost::array< boost::multi_array< double, 2 >::index, 2 >  _index;
    for ( unsigned int dir = 0; dir < NumberOfDimensions; dir++ )
    {
        _index[ dir ] = getIndex< NumberOfDimensions >( m, requestedElement, dir );
    }

    return _index;
}

template< int NumberOfDimensions >
boost::multi_array< Eigen::Vector3d, NumberOfDimensions > mergeNDimensionalCoefficients(
        boost::multi_array< double, NumberOfDimensions > xComponents,
        boost::multi_array< double, NumberOfDimensions > yComponents,
        boost::multi_array< double, NumberOfDimensions > zComponents )
{
    boost::multi_array< Eigen::Vector3d, NumberOfDimensions > vectorArray;

    // Check input consistency
    for( unsigned int i = 0; i < NumberOfDimensions; i++ )
    {
        if( !( xComponents.shape( )[ i ] == yComponents.shape( )[ i ] ) ||
                !( xComponents.shape( )[ i ] == zComponents.shape( )[ i ] ) )
        {
            throw std::runtime_error( "Error when creating N-D merged multi-array, input sizes are inconsistent" );
        }
    }

    std::vector< size_t > sizeVector;
    const size_t* arrayShape = xComponents.shape( );
    sizeVector.assign( arrayShape, arrayShape+ xComponents.num_dimensions( ) );

    vectorArray.resize( sizeVector );


    int numberOfEntries = xComponents.num_elements( );
    Eigen::Vector3d* vectorVector = new Eigen::Vector3d[ numberOfEntries ] ;

    typedef double tValue;
    typedef boost::multi_array<tValue,NumberOfDimensions> tArray;
    typedef typename tArray::index tIndex;
    typedef boost::array<tIndex, NumberOfDimensions> tIndexArray;


    tValue* p = xComponents.data();
    tIndexArray index;
    for( unsigned int i = 0; i < numberOfEntries; i++ )
    {
        index = getMultiArrayIndexArray< NumberOfDimensions >( xComponents, p );

        vectorVector[ i ] = ( Eigen::Vector3d( )<<xComponents( index ), yComponents( index ), zComponents( index ) ).finished( );
        ++p;
    }

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

//! Interface class for reading aerodynamic coefficients from files
/*!
 *  Interface class for reading aerodynamic coefficients from files.
 *  This class is used instead of a single templated free function to allow multi-arrays of different sizes to be created
 *  using the same interface. NOTE: The possibility of using a single templated implementation for arbitrary multi-array
 * size should be investigated in the future.
 */
template< int NumberOfDimensions >
class AerodynamicCoefficientReader
{
public:

    //! Function to read a list of aerodynamic coefficients and associated independent variables from a set of files
    /*!
     *  Function to read a list of aerodynamic coefficients and associated independent variables from a set of files.
     *  \param fileNames Vector of size 3, containing the file names for the x-, y- and z- components of the aerodynamic
     *  coefficients. Note that the independent variables for each components must be identical.
     *  \return  Pair: first entry containing multi-array of aerodynamic coefficients, second containing list of independent
     *  variables at which coefficients are defined.
     */
    static std::pair< boost::multi_array< Eigen::Vector3d, NumberOfDimensions >, std::vector< std::vector< double > > >
    readAerodynamicCoefficients( const std::vector< std::string >& fileNames );

    static std::pair< boost::multi_array< Eigen::Vector3d, NumberOfDimensions >, std::vector< std::vector< double > > >
    readAerodynamicCoefficients( const std::map< int, std::string >& fileNames );
};

//! Interface class for reading aerodynamic coefficients of 1 independent variables from files
template< >
class AerodynamicCoefficientReader< 1 >
{
public:

    //! Function to read a list of aerodynamic coefficients and associated independent variables from a set of files
    /*!
     *  Function to read a list of aerodynamic coefficients of 1 independent variables and associated independent variables
     *  from a set of files.
     *  \param fileNames Vector of size 3, containing the file names for the x-, y- and z- components of the aerodynamic
     *  coefficients. Note that the independent variables for each components must be identical.
     *  \return  Pair: first entry containing multi-array of aerodynamic coefficients, second containing list of independent
     *  variables at which coefficients are defined.
     */
    static std::pair< boost::multi_array< Eigen::Vector3d, 1 >, std::vector< std::vector< double > > >
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

        return readAerodynamicCoefficients( fileNameMap );
    }

    static std::pair< boost::multi_array< Eigen::Vector3d, 1 >, std::vector< std::vector< double > > >
    readAerodynamicCoefficients( const std::map< int, std::string >& fileNames )
    {
        // Read data from files.
        std::map< int, boost::multi_array< double, 1 > > rawCoefficientArrays;

        std::vector< boost::multi_array< double, 1 > > coefficientArrays;
        std::vector< std::vector< double > > independentVariables;

        for( std::map< int, std::string >::const_iterator fileIterator = fileNames.begin( ); fileIterator != fileNames.end( );
             fileIterator++ )
        {
            std::pair< boost::multi_array< double, 1 >, std::vector< std::vector< double > > > currentCoefficients =
                    MultiArrayFileReader< 1 >::readMultiArrayAndIndependentVariables( fileIterator->second );
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
                    throw std::runtime_error( "Error when reading 1-Dimensional aeroynamic coefficients, inconsistent aerodynamic coefficients" );
                }
            }

            coefficientArrays[ fileIterator->first ] = currentCoefficients.first;
        }

        if( coefficientArrays.size( ) == 0 )
        {
            throw std::runtime_error( "Error when reading aerodynamic coefficients, no files read" );
        }
        else
        {
            coefficientArrays.resize( 3 );
            boost::multi_array< double, 1 > firstMultiArray = rawCoefficientArrays.begin( )->second;
            for( unsigned int i = 0; i < 3; i++ )
            {
                if( rawCoefficientArrays.count( i ) != 0 )
                {
                    coefficientArrays[ i ] = rawCoefficientArrays.at( i );
                }
                else
                {
                    std::vector< size_t > sizeVector;
                    const size_t* arrayShape = firstMultiArray.shape( );
                    sizeVector.assign( arrayShape, arrayShape+ firstMultiArray.num_dimensions( ) );

                    coefficientArrays[ i ].resize( sizeVector );

                    std::fill( coefficientArrays[ i ].data(), coefficientArrays[ i ].data() + coefficientArrays[ i ].num_elements( ), 0.0 );
                }
            }
        }

        // Merge coefficient entries
        return std::make_pair(
                    mergeOneDimensionalCoefficients(
                        coefficientArrays.at( 0 ), coefficientArrays.at( 1 ), coefficientArrays.at( 2 ) ), independentVariables );
    }
};

//! Interface class for reading aerodynamic coefficients of 2 independent variables from files
template< >
class AerodynamicCoefficientReader< 2 >
{
public:

    //! Function to read a list of aerodynamic coefficients and associated independent variables from a set of files
    /*!
     *  Function to read a list of aerodynamic coefficients of 2 independent variables and associated independent variables
     *  from a set of files.
     *  \param fileNames Vector of size 3, containing the file names for the x-, y- and z- components of the aerodynamic
     *  coefficients. Note that the independent variables for each components must be identical.
     *  \return  Pair: first entry containing multi-array of aerodynamic coefficients, second containing list of independent
     *  variables at which coefficients are defined.
     */
    static std::pair< boost::multi_array< Eigen::Vector3d, 2 >, std::vector< std::vector< double > > >
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

        return readAerodynamicCoefficients( fileNameMap );
    }

    static std::pair< boost::multi_array< Eigen::Vector3d, 2 >, std::vector< std::vector< double > > >
    readAerodynamicCoefficients( const std::map< int, std::string >& fileNames )
    {
        // Read data from files.
        std::map< int, boost::multi_array< double, 2 > > rawCoefficientArrays;

        std::vector< boost::multi_array< double, 2 > > coefficientArrays;
        std::vector< std::vector< double > > independentVariables;

        for( std::map< int, std::string >::const_iterator fileIterator = fileNames.begin( ); fileIterator != fileNames.end( );
             fileIterator++ )
        {
            std::pair< boost::multi_array< double, 2 >, std::vector< std::vector< double > > > currentCoefficients =
                    MultiArrayFileReader< 2 >::readMultiArrayAndIndependentVariables( fileIterator->second );
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
                    throw std::runtime_error( "Error when reading 1-Dimensional aeroynamic coefficients, inconsistent aerodynamic coefficients" );
                }
            }
            utilities::copyMultiArray< double, 2 >(
                        currentCoefficients.first, rawCoefficientArrays[ fileIterator->first ] );

        }

        if( rawCoefficientArrays.size( ) == 0 )
        {
            throw std::runtime_error( "Error when reading aerodynamic coefficients, no files read" );
        }
        else
        {
            coefficientArrays.resize( 3 );
            boost::multi_array< double, 2 > firstMultiArray = rawCoefficientArrays.begin( )->second;

            for( unsigned int i = 0; i < 3; i++ )
            {
                if( rawCoefficientArrays.count( i ) != 0 )
                {
                    utilities::copyMultiArray< double, 2 >(
                                rawCoefficientArrays.at( i ), coefficientArrays[ i ] );
                }
                else
                {
                    std::vector< size_t > sizeVector;
                    const size_t* arrayShape = firstMultiArray.shape( );
                    sizeVector.assign( arrayShape, arrayShape+ firstMultiArray.num_dimensions( ) );

                    coefficientArrays[ i ].resize( sizeVector );

                    std::fill( coefficientArrays[ i ].data(), coefficientArrays[ i ].data() + coefficientArrays[ i ].num_elements( ), 0.0 );
                }
            }
        }

        // Merge coefficient entries
        return std::make_pair(
                    mergeNDimensionalCoefficients< 2 >(
                        coefficientArrays.at( 0 ), coefficientArrays.at( 1 ), coefficientArrays.at( 2 ) ), independentVariables );
    }

};

//! Interface class for reading aerodynamic coefficients of 3 independent variables from files
template< >
class AerodynamicCoefficientReader< 3 >
{
public:

    //! Function to read a list of aerodynamic coefficients and associated independent variables from a set of files
    /*!
     *  Function to read a list of aerodynamic coefficients of 3 independent variables and associated independent variables
     *  from a set of files.
     *  \param fileNames Vector of size 3, containing the file names for the x-, y- and z- components of the aerodynamic
     *  coefficients. Note that the independent variables for each components must be identical.
     *  \return  Pair: first entry containing multi-array of aerodynamic coefficients, second containing list of independent
     *  variables at which coefficients are defined.
     */
    static std::pair< boost::multi_array< Eigen::Vector3d, 3 >, std::vector< std::vector< double > > >
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

        return readAerodynamicCoefficients( fileNameMap );
    }

    static std::pair< boost::multi_array< Eigen::Vector3d, 3 >, std::vector< std::vector< double > > >
    readAerodynamicCoefficients( const std::map< int, std::string >& fileNames )
    {
        // Read data from files.
        std::map< int, boost::multi_array< double, 3 > > rawCoefficientArrays;

        std::vector< boost::multi_array< double, 3 > > coefficientArrays;
        std::vector< std::vector< double > > independentVariables;

        for( std::map< int, std::string >::const_iterator fileIterator = fileNames.begin( ); fileIterator != fileNames.end( );
             fileIterator++ )
        {
            std::pair< boost::multi_array< double, 3 >, std::vector< std::vector< double > > > currentCoefficients =
                    MultiArrayFileReader< 3 >::readMultiArrayAndIndependentVariables( fileIterator->second );
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
                    throw std::runtime_error( "Error when reading 1-Dimensional aeroynamic coefficients, inconsistent aerodynamic coefficients" );
                }
            }

            utilities::copyMultiArray< double, 3 >(
                        currentCoefficients.first, rawCoefficientArrays[ fileIterator->first ] );
        }

        if( rawCoefficientArrays.size( ) == 0 )
        {
            throw std::runtime_error( "Error when reading aerodynamic coefficients, no files read" );
        }
        else
        {
            coefficientArrays.resize( 3 );
            boost::multi_array< double, 3 > firstMultiArray = rawCoefficientArrays.begin( )->second;
            for( unsigned int i = 0; i < 3; i++ )
            {
                if( rawCoefficientArrays.count( i ) != 0 )
                {
                    utilities::copyMultiArray< double, 3 >(
                                rawCoefficientArrays.at( i ), coefficientArrays[ i ] );
                }
                else
                {
                    std::vector< size_t > sizeVector;
                    const size_t* arrayShape = firstMultiArray.shape( );
                    sizeVector.assign( arrayShape, arrayShape+ firstMultiArray.num_dimensions( ) );

                    coefficientArrays[ i ].resize( sizeVector );

                    std::fill( coefficientArrays[ i ].data( ), coefficientArrays[ i ].data( ) + coefficientArrays[ i ].num_elements( ), 0.0 );
                }
            }
        }

        // Merge coefficient entries
        return std::make_pair(
                    mergeThreeDimensionalCoefficients(
                        coefficientArrays.at( 0 ), coefficientArrays.at( 1 ), coefficientArrays.at( 2 ) ), independentVariables );
    }

};

}

}
