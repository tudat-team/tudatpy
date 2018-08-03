/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TABULATED_ATMOSPHERE_READER_H
#define TUDAT_TABULATED_ATMOSPHERE_READER_H

#include <map>
#include "Tudat/Basics/utilities.h"
#include "Tudat/Basics/basicTypedefs.h"

#include "Tudat/InputOutput/multiDimensionalArrayReader.h"

namespace tudat
{

namespace input_output
{

//! Function to compare if two lists of aerodynamic coefficient independent variables are equal
/*!
 * Function to compare if two lists of aerodynamic coefficient independent variables (vector of vector of doubles) are equal
 * \param list1 First list that is to be compared.
 * \param list2 Second list that is to be compared.
 * \return True of the two lists are completely equal in size and contents, false otherwise.
 */
bool compareIndependentVariables( const std::vector< std::vector< double > >& list1,
                                  const std::vector< std::vector< double > >& list2 );

//! Function to read a list of atmosphere parameters and associated independent variables from a set of files
/*!
 *  Function to read a list of atmosphere parameters of 2 independent variables and associated independent variables
 *  from a set of files.
 *  \param fileNames Map of file names, with the key  required to be 0, 1 and/or 2. These indices denote the  x-, y- and z-
 *  components of the atmosphere parameters. All indices that are not provided in this map are assumed to have associated
 *  coefficients equal to zero for all values of the independent variables
 *  Note that the independent variables for each components must be identical.
 *  \return  Pair: first entry containing multi-array of atmosphere parameters, second containing list of independent
 *  variables at which coefficients are defined.
 */
template< unsigned int NumberOfDimensions >
std::pair< std::vector< boost::multi_array< double, static_cast< size_t >( NumberOfDimensions ) > >,
std::vector< std::vector< double > > >
readTabulatedAtmosphere( const std::map< int, std::string >& fileNames )
{
    // Assing number of dependent variables
    unsigned int numberOfFiles = fileNames.size( );

    // Initialize output
    std::map< int, boost::multi_array< double, static_cast< size_t >( NumberOfDimensions ) > > rawAtmosphereArrays;
    std::vector< boost::multi_array< double, static_cast< size_t >( NumberOfDimensions ) > > atmosphereArrays;
    std::vector< std::vector< double > > independentVariables;

    // Iterate over files and read the contents into rawAtmosphereArrays/independentVariables.
    for( std::map< int, std::string >::const_iterator fileIterator = fileNames.begin( ); fileIterator != fileNames.end( );
         fileIterator++ )
    {
        // Read current coefficients/independent variables
        std::pair< boost::multi_array< double, static_cast< size_t >( NumberOfDimensions ) >,
                std::vector< std::vector< double > > > currentCoefficients =
                MultiArrayFileReader< NumberOfDimensions >::readMultiArrayAndIndependentVariables( fileIterator->second );

        // Save/check consistency of independent variables
        if( rawAtmosphereArrays.size( ) == 0 )
        {
            independentVariables = currentCoefficients.second;
        }
        else
        {
            bool areIndependentVariablesEqual = compareIndependentVariables(
                        independentVariables, currentCoefficients.second );

            if( !areIndependentVariablesEqual )
            {
                throw std::runtime_error( "Error when reading 1-Dimensional atmosphere parameters, inconsistent independent variables." );
            }
        }

        // Save file contents into rawAtmosphereArrays
        utilities::copyMultiArray< double, NumberOfDimensions >(
                    currentCoefficients.first, rawAtmosphereArrays[ fileIterator->first ] );
    }

    // Check if anything has been read from files.
    if( rawAtmosphereArrays.size( ) == 0 )
    {
        throw std::runtime_error( "Error when reading atmosphere parameters, no files read." );
    }
    else
    {
        atmosphereArrays.resize( numberOfFiles );
        boost::multi_array< double, static_cast< size_t >( NumberOfDimensions ) > firstMultiArray =
                rawAtmosphereArrays.begin( )->second;

        // Iterate over all entries
        for( unsigned int i = 0; i < numberOfFiles; i++ )
        {
            // Copy contents for current index into atmosphereArrays read from file.
            if( rawAtmosphereArrays.count( i ) != 0 )
            {
                utilities::copyMultiArray< double, NumberOfDimensions >(
                            rawAtmosphereArrays.at( i ), atmosphereArrays[ i ] );
            }
            // Set zero multi-array for current index.
            else
            {
                std::vector< size_t > sizeVector;
                const size_t* arrayShape = firstMultiArray.shape( );
                sizeVector.assign( arrayShape, arrayShape + firstMultiArray.num_dimensions( ) );

                atmosphereArrays[ i ].resize( sizeVector );

                std::fill( atmosphereArrays[ i ].data( ),
                           atmosphereArrays[ i ].data() + atmosphereArrays[ i ].num_elements( ), 0.0 );
            }
        }
    }

    // Merge coefficient entries
    return std::make_pair( atmosphereArrays, independentVariables );
}

} // namespace input_output

} // namespace tudat

#endif // TUDAT_TABULATED_ATMOSPHERE_READER_H
