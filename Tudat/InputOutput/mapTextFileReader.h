/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      If tabs are used as spaces, it doesn't work. The seperator should also be tabs then.
 *
 */

#ifndef TUDAT_MAP_TEXT_FILEREADER_H
#define TUDAT_MAP_TEXT_FILEREADER_H

#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/lexical_cast.hpp>

#include <Eigen/Core>

#include "Tudat/InputOutput/matrixTextFileReader.h"
#include "Tudat/InputOutput/streamFilters.h"


namespace tudat
{

namespace input_output
{

//! Read a map of Eigen vectors from a text file.
/*!
 * Read a textfile whith separated (space, tab, comma etc...) values.
 * The first element of each line is a key.
 * Successive elements in the current line are the associated values, converted to Eigen column vector.
 * \remarks The scalar type of `EigenVectorType` is used for the keys of the returned map.
 * \param relativePath Relative path to file.
 * \param separators Separators used, every character in the string will be used as separators.
 *         (multiple seperators possible).
 * \param skipLinesCharacter Skip lines starting with this character.
 * \return The map with data as read from the file.
 */
template< typename EigenVectorType = Eigen::VectorXd >
std::map< typename EigenVectorType::Scalar, EigenVectorType > readEigenVectorMapFromFile(
        const std::string& relativePath,
        const std::string& separators = "\t ;,",
        const std::string& skipLinesCharacter = "%" )
{
    // Load data into matrix
    Eigen::Matrix< typename EigenVectorType::Scalar, Eigen::Dynamic, Eigen::Dynamic > matrix =
            readMatrixFromFile< typename EigenVectorType::Scalar >( relativePath, separators, skipLinesCharacter );

    // Fill map using Eigen matrix
    std::map< typename EigenVectorType::Scalar, EigenVectorType > map;
    for ( int i = 0; i < matrix.rows( ); i++ )
    {
        const Eigen::Matrix< typename EigenVectorType::Scalar, Eigen::Dynamic, 1 > dynamicSizeColumnVector =
                matrix.row( i ).rightCols( matrix.cols( ) - 1 ).transpose( );
        EigenVectorType columnVector;
        if ( columnVector.rows( ) == 0 )
        {
            columnVector = EigenVectorType( dynamicSizeColumnVector.rows( ) );
        }
        else if ( dynamicSizeColumnVector.rows( ) != columnVector.rows( ) )
        {
            std::cerr << "The size of the vectors in \"" << relativePath << "\" (" << dynamicSizeColumnVector.rows( )
                      << ") is not consistent with the requested size (" << columnVector.rows( ) << ")." << std::endl;
            throw;
        }
        for ( int j = 0; j < columnVector.rows( ); j++ )
        {
            columnVector( j ) = dynamicSizeColumnVector( j );
        }
        map[ matrix( i, 0 ) ] = columnVector;
    }
    return map;
}


//! Read a map of std vector from a text file.
/*!
 * Read a textfile whith separated (space, tab, comma etc...) values.
 * The first element of each line is a key.
 * Successive elements in the current line are the associated values, stored as an STL vector.
 * \param relativePath Relative path to file.
 * \param separators Separators used, every character in the string will be used as separators.
 *         (multiple seperators possible).
 * \param skipLinesCharacter Skip lines starting with this character.
 * \return The map with data as read from the file.
 */
template< typename KeyType, typename ScalarValueType >
std::map< KeyType, std::vector< ScalarValueType > > readStlVectorMapFromFile(
        const std::string& relativePath,  const std::string& separators = "\t ;,",
        const std::string& skipLinesCharacter = "%" )
{
    // Open input and output.
    std::fstream file( relativePath.c_str( ), std::ios::in );
    if ( file.fail( ) )
    {
        throw std::runtime_error( "Data file could not be opened: " + relativePath );
    }

    std::stringstream filteredStream( std::ios::in | std::ios::out );
    {
        // Filter the file stream. This needs to be in its own scope, because filtering_stream::
        // flush( ) does not work if the underlying end point is a stringstream, so the flush has
        // to be forced by letting the filtering_stream go out of scope.
        boost::iostreams::filtering_ostream filterProcessor;
        for ( unsigned int i = 0; i < skipLinesCharacter.size( ); i++ )
        {
            // Remove all comments from the stream.
            filterProcessor.push( input_output::stream_filters::RemoveComment(
                                      skipLinesCharacter[ i ], true ) );
        }

        // Add the output to the filter chain.
        filterProcessor.push( filteredStream );

        // Copy the input to the filter.
        boost::iostreams::copy( file, filterProcessor );
    }

    // Seek stream back to start.
    filteredStream.seekg( 0, std::ios::beg );

    // Read the filtered stream into lines.
    std::vector< std::string > lines_;
    while ( !filteredStream.eof( ) )
    {
        std::string line_;
        getline( filteredStream, line_ );
        if ( !line_.empty( ) )
        {
            boost::trim_all( line_ );
            lines_.push_back( line_ );
        }
    }

    // Initialize the map.
    std::map< KeyType, std::vector< ScalarValueType > > dataMap_;

    // If there are no lines, return an empty matrix.
    if ( lines_.empty( ) )
    {
        return dataMap_;
    }

    const std::string realSeparators = std::string( separators ) + " ";

    for ( unsigned int rowIndex = 0; rowIndex < lines_.size( ); rowIndex++ )
    {
        // Determine the number of columns from.
        std::vector< std::string > lineSplit_;
        boost::algorithm::split( lineSplit_, lines_[ rowIndex ], boost::is_any_of( realSeparators ),
                                 boost::algorithm::token_compress_on );

        // Determine key and put single line entries into vector.
        KeyType key;
        std::vector< ScalarValueType > values;
        for ( unsigned int columnIndex = 0; columnIndex < lineSplit_.size( ); columnIndex++ )
        {
            boost::trim( lineSplit_.at( columnIndex ) );
            if ( columnIndex == 0 )
            {
                key = boost::lexical_cast< KeyType >( lineSplit_.at( columnIndex ) );
            }
            else
            {
                values.push_back( boost::lexical_cast< ScalarValueType >( lineSplit_.at( columnIndex ) ) );
            }
        }

        dataMap_[ key ] = values;
    }

    return dataMap_;
}

} // namespace input_output

} // namespace tudat

#endif // TUDAT_MAP_TEXT_FILEREADER_H
