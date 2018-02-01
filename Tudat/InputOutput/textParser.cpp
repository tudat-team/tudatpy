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

#include "Tudat/InputOutput/textParser.h"

namespace tudat
{
namespace input_output
{

// Using declaration.
using namespace parsed_data_vector_utilities;

// Parse string.
ParsedDataVectorPtr TextParser::parse( std::string& string )
{
    // Create stringstream object.
    std::stringstream stringStream( string );

    // Return parsed data from stringstream.
    return parse( stringStream );
}

// Parse stream.
ParsedDataVectorPtr TextParser::parse( std::istream& stream )
{
    // Clear parsedData variable.
    parsedData->clear( );

    if ( parseAsStream )
    {
        // Stream based parsing.
        parseStream( stream );
    }
    else
    {
        // Line based parsing
        while ( !stream.fail( ) && !stream.eof( ) )
        {
            std::string line;
            std::getline( stream, line );
            parseLine( line );
        }
    }

    // Return parsed data.
    return parsedData;
}

} // namespace input_output
} // namespace tudat
