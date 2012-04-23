/*    Copyright (c) 2010-2011 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      111206    S. Billemont      First creation of code.
 *      120326    D. Dirkx          Code checked, minor layout changes.
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
