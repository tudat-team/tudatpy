/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      111103    S. Billemont      Creation of code.
 *      120217    D.J. Gondelach    Code check.
 *
 *    References
 *
 */

#include <map>
#include <vector>
#include <utility>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <cmath>
#include <iostream>
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"
#include "Tudat/InputOutput/whiteSpaceParser.h"

namespace tudat
{
namespace input_output
{

//! Create a parser that parses based on white spaces given a field type list.
WhiteSpaceParser::WhiteSpaceParser( int numberOfFields, ... ) : TextParser ( false )
{
    // Copy number of fields.
    numberOfFields_ = numberOfFields;

    // Create a fancy vector (list) of all the fields:
    // Define argument list variable.
    va_list listOfArguments;

    // Initialize list. Point to last defined argument.
    va_start( listOfArguments, numberOfFields );

    for ( unsigned int i=0; i < numberOfFields_; i++ )
    {
        // Populate typeList with arguments from constructor.
        typeList.push_back( va_arg( listOfArguments, FieldType ) );
    }

    // Clean up the system stack.
    va_end( listOfArguments );
}

//! Parse line of data.
void WhiteSpaceParser::parseLine( std::string& line )
{
    // Using declaration.
    using namespace tudat::input_output::parsed_data_vector_utilities;

    // Define a new type: pair of field type and pointer to value.
    typedef std::pair< FieldType, FieldValuePtr > FieldDataPair;

    // Define a new type: vector of strings.
    typedef std::vector< std::string > vectorOfIndividualStrings;

    // Create a vector of individual strings.
    vectorOfIndividualStrings vectorOfIndividualStrings_;

    // Trim input string (removes all leading and trailing whitespaces).
    boost::algorithm::trim( line );

    // Split string into multiple strings, each containing one element from a line from the
    // data file.
    boost::algorithm::split( vectorOfIndividualStrings_,
                             line,
                             boost::algorithm::is_any_of( " " ),
                             boost::algorithm::token_compress_on );

    // Verify that number of individual vectors corresponds to the specified number of fields.
    if ( vectorOfIndividualStrings_.size( ) != numberOfFields_ )
    {
        std::cerr << "Number of elements in the line (" << vectorOfIndividualStrings_.size( )
        << ") does not match the specified number of fields (" << numberOfFields_ << ")"
        << std::endl;
    }

    // Create a new pointer to map of line data.
    ParsedDataLineMapPtr currentLineData = boost::make_shared< ParsedDataLineMap >(
                std::map< FieldType, FieldValuePtr >( ) );

    // Register the data line with the global current parsed data vector.
    parsedData->push_back( currentLineData );

    //Loop over all field type and field value pairs.
    for ( int unsigned currentFieldNumber = 0;
          currentFieldNumber < numberOfFields_;
          currentFieldNumber++ )
    {
        // Get the corresponding field type.
        FieldType fieldType( typeList.at( currentFieldNumber ) );

        // Define unit transformer.
        boost::shared_ptr< FieldTransform > transformer;

        // If type corresponds to one of the entries of the unit transformation map.
        if ( unitTransformationMap_.find( fieldType ) != unitTransformationMap_.end( ) )
        {
            // Set corresponding transformer.
            transformer = unitTransformationMap_.find( fieldType )->second;
        }

        // Else, do nothing.
        else
        {
            transformer = boost::shared_ptr< FieldTransform >( );
        }

        // Store the resulting field-value string.
        FieldValuePtr value( new FieldValue( fieldType,
                                             vectorOfIndividualStrings_.at( currentFieldNumber ),
                                             transformer ) );

        // Store the type and value in the current line data.
        currentLineData->insert( FieldDataPair( fieldType, value ) );
    }
}

} // namespace input_output
} // namespace tudat
