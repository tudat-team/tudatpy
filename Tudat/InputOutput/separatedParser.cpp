/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      120706    D.J. Gondelach    Code check.
 *      120718    A. Ronse          Replaced custom split functions by
 *                                  boost::algorithm::make_split_iterator method.
 *
 *    References
 *
 *    Notes
 *
 */

#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/range/algorithm/copy.hpp>

#include "Tudat/InputOutput/separatedParser.h"

namespace tudat
{
namespace input_output
{

//! Create a parser that parses based on a separator and given field type list.
SeparatedParser::SeparatedParser( std::string separator, int numberOfFields, ... )
    : TextParser ( false ), doTrim( true )
{
    // Copy number of fields.
    numberOfFields_ = numberOfFields;

    separator_ = separator;

    // Create a fancy vector (list) of all the fields.
    va_list listOfArguments;	        // Define argument list variable.
    va_start( listOfArguments,numberOfFields );  // init list; point to last defined argument.

    for ( int i = 0; i < numberOfFields; i++ )
    {
        typeList.push_back( va_arg( listOfArguments, FieldType ) ); // get next argument.
    }

    va_end ( listOfArguments );  // clean up the system stack.
}

//! Parses one line of text.
void SeparatedParser::parseLine( std::string& line )
{
    // Some short hand notations.
    using namespace std;
    using namespace parsed_data_vector_utilities;
    typedef std::pair< FieldType, FieldValuePtr > FieldDataPair;

    // Define a new type: vector of strings.
    typedef std::vector< std::string > vectorOfIndividualStrings;

    // Create a vector of individual strings.
    vectorOfIndividualStrings vectorOfIndividualStrings_;

    // Split string into multiple strings based on provided separator and place in vector.
    typedef boost::algorithm::split_iterator< std::string::iterator > string_split_iterator;
    for( string_split_iterator It=
         boost::algorithm::make_split_iterator( line, boost::algorithm::first_finder(
                                                    separator_, boost::algorithm::is_iequal( ) ) );
         It!=string_split_iterator( );
         ++It )
    {
        // Prevent empty rows due to a double separator.
        if ( !boost::copy_range< std::string >( *It ).empty( ) )
        {
            vectorOfIndividualStrings_.push_back( boost::copy_range< std::string >( *It ) );
        }
    }

    // Verify that number of individual vectors corresponds to the specified number of fields.
    if ( vectorOfIndividualStrings_.size( ) != numberOfFields_ )
    {
        std::cerr << "Number of elements in the line (" << vectorOfIndividualStrings_.size( )
        << ") does not match the specified number of fields (" << numberOfFields_ << ")"
        << std::endl;
    }

    // Create a new data line
    ParsedDataLineMapPtr currentLineData = boost::make_shared< ParsedDataLineMap >(
                std::map< FieldType, FieldValuePtr >( ) );

    // Register the data line with the global current parsed data vector
    parsedData->push_back( currentLineData );

    // Loop over all fields, until the end position pointer is located at the end of the string
    for ( int unsigned currentFieldNumber = 0; currentFieldNumber < numberOfFields_;
          currentFieldNumber++ )
    {
        // If we need to trim whitespace, do so.
        if ( doTrim )
        {
                boost::trim( vectorOfIndividualStrings_.at( currentFieldNumber ) );
        }

        // Get the corresponding field type.
        FieldType type( typeList.at( currentFieldNumber ) );

        // Define unit transformer.
        boost::shared_ptr< FieldTransform > transformer;

        // If type corresponds to one of the entries of the unit transformation map.
        if ( unitTransformationMap_.find( type ) != unitTransformationMap_.end( ) )
        {
            // Set corresponding transformer.
            transformer = unitTransformationMap_.find( type )->second;
        }
        // Else, do nothing.
        else
        {
            transformer = boost::shared_ptr< FieldTransform >( );
        }

        // Store the resulting string.
        FieldValuePtr value(
                    new FieldValue(type, vectorOfIndividualStrings_.at( currentFieldNumber ),
                                   transformer ) );

        // Store the type and value in the current line data
        currentLineData->insert( FieldDataPair( type, value ) );
    }
}

} // namespace input_output
} // namespace tudat
