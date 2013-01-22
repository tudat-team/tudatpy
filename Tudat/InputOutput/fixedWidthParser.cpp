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
 *      111209    D.J. Gondelach    File created.
 *      120718    A. Ronse          Code check. Implemented optional trim.
 *
 *    References
 *
 *    Notes
 *
 */

#include "Tudat/InputOutput/fixedWidthParser.h"

#include <boost/algorithm/string.hpp>

namespace tudat
{
namespace input_output
{

//! Default constructor.
FixedWidthParser::FixedWidthParser( int numberOfFields, ... ) : TextParser( false ), doTrim( true )
{
    numberOfFields_ = numberOfFields;

    // Create a fancy vector (list) of all the fields:
    // Define argument list variable.
    va_list	listOfArguments;

    // Initialize list. Point to last defined argument.
    va_start( listOfArguments, numberOfFields );

    for ( unsigned int i=0; i < numberOfFields_; i++ )
    {
        // Populate typeList with arguments from constructor.
        typeList.push_back( va_arg(listOfArguments, FieldType ) );
    }

    for ( unsigned int i=0; i < numberOfFields_; i++ )
    {
        // Populate sizeList with arguments from constructor.
        sizeList.push_back( va_arg( listOfArguments, int ) );
    }

    // Clean up the system stack.
    va_end( listOfArguments );
}

//! Parses one line of text.
void FixedWidthParser::parseLine(std::string& line)
{
    using namespace parsed_data_vector_utilities;

    // Define a new type: field type and pointer to value pair
    typedef std::pair< FieldType, FieldValuePtr > FieldDataPair;

    // Define a new type: vector of strings.
    typedef std::vector< std::string > vectorOfIndividualStrings;

    // Create a vector of individual strings.
    vectorOfIndividualStrings vectorOfIndividualStrings_;

    // Size of vector of individual strings equals number of fields.
    vectorOfIndividualStrings_.resize( numberOfFields_ );

    // Temporary string.
    std::string temp;

    // Start index of current field in the line.
    int currentFieldIndex = 0;

    for ( int unsigned currentFieldNumber = 0; currentFieldNumber < numberOfFields_;
          currentFieldNumber++ )
    {
        // Generate string of the field.
        temp = line.substr( currentFieldIndex, sizeList[currentFieldNumber] );

        // If we need to trim whitespace, do so
        if ( doTrim )
        {
            boost::trim( temp );
        }

        // Copy field value.
        vectorOfIndividualStrings_[currentFieldNumber] = temp;

        // Increase current field index to index of next field.
        currentFieldIndex += sizeList[currentFieldNumber];

        // Clear temporary string.
        temp.clear( );
    }

    // Verify that number of individual vectors corresponds to the specified number of fields.
    if ( vectorOfIndividualStrings_.size() != numberOfFields_ )
    {
        std::cerr << "Number of elements in the line (" << vectorOfIndividualStrings_.size( )
        << ") does not match the specified number of fields (" << numberOfFields_ << ")"
        << std::endl;
    }

    // Create a new data line.
    ParsedDataLineMapPtr currentLineData( new std::map< FieldType, FieldValuePtr >( ) );

    // Register the data line with the global current parsed data vector.
    parsedData->push_back( currentLineData );

    //Loop over all field type and field value pairs
    for ( int unsigned currentFieldNumber = 0; currentFieldNumber < numberOfFields_;
          currentFieldNumber++ )
    {
        // Get the corresponding field type
        FieldType type ( typeList.at( currentFieldNumber ) );

        // Define unit transformer
        boost::shared_ptr< FieldTransform > transformer;

        // If type corresponds to one of the entries of the unit transformation map
        if ( unitTransformationMap_.find( type ) != unitTransformationMap_.end( ) )
        {
            // Set corresponding transformer
            transformer = unitTransformationMap_.find( type )->second;
        }

        else
        {
            // Else, do nothing.
            transformer = boost::shared_ptr< FieldTransform >( );
        }

        // Store the resulting string.
        FieldValuePtr value( new FieldValue( type,
                                             vectorOfIndividualStrings_.at( currentFieldNumber ),
                                             transformer ) );

        // Store the type and value in the current line data.
        currentLineData->insert( FieldDataPair( type, value ) );
    }
}

} // namespace input_output
} // namespace tudat
