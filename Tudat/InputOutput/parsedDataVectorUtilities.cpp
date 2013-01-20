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
 *      120326    D. Dirkx          Code checked, minor layout changes, implementation moved
 *                                  to cpp file and static identifier removed to prevent compile
 *                                  warning.
 *      120424    T. Secretin       Code-check, layout changes.
 *
 *    References
 *
 *    Notes
 *
 */

#include <map>
#include <vector>

#include <boost/make_shared.hpp>

#include "Tudat/InputOutput/parsedDataVectorUtilities.h"

namespace tudat
{
namespace input_output
{
namespace parsed_data_vector_utilities
{

//! Filter the data vector for entries containing a given FieldType.
ParsedDataVectorPtr filterMapKey( ParsedDataVectorPtr datavector, int nrFields, ...)
{
    // Create a fancy vector (list) of all the fields:
    va_list	argumentList;                 // Define argument list variable.
    va_start( argumentList, nrFields );   // Initialize list; point to last defined argument.

    // Create a new datavector for the filtered data.
    ParsedDataVectorPtr newdatavector = boost::make_shared< ParsedDataVector >( );

    // Make a simple list to iterate over from all the FieldType arguments.
    std::vector< FieldType > checkForFieldTypes;
    for ( int i = 0; i < nrFields; i++ )
    {
        checkForFieldTypes.push_back( va_arg( argumentList, FieldType ) );
    }

    // Clean up the system stack.
    va_end( argumentList );

    // Go over every dataline in the current datavector.
    for ( ParsedDataVector::iterator currentDataLine = datavector->begin( );
          currentDataLine != datavector->end( );
          currentDataLine++ )
    {
        // Flag to indicate that all the FieldTypes from checkForFieldTypes are present in this
        // line.
        bool found = true;

        // Loop over each FieldType to check if it exists.
        for( std::vector< FieldType >::iterator currentFieldCheck
             = checkForFieldTypes.begin( );
             currentFieldCheck != checkForFieldTypes.end( );
             currentFieldCheck++ )
        {
            // Check if the FieldType is in the current data line.
            if ( ( *currentDataLine )->find( *currentFieldCheck )
                 == ( *currentDataLine )->end( ) )
            {
                // If not, mark that not all entries are present, and stop searching for others.
                found = false;
                break;
            }
        }

        // If all the fields are present, add the current line to the new (filtered) vector.
        if ( found )
        {
            newdatavector->push_back( *currentDataLine );
        }
    }

    // Return the filtered data vector.
    return newdatavector;
}

//! Filter the data vector vector for entries containing a given FieldType and a matching
//! FieldValue regex.
ParsedDataVectorPtr filterMapKeyValue( ParsedDataVectorPtr datavector, int nrFields, ... )
{
    // Create a fancy vector (list) of all the fields:
    va_list	argumentList;                 // Define argument list variable.
    va_start( argumentList, nrFields );   // Initialize list; point to last defined argument.

    // Create a new data vector for the filtered data.
    ParsedDataVectorPtr newdatavector = boost::make_shared< ParsedDataVector>( );

    // Make a simple list to iterate over with the FieldType arguments and respective regex
    // expressions.
    std::map<FieldType, boost::regex> checkForFieldTypes;
    for ( int i=0; i < nrFields; i++ )
    {
        FieldType type = va_arg( argumentList, FieldType );
        boost::regex regex = boost::regex( va_arg( argumentList, char* ) );
        checkForFieldTypes.insert( std::pair< FieldType, boost::regex >( type, regex ) );
    }

    // Clean up the system stack.
    va_end ( argumentList );

    // Go over every dataline in the current data vector.
    for (ParsedDataVector::iterator currentDataLine = datavector->begin( );
         currentDataLine != datavector->end( );
         currentDataLine++)
    {
        // Flag to indicate that all the FieldTypes from checkForFieldTypes are present in this
        // line.
        bool found = true;

        // Loop over each FieldType to check if it exists.
        for(std::map<FieldType, boost::regex>::iterator currentFieldCheck
            = checkForFieldTypes.begin( );
            currentFieldCheck != checkForFieldTypes.end( );
            currentFieldCheck++)
        {
            // Check if the FieldType is in the current dataline
            ParsedDataLineMap::iterator entry =
                    ( *currentDataLine )->find( currentFieldCheck->first );
            if ( entry == ( *currentDataLine )->end( ) )
            {
                // If not, mark that not all entry pairs are present, and stop searching for
                // others.
                found = false;
                break;
            }

            // Get the field value string.
            boost::shared_ptr<std::string> str = entry->second->get( );

            // Check if the field value regex matches.
            if ( !boost::regex_search( str->begin( ), str->end( ),
                                       currentFieldCheck->second ) )
            {
                // If not, mark this, and stop searching for others.
                found = false;
                break;
            }
        }

        // If all the fields are present, add the current line to the new (filtered) vector.
        if ( found )
        {
            newdatavector->push_back( *currentDataLine );
        }
    }

    // Return the filtered data vector.
    return newdatavector;
}

//! Dump the content of a data map to an ostream.
std::ostream& dump( std::ostream& stream, ParsedDataLineMapPtr data, bool showTransformed )
{
    // Initial character to separate elements.
    stream << "|";

    // Loop over data map.
    for( ParsedDataLineMap::iterator element = data->begin( );
         element != data->end( ); element++ )
    {
        // If transformed flag is true, dump transformed values.
        if ( showTransformed )
        {
            stream << element->second->get( )->c_str( );
        }

        // If not, dump raw values.
        else
        {
            stream << element->second->getRaw( )->c_str( );
        }

        // Final character to separate elements.
        stream << "\t|";
    }

    stream << "\n";
    return stream;
}

//! Dump the content of a data vector to an ostream (eg std::cout).
std::ostream& dump( std::ostream& stream, ParsedDataVectorPtr data, bool showTransformed )
{
    // Loop over vector.
    for( std::size_t i=0; i < data->size( ); i++ )
    {
        // Get data line.
        ParsedDataLineMapPtr line = data->at( i );

        // Call dump function for single line (datamap).
        dump ( stream, line, showTransformed );
    }

    return stream;
}

} // namespace parsed_data_vector_utilities
} // namespace input_output
} // namespace tudat
