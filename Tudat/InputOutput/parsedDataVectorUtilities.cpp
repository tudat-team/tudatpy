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
    // Define argument list variable.
    va_list argumentList;
    // Initialize list; point to last defined argument.
    va_start( argumentList, nrFields );   

    // Create a new datavector for the filtered data.
    ParsedDataVectorPtr newDataVector = std::make_shared< ParsedDataVector >( );

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
        for ( std::vector< FieldType >::iterator currentFieldCheck
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
            newDataVector->push_back( *currentDataLine );
        }
    }

    // Return the filtered data vector.
    return newDataVector;
}

//! Filter the data vector vector for entries containing a given FieldType and a matching
//! FieldValue regex.
ParsedDataVectorPtr filterMapKeyValue( ParsedDataVectorPtr datavector, int nrFields, ... )
{
    // Create a fancy vector (list) of all the fields:
    // Define argument list variable.
    va_list argumentList;
    // Initialize list; point to last defined argument.
    va_start( argumentList, nrFields );

    // Create a new data vector for the filtered data.
    ParsedDataVectorPtr newDataVector = std::make_shared< ParsedDataVector>( );

    // Make a simple list to iterate over with the FieldType arguments and respective regex
    // expressions.
    std::map< FieldType, boost::regex > checkForFieldTypes;
    for ( int i = 0; i < nrFields; i++ )
    {
        FieldType type = va_arg( argumentList, FieldType );
        boost::regex regex = boost::regex( va_arg( argumentList, char* ) );
        checkForFieldTypes.insert( std::pair< FieldType, boost::regex >( type, regex ) );
    }

    // Clean up the system stack.
    va_end ( argumentList );

    // Go over every dataline in the current data vector.
    for ( ParsedDataVector::iterator currentDataLine = datavector->begin( );
          currentDataLine != datavector->end( );
          currentDataLine++ )
    {
        // Flag to indicate that all the FieldTypes from checkForFieldTypes are present in this
        // line.
        bool found = true;

        // Loop over each FieldType to check if it exists.
        for( std::map< FieldType, boost::regex >::iterator currentFieldCheck
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
            const std::string str = entry->second->getTransformed( );

            // Check if the field value regex matches.
            if ( !boost::regex_search( str.begin( ), str.end( ), currentFieldCheck->second ) )
            {
                // If not, mark this, and stop searching for others.
                found = false;
                break;
            }
        }

        // If all the fields are present, add the current line to the new (filtered) vector.
        if ( found )
        {
            newDataVector->push_back( *currentDataLine );
        }
    }

    // Return the filtered data vector.
    return newDataVector;
}

//! Dump the content of a data map to an ostream.
std::ostream& dump( std::ostream& stream, ParsedDataLineMapPtr data, bool showTransformed )
{
    // Initial character to separate elements.
    stream << "|";

    // Loop over data map.
    for ( ParsedDataLineMap::iterator element = data->begin( );
          element != data->end( ); element++ )
    {
        // If transformed flag is true, dump transformed values.
        if ( showTransformed )
        {
            stream << element->second->getTransformed( );
        }

        // If not, dump raw values.
        else
        {
            stream << element->second->getRaw( );
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
    for ( std::size_t i=0; i < data->size( ); i++ )
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
