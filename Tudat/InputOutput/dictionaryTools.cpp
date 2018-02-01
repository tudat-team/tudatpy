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

#include <sstream>

#include <boost/make_shared.hpp>

#include "Tudat/InputOutput/dictionaryTools.h"

namespace tudat
{
namespace input_output
{
namespace dictionary
{

//! Check required parameters.
void checkRequiredParameters( const DictionaryPointer& aDictionary )
{
    // Create new list of missing, required parameters.
    RequiredParametersList missingRequiredParameters;

    // Loop through all dictionary entries to check if any required parameters have not been
    // extracted. Add any missing, required parameters to the list.
    for ( Dictionary::const_iterator iteratorDictionary = aDictionary->begin( );
          iteratorDictionary != aDictionary->end( ); iteratorDictionary++ )
    {
        if ( !( *iteratorDictionary )->isExtracted && ( *iteratorDictionary )->isRequired )
        {
            missingRequiredParameters.insert( *iteratorDictionary );
        }
    }

    // Check if list is non-empty. If it is non-empty, throw an error, indicating which required
    // parameters are missing.
    if ( missingRequiredParameters.size( ) > 0 )
    {
        // Create error message.
        std::stringstream missingRequiredParametersMessage;
        missingRequiredParametersMessage << "The following required parameters are missing: ";

        // Set iterator to last element in list.
        RequiredParametersList::const_iterator lastIterator = missingRequiredParameters.end( );
        lastIterator--;

        // Loop through the list and add the name of the missing parameter to the error message.
        for ( RequiredParametersList::const_iterator missingParametersIterator
              = missingRequiredParameters.begin( );
              missingParametersIterator != lastIterator; missingParametersIterator++ )
        {
            missingRequiredParametersMessage
                    << "\"" << ( *missingParametersIterator )->parameterName << "\", ";
        }

        missingRequiredParametersMessage << "\"" << ( *lastIterator )->parameterName << "\".";

        // Thrown exception based on constructed error message.
        throw std::runtime_error( missingRequiredParametersMessage.str( ) );
    }
}

//! Add entry.
void addEntry( const DictionaryPointer& dictionary, const std::string& parameterName,
               const bool isRequired, const bool isCaseSensitive, const StringSet& someSynonyms )
{
    dictionary->insert( boost::make_shared< DictionaryEntry >(
                            parameterName, isRequired, isCaseSensitive, someSynonyms ) );
}

//! Find entry.
DictionaryIterator findEntry( const DictionaryPointer dictionary,
                              const std::string& parameterName )
{
    DictionaryIterator iteratorDictionary
            = std::find_if( dictionary->begin( ), dictionary->end( ),
                            DictionaryComparer( parameterName ) );

    std::stringstream errorMessage;
    errorMessage << "Dictionary entry \"" << parameterName << "\" not found!" << std::endl;

    if ( iteratorDictionary == dictionary->end( ) )
    {
        throw std::runtime_error( errorMessage.str( ) );
    }

    return iteratorDictionary;
}

} // namespace dictionary
} // namespace input_output
} // namespace tudat
