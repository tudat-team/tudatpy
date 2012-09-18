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
 *      120806    K. Kumar          File created.
 *      120815    K. Kumar          Added addEntry() function.
 *      120910    K. Kumar          Added findEntry() function.
 *
 *    References
 *
 *    Notes
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
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( missingRequiredParametersMessage.str( ) ) ) );
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
        boost::throw_exception( boost::enable_error_info(
                                    std::runtime_error( errorMessage.str( ) ) ) );
    }

    return iteratorDictionary;
}

} // namespace dictionary
} // namespace input_output
} // namespace tudat
