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

#include "Tudat/InputOutput/dictionaryComparer.h"

namespace tudat
{
namespace input_output
{
namespace dictionary
{

//! Overload ()-operator to compare data line with dictionary entry.
bool DictionaryComparer::operator( )(
        const input_output::parsed_data_vector_utilities::ParsedDataLineMapPtr& dataLine ) const
{
    // Get parameter name from data line map.
    std::string parameterName
            = input_output::parsed_data_vector_utilities::getField< std::string >(
                dataLine, input_output::field_types::general::parameterName );

    // Store case-sensitive flag locally.
    bool isCaseSensitive = dictionaryEntry->isCaseSensitive;

    // Perform comparison check, taking state of isCaseSensitive in account.
    return ( !dictionaryEntry->parameterName.compare( parameterName )
             || !( dictionaryEntry->synonyms.find( parameterName )
                   == dictionaryEntry->synonyms.end( ) ) )
            || ( !isCaseSensitive &&
                 ( boost::iequals( dictionaryEntry->parameterName, parameterName )
                   || !( std::find_if( dictionaryEntry->synonyms.begin( ),
                                       dictionaryEntry->synonyms.end( ),
                                       DictionaryComparer( parameterName ) )
                         == dictionaryEntry->synonyms.end( ) ) ) );
}

} // namespace dictionary
} // namespace input_output
} // namespace tudat
