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

#ifndef TUDAT_DICTIONARY_TOOLS_H
#define TUDAT_DICTIONARY_TOOLS_H

#include <algorithm>
#include <set>
#include <stdexcept>

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/InputOutput/dictionaryComparer.h"
#include "Tudat/InputOutput/dictionaryEntry.h"

namespace tudat
{
namespace input_output
{
namespace dictionary
{

//! Typedef for dictionary containing dictionary entries.
typedef std::set< DictionaryEntryPointer, DictionaryComparer > Dictionary;

//! Typedef iterator pointing to elements of Dictionary.
typedef Dictionary::const_iterator DictionaryIterator;

//! Typedef for pointer to dictionary containing dictionary entries.
typedef boost::shared_ptr< Dictionary > DictionaryPointer;

//! Typedef for list of required parameters.
typedef std::set< DictionaryEntryPointer, DictionaryComparer > RequiredParametersList;

//! Typedef iterator pointing to elements of ParsedDataVector.
typedef parsed_data_vector_utilities::ParsedDataVector::const_iterator DataLineIterator;

//! Typedef for set of strings.
typedef std::set< std::string > StringSet;

//! Check required parameters.
/*!
 * Checks that all required parameters have been defined. This function checks that all required
 * parameters have been extracted by looping through the given dictionary and checking the
 * isExtracted flag for each entry. If any required parameters are missing, an error is thrown,
 * listing all the missing parameters.
 * \param aDictionary Dictionary to search through.
 */
void checkRequiredParameters( const DictionaryPointer& aDictionary );

//! Add entry.
/*!
 * Adds a new dictionary entry in a given dictionary as a shared-pointer.
 * \param dictionary Parameter dictionary.
 * \param parameterName Parameter name.
 * \param isRequired Flag to indicate if parameter is required.
 * \param isCaseSensitive Flag to indicate if parameter name (and any synonyms) are case-sensitive.
 * \param someSynonyms Set of synonyms.
 */
void addEntry( const DictionaryPointer& dictionary, const std::string& parameterName,
               const bool isRequired, const bool isCaseSensitive,
               const StringSet& someSynonyms = StringSet( ) );

//! Find entry.
/*!
 * Searches through the specified dictionary for a given parameter name and tries to find a match.
 * The string matching is EXACT (case-sensitive). If the parameter name cannot be found in the
 * dictionary, a runtime error is thrown indicating this.
 * \param dictionary Parameter dictionary.
 * \param parameterName Parameter name as string.
 * \return Iterator to matched dictionary entry (this will not be returned if the runtime error is
 *          triggered, in case the entry cannot be found).
 */
DictionaryIterator findEntry( const  DictionaryPointer dictionary,
                              const std::string& parameterName );

//! Convert dummy value.
/*!
 * Converts dummy value. This function simply returns the input value.
 * \tparam DataType Data type of input and output values.
 * \param value Input value.
 * \return Unaffected input value.
 */
template< typename DataType > DataType convertDummy( DataType value ) { return value; }

//! Extract parameter value.
/*!
 * Extracts the value of a parameter, given a range of data lines to search through and a
 * dictionary entry containing the characteristics of the parameter. This function searches
 * through the given set of data lines for a match with the required parameter, specified by an
 * iterator to a DictionaryEntry object. The search takes synonyms and case-(in)sensitivity into
 * account using the DictionaryComparer functor. A conversion function can also be specified to
 * convert the parameter value on the fly. Also, in case the parameter is not found, and is not
 * required, a default value can be specified.
 * \tparam DataType Data type of parameter value.
 * \param firstDataLine Iterator to first ParsedDataVector element to search through in range.
 * \param lastDataLine Iterator to last ParsedDataVector element to search through in range.
 * \param dictionaryEntry Iterator to a DictionaryEntry object for the requested parameter.
 * \param defaultValue Default value that is used in case the requested parameter is optional and
 *          can't be found.
 * \param convert Pointer to conversion function, to convert parameter value on the fly. The
 *          default value for this pointer is a dummy function that does nothing to the parameter
 *          value.
 * \sa DictionaryComparer, DictionaryEntry, ParsedDataVector, convertDummy().
 */
template< typename DataType >
DataType extractParameterValue( const DataLineIterator& firstDataLine,
                                const DataLineIterator& lastDataLine,
                                const DictionaryIterator& dictionaryEntry,
                                const DataType& defaultValue = DataType( ),
                                const boost::function< DataType( DataType ) >& convert
                                = &convertDummy< DataType > )
{
    // Attempt to match dictionary entry with any data lines.
    DataLineIterator parsedDataVectorIterator
            = std::find_if( firstDataLine, lastDataLine, DictionaryComparer( *dictionaryEntry ) );

    // If the data line found is the last data line (typically iterator to the end of the
    // container), then check if the dictionary entry is required.
    if ( parsedDataVectorIterator == lastDataLine )
    {
        // If the entry is required, throw an error indicating that the parameter cannot be found
        // in the input stream.
        if ( ( *dictionaryEntry )->isRequired )
        {
           throw std::runtime_error(
                                "Required parameter " + ( *dictionaryEntry )->parameterName
                                + " not found in input stream! "  );
        }

        // Else, return the default value specified.
        else
        {
            return defaultValue;
        }
    }

    // Else, set the flag in the dictionary entry indicating that it is extracted to true and
    // return the extract value using the getField() function.
    else
    {
        ( *dictionaryEntry )->isExtracted = true;

        return convert( parsed_data_vector_utilities::getField< DataType >(
                    *parsedDataVectorIterator, field_types::general::parameterValue ) );
    }

    // Dummy return to ensure no warnings are output.
    return DataType( );
}

} // namespace dictionary
} // namespace input_output
} // namespace tudat

#endif // TUDAT_DICTIONARY_TOOLS_H
