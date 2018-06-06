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

#ifndef TUDAT_DICTIONARY_ENTRY_H
#define TUDAT_DICTIONARY_ENTRY_H

#include <set>
#include <string>

#include <memory>

namespace tudat
{
namespace input_output
{
namespace dictionary
{

//! Dictionary entry.
/*!
 * This struct contains all of the information pertaining to a single dictionary entry:
 * a parameter name, a flag to indicate if the parameter is required or optional, a flag to
 * indicate if the parameter name is case-sensitive or not, a flag to indicate if the parameter has
 * been extracted, and a list of synonyms for the given parameter.
 */
struct DictionaryEntry
{
public:

    //! Typedef for set of strings.
    typedef std::set< std::string > StringSet;

    //! Zombie constructor.
    /*!
     * Zombie constructor, initializing all of the class members to zombie state. This is needed to
     * create a dictionary, which is an STL set of shared-pointers to DictionaryEntry objects.
     * \sa Dictionary.
     */
    DictionaryEntry( )
        : parameterName( "" ),
          isRequired( true ),
          isCaseSensitive( false ),
          isExtracted( false ),
          synonyms( )
    { }

    //! Default constructor.
    /*!
     * Default constructor, initializing all of the class members.
     * \param aParameterName Name of parameter for dictionary entry.
     * \param required Boolean indicating if parameter is required. Default set to true.
     * \param caseSensitive Boolean indicating if parameter is caseSensitive. Default set to false.
     * \param someSynonyms Set of synonyms.
     */
    DictionaryEntry( std::string aParameterName,
                     bool required = true,
                     bool caseSensitive = false,
                     StringSet someSynonyms = StringSet( ) )
        : parameterName( aParameterName ),
          isRequired( required ),
          isCaseSensitive( caseSensitive ),
          isExtracted( false ),
          synonyms( someSynonyms )
    { }

    //! Parameter name.
    /*!
     * String name of parameter.
     */
    std::string parameterName;

    //! Flag that sets if parameter is required or not.
    /*!
     * This flag indicates if the given dictionary entry must be provided as input. If the
     * required parameter name is not provided, an error is thrown. Default is set to true, i.e.,
     * parameter IS required.
     */
    bool isRequired;

    //! Flag that sets if parameter is case-sensitive.
    /*!
     * This flag indicates if the given dictionary entry is case-sensitive. If it is not
     * case-senstive, both the input parameter and the dictionary entry are compared after
     * transforming them to upper-case. Default is set to false, i.e., parameter IS NOT
     * case-sensitive.
     */
    bool isCaseSensitive;

    //! Flag that sets if parameter has been extracted already.
    /*!
     * This flag indicates if the given dictionary entry has already been extracted. The
     * default is set to false.
     */
    bool isExtracted;

    //! List of synonyms.
    /*!
     * List of synonyms stored as a set of strings.
     */
    StringSet synonyms;
};

//! Typedef for shared-pointer to dictionary entry.
typedef std::shared_ptr< DictionaryEntry > DictionaryEntryPointer;

} // namespace dictionary
} // namespace input_output
} // namespace tudat

#endif // TUDAT_DICTIONARY_ENTRY_H
