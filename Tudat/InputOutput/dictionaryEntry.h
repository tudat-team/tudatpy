/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      120712    K. Kumar          File created.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_DICTIONARY_ENTRY_H
#define TUDAT_DICTIONARY_ENTRY_H

#include <set>
#include <string>

#include <boost/shared_ptr.hpp>

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
typedef boost::shared_ptr< DictionaryEntry > DictionaryEntryPointer;

} // namespace dictionary
} // namespace input_output
} // namespace tudat

#endif // TUDAT_DICTIONARY_ENTRY_H
