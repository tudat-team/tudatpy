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
 *      111103    S. Billemont      File created.
 *      120521    T. Secretin       Layout changes.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_EXTRACTOR_H
#define TUDAT_EXTRACTOR_H

#include <boost/exception/all.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <map>
#include <stdarg.h>

#include "Tudat/InputOutput/parsedDataVectorUtilities.h"

namespace tudat
{
namespace input_output
{

template < class T >
class Extractor
{
public:

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~Extractor( ) { }

    //! Extract function.
    /*!
     * Extracts parsed data into a given data type and returns a shared pointer to it.
     */
    virtual boost::shared_ptr< T > extract(
            tudat::input_output::parsed_data_vector_utilities::ParsedDataLineMapPtr data ) = 0;

protected:

    //! Short-hand notations.
    typedef parsed_data_vector_utilities::FieldValuePtr FieldValuePtr;
    typedef parsed_data_vector_utilities::ParsedDataLineMapPtr ParsedDataLineMapPtr;
    typedef parsed_data_vector_utilities::ParsedDataVectorPtr ParsedDataVectorPtr;

    //! Check if the given data map contains all the passed FieldTypes. Returns false on fail.
    /*!
     * When this method encounters a FieldType that is not present in the passed data map, it
     * returns false otherwise true.
     * \param dataLineMap Check this datamap for all the requested FieldTypes.
     * \param numberOfFields Number of FieldTypes passed along to the function (after the
     * numberOfFields field).
     * \param ... Variable number of Fieldtypes to search for.
     * \return True if all the passed FieldTypes are present in the passed datamap.
     */
    bool checkOptionalFieldType( ParsedDataLineMapPtr dataLineMap, int numberOfFields, ... )
    {
        // Create a fancy vector (list) of all the fields:
        // Define argument list variable.
        va_list listOfArguments;

        // Initialize list. Point to last defined argument.
        va_start( listOfArguments, numberOfFields );

        for ( int i = 0; i < numberOfFields; i++ )
        {
            // Assign current field type being tested.
            tudat::input_output::FieldType testForType = va_arg( listOfArguments,
                                                                 tudat::input_output::FieldType );

            // Check if field type is found in the data map.
            if ( dataLineMap->find( testForType ) == dataLineMap->end( ) )
            {
                // Clean up the system stack.
                va_end( listOfArguments );

                // Not found in map, return false.
                return false;
            }
        }

        // Clean up the system stack.
        va_end ( listOfArguments );

        // All field types are present, return true.
        return true;
    }

    //! Check if the given data map contains all the passed FieldTypes. Throws an exception on fail.
    /*!
     * When this method encounters a FieldType that is not present in the passed data map, it throws
     * a boost::enable_error_info (std::runtime_error).
     * \param dataLineMap Check this datamap for all the requested FieldTypes.
     * \param numberOfFields Number of FieldTypes passed along to the function (after the
     * numberOfFields field).
     * \param ... Variable number of Fieldtypes to search for.
     */
    void checkRequiredFieldType( ParsedDataLineMapPtr dataLineMap, int numberOfFields, ... )
    {
        // Create a fancy vector (list) of all the fields:
        // Define argument list variable.
        va_list listOfArguments;

        // Initialize list. Point to last defined argument.
        va_start( listOfArguments, numberOfFields );

        for ( int i = 0; i < numberOfFields; i++ )
        {
            // Assign current field type being tested.
            tudat::input_output::FieldType testForType = va_arg( listOfArguments,
                                                                 tudat::input_output::FieldType );

            // Check if field type is found in the data map.
            if ( dataLineMap->find( testForType ) == dataLineMap->end( ) )
            {
                // Clean up the system stack.
                va_end( listOfArguments );

                // Not found in map, throw exception.
                boost::throw_exception( boost::enable_error_info( std::runtime_error(
                    "One of the types required for extracting this dataLineMap type is not present"
                                                                      ) ) );
            }
        }

        // Clean up the system stack.
        va_end ( listOfArguments );
    }

private:
};

} // namespace input_output
} // namespace tudat

#endif // TUDAT_EXTRACTOR_H
