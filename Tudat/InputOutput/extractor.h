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

#ifndef TUDAT_EXTRACTOR_H
#define TUDAT_EXTRACTOR_H

#include <memory>
#include <map>
#include <cstdarg>

#include "Tudat/InputOutput/parsedDataVectorUtilities.h"

namespace tudat
{
namespace input_output
{

//! Extractor base class.
/*!
 * This abstract class belongs to the parser-extractor architecture implemented in Tudat. This base
 * class can be used to derive specific extractors that take data in the form of a
 * ParsedDataLineMapPtr object and return extracted values.
 * \tparam T Generic-type returned by extract() function.
 * \sa ParsedDataLineMapPtr, Parser
 */
template< class T >
class Extractor
{
public:

    //! Default destructor.
    virtual ~Extractor( ) { }

    //! Extract function.
    /*!
     * Extracts parsed data into a given data type and returns a shared pointer to it.
     * \param data Raw (parsed) data to extract values from.
     * \return Extracted value.
     */
    virtual std::shared_ptr< T > extract(
        parsed_data_vector_utilities::ParsedDataLineMapPtr data ) = 0;

protected:

    //! Short-hand notations.
    //! Shortcut to field value shared-pointer.
    typedef parsed_data_vector_utilities::FieldValuePtr FieldValuePtr;
    //! Shortcut to parsed data line map shared-pointer.
    typedef parsed_data_vector_utilities::ParsedDataLineMapPtr ParsedDataLineMapPtr;
    //! Shortcut to parsed data vector shared-pointer.
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
            input_output::FieldType testForType = va_arg( listOfArguments,
                                                                 input_output::FieldType );

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
        va_end( listOfArguments );

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
            input_output::FieldType testForType = va_arg( listOfArguments,
                                                                 input_output::FieldType );

            // Check if field type is found in the data map.
            if ( dataLineMap->find( testForType ) == dataLineMap->end( ) )
            {
                // Clean up the system stack.
                va_end( listOfArguments );

                // Not found in map, throw exception.
                throw std::runtime_error(
                    "One of the types required for extracting this dataLineMap type is not present" );
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
