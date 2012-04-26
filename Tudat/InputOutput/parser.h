/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      111103    S. Billemont      First creation of code.
 *      120326    D. Dirkx          Code checked, minor layout changes.
 *
 *    References
 *
 */

#ifndef TUDAT_PARSER_H
#define TUDAT_PARSER_H

#include "Tudat/InputOutput/parsedDataVectorUtilities.h"

namespace tudat
{
namespace input_output
{

//! Parser interface, indicates that the class can parse passed string and stream data.
/*!
 * A parser is an interface to parse a specific file or filegroup. This can be either plain
 * text (e.g. ascii) or a binary file. When the data has been parsed, it returns a vector of
 * parsed lines or equivalent. Each line contains a map of (FieldType, FieldValue). This key
 * of the map is the type of the specific field (e.g. date) and the entry is the parsed value
 * for the specific line (e.g. 2000-1-1)
 *
 * The data vector can then be passed to an extractor to extract, for example, the date in Julian
 * date form
 *
 * NOTE: This TextParser works with the FieldValue/FieldType architecture.
 * For simpler file reading, use, for instance, matrixTextFileReader.
 */
class Parser
{
public:

    //! Default constructor, interface therefore empty.
    Parser( ) { }

    //! Default destructor, interface therefore empty.
    virtual ~Parser( ) { }

    //! Parse data from a string.
    /*!
     * \param string String reference to parse.
     * \return  Parsed data from the string as ParsedDataVector.
     */
    virtual parsed_data_vector_utilities::ParsedDataVectorPtr parse( std::string& string ) = 0;

    //! Parse data from a stream.
    /*!
     * \param stream Stream reference to parse.
     * \return  Parsed data from the stream as ParsedDataVector.
     */
    virtual parsed_data_vector_utilities::ParsedDataVectorPtr parse( std::istream& stream ) = 0;

protected:

private:
};

} // namespace input_output
} // namespace tudat

#endif // TUDAT_PARSER_H
