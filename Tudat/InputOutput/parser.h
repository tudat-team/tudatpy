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

#ifndef TUDAT_PARSER_H
#define TUDAT_PARSER_H

#include <boost/shared_ptr.hpp>

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

//! Typedef for shared-pointer to Parser object.
typedef boost::shared_ptr< Parser > ParserPointer;

} // namespace input_output
} // namespace tudat

#endif // TUDAT_PARSER_H
