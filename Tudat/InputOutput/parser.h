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
 *      111103    S. Billemont      Creation of code.
 *      120326    D. Dirkx          Code checked, minor layout changes.
 *
 *    References
 *
 *    Notes
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
