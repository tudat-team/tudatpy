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

#ifndef TUDAT_TEXT_PARSER_H
#define TUDAT_TEXT_PARSER_H

#include <iostream>

#include <boost/make_shared.hpp>

#include "Tudat/InputOutput/parser.h"

namespace tudat
{
namespace input_output
{

//! A TextParser is a Parser that specialized in parsing text files.
/*!
 * A TextParser provides a set of functions to streamline parsing to simplify the parsing
 * process. The inheriting parsers can chose there prefered way of processing by calling 
 * TextParser(bool).
 *
 * NOTE: This TextParser works with the FieldValue/FieldType architecture.
 * For simpler file reading, use, for instance, matrixTextFileReader.
 */
class TextParser : public Parser
{
public:

    //! Create the default TextParser with processAsStream == false.
    /*!
     * The default constructor for TextParser causes the parser to behave as a line based parser.
     */
    TextParser( )
        : parsedData( boost::make_shared< parsed_data_vector_utilities::ParsedDataVector >( ) ),
          parseAsStream( false )
    { }

    //! Create the TextParser in the given process mode.
    /*!
     * Create a Textparser that deligates parsing either as stream or on a line-by-line basis.
     * \param processAsStream Boolean that determines whether the parsing should be done as a
     * stream (string is default).
     */
    TextParser( bool processAsStream )
        : parsedData( boost::make_shared< parsed_data_vector_utilities::ParsedDataVector >( ) ),
          parseAsStream( processAsStream )
    { }

    //! Default destructor, no new objects besides smart ones.
    virtual ~TextParser( ) { }

    //! Implementation of parse for strings.
    /*!
     * \param string String that is to be parsed.
     * \see Parser::parse(std::string& string).
     */
    parsed_data_vector_utilities::ParsedDataVectorPtr parse( std::string& string );

    //! Implementation of parse for an istream.
    /*!
     * \param stream Stream that is to be parsed.
     * \see Parser::parse(std::istream& stream).
     */
    parsed_data_vector_utilities::ParsedDataVectorPtr parse( std::istream& stream );
	
protected:

    //! Data container of the parsed data.
    /*!
     * Cleared and refilled on every call to parse(string) or parse(stream).
     */
    parsed_data_vector_utilities::ParsedDataVectorPtr parsedData;

    //! Flag to identify either stream parsing or line based parsing.
    /*!
     * Clients doing stream based parsing must override parseStream!
     * Clients doing line   based parsing must override parseLine!
     * Clients doing both parsing techniques must override parseLine and parseStream.
     */
    bool parseAsStream;
	
    //! Parse the given line content and append the resulting data lines to parsedData.
    /*!
     * Parse all lines/fields from the passed string and store (append) them to parsedData.
     * 
     * Needs to be overwritten if parseAsStream==false otherwise an exception will be thrown!
     * 
     * \param line String to parse
     */
    virtual void parseLine( std::string& line )
    {
        boost::throw_exception( boost::enable_error_info( std::runtime_error
                                                         ( "Must be overriden to be used" ) ) );
    }

    //! Parse the given stream content and append the resulting data lines to parsedData.
    /*!
     * Parse all lines/fields from the passed stream and store (append) them to parsedData.
     *
     * Needs to be overwritten if parseAsStream==true otherwise an exception will be thrown!
     *
     * \param stream Stream to parse
     */
    virtual void parseStream( std::istream& stream )
    {
        boost::throw_exception( boost::enable_error_info( std::runtime_error
                                                         ( "Must be overriden to be used" ) ) );
    }

private:
};

} // namespace input_output
} // namespace tudat

#endif // TUDAT_TEXT_PARSER_H
