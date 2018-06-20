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

#ifndef TUDAT_TEXT_PARSER_H
#define TUDAT_TEXT_PARSER_H

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
        : parsedData( std::make_shared< parsed_data_vector_utilities::ParsedDataVector >( ) ),
          parseAsStream( false )
    { }

    //! Create the TextParser in the given process mode.
    /*!
     * Create a Textparser that deligates parsing either as stream or on a line-by-line basis.
     * \param processAsStream Boolean that determines whether the parsing should be done as a
     * stream (string is default).
     */
    TextParser( bool processAsStream )
        : parsedData( std::make_shared< parsed_data_vector_utilities::ParsedDataVector >( ) ),
          parseAsStream( processAsStream )
    { }

    //! Default destructor, no new objects besides smart ones.
    virtual ~TextParser( ) { }

    //! Implementation of parse for strings.
    /*!
     * \param string String that is to be parsed.
     * \return Parsed line data, stored in a ParsedDataVectorPtr object.
     * \see Parser::parse(std::string& string).
     */
    parsed_data_vector_utilities::ParsedDataVectorPtr parse( std::string& string );

    //! Implementation of parse for an istream.
    /*!
     * \param stream Stream that is to be parsed.
     * \return Parsed line data, stored in a ParsedDataVectorPtr object.
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
        throw std::runtime_error( "Function parseLine must be overriden to be used" );
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
        throw std::runtime_error( "Function parseStream must be overriden to be used" );

    }

private:
};

//! Typedef for shared-pointer to TextParser object.
typedef std::shared_ptr< TextParser > TextParserPointer;

} // namespace input_output
} // namespace tudat

#endif // TUDAT_TEXT_PARSER_H
