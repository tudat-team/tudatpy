/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TEXT_PARSING_H
#define TUDAT_TEXT_PARSING_H

#include <algorithm>
#include <cmath>
#include <locale>
#include <optional>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

namespace tudat
{
namespace input_output
{
namespace text_parsing
{

//! Remove leading and trailing characters from the given set (default: spaces, tabs and line ends).
inline std::string trim( const std::string& text, const char* whitespace = " \t\r\n" )
{
    const auto first = text.find_first_not_of( whitespace );
    return first == std::string::npos ? "" : text.substr( first, text.find_last_not_of( whitespace ) - first + 1 );
}

//! Split a line into whitespace-separated tokens, independently of the global locale.
inline std::vector< std::string > splitIntoTokens( const std::string& text )
{
    std::istringstream stream( text );
    stream.imbue( std::locale::classic( ) );
    std::vector< std::string > tokens;
    std::string token;
    while( stream >> token )
    {
        tokens.push_back( token );
    }
    return tokens;
}

//! Parse a complete token as a finite number, accepting Fortran D exponents; empty if the token is invalid.
inline std::optional< double > parseFiniteNumber( std::string token )
{
    std::replace( token.begin( ), token.end( ), 'D', 'E' );
    std::replace( token.begin( ), token.end( ), 'd', 'e' );
    try
    {
        std::size_t consumed = 0;
        const double value = std::stod( token, &consumed );
        if( consumed == token.size( ) && std::isfinite( value ) )
        {
            return value;
        }
    }
    catch( const std::exception& )
    {}
    return std::nullopt;
}

//! Parse a token of decimal digits as a non-negative integer; empty if the token is invalid or out of range.
inline std::optional< int > parseNonNegativeInteger( const std::string& token )
{
    if( token.empty( ) || !std::all_of( token.begin( ), token.end( ), []( const char c ) { return c >= '0' && c <= '9'; } ) )
    {
        return std::nullopt;
    }
    try
    {
        return std::stoi( token );
    }
    catch( const std::exception& )
    {
        return std::nullopt;
    }
}

}  // namespace text_parsing
}  // namespace input_output
}  // namespace tudat

#endif  // TUDAT_TEXT_PARSING_H
