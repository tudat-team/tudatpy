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

#ifndef TUDAT_STREAM_FILTERS_H
#define TUDAT_STREAM_FILTERS_H

#include <string>

#include <boost/iostreams/filter/line.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/regex.hpp>

namespace tudat
{
namespace input_output
{
namespace stream_filters
{

//! Filter that removes comments from a stream.
/*!
 * This class allows for a data stream to be filtered to remove comments starting with a
 * given character, e.g., "#".
 * When this character is the first of a line, the entire line is considered comment, and is removed.
 * When the comment appears in the middle of a line, only the text after this line is removed.
 * NOTE: The start of a comment is always a single character!
 * NOTE: You cannot escape the comment character!
 * NOTE: A comment always runs until the end of the same line!
 */
class RemoveComment : public boost::iostreams::line_filter
{
public:

    //! Create a comment filter for a given start comment character.
    /*!
     * Creates a comment filter for a given start comment character.
     * Note that the variable "line_filter" is set to true, to indicate to Boost not to append
     * newline characters.
     * \param skipCharacter Character that initiates a comment.
     * \param isOmitIfEmpty True if a line should be removed when a complete line is filtered away.
     *                      (When the 'skipCharacter' is the first character on the line)
     */
    RemoveComment( char skipCharacter = '#', bool isOmitIfEmpty = true )
        : boost::iostreams::line_filter( true ), 
          skipCharacter_( skipCharacter ),
          isOmitIfEmpty_( isOmitIfEmpty ) 
    { }

    //! Remove comments in a single line.
    /*!
     * Executes filter to remove comments from input line of data. 
     * For more information see boost::iostreams::line_filter (overrides definition of line_filter
     * to enable chaining of filters).
     * \param line Unfiltered line.
     * \return Filtered line.
     */
    std::string do_filter( const std::string& line );

private:

    //! Skip character.
    /*!
     * Character that initiates a comment.
     */
    char skipCharacter_;

    //! Omit state for filtered line.
    /*!
     * If filtered line is empty, and omit flag is set to true, no newline empty line is returned.
     */
    bool isOmitIfEmpty_;
};

//! Filter for skipping the first several lines in a stream.
/*!
 * This class allows for a data stream to be filtered by skipping a given number of lines from the
 * start of the stream.
 */
class SkipFirstLines : public boost::iostreams::line_filter
{
public:

    //! Create a filter to skip a given amount of lines.
    /*!
     * Creates a filter to skip a given amount of lines by passing in the number of lines to skip
     * and the action for empty lines. 
     * Note that the variable "line_filter" is set to true, to indicate to Boost not to append
     * newline characters.
     * \param linesToSkip Number of lines (counted sequentially from the stream head) to discard.
     * \param isOmitIfEmpty True if a line should be removed when a complete line is filtered away.
     *       (If false, there will be 'linesToSkip' newline characters is the start of the stream.)
     */
    SkipFirstLines( unsigned int linesToSkip = 0, bool isOmitIfEmpty = true )
        : boost::iostreams::line_filter( true ), 
          linesToSkip_( linesToSkip ),
          numberOfSkippedLines_( 0 ), 
          isOmitIfEmpty_( isOmitIfEmpty ) 
    { }

    //! Execute filter on the input line by skipping if required.
    /*!
     * This object will ignore the first 'linesToSkip_' number of passed in lines before returning 
     * the same line.
     * For more information see boost::iostreams::line_filter (overrides definition of line_filter
     * to enable chaining of filters).
     * \param line Unfiltered line.
     * \return Filtered line.
     */
    std::string do_filter( const std::string& line );

private:

    //! Lines to skip.
    /*!
     * Lines to skip, counted sequentially from the head of a stream.
     */
    unsigned int linesToSkip_;

    //! Number of skipped lines.
    /*!
     * Counter of the number of lines already skipped in the stream.
     */
    unsigned int numberOfSkippedLines_;

    //! Omit state for filtered line.
    /*!
     * If filtered line is empty, and omit flag is set to true, no newline empty line is returned.
     */
    bool isOmitIfEmpty_;
};

//! Filter for searching and replacing text in a stream
/*!
 * This filter allows the user to perform a search&replace operation on a stream (line based).
 * If a match is found, it is replaced by the given replace string.
 *
 * The search pattern is created using a regex expression (boost::regex). The following are example 
 * search queries:
 * 
 * Search for keyword:                      'myKeyWord'
 *      This will match all exact matches of 'myKeyWord'.
 *
 * Search for entire line with keyword:     '^.*myKeyWord.*$'
 *      This will match the line (from start to finish) if it contains the string "myKeyWord".
 *
 * Search for numbers with two decimals:    '[\+\-]*[0-9]*\.[0-9]{2}'
 *      This will match numbers up to two decimals (including preceding + and -).
 *
 * More examples: http://www.boost.org/libs/regex/
 *                http://www.regular-expressions.info/examples.html
 */
class ReplaceElements : public boost::iostreams::line_filter
{
public:

    //! Create filter with a regex object and replace string.
    /*!
     * Creates a ReplaceElements filter with a provided regex search pattern object and a
     * replacement for any matches. If you want to remove any match leave replaceString empty.
     * Note that "line_filter" is set to true, to indicate to Boost not to append newline 
     * characters.
     * \param searchFilter Regex search pattern to search for.
     * \param replaceString String to replace any search matches with.
     * \param isOmitIfEmpty True if a line should be removed when a complete line is filtered away.
     *              (If true, if an entire line is removed, no newline character will be returned.)
     */
    ReplaceElements( boost::regex searchFilter, std::string replaceString="",
                     bool isOmitIfEmpty = true )
        : boost::iostreams::line_filter( true ), 
          searchFilter_( searchFilter ),
          replaceString_( replaceString ), 
          isOmitIfEmpty_( isOmitIfEmpty ) 
    { }
    
    //! Create filter with a basic search and replace string.
    /*!
     * Creates a ReplaceElements filter with a provided search string object and a
     * replacement for any matches. Note this is for exact matches only (creates a regex object and
     * escapes the following charachters: .[]{}()\*+?|^$ ).
     * If you want to remove any match leave replaceString empty.
     * \param searchFilter Regex search pattern to search for
     * \param replaceString String to replace any search matches with
     * \param isOmitIfEmpty True if a line should be removed when a complete line is filtered away.
     *              (If true, if an entire line is removed, no newline character will be returned.)
     */
    ReplaceElements( std::string searchFilter, std::string replaceString = "",
                     bool isOmitIfEmpty = true );

    //! Execute filter on a single line to replace matched elements with replace string.
    /*!
     * Find search pattern matches and replace them with a predefined string.
     * For more information see boost::iostreams::line_filter (overrides definition of line_filter
     * to enable chaining of filters).
     * \param line Unfiltered line.
     * \return Filtered line.
     */
    std::string do_filter( const std::string& line );

private:

    //! Search pattern.
    /*!
     * Search pattern used to find matches to replace.
     * This is an encapsulated regular expression, see http://www.regular-expressions.info/
     */
    boost::regex searchFilter_;

    //! Replace string.
    /*!
     * Search matches are replaced with this string.
     */
    std::string replaceString_;

    //! Omit state for filtered line.
    /*!
     * If filtered line is empty, and omit flag is set to true, no newline empty line is returned.
     */
    bool isOmitIfEmpty_;
};

} // namespace stream_filters
} // namespace input_output
} // namespace tudat

#endif // TUDAT_STREAM_FILTERS_H
