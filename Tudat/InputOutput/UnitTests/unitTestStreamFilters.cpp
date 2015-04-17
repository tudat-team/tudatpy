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
 *      120127    S. Billemont      File created.
 *      120127    K. Kumar          Added missing comments and clarified variable-naming.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <iostream>
#include <sstream>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/InputOutput/streamFilters.h"

namespace tudat
{
namespace unit_tests
{

//! Test fixture for stream filters.
/*!
 * This defines all the predefined memebers and functions that are available to each test.
 */
struct StreamFilterFixture
{
public:

    //! Constructor.
    /*!
     * Constructor. DO NOT INDENT testString initialization all spaces will end up in the test
     * string!
     */
    StreamFilterFixture( ) :
    testString("\
The first header line\n\
A line with partial comment # this is comment\n\
# Complete comment line\n\
-> odd data <-\n\0") { }

    //! Container for the test data.
    /*!
     * Container for the test data; reinitalized for each test.
     */
    std::string testString;

    //! Container for the resulting filtered data.
    /*!
     * Container for the resulting filtered data.
     */
    std::string filteredData;

protected:

private:
};

//! Create Boost fixture test suite for all the stream filter unit tests.
BOOST_FIXTURE_TEST_SUITE( test_suite_streamFilters, StreamFilterFixture )

//! Test with remove line endings on the skipped lines.
BOOST_AUTO_TEST_CASE( test_skipFirstLines_noEmptyLines )
{  
    // Create a filter chain, attach the test filter and push in the testStream.
    // Filter chain object.
    boost::iostreams::filtering_ostream filterProcessor;

    // Add skip first lines filter.
    filterProcessor.push( input_output::stream_filters::SkipFirstLines( 2, true ) );

    // Last step in the chain; store the resulting string in result.
    filterProcessor.push(boost::iostreams::back_inserter( filteredData ) );

    // Push in test data.
    filterProcessor << testString;

    // Make sure the data is processed by the chain.
    filterProcessor.flush( );

    // Check if the results match the expected result.
    BOOST_CHECK( filteredData.compare(
"# Complete comment line\n-> odd data <-\n"
                     ) == 0 );
}

//! Test with don't remove line endings on the skipped lines.
BOOST_AUTO_TEST_CASE( test_skipFirstLines_emptyLines )
{
    // Create a filter chain, attach the test filter and push in the testStream.
    // Filter chain object.
    boost::iostreams::filtering_ostream filterProcessor;

    // Add skip first lines filter.
    filterProcessor.push( input_output::stream_filters::SkipFirstLines( 2, false ) );

    // Last step in the chain; store the resulting string in result.
    filterProcessor.push( boost::iostreams::back_inserter( filteredData ) );

    // Push in test data.
    filterProcessor << testString;

    // Make sure the data is processed by the chain
    filterProcessor.flush( );

    // Check if the results match the expected result.
    BOOST_CHECK( filteredData.compare(
"\n\n# Complete comment line\n-> odd data <-\n"
                     ) == 0 );
}

//! Test with remove line endings on the tested lines.
BOOST_AUTO_TEST_CASE( test_removeComment_noEmptyLines )
{
    // Create a filter chain, attach the test filter and push in the testStream.
    // Filter chain object
    boost::iostreams::filtering_ostream filterProcessor;

    // Add skip first lines filter.
    filterProcessor.push( input_output::stream_filters::RemoveComment( '#', true ) );

    // Last step in the chain; store the resulting string in result.
    filterProcessor.push( boost::iostreams::back_inserter( filteredData ) );

    // Push in test data.
    filterProcessor << testString;

    // Make sure the data is processed by the chain.
    filterProcessor.flush( );

    // Check if the results match the expected result
    BOOST_CHECK( filteredData.compare(
"The first header line\nA line with partial comment \n-> odd data <-\n"
                     ) == 0 );
}

//! Test with dont remove line endings on the tested lines.
BOOST_AUTO_TEST_CASE( test_removeComment_emptyLines )
{
    // Create a filter chain, attach the test filter and push in the testStream.
    // Filter chain object.
    boost::iostreams::filtering_ostream filterProcessor;

    // Add skip first lines filter.
    filterProcessor.push( input_output::stream_filters::RemoveComment( '#', false ) );

    // Last step in the chain; store the resulting string in result
    filterProcessor.push(boost::iostreams::back_inserter( filteredData ) );

    // Push in test data.
    filterProcessor << testString;

    // Make sure the data is processed by the chain.
    filterProcessor.flush( );

    // Check if the results match the expected result.
    BOOST_CHECK( filteredData.compare(
"The first header line\nA line with partial comment \n\n-> odd data <-\n"
                    ) == 0 );
}

//! Test with remove line endings on the tested lines.
BOOST_AUTO_TEST_CASE( test_replaceElements_delete_noEmptyLines )
{
    // Create a filter chain, attach the test filter and push in the testStream.
    // Filter chain object.
    boost::iostreams::filtering_ostream filterProcessor;

    // Regex '-.*-' is everything between two dashes and replace with nothing.
    filterProcessor.push( input_output::stream_filters::ReplaceElements(
                    boost::regex("-.*-"), "", false) );

    // Last step in the chain; store the resulting string in result.
    filterProcessor.push( boost::iostreams::back_inserter( filteredData ) );

    // Push in test data.
    filterProcessor << testString;

    // Make sure the data is processed by the chain.
    filterProcessor.flush( );

    // Check if the results match the expected result.
    BOOST_CHECK( filteredData.compare(
// DONT INDENT! all spaces will end up in the test line
"The first header line\nA line with partial comment # this is comment\n\
# Complete comment line\n\n"
        ) == 0 );
}

//! Test with don't remove line endings on the tested lines.
BOOST_AUTO_TEST_CASE( test_replaceElements_delete_emptyLines )
{
    // Create a filter chain, attach the test filter and push in the testStream.
    // Filter chain object.
    boost::iostreams::filtering_ostream filterProcessor;

    // Regex '-.*-' is everything between two dashes and replace with nothing.
    filterProcessor.push( input_output::stream_filters::ReplaceElements(
                     boost::regex( "-.*-" ), "", false ) );

    // Last step in the chain; store the resulting string in result.
    filterProcessor.push( boost::iostreams::back_inserter( filteredData ) );

    // Push in test data.
    filterProcessor << testString;

    // Make sure the data is processed by the chain.
    filterProcessor.flush( );

    // Check if the results match the expected result
    BOOST_CHECK( filteredData.compare(
// DONT INDENT! all spaces will end up in the test line
"The first header line\nA line with partial comment # this is comment\n\
# Complete comment line\n\n"
        ) == 0 );
}

//! Test with remove line endings on the tested lines.
BOOST_AUTO_TEST_CASE(test_replaceElements_replace)
{
    // Create a filter chain, attach the test filter and push in the testStream.
    // Filter chain object.
    boost::iostreams::filtering_ostream filterProcessor;

    // Regex '>.*<' is everything between two angle brackets and replace with 'foobar'.
    filterProcessor.push( input_output::stream_filters::ReplaceElements(
                    boost::regex(">.*<"), "foobar", true) );

    // Last step in the chain; store the resulting string in result.
    filterProcessor.push( boost::iostreams::back_inserter( filteredData ) );

    // Push in test data.
    filterProcessor << testString;

    // Make sure the data is processed by the chain.
    filterProcessor.flush( );

    // Check if the results match the expected result.
    BOOST_CHECK( filteredData.compare(
// DONT INDENT! all spaces will end up in the test line
"The first header line\nA line with partial comment # this is comment\n\
# Complete comment line\n-foobar-\n"
        ) == 0 );
}

//! Test with creating a literal search string (so no regex object but just literal string)
BOOST_AUTO_TEST_CASE(test_replaceElements_literalConstuctor)
{
    // Create a filter chain, attach the test filter and push in the testStream.
    // Filter chain object.
    boost::iostreams::filtering_ostream filterProcessor;

    // Search for all the special regex characters and replace with foobar
    filterProcessor.push( input_output::stream_filters::ReplaceElements(
        ".[]{}()\\*+?|^$", "foobar", true) );

    // Last step in the chain; store the resulting string in result.
    filterProcessor.push( boost::iostreams::back_inserter( filteredData ) );

    // Push in test data.
    filterProcessor << "->.[]{}()\\*+?|^$<-\n";

    // Make sure the data is processed by the chain.
    filterProcessor.flush( );

    // Check if the results match the expected result.
    BOOST_CHECK( filteredData.compare("->foobar<-\n") == 0 );
}

BOOST_AUTO_TEST_SUITE_END( );

} // namespace unit_tests
} // namespace tudat
