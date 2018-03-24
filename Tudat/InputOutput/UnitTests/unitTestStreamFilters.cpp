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

#define BOOST_TEST_MAIN

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
