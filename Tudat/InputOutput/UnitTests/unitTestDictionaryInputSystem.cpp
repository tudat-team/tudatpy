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

#include <fstream>
#include <stdexcept>
#include <string>

#include <boost/algorithm/string.hpp>
#include <boost/assign.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/InputOutput/streamFilters.h"

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/dictionaryComparer.h"
#include "Tudat/InputOutput/dictionaryEntry.h"
#include "Tudat/InputOutput/dictionaryTools.h"
#include "Tudat/InputOutput/fieldType.h"
#include "Tudat/InputOutput/parsedDataVectorUtilities.h"
#include "Tudat/InputOutput/separatedParser.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_dictionary_input_system )

//! Get example dictionary.
/*!
 * Returns an example dictionary with entries used by the unit tests.
 * \return Shared pointer to an example dictionary.
 */
input_output::dictionary::DictionaryPointer getExampleDictionary( )
{
    using namespace input_output::dictionary;
    using boost::assign::list_of;
    using boost::make_shared;

    // Create new dictionary.
    DictionaryPointer dictionary = boost::make_shared< Dictionary >( );

    // Insert example entries.
    addEntry( dictionary, "stringParameter",             true,  false                       );
    addEntry( dictionary, "stringParameterWithSynonym",  false, true, list_of( "sParam" )   );
    addEntry( dictionary, "doubleParameter",             false, false                       );
    addEntry( dictionary, "integerParameterWithSynonym", true,  true, list_of( "intParam" ) );
    addEntry( dictionary, "missingOptionalParameter",    false, false                       );
    addEntry( dictionary, "missingRequiredParameter",    true,  false                       );

    return dictionary;
}

//! Read in input file and filters out comment lines.
/*!
 * Reads in an input file (ASCII) and filters out comment lines.
 * \param inputFileName input file name.
 * \param commentCharacter Comment character used to denote comment lines
 *          (default is taken as '#').
 * \return filteredData Filtered data as string.
 */
std::string readAndFilterInputFile( const std::string& inputFileName,
                                    const char commentCharacter = '#' )
{
    // Create input file stream.
    std::ifstream inputFileStream( inputFileName.c_str( ) );

    // Declare unfiltered data string.
    std::string unfilteredData;

    // Declare line data string.
    std::string line;

    // Read input file line-by-line and store in unfilteredData.
    if ( inputFileStream.is_open( ) )
    {
        while ( std::getline( inputFileStream, line ) )
        {
            if ( !line.empty( ) ) { unfilteredData += line + "\n"; }
        }
    }

    // Close input file.
    inputFileStream.close( );

    // Create filter processor to filter out comment lines.
    boost::iostreams::filtering_ostream filterProcessor;

    // Add remove comment lines filter.
    filterProcessor.push( input_output::stream_filters::RemoveComment( commentCharacter ) );

    // Create filtered data string.
    std::string filteredData;

    // Last step in the chain; store the resulting string in filteredData.
    filterProcessor.push( boost::iostreams::back_inserter( filteredData ) );

    // Push in unfiltered data.
    filterProcessor << unfilteredData;

    // Make sure the data is processed by the chain.
    filterProcessor.flush( );

    // Trim all stray characters.
    boost::trim( filteredData );

    // Return filtered data string.
    return filteredData;
}

//! Test if parsing input file using dictionary works.
BOOST_AUTO_TEST_CASE( testInputFileParsingUsingDictionary )
{
    using namespace input_output::field_types;
    using namespace input_output::parsed_data_vector_utilities;
    using namespace input_output::dictionary;
    using unit_conversions::convertKilometersToMeters;
    using boost::assign::list_of;

    // Set input file.
    const std::string inputFile = input_output::getTudatRootPath( )
            + "InputOutput/UnitTests/exampleInputFile.in";

    // Read in input file and filter out comment lines.
    std::string filteredData = readAndFilterInputFile( inputFile );

    // Declare a separated parser.
    input_output::SeparatedParser parser( std::string( ": " ), 2,
                                                 general::parameterName,
                                                 general::parameterValue );

    // Parse filtered data.
    const ParsedDataVectorPtr parsedData = parser.parse( filteredData );

    // Get example dictionary.
    const DictionaryPointer dictionary = getExampleDictionary( );

    // Test 1: Check if any required parameters are missing, and catch the error thrown.
    {
        // Set flag indicating that required parameters are missing to false. True means
        // "missing", false means "not missing".
        bool areRequiredParametersMissingBeforeExtraction = false;

        // Check if required parameters are missing.
        try
        {
            checkRequiredParameters( dictionary );
        }

        // Catch the error thrown.
        catch ( std::runtime_error& requiredParametersError )
        {
            areRequiredParametersMissingBeforeExtraction = true;
        }

        // Check that flag is set to true, i.e. parameters are missing, as expected.
        BOOST_CHECK( areRequiredParametersMissingBeforeExtraction );
    }

    // Test 2: test extraction of "stringParameter" (required, not case-sensitive, no conversion,
    //         not synonym).
    {
        std::string stringParameter = extractParameterValue< std::string >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "stringParameter" ) );

        BOOST_CHECK_EQUAL( stringParameter, "this is a string parameter" );
    }

    // Test 3: test extraction of "sParam" (optional, case-sensitive, no conversion, synonym for
    //         "stringParameterWithSynonym").
    {
        std::string stringParameterWithSynonym = extractParameterValue< std::string >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "stringParameterWithSynonym" ) );

        BOOST_CHECK_EQUAL( stringParameterWithSynonym,
                           "this string was extracted using a synonym" );
    }

    // Test 4: test extraction of "intParam" (required, case-sensitive, no conversion, synonym for
    //         "integerParameterWithSynonym").
    {
        int integerParameterWithSynonym = extractParameterValue< int >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "integerParameterWithSynonym" ) );

        BOOST_CHECK_EQUAL( integerParameterWithSynonym, 12 );
    }

    // Test 5: test extraction of "doubleParameter" (required, not case-sensitive,
    //         kilometers-to-meters conversion, not synonym).
    {
        double doubleParameter = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "doubleParameter" ), 0.0,
                    &convertKilometersToMeters< double > );

        BOOST_CHECK_EQUAL( doubleParameter, 987650.0 );
    }

    // Test 6: test extraction of missing, optional parameter in input stream.
    {
        double missingOptionalParameter = extractParameterValue< double >(
                    parsedData->begin( ), parsedData->end( ),
                    findEntry( dictionary, "missingOptionalParameter" ), 123.45 );

        BOOST_CHECK_EQUAL( missingOptionalParameter, 123.45 );
    }

    // Test 7: test extraction of missing, required parameter in input stream.
    {
        // Set flag indicating that required parameter is not missing in input stream.
        bool isRequiredParameterMissing = false;

        // Try to extract missing parameter.
        try
        {
            std::string integerParameterWithSynonym = extractParameterValue< std::string >(
                        parsedData->begin( ), parsedData->end( ),
                        findEntry( dictionary, "missingRequiredParameter" ) );
        }

        // Catch error thrown.
        catch ( std::runtime_error& nonExistentRequireParameterError )
        {
            isRequiredParameterMissing = true;
        }

        // Check that flag is set to true, i.e. required parameter is missing, as expected.
        BOOST_CHECK( isRequiredParameterMissing );
    }

    // Test 8: Check if any required parameters are missing.
    {
        // Set flag indicating that required parameters are missing to true. True means
        // "missing", false means "not missing".
        bool areRequiredParametersMissingAfterFullExtraction = false;

        // Set "missingParameter" dictionary entry to extracted.
        Dictionary::const_iterator missingParameterIterator
                = findEntry( dictionary, "missingRequiredParameter" );
        ( *missingParameterIterator )->isExtracted = true;

        // Check if required parameters are missing.
        try
        {
            checkRequiredParameters( dictionary );
        }

        // Catch the error thrown.
        catch ( std::runtime_error& requiredParametersError )
        {
            areRequiredParametersMissingAfterFullExtraction = true;
        }

        // Check that flag is set to false, i.e. no parameters are missing, as expected.
        BOOST_CHECK( !areRequiredParametersMissingAfterFullExtraction );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
