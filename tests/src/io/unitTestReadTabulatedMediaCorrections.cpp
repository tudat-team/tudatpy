/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include "tudat/basics/testMacros.h"

#include "tudat/io/readTabulatedMediaCorrections.h"

namespace tudat
{
namespace unit_tests
{

using namespace input_output;

BOOST_AUTO_TEST_SUITE( test_read_tabulated_media_corrections )

BOOST_AUTO_TEST_CASE( testCspCommandsExtraction )
{

    // Reading troposphere file
    {
        std::shared_ptr< CspRawFile > cspFile = std::make_shared< CspRawFile >(
                "/Users/pipas/Documents/mro-data/tro/mromagr2017_091_2017_121.tro.txt" );
        std::vector< std::shared_ptr< CspCommand > > extractedCspCommands = cspFile->getCspCommands( );

        std::vector< unsigned long > cspCommandsIds = { 0, 1, extractedCspCommands.size( ) - 1 };

        for ( unsigned int cspCommandId: cspCommandsIds )
        {
            std::shared_ptr< AtmosphericCorrectionCspCommand > atmosphericCommand = std::dynamic_pointer_cast< AtmosphericCorrectionCspCommand >(
                    extractedCspCommands.at( cspCommandId ) );

            // Expected start and final time obtained from: https://nsidc.org/data/icesat/glas-date-conversion-tool/date_convert/
            if ( cspCommandId == 0 )
            {
                std::vector< double > expectedCoefficients = {
                        -.0048, .0077, .0238, -.0872, -.0445, .1499, .0303, -.0995, -.0066, .0228 };
                BOOST_CHECK_EQUAL ( atmosphericCommand->modelIdentifier_, "WET NUPART" );
                BOOST_CHECK_EQUAL ( atmosphericCommand->dataTypesIdentifier_, "ALL" );
                BOOST_CHECK_EQUAL ( atmosphericCommand->computationSpecifier_, "BY NRMPOW" );
                for ( unsigned int i = 0; i < expectedCoefficients.size( ); ++i )
                {
                    BOOST_CHECK_EQUAL ( atmosphericCommand->computationCoefficients_.at( i ),
                                        expectedCoefficients.at( i ) );
                }
                BOOST_CHECK_EQUAL ( atmosphericCommand->groundStationsId_, "C10" );
                BOOST_CHECK_EQUAL ( atmosphericCommand->sourceSpecifier_, "" );
                BOOST_CHECK_EQUAL ( atmosphericCommand->sourceId_, 0 );
                BOOST_CHECK_EQUAL ( atmosphericCommand->startTime_, 544287600.00100 );
                BOOST_CHECK_EQUAL ( atmosphericCommand->endTime_, 544309200.000000 );
            }
            else if ( cspCommandId == 1 )
            {
                std::vector< double > expectedCoefficients = { .0018, .0039, -.0053, -.0002, .0015 };
                BOOST_CHECK_EQUAL ( atmosphericCommand->modelIdentifier_, "DRY NUPART" );
                BOOST_CHECK_EQUAL ( atmosphericCommand->dataTypesIdentifier_, "ALL" );
                BOOST_CHECK_EQUAL ( atmosphericCommand->computationSpecifier_, "BY NRMPOW" );
                for ( unsigned int i = 0; i < expectedCoefficients.size( ); ++i )
                {
                    BOOST_CHECK_EQUAL ( atmosphericCommand->computationCoefficients_.at( i ),
                                        expectedCoefficients.at( i ) );
                }
                BOOST_CHECK_EQUAL ( atmosphericCommand->groundStationsId_, "C10" );
                BOOST_CHECK_EQUAL ( atmosphericCommand->sourceSpecifier_, "" );
                BOOST_CHECK_EQUAL ( atmosphericCommand->sourceId_, 0 );
                BOOST_CHECK_EQUAL ( atmosphericCommand->startTime_, 544287600.00100 );
                BOOST_CHECK_EQUAL ( atmosphericCommand->endTime_, 544309200.000000 );
            }
            else if ( cspCommandId == extractedCspCommands.size( ) - 1 )
            {
                std::vector< double > expectedCoefficients = { .0031, .0043, -.0006, -.0054, -.0015, .0069, .0008,
                                                               -.0024 };
                BOOST_CHECK_EQUAL ( atmosphericCommand->modelIdentifier_, "DRY NUPART" );
                BOOST_CHECK_EQUAL ( atmosphericCommand->dataTypesIdentifier_, "ALL" );
                BOOST_CHECK_EQUAL ( atmosphericCommand->computationSpecifier_, "BY NRMPOW" );
                for ( unsigned int i = 0; i < expectedCoefficients.size( ); ++i )
                {
                    BOOST_CHECK_EQUAL ( atmosphericCommand->computationCoefficients_.at( i ),
                                        expectedCoefficients.at( i ) );
                }
                BOOST_CHECK_EQUAL ( atmosphericCommand->groundStationsId_, "C60" );
                BOOST_CHECK_EQUAL ( atmosphericCommand->sourceSpecifier_, "" );
                BOOST_CHECK_EQUAL ( atmosphericCommand->sourceId_, 0 );
                BOOST_CHECK_EQUAL ( atmosphericCommand->startTime_, 546858000.00100 );
                BOOST_CHECK_EQUAL ( atmosphericCommand->endTime_, 546879600.000000 );
            }
        }
    }

    // Reading ionosphere file
    {
        std::shared_ptr< CspRawFile > cspFile = std::make_shared< CspRawFile >(
                "/Users/pipas/Documents/messenger-data/ion/mess_rs_2007335_001_dp_ion.txt" );
        std::vector< std::shared_ptr< CspCommand > > extractedCspCommands = cspFile->getCspCommands( );

        std::shared_ptr< AtmosphericCorrectionCspCommand > atmosphericCommand = std::dynamic_pointer_cast< AtmosphericCorrectionCspCommand >(
                    extractedCspCommands.at( 0 ) );

        std::vector< double > expectedCoefficients = { 1.4304,  -0.2192,   3.2505,   3.7560, -11.1495,
                                                       -8.3949,  17.8637,   7.6619,  -8.8524,  -2.3772 };

        BOOST_CHECK_EQUAL ( atmosphericCommand->modelIdentifier_, "CHPART" );
        BOOST_CHECK_EQUAL ( atmosphericCommand->dataTypesIdentifier_, "DOPRNG" );
        BOOST_CHECK_EQUAL ( atmosphericCommand->computationSpecifier_, "BY NRMPOW" );
        for ( unsigned int i = 0; i < expectedCoefficients.size( ); ++i )
        {
            BOOST_CHECK_EQUAL ( atmosphericCommand->computationCoefficients_.at( i ),
                                expectedCoefficients.at( i ) );
        }
        BOOST_CHECK_EQUAL ( atmosphericCommand->groundStationsId_, "C60" );
        BOOST_CHECK_EQUAL ( atmosphericCommand->sourceSpecifier_, "SCID" );
        BOOST_CHECK_EQUAL ( atmosphericCommand->sourceId_, 236 );
        BOOST_CHECK_EQUAL ( atmosphericCommand->startTime_, 249766560.000000 );
        BOOST_CHECK_EQUAL ( atmosphericCommand->endTime_, 249799620.000000 );

    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat