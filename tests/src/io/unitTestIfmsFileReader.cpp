/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      If tabs are used as spaces, it doesn't work. The seperator should also be tabs then.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <iostream>
#include <utility>
#include "tudat/basics/testMacros.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/estimation_setup/observations.h"

#include "tudat/io/readTrackingTxtFile.h"
#include "tudat/simulation/estimation_setup/processTrackingTxtFile.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"

// Some simplifications for shorter lines
using namespace tudat::input_output;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat;

namespace tudat
{
namespace unit_tests
{



//! Starting the entire test suite
BOOST_AUTO_TEST_SUITE(test_ifms_file_reader);

//! Test reading of mars express IFMS files.
BOOST_AUTO_TEST_CASE(testIfmsFileReader)
{

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create bodies.

    simulation_setup::BodyListSettings bodySettings =
        simulation_setup::getDefaultBodySettings( { "Earth"}, "SSB", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->groundStationSettings.push_back( std::make_shared< GroundStationSettings >(
        "NWNORCIA", getCombinedApproximateGroundStationPositions( ).at( "NWNORCIA" ) ) );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "MeX" );

    std::vector< std::shared_ptr< TrackingTxtFileContents > > rawIfmsFiles;

    rawIfmsFiles.push_back( readIfmsFile( "/home/dominic/Downloads/M32ICL2L02_D2X_133630120_00.TAB.txt" ) );
    rawIfmsFiles.push_back( readIfmsFile( "/home/dominic/Downloads/M32ICL2L02_D2X_133630203_00.TAB.txt" ) );
    rawIfmsFiles.push_back( readIfmsFile( "/home/dominic/Downloads/M32ICL2L02_D2X_133631902_00.TAB.txt" ) );
    rawIfmsFiles.push_back( readIfmsFile( "/home/dominic/Downloads/M32ICL2L02_D2X_133632221_00.TAB.txt" ) );
    rawIfmsFiles.push_back( readIfmsFile( "/home/dominic/Downloads/M32ICL2L02_D2X_133632301_00.TAB.txt" ) );

    for( unsigned int i = 0; i < rawIfmsFiles.size( ) + 1; i++ )
    {
        std::vector< std::shared_ptr< ProcessedTrackingTxtFileContents< > > > processedIfmsFiles;
        std::cout<<i<<std::endl;
        if( i < rawIfmsFiles.size( ) )
        {
            rawIfmsFiles.at( i )->addMetaData( TrackingDataType::receiving_station_name, "NWNORCIA" );
            rawIfmsFiles.at( i )->addMetaData( TrackingDataType::transmitting_station_name, "NWNORCIA" );
            processedIfmsFiles.push_back( std::make_shared<observation_models::ProcessedTrackingTxtFileContents<> >(
                rawIfmsFiles.at( i ), "MeX", simulation_setup::getCombinedApproximateGroundStationPositions( )));
        }
        else
        {
            for( unsigned int j = 0; j < rawIfmsFiles.size( ); j++ )
            {
                processedIfmsFiles.push_back( std::make_shared<observation_models::ProcessedTrackingTxtFileContents<> >(
                    rawIfmsFiles.at( j ), "MeX", simulation_setup::getCombinedApproximateGroundStationPositions( )));
            }
        }



        setTrackingDataInformationInBodies( processedIfmsFiles, bodies, { dsn_n_way_averaged_doppler } );

        auto observationCollection = observation_models::createTrackingTxtFilesObservationCollection<double, double>(
            processedIfmsFiles, { dsn_n_way_averaged_doppler } );


        Eigen::VectorXd observationTimes =
            utilities::convertStlVectorToEigenVector( observationCollection->getConcatenatedTimeVector( ));
        Eigen::VectorXd observationValues = observationCollection->getObservationVector( );

        std::cout << "TEST: " << observationTimes( 0 ) << " " << observationValues( 1 ) << std::endl;
//    input_output::writeMatrixToFile( observationValues, "ifms_doppler.dat", 16 );
//    input_output::writeMatrixToFile( observationTimes, "ifms_times.dat", 16 );

    }

}

// End test suite
BOOST_AUTO_TEST_SUITE_END();

}// namespace unit_tests

}// namespace tudat
