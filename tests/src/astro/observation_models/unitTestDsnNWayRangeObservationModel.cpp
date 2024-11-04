/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"
#include "tudat/simulation/estimation_setup.h"

#include "tudat/io/readOdfFile.h"
#include "tudat/io/readTabulatedMediaCorrections.h"
#include "tudat/io/readTabulatedWeatherData.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include <boost/date_time/gregorian/gregorian.hpp>

#include "tudat/astro/ground_stations/transmittingFrequencies.h"

//namespace tudat
//{
//namespace unit_tests
//{

using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::input_output;
using namespace tudat::simulation_setup;
using namespace tudat::unit_tests;
using namespace tudat;


//BOOST_AUTO_TEST_SUITE( test_dsn_n_way_range_observation_model )
//
//BOOST_AUTO_TEST_CASE( testDsnNWayRangeModel )
int main( )
{

    spice_interface::loadStandardSpiceKernels( );
    // Verma (2022) uses DE438, but here we use the standard Tudat SPICE kernels as the difference produced by the kernels
    // is way below the level of the current residuals.
    // spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/planet-spice/de438.bsp" );
    spice_interface::loadSpiceKernelInTudat( tudat::paths::getTudatTestDataPath( ) + "dsn_n_way_doppler_observation_model/mgs_map1_ipng_mgs95j.bsp" );

    std::vector< double > testResiduals;

    int counter = 0;
    std::string spacecraftName = "MGS";
    std::string ephemeridesOrigin = "SSB";

    std::vector< std::string > odfFiles = {
            tudat::paths::getTudatTestDataPath( ) + "dsn_n_way_doppler_observation_model/9068068a.odf",
            tudat::paths::getTudatTestDataPath( ) + "dsn_n_way_doppler_observation_model/9068071a.odf" };

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Mars" };

    // Create bodies settings needed in simulation
    BodyListSettings bodySettings;
    bodySettings = getDefaultBodySettings( bodiesToCreate, ephemeridesOrigin, "J2000" );


    bodySettings.at( "Earth" )->shapeModelSettings = fromSpiceOblateSphericalBodyShapeSettings( );
    bodySettings.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings(
            basic_astrodynamics::iau_2006, "J2000" );
    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    // Create spacecraft
    bodySettings.addSettings( spacecraftName );
    bodySettings.at( spacecraftName )->ephemerisSettings =
            std::make_shared< DirectSpiceEphemerisSettings >( ephemeridesOrigin, "J2000" );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );

    // Read and process ODF file data
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > observedObservationCollection;


    std::vector<std::shared_ptr<input_output::OdfRawFileContents> > rawOdfDataVector;
    for ( std::string odfFile: odfFiles )
        rawOdfDataVector.push_back( std::make_shared<OdfRawFileContents>( odfFile ));

    std::shared_ptr<ProcessedOdfFileContents<Time> > processedOdfFileContents =
        std::make_shared<ProcessedOdfFileContents<Time> >( rawOdfDataVector, spacecraftName );

    // Create observed observation collection
    observedObservationCollection =
        observation_models::createOdfObservedObservationCollection<long double, Time>(
            processedOdfFileContents, { dsn_n_way_range } );

    observation_models::setOdfInformationInBodies<Time>( processedOdfFileContents, bodies );

    // Create computed observation collection
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;

    std::map < observation_models::ObservableType, std::vector< observation_models::LinkEnds > > linkEndsPerObservable =
            observedObservationCollection->getLinkEndsPerObservableType( );
    for ( auto it = linkEndsPerObservable.begin(); it != linkEndsPerObservable.end(); ++it )
    {
        for ( unsigned int i = 0; i < it->second.size(); ++i )
        {
            if ( it->first == observation_models::dsn_n_way_range )
            {
                const std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections =
                    std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( std::vector< std::string >( { "Sun" } ) );
                observationModelSettingsList.push_back(
                    std::make_shared< observation_models::ObservationModelSettings >( dsn_n_way_range, it->second.at( i ), lightTimeCorrections ) );
            }
        }
    }

    std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< long double, Time > > >
            observationSimulators = observation_models::createObservationSimulators< long double, Time >(
                    observationModelSettingsList, bodies );

    std::vector< std::shared_ptr< ObservationSimulationSettings< Time > > > observationSimulationSettings =
            getObservationSimulationSettingsFromObservations< long double, Time >( observedObservationCollection, bodies );

    std::cout<<"Pre-simulation "<<observationSimulationSettings.size( )<<std::endl;

    std::shared_ptr< observation_models::ObservationCollection< long double, Time > >
            simulatedObservationCollection = simulation_setup::simulateObservations< long double, Time >(
                    observationSimulationSettings, observationSimulators, bodies );
    std::cout<<"Post-simulation"<<std::endl;

    std::cout<<simulatedObservationCollection->getConcatenatedObservations( ).transpose( ) -
    observedObservationCollection->getConcatenatedObservations( ).transpose( )<<std::endl;

    LinkEnds dss45MgsLinkEnds;
    dss45MgsLinkEnds[ transmitter ] = LinkEndId( "Earth", "DSS-45" );
    dss45MgsLinkEnds[ retransmitter ] = LinkEndId( "MGS" );
    dss45MgsLinkEnds[ receiver ] = LinkEndId( "Earth", "DSS-45" );

     std::vector< std::shared_ptr< SingleObservationSet< long double, Time > > > singleLinkSimulatedObservations =
             simulatedObservationCollection->getSingleLinkAndTypeObservationSets( dsn_n_way_range, dss45MgsLinkEnds );
     std::vector< std::shared_ptr< SingleObservationSet< long double, Time > > > singleLinkObservedObservations =
             observedObservationCollection->getSingleLinkAndTypeObservationSets( dsn_n_way_range, dss45MgsLinkEnds );

}
//
//BOOST_AUTO_TEST_SUITE_END( )
//
//}
//
//}