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

#include <boost/date_time/gregorian/gregorian.hpp>

#include "tudat/astro/ground_stations/transmittingFrequencies.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::input_output;
using namespace tudat::simulation_setup;
using namespace tudat;

void runSimulation(
        std::string saveDirectory,
        std::string fileTag,
        bool useInterpolatedEphemerides,
        double epehemeridesTimeStep,
        std::vector< LightTimeCorrectionType > lightTimeCorrectionTypes )
{

//    std::make_shared< OdfRawFileContents >( "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/odf/5327332a.odf" )->writeOdfToTextFile(
//            saveDirectory + "5327332a.txt");
//    std::make_shared< OdfRawFileContents >( "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/odf/5332333a.odf" )->writeOdfToTextFile(
//            saveDirectory + "5332333a.txt");
//    std::make_shared< OdfRawFileContents >( "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/odf/5333334a.odf" )->writeOdfToTextFile(
//            saveDirectory + "5333334a.txt");
//    return;

    // Official mission trajectory kernels V2 (already includes planetary kernels)
//    spice_interface::loadStandardSpiceKernels(
//            { "/Users/pipas/Documents/messenger-spice/msgr_040803_150430_150430_od431sc_2.bsp" } );

    // Official mission trajectory kernels V0
//    spice_interface::loadStandardSpiceKernels( );
//    spice_interface::loadSpiceKernelInTudat(
//            { "/Users/pipas/Documents/messenger-spice/msgr_040803_150430_150430_od431sc_0.bsp" } );

    // GSFC kernels (more accurate)
//    spice_interface::loadStandardSpiceKernels( );
//    spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/messenger-spice/msgr_110318_110714_recon_gsfc_1.bsp" );
//    spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/messenger-spice/msgr_110716_120430_recon_gsfc_1.bsp" );
//    spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/messenger-spice/msgr_120501_130430_recon_gsfc_1.bsp" );
//    spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/messenger-spice/msgr_130501_140430_recon_gsfc_1.bsp" );
//    spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/messenger-spice/msgr_140501_150429_recon_gsfc_1.bsp" );

    // MRO SSL kernels
//    spice_interface::loadStandardSpiceKernels( );
//    spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/mro-spice/mro_psp43_ssd_mro95a.bsp" );

//    spice_interface::loadStandardSpiceKernels(
//            { "/Users/pipas/Documents/mro-spice/de414.bsp",
//              "/Users/pipas/Documents/mro-spice/mar063.bsp",
//              "/Users/pipas/Documents/mro-spice/mro_psp43_ssd_mro95a.bsp"} );

    // MRO NAV kernerls
//    spice_interface::loadStandardSpiceKernels(
//            { "/Users/pipas/Documents/mro-spice/de414.bsp",
//              "/Users/pipas/Documents/mro-spice/mar063.bsp",
//              "/Users/pipas/Documents/mro-spice/mro_psp43.bsp"} );

//    spice_interface::loadStandardSpiceKernels( );
//    spice_interface::loadSpiceKernelInTudat( { "/Users/pipas/Documents/mro-spice/mro_psp43.bsp" } );

    // MGS kernels
//    spice_interface::loadStandardSpiceKernels( {
//        "/Users/pipas/Documents/mro-spice/de414.bsp",
//        "/Users/pipas/Documents/mgs-spice/mar063.bsp",
//        "/Users/pipas/Documents/mgs-spice/mgs_ext5_ipng_mgs95j.bsp",
//        "/Users/pipas/Documents/mgs-spice/mgs_ext6_ipng_mgs95j.bsp" } );
//    spice_interface::loadStandardSpiceKernels( {
//        "/Users/pipas/Documents/mro-spice/de414.bsp",
//        "/Users/pipas/Documents/mgs-spice/mar063.bsp",
//        "/Users/pipas/Documents/mgs-spice/mgs_ext22_ipng_mgs95j.bsp" } );
    spice_interface::loadStandardSpiceKernels( );
    spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/mgs-spice/mgs_ext22_ipng_mgs95j.bsp" );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Mercury", "Venus", "Mars" };

    // Define light-time perturbing bodies
    std::vector< std::string > lightTimePerturbingBodies = bodiesToCreate;

    // Specify initial time
//    Time initialEphemerisTime = Time( 234224663.2 - 10.0 * 86400.0 ); // 2007
//    Time finalEphemerisTime = Time( 234349306.2 + 10.0 * 86400.0 ); // 2007
//    Time initialEphemerisTime = Time( 329140800 - 100.0 * 86400.0 ); // 2010
//    Time finalEphemerisTime = Time( 329400000 + 100.0 * 86400.0 ); // 2010
//    Time initialEphemerisTime = Time( 294455561.18548876 - 10.0 * 86400.0 ); // 2009
//    Time finalEphemerisTime = Time( 294489941.185480893 + 10.0 * 86400.0 ); // 2009
//    Time initialEphemerisTime = Time( 406630870.68285096 - 10.0 * 86400.0 ); // 2011
//    Time finalEphemerisTime = Time( 406717262.6828746 + 10.0 * 86400.0 ); // 2011
//    Time initialEphemerisTime = Time( 544795200 - 5.0 * 86400.0 ); // 2017
//    Time finalEphemerisTime = Time( 544881600 + 5.0 * 86400.0 ); // 2017
//    Time initialEphemerisTime = Time( 70545600 - 5.0 * 86400.0 ); // End of March 2002
//    Time finalEphemerisTime = Time( 70891200 + 1.0 * 86400.0 ); // End of March 2002
    Time initialEphemerisTime = Time( 185976000 - 2.0 * 86400.0 ); // End of November 2005
    Time finalEphemerisTime = Time( 186580800 + 5.0 * 86400.0 ); // End of November 2005
    Time ephemerisTimeStepPlanets = Time( epehemeridesTimeStep );
    Time bufferPlanets = Time( 10.0 * ephemerisTimeStepPlanets );
    Time ephemerisTimeStepSpacecraft = Time( epehemeridesTimeStep );
    Time bufferSpacecraft = Time( 10.0 * ephemerisTimeStepSpacecraft );

    // Create bodies settings needed in simulation
    BodyListSettings bodySettings;
    if ( useInterpolatedEphemerides )
    {
        bodySettings = getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - bufferPlanets, finalEphemerisTime + bufferPlanets,
                    "SSB", "J2000", ephemerisTimeStepPlanets );
    }
    else
    {
        bodySettings = getDefaultBodySettings( bodiesToCreate, "SSB", "J2000" );
    }


    bodySettings.at( "Earth" )->shapeModelSettings = fromSpiceOblateSphericalBodyShapeSettings( );
    bodySettings.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings(
            basic_astrodynamics::iau_2006, "J2000" );
    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    // Create spacecraft
//    std::string spacecraftName = "Messenger";
    std::string spacecraftName = "MGS";
    bodySettings.addSettings( spacecraftName );
    if ( useInterpolatedEphemerides )
    {
        bodySettings.at( spacecraftName )->ephemerisSettings =
                std::make_shared< InterpolatedSpiceEphemerisSettings >(
                        initialEphemerisTime - bufferSpacecraft, finalEphemerisTime + bufferSpacecraft,
                        ephemerisTimeStepSpacecraft, "SSB", "J2000",
                        std::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 ), spacecraftName );
    }
    else
    {
        bodySettings.at( spacecraftName )->ephemerisSettings =
                std::make_shared< DirectSpiceEphemerisSettings >( "SSB", "J2000" );
    }

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );

    // Print initial state
//    Eigen::Matrix<double, 6, 1> state = bodies.at( spacecraftName )->getStateInBaseFrameFromEphemeris( initialEphemerisTime ) -
//                bodies.at( "Mars" )->getStateInBaseFrameFromEphemeris( initialEphemerisTime );
//    std::cout << "State: " << convertCartesianToKeplerianElements(
//            state, bodies.at( "Mars" )->getGravitationalParameter() ).transpose() << std::endl;
//    std::cout << "Orbital period: " << 2 * mathematical_constants::PI * std::sqrt(
//            std::pow( convertCartesianToKeplerianElements( state, bodies.at( "Mars" )->getGravitationalParameter() )( 0 ), 3.0 ) /
//            bodies.at( "Mars" )->getGravitationalParameter() ) << std::endl;


    std::shared_ptr< OdfRawFileContents > rawOdfFileContents =
//            std::make_shared< OdfRawFileContents >( "/Users/pipas/Documents/dsn_trk-2-18/odf07155.dat" );
//            std::make_shared< OdfRawFileContents >( "/Users/pipas/Documents/messenger-rawdata-odf/mess_rs_10162_163_odf.dat" );
//            std::make_shared< OdfRawFileContents >( "/Users/pipas/Documents/messenger-rawdata-odf/mess_rs_09121_121_odf.dat" );
//            std::make_shared< OdfRawFileContents >( "/Users/pipas/Documents/messenger-rawdata-odf/mess_rs_11336_2100_odf.dat" );
//            std::make_shared< OdfRawFileContents >( "/Users/pipas/Documents/mro-rawdata-odf/mromagr2017_097_1335xmmmv1.odf" );
//            std::make_shared< OdfRawFileContents >( "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_0720/odf/2088089a.odf" );
            std::make_shared< OdfRawFileContents >( "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/odf/5332333a.odf" );
//            std::make_shared< OdfRawFileContents >( "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/odf/5327332a.odf" );

//    rawOdfFileContents->writeOdfToTextFile( saveDirectory + "mess_rs_09121_121_odf.txt");

//    std::cout << std::setprecision( 18 );
//    std::cout << "Start time from 1950 UTC: " << rawOdfFileContents->getDataBlocks( ).front()->getCommonDataBlock( )->getObservableTime() << std::endl;
//    std::cout << "End time from 1950 UTC: " << rawOdfFileContents->getDataBlocks( ).back()->getCommonDataBlock( )->getObservableTime() << std::endl;

//    for ( unsigned int i = 0; i < rawOdfFileContents->getDataBlocks( ).size(); ++i )
//    {
//        if ( rawOdfFileContents->getDataBlocks( ).at( i )->getObservableSpecificDataBlock( )->dataType_ == 12 )
//        {
//            std::cout << "Start time from 1950 UTC (12): " << rawOdfFileContents->getDataBlocks( ).at(i)->getCommonDataBlock( )->getObservableTime() << std::endl;
//            break;
//        }
//    }
//    for ( unsigned int i = rawOdfFileContents->getDataBlocks( ).size() - 1; i >= 0 ; --i )
//    {
//        if ( rawOdfFileContents->getDataBlocks( ).at( i )->getObservableSpecificDataBlock( )->dataType_ == 12 )
//        {
//            std::cout << "End time from 1950 UTC (12): " << rawOdfFileContents->getDataBlocks( ).at(i)->getCommonDataBlock( )->getObservableTime() << std::endl;
//            break;
//        }
//    }

//    std::shared_ptr< OdfRawFileContents > rawOdfFileContents2 =
//            std::make_shared< OdfRawFileContents >( "/Users/pipas/Documents/mro-rawdata-odf/mromagr2017_098_1555xmmmv1.odf" );
//    std::shared_ptr< OdfRawFileContents > rawOdfFileContents0 =
//            std::make_shared< OdfRawFileContents >( "/Users/pipas/Documents/mro-rawdata-odf/mromagr2017_096_0820xmmmv1.odf" );
//    std::shared_ptr< OdfRawFileContents > rawOdfFileContentsNeg1 =
//            std::make_shared< OdfRawFileContents >( "/Users/pipas/Documents/mro-rawdata-odf/mromagr2017_095_0825xmmmv1.odf" );

    // Read and process ODF file data
    std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector = { rawOdfFileContents };
    std::shared_ptr< ProcessedOdfFileContents > processedOdfFileContents =
            std::make_shared< ProcessedOdfFileContents >(
                    rawOdfDataVector, bodies.getBody( "Earth" ), true, spacecraftName );

//    std::cout << std::endl << "Spacecraft name: " << processedOdfFileContents->getSpacecraftName( ) << std::endl;
//    std::vector< std::string > groundStationsNames = processedOdfFileContents->getGroundStationsNames( );
//    std::cout << "Ground stations: ";
//    for ( unsigned int i = 0; i < groundStationsNames.size( ); ++i )
//    {
//        std::cout << groundStationsNames.at( i ) << " ";
//    }
//    std::cout << std::endl;
//
//    std::cout << std::setprecision( 18 );
//    std::cout << "Start time: " << processedOdfFileContents->getStartAndEndTime( ).first << std::endl;
//    std::cout << "End time: " << processedOdfFileContents->getStartAndEndTime( ).second << std::endl;
//
//    std::vector< observation_models::ObservableType > observableTypes = processedOdfFileContents->getProcessedObservableTypes( );
//    std::cout << "Observable types: ";
//    for ( unsigned int i = 0; i < observableTypes.size( ); ++i )
//    {
//        std::cout << observableTypes.at( i ) << " ";
//    }
//    std::cout << std::endl;

    // Create observed observation collection
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > observedObservationCollection =
            observation_models::createOdfObservedObservationCollection< long double, Time >(
                    processedOdfFileContents, { dsn_n_way_averaged_doppler } );

//    std::cout << std::endl << "Observation type start and size:" << std::endl;
//    std::map< observation_models::ObservableType, std::pair< int, int > > observationTypeStartAndSize =
//            observedObservationCollection->getObservationTypeStartAndSize( );
//    for ( auto it = observationTypeStartAndSize.begin(); it != observationTypeStartAndSize.end(); ++it )
//    {
//        std::cout << it->first << " " << std::get<0>(it->second) << " " << std::get<1>(it->second) << std::endl;
//    }

//    std::cout << "Observed observables: "<< observedObservationCollection->getObservationVector( ) << std::endl;

    // Create computed observation collection
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;

    std::shared_ptr< input_output::CspRawFile > troposphericCspFile = std::make_shared< input_output::CspRawFile >(
//            "/Users/pipas/Documents/mro-data/tro/mromagr2017_091_2017_121.tro.txt"
//            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_0720/tro/2060091a.tro"
            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/tro/5305337a.tro"
            );
    std::shared_ptr< input_output::CspRawFile > ionosphericCspFile1 = std::make_shared< input_output::CspRawFile >(
//            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_0720/ion/2060121b.ion"
            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/ion/5305335g.ion"
            );
    std::shared_ptr< input_output::CspRawFile > ionosphericCspFile2 = std::make_shared< input_output::CspRawFile >(
//            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_0720/ion/2060121b.ion"
            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/ion/5335001a.ion"
            );

    std::map< int, std::string > spacecraftNamePerSpacecraftId;
    spacecraftNamePerSpacecraftId[ 74 ] = "MRO";
    spacecraftNamePerSpacecraftId[ 94 ] = "MGS";

    std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > > lightTimeCorrectionSettings;

    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), first_order_relativistic ) )
    {
        lightTimeCorrectionSettings.push_back(
                std::make_shared< observation_models::FirstOrderRelativisticLightTimeCorrectionSettings >(
                    lightTimePerturbingBodies ) );
    }
    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), tabulated_tropospheric ) )
    {
        lightTimeCorrectionSettings.push_back(
                std::make_shared< observation_models::TabulatedTroposphericCorrectionSettings >(
                        createTroposphericDryCorrectionAdjustment( { troposphericCspFile } ),
                        createTroposphericDryCorrectionAdjustment( { troposphericCspFile } ) ) );
    }
    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), tabulated_ionospheric ) )
    {
        lightTimeCorrectionSettings.push_back(
                std::make_shared< observation_models::TabulatedIonosphericCorrectionSettings >(
                      createIonosphericCorrection( { ionosphericCspFile1, ionosphericCspFile2 },
                                                   spacecraftNamePerSpacecraftId ) ) );
    }
    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), saastamoinen_tropospheric ) )
    {
        input_output::setDsnWeatherDataInGroundStations(
            bodies,
            std::vector< std::string >
                    { "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/wea/50013321.wea",
                      "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/wea/50013324.wea",
                      "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/wea/50013326.wea" } );

        lightTimeCorrectionSettings.push_back( std::make_shared< SaastamoinenTroposphericCorrectionSettings >( ) );
    }
    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), jakowski_vtec_ionospheric ) )
    {
        lightTimeCorrectionSettings.push_back( std::make_shared< JakowskiIonosphericCorrectionSettings >( ) );
    }

    std::map < observation_models::ObservableType, std::vector< observation_models::LinkEnds > > linkEndsPerObservable =
            observedObservationCollection->getLinkEndsPerObservableType( );
    std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
            std::make_shared< LightTimeConvergenceCriteria >( true );
    for ( auto it = linkEndsPerObservable.begin(); it != linkEndsPerObservable.end(); ++it )
    {
        for ( unsigned int i = 0; i < it->second.size(); ++i )
        {
//            observationModelSettingsList.push_back(
//                    std::make_shared< observation_models::ObservationModelSettings >(
//                            it->first, it->second.at( i ), lightTimeCorrectionSettings, nullptr, nullptr ) );
            if ( it->first == observation_models::dsn_n_way_averaged_doppler )
            {
                observationModelSettingsList.push_back(
                    std::make_shared< observation_models::DsnNWayAveragedDopplerObservationSettings >(
                            it->second.at( i ), lightTimeCorrectionSettings, nullptr,
                            lightTimeConvergenceCriteria ) );
//                observationModelSettingsList.push_back(
//                    std::make_shared< observation_models::NWayDifferencedRangeObservationSettings >(
//                            it->second.at( i ), lightTimeCorrectionSettings, nullptr,
//                            lightTimeConvergenceCriteria ) );
//                observationModelSettingsList.push_back(
//                    std::make_shared< observation_models::NWayRangeObservationSettings >(
//                            it->second.at( i ),
//                            nullptr, 3, nullptr,
//                            lightTimeConvergenceCriteria ) );
            }
        }
    }

    observation_models::setOdfInformationInBodies( processedOdfFileContents, bodies );
    
    std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< long double, Time > > >
            observationSimulators = observation_models::createObservationSimulators< long double, Time >(
                    observationModelSettingsList, bodies );

    std::vector< std::shared_ptr< ObservationSimulationSettings< Time > > > observationSimulationSettings =
            observation_models::createOdfObservationSimulationSettingsList< long double, Time >(
                    observedObservationCollection );


    std::shared_ptr< observation_models::ObservationCollection< long double, Time > >
            simulatedObservationCollection = simulation_setup::simulateObservations< long double, Time >(
                    observationSimulationSettings, observationSimulators, bodies );

//    std::cout << std::endl << "Observation type start and size:" << std::endl;
//    observationTypeStartAndSize = simulatedObservationCollection->getObservationTypeStartAndSize( );
//    for ( auto it = observationTypeStartAndSize.begin(); it != observationTypeStartAndSize.end(); ++it )
//    {
//        std::cout << it->first << " " << std::get<0>(it->second) << " " << std::get<1>(it->second) << std::endl;
//    }

//    std::cout << simulatedObservationCollection->getObservationVector( ) << std::endl;

//    if ( fileTag == "2007" || fileTag == "2017_nav" )
//    {
//        std::vector< std::string > stations;
//
//        if ( fileTag == "2007" )
//            stations = { "DSS-63", "DSS-14", "DSS-43" };
//        else
//            stations = { "DSS-55", "DSS-14", "DSS-35" };
//
//        for ( unsigned int i = 0; i < 3; ++i )
//        {
//            std::shared_ptr< ground_stations::PiecewiseLinearFrequencyInterpolator > rampInterpolator =
//                    std::dynamic_pointer_cast< ground_stations::PiecewiseLinearFrequencyInterpolator >(
//                bodies.at("Earth")->getGroundStation( stations.at( i ) )->getTransmittingFrequencyCalculator( ) );
//
//            std::vector< double > startTimes = rampInterpolator->getStartTimes( );
//            std::vector< double > endTimes = rampInterpolator->getEndTimes( );
//            std::vector< double > rampRates = rampInterpolator->getRampRates( );
//            std::vector< double > startFrequency = rampInterpolator->getStartFrequencies( );
//
//            Eigen::MatrixXd rampData ( startTimes.size( ), 4 );
//            for ( unsigned int j = 0; j < startTimes.size( ); ++j )
//            {
//                rampData( j, 0 ) = startTimes.at( j );
//                rampData( j, 1 ) = endTimes.at( j );
//                rampData( j, 2 ) = rampRates.at( j );
//                rampData( j, 3 ) = startFrequency.at( j );
//            }
//
//            std::ofstream file("/Users/pipas/tudatpy-testing/rampData_" + fileTag + stations.at( i ) + ".txt");
//            file << std::setprecision( 15 ) << rampData;
//            file.close();
//        }
//    }

    Eigen::Matrix< long double, Eigen::Dynamic, 1 > simulatedObservations =
            simulatedObservationCollection->getObservationVector( );

    Eigen::Matrix< long double, Eigen::Dynamic, 4 > observations;
    observations.resize( simulatedObservationCollection->getObservationVector( ).size( ), 4 );

    observations.col( 1 ) = simulatedObservationCollection->getObservationVector( );
    observations.col( 3 ) = observedObservationCollection->getObservationVector( );

    for ( unsigned int i = 0; i < simulatedObservationCollection->getObservationVector( ).size( ); ++i )
    {
        observations( i, 0 ) = static_cast< Time >( simulatedObservationCollection->getConcatenatedTimeVector( ).at( i ) ).getSeconds< long double >();
        observations( i, 2 ) = static_cast< Time >( observedObservationCollection->getConcatenatedTimeVector( ).at( i ) ).getSeconds< long double >();
    }

    std::ofstream file(saveDirectory + "observations_" + fileTag + ".txt");
    file << std::setprecision( 15 ) << observations;
    file.close();

    std::ofstream file2(saveDirectory + "observationsStartAndSize_" + fileTag + ".txt");
    std::map< ObservableType, std::map< int, std::vector< std::pair< int, int > > > > observationSetStartAndSize =
            observedObservationCollection->getObservationSetStartAndSizePerLinkEndIndex();
    for ( auto it = observationSetStartAndSize.begin(); it != observationSetStartAndSize.end(); ++it )
    {
        ObservableType observable = it->first;
        for ( auto it2 = it->second.begin(); it2 != it->second.end(); ++it2 )
        {
            int linkEnd = it2->first;
            for ( unsigned int i = 0; i < it2->second.size(); ++i )
            {
                file2 << std::setprecision( 15 ) << observable << " " << linkEnd << " " << it2->second.at( i ).first <<
                    " " << it2->second.at( i ).second << std::endl;
            }
        }
    }
    file2.close();

}

BOOST_AUTO_TEST_SUITE( test_dsn_odf_observation_models )

BOOST_AUTO_TEST_CASE( testDsnNWayAveragedDopplerModel )
{
    std::string saveDirectory = "/Users/pipas/tudatpy-testing/mgs/mors_2190/";

    int testCase = 3;

//    std::string fileTag = "2007";
//    std::string fileTag = "2009";
//    std::string fileTag = "2011";
//    std::string fileTag = "2017_ssd";
//    std::string fileTag = "2017_096_nav";


    // Default
    if ( testCase == 0 )
    {
        std::string fileTag = "5332333aOdf";
        double ephemeridesTimeStep = 100.0;
        runSimulation( saveDirectory, fileTag + "_troCorr", true, ephemeridesTimeStep,
                       { tabulated_tropospheric } );
    }
    else if ( testCase == 2 )
    {
        std::string fileTag = "5332333aOdf_iau2000b";
        runSimulation( saveDirectory, fileTag + "_noCorr", false, TUDAT_NAN, { } );
    }
    else if ( testCase == 3 )
    {
//        std::string fileTag = "5332333aOdf_interpState50";
        std::string fileTag = "5332333aOdf_interpState50";
        double ephemeridesTimeStep = 50.0;

        runSimulation( saveDirectory, fileTag + "_noCorr", true, ephemeridesTimeStep, { } );

        runSimulation( saveDirectory, fileTag + "_relCorr", true, ephemeridesTimeStep, { first_order_relativistic } );

        runSimulation( saveDirectory, fileTag + "_troCorr", true, ephemeridesTimeStep, { tabulated_tropospheric } );

//        runSimulation( saveDirectory, fileTag + "_ionCorr", true, ephemeridesTimeStep, { tabulated_ionospheric } );

//        runSimulation( saveDirectory, fileTag + "_troGdCorr", true, ephemeridesTimeStep, { saastamoinen_tropospheric } );

//        runSimulation( saveDirectory, fileTag + "_ionGdCorr", true, ephemeridesTimeStep, { jakowski_vtec_ionospheric } );
    }
    else if ( testCase == 1 )
    {
        // Test of SPICE interpolation step size
        std::string fileTag = "5332333aOdf";
        for ( unsigned int i = 0; i < 2; ++i )
        {
            std::vector< LightTimeCorrectionType > lightTimeCorrections;
            std::string typeTag;
            if ( i == 0 )
            {
                lightTimeCorrections = { };
                typeTag = "noCorr";
            }
            else
            {
                lightTimeCorrections = { first_order_relativistic };
                typeTag = "relCorr";
            }

            runSimulation( saveDirectory + "data_spice_test/", fileTag + "_" + typeTag + "_SpiceDirect", false,
                           TUDAT_NAN,
                           lightTimeCorrections );

            std::vector< int > ephemeridesTimeSteps = { 1, 2, 5, 10, 20, 50, 100, 200, 300, 400 };

            for ( int step: ephemeridesTimeSteps )
            {
                runSimulation(
                        saveDirectory + "data_spice_test/",
                        fileTag + "_" + typeTag + "_SpiceInterp" + std::to_string( step ),
                        true, step, lightTimeCorrections );
            }
        }
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}