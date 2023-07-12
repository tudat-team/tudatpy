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
        std::string spacecraftName,
        std::vector< std::string > odfFiles,
        std::vector< std::string > troposphericCorrectionFileNames,
        std::vector< std::string > ionosphericCorrectionFileNames,
        std::vector< std::string > weatherFileNames,
        std::string saveDirectory,
        std::string fileTag,
        Time initialEphemerisTime,
        Time finalEphemerisTime,
        bool useInterpolatedEphemerides,
        double epehemeridesTimeStep,
        std::vector< LightTimeCorrectionType > lightTimeCorrectionTypes,
        std::string spacecraftEphemeridesOrigin,
        std::string planetaryEphemeridesOrigin,
        bool writeOdfToTxtFile,
        std::pair< Time, Time > startAndEndTimesToProcess = std::make_pair< Time, Time >( TUDAT_NAN, TUDAT_NAN ) )
{

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Mercury", "Venus", "Mars" };

    // Define light-time perturbing bodies
    std::vector< std::string > lightTimePerturbingBodies = bodiesToCreate;

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
                    planetaryEphemeridesOrigin, "J2000", ephemerisTimeStepPlanets );
    }
    else
    {
        bodySettings = getDefaultBodySettings( bodiesToCreate, planetaryEphemeridesOrigin, "J2000" );
    }


    bodySettings.at( "Earth" )->shapeModelSettings = fromSpiceOblateSphericalBodyShapeSettings( );
    bodySettings.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings(
            basic_astrodynamics::iau_2006, "J2000" );
    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    // Create spacecraft
    bodySettings.addSettings( spacecraftName );
    if ( useInterpolatedEphemerides )
    {
        bodySettings.at( spacecraftName )->ephemerisSettings =
                std::make_shared< InterpolatedSpiceEphemerisSettings >(
                        initialEphemerisTime - bufferSpacecraft, finalEphemerisTime + bufferSpacecraft,
                        ephemerisTimeStepSpacecraft, spacecraftEphemeridesOrigin, "J2000",
                        std::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 ), spacecraftName );
    }
    else
    {
        bodySettings.at( spacecraftName )->ephemerisSettings =
                std::make_shared< DirectSpiceEphemerisSettings >( spacecraftEphemeridesOrigin, "J2000" );
    }

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );

    // Read and process ODF file data
    std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector;
    for ( std::string odfFile : odfFiles )
        rawOdfDataVector.push_back( std::make_shared< OdfRawFileContents >( odfFile ) );

    if ( writeOdfToTxtFile )
        for ( std::shared_ptr< input_output::OdfRawFileContents > rawOdfData : rawOdfDataVector )
            rawOdfData->writeOdfToTextFile( saveDirectory + rawOdfData->fileName_.substr( rawOdfData->fileName_.find_last_of( "/" ) + 1 ) + ".txt" );

    std::shared_ptr< ProcessedOdfFileContents > processedOdfFileContents =
            std::make_shared< ProcessedOdfFileContents >(
                    rawOdfDataVector, spacecraftName, true );

    observation_models::setOdfInformationInBodies( processedOdfFileContents, bodies );

    // Create observed observation collection
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > observedObservationCollection =
            observation_models::createOdfObservedObservationCollection< long double, Time >(
                    processedOdfFileContents, { dsn_n_way_averaged_doppler },
                    startAndEndTimesToProcess );

    // Create computed observation collection
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;

    std::map< int, std::string > spacecraftNamePerSpacecraftId;
    spacecraftNamePerSpacecraftId[ 74 ] = "MRO";
    spacecraftNamePerSpacecraftId[ 94 ] = "MGS";

    std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > > lightTimeCorrectionSettings;

    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), first_order_relativistic ) )
    {
        lightTimeCorrectionSettings.push_back( firstOrderRelativisticLightTimeCorrectionSettings( lightTimePerturbingBodies ) );
    }
    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), tabulated_tropospheric ) )
    {
        lightTimeCorrectionSettings.push_back( tabulatedTroposphericCorrectionSettings( troposphericCorrectionFileNames ) );
    }
    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), tabulated_ionospheric ) )
    {
        lightTimeCorrectionSettings.push_back( tabulatedIonosphericCorrectionSettings( ionosphericCorrectionFileNames, spacecraftNamePerSpacecraftId ) );
    }
    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), saastamoinen_tropospheric ) )
    {
        input_output::setDsnWeatherDataInGroundStations(bodies, weatherFileNames  );

        lightTimeCorrectionSettings.push_back( saastamoinenTroposphericCorrectionSettings( ) );
    }
    if ( std::count( lightTimeCorrectionTypes.begin(), lightTimeCorrectionTypes.end(), jakowski_vtec_ionospheric ) )
    {
        lightTimeCorrectionSettings.push_back( jakowskiIonosphericCorrectionSettings( ) );
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
                    observation_models::dsnNWayAveragedDopplerObservationSettings(
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
    
    std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< long double, Time > > >
            observationSimulators = observation_models::createObservationSimulators< long double, Time >(
                    observationModelSettingsList, bodies );

    std::vector< std::shared_ptr< ObservationSimulationSettings< Time > > > observationSimulationSettings =
            observation_models::createOdfObservationSimulationSettingsList< long double, Time >(
                    observedObservationCollection );


    std::shared_ptr< observation_models::ObservationCollection< long double, Time > >
            simulatedObservationCollection = simulation_setup::simulateObservations< long double, Time >(
                    observationSimulationSettings, observationSimulators, bodies );

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

    // MESSENGER: Official mission trajectory kernels V2 (already includes planetary kernels)
//    spice_interface::loadStandardSpiceKernels(
//            { "/Users/pipas/Documents/messenger-spice/msgr_040803_150430_150430_od431sc_2.bsp" } );

    // MESSENGER: GSFC kernels (more accurate)
//    spice_interface::loadStandardSpiceKernels( );
//    spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/messenger-spice/msgr_110318_110714_recon_gsfc_1.bsp" );
//    spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/messenger-spice/msgr_110716_120430_recon_gsfc_1.bsp" );
//    spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/messenger-spice/msgr_120501_130430_recon_gsfc_1.bsp" );
//    spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/messenger-spice/msgr_130501_140430_recon_gsfc_1.bsp" );
//    spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/messenger-spice/msgr_140501_150429_recon_gsfc_1.bsp" );

    // MRO SSD kernels
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

    std::vector< std::string > odfFiles = { "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/odf/5332333a.odf" };
//    std::vector< std::string > odfFiles = { "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/odf/5327332a.odf" };

    std::vector< std::string > tropCorrectionFiles = { "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/tro/5305337a.tro" };
    std::vector< std::string > ionCorrectionFiles = {
            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/ion/5305335g.ion",
            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/ion/5335001a.ion" };
    std::vector< std::string > weatherFiles = {
            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/wea/50013321.wea",
            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/wea/50013324.wea",
            "/Users/pipas/Documents/mgs-m-rss-1-ext-v1/mors_2190/wea/50013326.wea" };

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
    Time initialEphemerisTime = Time( 185976000 - 2.0 * 86400.0 ); // End of November 2005
    Time finalEphemerisTime = Time( 186580800 + 5.0 * 86400.0 ); // End of November 2005

    std::string spacecraftName = "MGS";

    int testCase = 3;


    if ( testCase == 0 )
    {
        std::string saveDirectory = "/Users/pipas/tudatpy-testing/mgs/mors_2190/";
        std::string fileTag = "5332333aOdf";
        std::string ephemeridesOrigin = "SSB";
        double ephemeridesTimeStep = 100.0;
        runSimulation( spacecraftName, odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_troCorr",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep,
                       { tabulated_tropospheric }, ephemeridesOrigin, ephemeridesOrigin, false );
    }
    else if ( testCase == 2 )
    {
        std::string saveDirectory = "/Users/pipas/tudatpy-testing/mgs/mors_2190/";
        std::string fileTag = "5332333aOdf_iau2000b";
        std::string ephemeridesOrigin = "SSB";
        runSimulation( spacecraftName, odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_noCorr",
                       initialEphemerisTime, finalEphemerisTime, false, TUDAT_NAN, { },
                       ephemeridesOrigin, ephemeridesOrigin, false );
    }
    else if ( testCase == 3 )
    {
        std::string saveDirectory = "/Users/pipas/tudatpy-testing/mgs/mors_2190/";
//        std::string fileTag = "5332333aOdf_interpState50";
        std::string ephemeridesOrigin = "SSB";
        std::string fileTag = "5332333aOdf_interpState50";
        double ephemeridesTimeStep = 50.0;

        runSimulation( spacecraftName, odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_noCorr",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep, { },
                       ephemeridesOrigin, ephemeridesOrigin, false );

//        runSimulation( spacecraftName, odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
//                       saveDirectory, fileTag + "_relCorr",
//                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep, { first_order_relativistic },
//                       ephemeridesOrigin, ephemeridesOrigin, false );
//
//        runSimulation( spacecraftName, odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
//                       saveDirectory, fileTag + "_troCorr",
//                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep, { tabulated_tropospheric },
//                       ephemeridesOrigin, ephemeridesOrigin, false );
//
//        runSimulation( spacecraftName, odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
//                       saveDirectory, fileTag + "_ionCorr",
//                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep, { tabulated_ionospheric },
//                       ephemeridesOrigin, ephemeridesOrigin, false );

//        runSimulation( spacecraftName, odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
//                       saveDirectory, fileTag + "_troGdCorr",
//                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep, { saastamoinen_tropospheric },
//                       ephemeridesOrigin, ephemeridesOrigin, false );
//
//        runSimulation( spacecraftName, odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
//                       saveDirectory, fileTag + "_ionGdCorr",
//                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep, { jakowski_vtec_ionospheric },
//                       ephemeridesOrigin, ephemeridesOrigin, false );
    }
    // Effect of ephemerides origin on residuals
    else if ( testCase == 4 )
    {
        std::string saveDirectory = "/Users/pipas/tudatpy-testing/mgs/mors_2190/data_origin_test/";
        std::string fileTag = "5332333aOdf_interpState50";
        double ephemeridesTimeStep = 50.0;

        runSimulation( spacecraftName, odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_noCorr_SsbPlanet_SsbSpacecraft",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep, { },
                       "SSB", "SSB", false );

        runSimulation( spacecraftName, odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_noCorr_SsbPlanet_MarsSpacecraft",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep, { },
                       "Mars", "SSB", false );

        runSimulation( spacecraftName, odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_noCorr_MarsPlanet_SsbSpacecraft",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep, { },
                       "SSB", "Mars", false );

        runSimulation( spacecraftName, odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_noCorr_MarsPlanet_MarsSpacecraft",
                       initialEphemerisTime, finalEphemerisTime, true, ephemeridesTimeStep, { },
                       "Mars", "Mars", false );

    }
    // Validation against Verma (2022), MGS
    else if ( testCase == 5 )
    {
        std::string saveDirectory = "/Users/pipas/tudatpy-testing/mgs/mors_0401/";
        std::string ephemeridesOrigin = "SSB";
        std::string fileTag = "mors0401";
        double ephemeridesTimeStep = 50.0;
        odfFiles = {
                "/Users/pipas/Documents/mgs-m-rss-1-map-v1/mors_0401/odf/9068068a.odf",
//                "/Users/pipas/Documents/mgs-m-rss-1-map-v1/mors_0401/odf/9068068b.odf",
                "/Users/pipas/Documents/mgs-m-rss-1-map-v1/mors_0401/odf/9068071a.odf"
//                "/Users/pipas/Documents/mgs-m-rss-1-map-v1/mors_0401/odf/9071074a.odf",
//                "/Users/pipas/Documents/mgs-m-rss-1-map-v1/mors_0401/odf/9072073a.odf"
        };
        initialEphemerisTime = Time( -25754400 );
        finalEphemerisTime = Time( -25272000 + 2.0 * 86400.0 );

        spice_interface::clearSpiceKernels( );
        spice_interface::loadStandardSpiceKernels( );
        spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/planet-spice/de438.bsp" );
        spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/mgs-spice/mgs_map1_ipng_mgs95j.bsp" );

        runSimulation( spacecraftName, odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_noCorr",
                       initialEphemerisTime, finalEphemerisTime, false, ephemeridesTimeStep, { },
                       ephemeridesOrigin, ephemeridesOrigin, true );
    }
    else if ( testCase == 6 )
    {
        std::string saveDirectory = "/Users/pipas/tudatpy-testing/magellan/mg_2601/";
        std::string ephemeridesOrigin = "SSB";
        std::string fileTag = "mg2601";
        double ephemeridesTimeStep = 50.0;
        spacecraftName = "Magellan";
        odfFiles = {
                "/Users/pipas/Documents/mgn-v-rss-1-tracking-v1/mg_2601/odf/3004003a.odf"
        };
        initialEphemerisTime = Time( TUDAT_NAN );
        finalEphemerisTime = Time( TUDAT_NAN );
        std::pair< Time, Time > startAndEndTimesToProcess = std::make_pair( Time( -212932800 ), Time( -212328000 ) );

        spice_interface::clearSpiceKernels( );
        spice_interface::loadStandardSpiceKernels( );
        spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/planet-spice/de438.bsp" );
        spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/mgn-spice/mgn_cycle4_grav.bsp" );

        runSimulation( spacecraftName, odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                       saveDirectory, fileTag + "_noCorr",
                       initialEphemerisTime, finalEphemerisTime, false, ephemeridesTimeStep, { },
                       ephemeridesOrigin, ephemeridesOrigin, true,
                       startAndEndTimesToProcess );
    }
    // Test of SPICE interpolation step size
    else if ( testCase == 1 )
    {
        std::string saveDirectory = "/Users/pipas/tudatpy-testing/mgs/mors_2190/";
        std::string fileTag = "5332333aOdf";
        std::string ephemeridesOrigin = "SSB";
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

            runSimulation( spacecraftName, odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                           saveDirectory + "data_spice_test/", fileTag + "_" + typeTag + "_SpiceDirect",
                           initialEphemerisTime, finalEphemerisTime, false, TUDAT_NAN,
                           lightTimeCorrections, ephemeridesOrigin, ephemeridesOrigin, false );

            std::vector< int > ephemeridesTimeSteps = { 1, 2, 5, 10, 20, 50, 100, 200, 300, 400 };

            for ( int step: ephemeridesTimeSteps )
            {
                runSimulation(
                        spacecraftName, odfFiles, tropCorrectionFiles, ionCorrectionFiles, weatherFiles,
                        saveDirectory + "data_spice_test/",
                        fileTag + "_" + typeTag + "_SpiceInterp" + std::to_string( step ),
                        initialEphemerisTime, finalEphemerisTime, true, step, lightTimeCorrections,
                        ephemeridesOrigin, ephemeridesOrigin, false );
            }
        }
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}