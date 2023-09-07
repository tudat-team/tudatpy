/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


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

using namespace tudat::propagators;
using namespace tudat::estimatable_parameters;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::input_output;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::numerical_integrators;
using namespace tudat::basic_astrodynamics;

using namespace tudat;

int main( )
{
    double initialTimeEnvironment = Time(107561, 2262.19) - 2.0 * 3600.0;
    double finalTimeEnvironment = Time(108258, 2771.19) + 2.0 * 3600.0;

    // Load spice kernels
    spice_interface::loadStandardSpiceKernels( );
    spice_interface::loadSpiceKernelInTudat( "/home/dominic/Tudat/Data/GRAIL_Spice/grail_120301_120529_sci_v02.bsp" );

    // Create settings for default bodies
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Moon" };
    std::string globalFrameOrigin = "Earth";
    std::string globalFrameOrientation = "J2000";
    BodyListSettings bodySettings = getDefaultBodySettings(
        bodiesToCreate, initialTimeEnvironment, finalTimeEnvironment, globalFrameOrigin, globalFrameOrientation, 120.0 );

    // Add high-accuracy Earth settings
    bodySettings.at( "Earth" )->shapeModelSettings = fromSpiceOblateSphericalBodyShapeSettings( );
    bodySettings.at( "Earth" )->rotationModelSettings = gcrsToItrsRotationModelSettings(
        basic_astrodynamics::iau_2006, globalFrameOrientation,
        std::make_shared< interpolators::InterpolatorGenerationSettings< double > >(
            interpolators::cubicSplineInterpolation( ), initialTimeEnvironment, finalTimeEnvironment, 3600.0 ),
        std::make_shared< interpolators::InterpolatorGenerationSettings< double > >(
            interpolators::cubicSplineInterpolation( ), initialTimeEnvironment, finalTimeEnvironment, 3600.0 ),
        std::make_shared< interpolators::InterpolatorGenerationSettings< double > >(
            interpolators::cubicSplineInterpolation( ), initialTimeEnvironment, finalTimeEnvironment, 60.0 ));
    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );
    bodySettings.at( "Moon" )->ephemerisSettings->resetFrameOrigin( "Earth" );

    // Add spacecraft settings
    std::string spacecraftName = "GRAIL-A";
    std::string spacecraftCentralBody = "Moon";
    bodySettings.addSettings( spacecraftName );
    bodySettings.at( spacecraftName )->ephemerisSettings =
        std::make_shared< InterpolatedSpiceEphemerisSettings >(
                initialTimeEnvironment, finalTimeEnvironment, 10.0, spacecraftCentralBody, globalFrameOrientation );
//    107560 * 3600.0 , 107579.0 * 3600.0, 10.0, spacecraftCentralBody, globalFrameOrientation );


//        std::make_shared< DirectSpiceEphemerisSettings >( spacecraftCentralBody, "J2000" );

    // Create spacecraft
    bodySettings.at( spacecraftName )->constantMass = 1200.0;

    // Create radiation pressure settings
    double referenceAreaRadiation = 10.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies = { "Mars" };
    bodySettings.at( spacecraftName )->radiationPressureSettings[ "Sun" ] =
        cannonBallRadiationPressureSettings( "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );

    // Define ODF data paths
    std::string dataDirectory = "/home/dominic/Tudat/Data/GRAIL_ODF/";
    std::vector< std::string > odfFiles = { "gralugf2012_100_0540smmmv1.odf", "gralugf2012_101_0235smmmv1.odf",
                                            "gralugf2012_102_0358smmmv1.odf", "gralugf2012_103_0145smmmv1.odf",
                                            "gralugf2012_105_0352smmmv1.odf", "gralugf2012_107_0405smmmv1.odf",
                                            "gralugf2012_108_0450smmmv1.odf", "gralugf2012_109_1227smmmv1.odf",
                                            "gralugf2012_111_1332smmmv1.odf", "gralugf2012_114_0900smmmv1.odf" };
//                                            "gralugf2012_117_0300smmmv1.odf", "gralugf2012_119_0545smmmv1.odf",
//                                            "gralugf2012_121_0155smmmv1.odf", "gralugf2012_122_1113smmmv1.odf",
//                                            "gralugf2012_123_0607smmmv1.odf", "gralugf2012_124_0253smmmv1.odf",
//                                            "gralugf2012_125_0400smmmv1.odf", "gralugf2012_126_0203smmmv1.odf",
//                                            "gralugf2012_128_0233smmmv1.odf"};

    // Laod raw ODF data
    std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector;
    for ( std::string odfFile : odfFiles )
    {
        rawOdfDataVector.push_back( std::make_shared< OdfRawFileContents >( dataDirectory + odfFile ) );
    }

    // Process ODF file data
    std::shared_ptr< ProcessedOdfFileContents > processedOdfFileContents =
        std::make_shared< ProcessedOdfFileContents >(
            rawOdfDataVector, spacecraftName, true );
    observation_models::setOdfInformationInBodies( processedOdfFileContents, bodies );

    // Create data structure that handles Observed Data in Tudat
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > observedObservationCollection =
        splitObservationSetsIntoArcs< long double, Time >( observation_models::createOdfObservedObservationCollection< long double, Time >(
            processedOdfFileContents, { dsn_n_way_averaged_doppler } ), 60.0, 10 );


        std::map< int, observation_models::LinkEnds > linkEndIds = observedObservationCollection->getInverseLinkEndIdentifierMap( );
    for( auto it : linkEndIds )
    {
        std::cout<<it.first<<", ("<<it.second[ transmitter ].bodyName_<<", "<<it.second[ transmitter ].stationName_<<"); "
            <<", ("<<it.second[ retransmitter ].bodyName_<<", "<<it.second[ retransmitter ].stationName_<<"); "
            <<", ("<<it.second[ receiver ].bodyName_<<", "<<it.second[ receiver ].stationName_<<")"<<std::endl;

    }

    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< double, double > > > > arcStartEndTimes;
    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > > arcStartEndIndices;

    std::pair< Time, Time > timeBounds = observedObservationCollection->getTimeBounds( );
    Time initialTime = timeBounds.first - 3600.0;
    Time finalTime = timeBounds.second + 3600.0;

    std::cout<<initialTime<<" "<<finalTime<<std::endl;

    // Define observation model settings
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;



    std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back( firstOrderRelativisticLightTimeCorrectionSettings( { "Sun" } ) );
    std::vector< std::string > troposphericCorrectionFileNames =
        {"/home/dominic/Tudat/Data/GRAIL_Ancilliary/grxlugf2012_092_2012_122.tro", "/home/dominic/Tudat/Data/GRAIL_Ancilliary/grxlugf2012_122_2012_153.tro"};
    std::vector< std::string > ionosphericCorrectionFileNames =
        {"/home/dominic/Tudat/Data/GRAIL_Ancilliary/gralugf2012_092_2012_122.ion", "/home/dominic/Tudat/Data/GRAIL_Ancilliary/gralugf2012_122_2012_153.ion"};
    std::map< int, std::string > spacecraftNamePerSpacecraftId;
    spacecraftNamePerSpacecraftId[ 177 ] = "GRAIL-A";
    spacecraftNamePerSpacecraftId[ 181 ] = "GRAIL-B";

    lightTimeCorrectionSettings.push_back( tabulatedTroposphericCorrectionSettings( troposphericCorrectionFileNames ) );
    lightTimeCorrectionSettings.push_back( tabulatedIonosphericCorrectionSettings( ionosphericCorrectionFileNames, spacecraftNamePerSpacecraftId ) );

    std::map < observation_models::ObservableType, std::vector< observation_models::LinkEnds > > linkEndsPerObservable =
        observedObservationCollection->getLinkEndsPerObservableType( );

    for ( auto it = linkEndsPerObservable.begin(); it != linkEndsPerObservable.end(); ++it )
    {
        for ( unsigned int i = 0; i < it->second.size( ); ++i )
        {
            if ( it->first == observation_models::dsn_n_way_averaged_doppler )
            {
                observationModelSettingsList.push_back(
                    observation_models::dsnNWayAveragedDopplerObservationSettings(
                        it->second.at( i ), lightTimeCorrectionSettings, nullptr,
                        std::make_shared< LightTimeConvergenceCriteria >( true ) ) );
            }
        }
    }

    std::vector< std::shared_ptr< ObservationSimulatorBase< long double, Time > > > observationSimulators =
        createObservationSimulators< long double, Time >( observationModelSettingsList, bodies );

    std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< Time > > > observationSimulationSettings =
        getObservationSimulationSettingsFromObservations( observedObservationCollection );
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > computedObservationCollection =
        simulateObservations( observationSimulationSettings, observationSimulators, bodies );

    std::cout<<"Create residuals"<<std::endl;
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > residualObservationCollection =
        createResidualCollection( observedObservationCollection, computedObservationCollection );

    std::cout<<"Filter observations"<<std::endl;
    std::map< ObservableType, double > residualCutoffValuePerObservable;
    residualCutoffValuePerObservable[ dsn_n_way_averaged_doppler ] = 0.05;

    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > filteredObservedObservationCollection =
        filterResidualOutliers( observedObservationCollection, residualObservationCollection, residualCutoffValuePerObservable );

    std::cout<<"Filtered size: "<<observedObservationCollection->getTotalObservableSize( )<<" "<<filteredObservedObservationCollection->getTotalObservableSize( )<<std::endl;


    std::cout<<"Computed filtered observations"<<std::endl;
    std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< Time > > > filteredObservationSimulationSettings =
        getObservationSimulationSettingsFromObservations( filteredObservedObservationCollection );
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > filteredComputedObservationCollection =
        simulateObservations( filteredObservationSimulationSettings, observationSimulators, bodies );

    std::cout<<"Create filtered residuals"<<std::endl;
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > filteredResidualObservationCollection =
        createResidualCollection( filteredObservedObservationCollection, filteredComputedObservationCollection );

    Eigen::VectorXd numericalTimeBiasPartials = getNumericalObservationTimePartial< long double, Time >(
        filteredObservationSimulationSettings, observationSimulators, bodies, 5.0 );


    input_output::writeMatrixToFile( numericalTimeBiasPartials, "grailTestTimeDerivative.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults/");

    std::vector< double > timeBiases;
    std::vector< Eigen::VectorXd > polynomialCoefficientsList;
    Eigen::VectorXd correctedResiduals;
    estimateTimeBiasAndPolynomialFitPerSet< long double, Time >(
        filteredResidualObservationCollection, numericalTimeBiasPartials, timeBiases, polynomialCoefficientsList, correctedResiduals );

    Eigen::VectorXd startTimes;
    Eigen::VectorXd durations;
    Eigen::VectorXd meanValues;
    Eigen::VectorXd rmsValues;

    getResidualStatistics(
        filteredResidualObservationCollection, correctedResiduals, startTimes, durations, meanValues, rmsValues );

//    for( unsigned int i = 0; i < timeBiases.size( ); i++ )
//    {
//        std::cout<<startTimes( i )<<" "<<durations( i )<<" "<<meanValues( i )<<" "<<rmsValues( i )<<" "<<timeBiases.at( i )<<" "<<polynomialCoefficientsList.at( i ).transpose( )<<std::endl;
//    }
    input_output::writeMatrixToFile( startTimes, "grailTestStartTimes.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults/");
    input_output::writeMatrixToFile( durations, "grailTestDurations.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults/");
    input_output::writeMatrixToFile( meanValues, "grailMeanValues.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults/");
    input_output::writeMatrixToFile( rmsValues, "grailRmsValues.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults/");
    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenVector( timeBiases ), "grailTimeBiases.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults/");
//    input_output::writeMatrixToFile( utilities::convertStlVectorToEigenMatrix( polynomialCoefficientsList ), "grailBiases.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults/");

    //
//    std::cout<<"Filter observations"<<std::endl;
//    std::map< ObservableType, double > residualCutoffValuePerObservable;
//    residualCutoffValuePerObservable[ dsn_n_way_averaged_doppler ] = 0.3;
//    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > filteredObservedObservationCollection =
//        filterResidualOutliers( observedObservationCollection, residualObservationCollection, residualCutoffValuePerObservable );
//
//
//    std::cout<<"Computed filtered observations"<<std::endl;
//    std::vector< std::shared_ptr< simulation_setup::ObservationSimulationSettings< Time > > > filteredObservationSimulationSettings =
//        getObservationSimulationSettingsFromObservations( filteredObservedObservationCollection );
//    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > filteredComputedObservationCollection =
//        simulateObservations( filteredObservationSimulationSettings, observationSimulators, bodies );
//
//    std::cout<<"Create filtered residuals"<<std::endl;
//    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > filteredResidualObservationCollection =
//        createResidualCollection( filteredObservedObservationCollection, filteredComputedObservationCollection );



    {
        Eigen::VectorXd residuals = filteredResidualObservationCollection->getObservationVector( ).template cast< double >( );
        input_output::writeMatrixToFile( residuals, "grailTestResiduals.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults/");

        input_output::writeMatrixToFile( correctedResiduals, "grailTestCorrectedResiduals.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults/");

        Eigen::VectorXd observationTimes = utilities::convertStlVectorToEigenVector(
            filteredResidualObservationCollection->getConcatenatedTimeVector( ) ).template cast< double >( );
        input_output::writeMatrixToFile( observationTimes, "grailTestTimes.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults/");

        Eigen::VectorXd observationLinkEndsIds = utilities::convertStlVectorToEigenVector(
            filteredResidualObservationCollection->getConcatenatedLinkEndIds( ) ).template cast< double >( );
        input_output::writeMatrixToFile(observationLinkEndsIds , "grailTestLinkEnds.dat", 16, "/home/dominic/Tudat/Data/GRAIL_TestResults/");
    }
}