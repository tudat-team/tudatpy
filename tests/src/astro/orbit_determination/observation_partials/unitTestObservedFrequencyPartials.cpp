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
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"
#include "tudat/io/readOdfFile.h"

#include "tudat/io/readOdfFile.h"
#include "tudat/io/readTabulatedMediaCorrections.h"
#include "tudat/io/readTabulatedWeatherData.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"
#include "tudat/simulation/estimation_setup/processTrackingTxtFile.h"

#include <boost/date_time/gregorian/gregorian.hpp>

#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/support/observationPartialTestFunctions.h"


//Using declarations.
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;
using namespace tudat::coordinate_conversions;
using namespace tudat::ground_stations;
using namespace tudat::observation_models;
using namespace tudat::ephemerides;
using namespace tudat::electromagnetism;
using namespace tudat::input_output;
using namespace tudat;

namespace tudat
{
namespace unit_tests
{



BOOST_AUTO_TEST_SUITE( test_frequency_observation_partials)

class LinearStateWrapper
{
public:
    LinearStateWrapper( const double referenceTime, Eigen::Vector6d referenceState  ):
        referenceState_( referenceState ), referenceTime_( referenceTime )
    { }

    Eigen::Vector6d getLinearState( const double currentTime )
    {
        Eigen::Vector6d currentState = referenceState_;
        currentState.segment( 0, 3 ) += ( currentTime - referenceTime_ ) * referenceState_.segment( 3, 3 );
        return currentState;
    }

    Eigen::Vector6d getState( )
    {
        return referenceState_;
    }

    void setState( const Eigen::Vector6d& state )
    {
        referenceState_ = state;
    }

    Eigen::Vector6d referenceState_;
    double referenceTime_;
};

void testPartials(
    const std::shared_ptr< observation_models::ObservationModel< 1, double, double > >  observationModel,
    const std::shared_ptr< ObservationPartial< 1 > > observationPartial,
    const std::shared_ptr< PositionPartialScaling > observationPartialScaling,
    LinearStateWrapper& stateWrapper,
    const double testObservationTime,
    const std::shared_ptr<ObservationAncilliarySimulationSettings> ancilliarySettings = nullptr )
{
    std::vector< Eigen::Vector6d > vectorOfStates;
    std::vector< double > vectorOfTimes;

    // Compute nominal observation
    Eigen::Matrix< double, 1, 1 > currentObservation = observationModel->computeObservationsWithLinkEndData(
        testObservationTime, receiver, vectorOfTimes, vectorOfStates, ancilliarySettings );

    // Update scaling
    observationPartialScaling->update( vectorOfStates, vectorOfTimes, receiver, currentObservation );

    // Compute partial
    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > singlePartialSet =
        observationPartial->calculatePartial( vectorOfStates, vectorOfTimes, receiver, ancilliarySettings, currentObservation.template cast< double >( ) );

    // Combine partial
    Eigen::Matrix< double, 1, 6 > totalPartial;
    totalPartial.setZero( );
    for( int i = 0; i < singlePartialSet.size( ); i++ )
    {
        totalPartial += singlePartialSet.at( i ).first;
    }

    // Initialize numerical aprtial
    Eigen::Matrix< double, 1, 6 > numericalPartial = Eigen::Matrix< double, 1, 6 >::Zero( );
    double positionPerturbation = 10000.0;
    double velocityPerturbation = 1.0;

    // Compute numerical state partials
    Eigen::Vector6d nominalReferenceState = stateWrapper.getState( );
    Eigen::Vector6d perturbedReferenceState = nominalReferenceState;
    for( int i = 0; i < 3; i++ )
    {
        // Perturb parameter up
        perturbedReferenceState = nominalReferenceState;
        perturbedReferenceState( i ) += positionPerturbation;
        stateWrapper.setState( perturbedReferenceState );
        Eigen::Matrix< double, 1, 1 > upperturbedObservation = observationModel->computeObservationsWithLinkEndData(
            testObservationTime, receiver, vectorOfTimes, vectorOfStates, ancilliarySettings );

        // Perturb parameter down
        perturbedReferenceState = nominalReferenceState;
        perturbedReferenceState( i ) -= positionPerturbation;
        stateWrapper.setState( perturbedReferenceState );
        Eigen::Matrix< double, 1, 1 > downperturbedObservation = observationModel->computeObservationsWithLinkEndData(
            testObservationTime, receiver, vectorOfTimes, vectorOfStates, ancilliarySettings );
        numericalPartial( i ) = ( upperturbedObservation( 0 ) - downperturbedObservation( 0 ) ) / ( 2.0 * positionPerturbation );


        perturbedReferenceState = nominalReferenceState;
        perturbedReferenceState( i + 3 ) += velocityPerturbation;
        stateWrapper.setState( perturbedReferenceState );
        upperturbedObservation = observationModel->computeObservationsWithLinkEndData(
            testObservationTime, receiver, vectorOfTimes, vectorOfStates, ancilliarySettings );

        perturbedReferenceState = nominalReferenceState;
        perturbedReferenceState( i + 3 ) -= velocityPerturbation;
        stateWrapper.setState( perturbedReferenceState );
        downperturbedObservation = observationModel->computeObservationsWithLinkEndData(
            testObservationTime, receiver, vectorOfTimes, vectorOfStates, ancilliarySettings );
        stateWrapper.referenceState_( i + 3 ) += velocityPerturbation;
        numericalPartial( i + 3 ) = ( upperturbedObservation( 0 ) - downperturbedObservation( 0 ) ) / ( 2.0 * velocityPerturbation );
    }
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( totalPartial, numericalPartial, 1.0E-4 );
}

BOOST_AUTO_TEST_CASE( testFrequencyDopplerPartialsDirect )
{
    // Define reference time
    double referenceTime = 777447060.682842;

    // Load Spice kernels
    std::string kernelsPath = paths::getSpiceKernelPath( );
    spice_interface::loadStandardSpiceKernels( );
    spice_interface::loadSpiceKernelInTudat( tudat::paths::getTudatTestDataPath( ) + "juice_orbc_000074_230414_310721_v01.bsp" );

    // Define bodies to use.
    std::vector<std::string> bodiesToCreate = { "Earth", "Moon", "Sun", "Jupiter" };
    std::string globalFrameOrigin = "SSB";
    std::string globalFrameOrientation = "J2000";

    // Create bodies settings needed in simulation
    BodyListSettings bodySettings = getDefaultBodySettings(
        bodiesToCreate, globalFrameOrigin, globalFrameOrientation );

    // Add ground stations
    std::shared_ptr<GroundStationSettings> nnorciaSettings = std::make_shared<GroundStationSettings>(
        "NWNORCIA", getCombinedApproximateGroundStationPositions( ).at( "NWNORCIA" ));
    nnorciaSettings->addStationMotionSettings(
        std::make_shared<LinearGroundStationMotionSettings>(
            ( Eigen::Vector3d( ) << -45.00, 10.00, 47.00 ).finished( ) / 1.0E3 / physical_constants::JULIAN_YEAR,
            0.0 ));
    std::shared_ptr<GroundStationSettings> yarragadeeSettings = std::make_shared<GroundStationSettings>(
        "YARRAGAD", getCombinedApproximateGroundStationPositions( ).at( "YARRAGAD" ));
    yarragadeeSettings->addStationMotionSettings(
        std::make_shared<LinearGroundStationMotionSettings>(
            ( Eigen::Vector3d( ) << -47.45, 9.12, 51.76 ).finished( ) / 1.0E3 / physical_constants::JULIAN_YEAR, 0.0 ));
    bodySettings.at( "Earth" )->groundStationSettings.push_back( nnorciaSettings );
    bodySettings.at( "Earth" )->groundStationSettings.push_back( yarragadeeSettings );

    // Create Spacecraft
    const std::string spacecraftName = "JUICE";
    bodiesToCreate.push_back( spacecraftName );
    bodySettings.addSettings( spacecraftName );

    // Create custom spacecraft ephemeris
    LinearStateWrapper stateWrapper( referenceTime, spice_interface::getBodyCartesianStateAtEpoch(
        "JUICE", "Earth", "J2000", "None", referenceTime ) );
    bodySettings.get( spacecraftName )->ephemerisSettings = customEphemerisSettings(
        std::bind( &LinearStateWrapper::getLinearState, &stateWrapper, std::placeholders::_1 ), "Earth", "J2000" );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Set turnaround ratios in spacecraft (ground station)
    std::shared_ptr<system_models::VehicleSystems> vehicleSystems = std::make_shared<system_models::VehicleSystems>( );
    vehicleSystems->setTransponderTurnaroundRatio( &getDsnDefaultTurnaroundRatios );
    double xBandRatio = getDsnDefaultTurnaroundRatios( x_band, x_band );
    bodies.at( spacecraftName )->setVehicleSystems( vehicleSystems );
    bodies.processBodyFrameDefinitions( );

    // Set station transmitting frequencies
    double transmissionFrequency = 7180127320;
    std::shared_ptr<ground_stations::StationFrequencyInterpolator> transmittingFrequencyCalculator =
        std::make_shared<ground_stations::ConstantFrequencyInterpolator>( transmissionFrequency );
    bodies.at( "Earth" )->getGroundStation( "NWNORCIA" )->setTransmittingFrequencyCalculator(
        transmittingFrequencyCalculator );

    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[ transmitter ] = std::make_pair<std::string, std::string>( "Earth", static_cast<std::string>( "NWNORCIA" ));
    linkEnds[ retransmitter ] = std::make_pair<std::string, std::string>( static_cast<std::string>(spacecraftName), "" );
    linkEnds[ receiver ] = std::make_pair<std::string, std::string>( "Earth", static_cast<std::string>( "YARRAGAD" ));

    // Create observation settings
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionsList;
//        lightTimeCorrectionsList.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( std::vector< std::string >( { "Sun", "Moon", "Earth" } ) ) );

    // Create Doppler model (m/s)
    std::shared_ptr<DopplerMeasuredFrequencyObservationModel<double, double> > dopplerFrequencyObservationModel =
        std::dynamic_pointer_cast<DopplerMeasuredFrequencyObservationModel<double, double>>(
            ObservationModelCreator<1, double, double>::createObservationModel(
                std::make_shared<ObservationModelSettings>( doppler_measured_frequency, linkEnds, lightTimeCorrectionsList   ), bodies ));

    // Create ancilliary settings
    std::shared_ptr<ObservationAncilliarySimulationSettings> ancillarySettings =
        std::make_shared<ObservationAncilliarySimulationSettings>( );
    ancillarySettings->setAncilliaryDoubleVectorData( frequency_bands, { x_band, x_band } );

    // Create Doppler model (frequency)
    std::shared_ptr<TwoWayDopplerObservationModel<double, double> > twoWayDopplerObservationModel =
        std::dynamic_pointer_cast<TwoWayDopplerObservationModel<double, double>>(
            ObservationModelCreator<1, double, double>::createObservationModel(
                std::make_shared<TwoWayDopplerObservationSettings>( linkEnds, lightTimeCorrectionsList  ), bodies ) );

    // Define link end
    LinkEndType referenceLinkEnd = receiver;

    // Create acceleration settings
    SelectedAccelerationMap accelerationMap;
    std::map<std::string, std::vector<std::shared_ptr<AccelerationSettings> > > accelerationsOfSpacecraft;
    accelerationsOfSpacecraft[ "Earth" ].push_back( std::make_shared<AccelerationSettings >( point_mass_gravity ));
    accelerationMap[ "JUICE" ] = accelerationsOfSpacecraft;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector<std::string> bodiesToEstimate = { "JUICE" };
    std::vector<std::string> centralBodies = { "Earth" };
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
        bodies, accelerationMap, bodiesToEstimate, centralBodies );

    // Get initial state from Spice kernel
    Eigen::VectorXd initialState = getInitialStateOfBody( "JUICE", "Earth", bodies, referenceTime );

    // Define propagator settings
    std::shared_ptr<TranslationalStatePropagatorSettings<double, double> > propagatorSettings =
        std::make_shared<TranslationalStatePropagatorSettings<double, double> >
            ( centralBodies, accelerationModelMap, bodiesToEstimate, initialState.template cast<double>( ),
              double( referenceTime ), numerical_integrators::rungeKuttaFixedStepSettings<double>( 10.0, numerical_integrators::rungeKuttaFehlberg78 ),
              std::make_shared<PropagationTimeTerminationSettings>( referenceTime + 3600.0 ));

    // Create parameters
    std::vector<std::shared_ptr<EstimatableParameterSettings> > parameterNames = getInitialStateParameterSettings<double, double>( propagatorSettings, bodies );
    std::shared_ptr<estimatable_parameters::EstimatableParameterSet<double> > parametersToEstimate =
        createParametersToEstimate<double, double>( parameterNames, bodies, propagatorSettings );

    // Manually create partials for doppler observable (m/s)
    std::map< LinkEnds, std::shared_ptr< observation_models::ObservationModel< 1, double, double > > > twoWayDopplerObservationModelList;
    twoWayDopplerObservationModelList[ linkEnds ] = twoWayDopplerObservationModel;
    std::map< LinkEnds, std::pair< std::map< std::pair< int, int >,
        std::shared_ptr< ObservationPartial< 1 > > >,
        std::shared_ptr< PositionPartialScaling > > > twoWayDopplerObservationPartialsAndScaler = createObservablePartialsList(
            twoWayDopplerObservationModelList, bodies, parametersToEstimate, false, false );
    std::shared_ptr< ObservationPartial< 1 > > twoWayDopplerPartial = twoWayDopplerObservationPartialsAndScaler.begin( )->second.first.begin( )->second;
    std::shared_ptr< PositionPartialScaling > twoWayDopplerScaling = twoWayDopplerObservationPartialsAndScaler.begin( )->second.second;

    // Manually create partials for doppler observable (frequency)
    std::map< LinkEnds, std::shared_ptr< observation_models::ObservationModel< 1, double, double > > > frequencyDopplerObservationModelList;
    frequencyDopplerObservationModelList[ linkEnds ] = dopplerFrequencyObservationModel;
    std::map< LinkEnds, std::pair< std::map< std::pair< int, int >,
        std::shared_ptr< ObservationPartial< 1 > > >,
        std::shared_ptr< PositionPartialScaling > > > frequencyDopplerObservationPartialsAndScaler = createObservablePartialsList(
        frequencyDopplerObservationModelList, bodies, parametersToEstimate, false, false );
    std::shared_ptr< ObservationPartial< 1 > > frequencyDopplerPartial = frequencyDopplerObservationPartialsAndScaler.begin( )->second.first.begin( )->second;
    std::shared_ptr< PositionPartialScaling > frequencyDopplerScaling = frequencyDopplerObservationPartialsAndScaler.begin( )->second.second;

    // Compute manual
    testPartials( twoWayDopplerObservationModel, twoWayDopplerPartial, twoWayDopplerScaling, stateWrapper, referenceTime );
    testPartials( dopplerFrequencyObservationModel, frequencyDopplerPartial, frequencyDopplerScaling, stateWrapper, referenceTime, ancillarySettings );

}

//! Test partial derivatives of observed frequency observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testFrequencyDopplerPartials )
{
    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 2 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );


    // Test partials with constant ephemerides (allows test of position partials)
    {
        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );

        std::shared_ptr< system_models::VehicleSystems > vehicleSystems = std::make_shared< system_models::VehicleSystems >( );
        vehicleSystems->setTransponderTurnaroundRatio( [=]( FrequencyBands uplinkBand, FrequencyBands downlinkBand ){ return 1.0; } );
        bodies.getBody( "Earth" )->getGroundStation( "Graz" )->setVehicleSystems( vehicleSystems );
        bodies.getBody( "Mars" )->getGroundStation( "MSL" )->setTransmittingFrequencyCalculator( std::make_shared< ConstantFrequencyInterpolator >( 1.0 ) );

        // Set link ends for observation model
        LinkDefinition linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ reflector1 ] = groundStations[ 0 ];
        linkEnds[ receiver ] = groundStations[ 1 ];


            std::cout << "Case A "  << std::endl;
            // Generate one-way doppler model
            std::shared_ptr< ObservationModel< 1 > > twoWayDopplerModel;
            std::vector< std::string > perturbingBodies;
            perturbingBodies.push_back( "Earth" );

            twoWayDopplerModel =
                observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                    std::make_shared< TwoWayDopplerObservationSettings >
                        (  std::make_shared< observation_models::OneWayDopplerObservationSettings >(
                               getUplinkFromTwoWayLinkEnds( linkEnds ),
                               std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                   perturbingBodies ),
                               std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" ),
                               std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                               nullptr, std::make_shared< LightTimeConvergenceCriteria>( ), true ),
                           std::make_shared< observation_models::OneWayDopplerObservationSettings >(
                               getDownlinkFromTwoWayLinkEnds( linkEnds ),
                               std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                   perturbingBodies ),
                               std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                               std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" ),
                               nullptr, std::make_shared< LightTimeConvergenceCriteria>( ), true ),
                           doppler_measured_frequency ), bodies );

            // Create parameter objects.
            std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet;
            Eigen::VectorXd parameterPerturbationMultipliers = Eigen::Vector4d::Constant( 1.0 );
            fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );

            printEstimatableParameterEntries( fullEstimatableParameterSet );

            testObservationPartials< 1 >(
                twoWayDopplerModel, bodies, fullEstimatableParameterSet, linkEnds, doppler_measured_frequency, 1.0E-5,
                true, true, 10.0, parameterPerturbationMultipliers, getDopplerMeasuredFrequencyAncilliarySettings(
                    std::vector< FrequencyBands >{ x_band, x_band } ), 1.1E7, 100.0 );
    }

    // Test partials with real ephemerides (without test of position partials)
    {

        // Create environment
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

        std::shared_ptr< system_models::VehicleSystems > vehicleSystems = std::make_shared< system_models::VehicleSystems >( );
        vehicleSystems->setTransponderTurnaroundRatio( [=]( FrequencyBands uplinkBand, FrequencyBands downlinkBand ){ return 1.0; } );
        bodies.getBody( "Earth" )->getGroundStation( "Graz" )->setVehicleSystems( vehicleSystems );
        bodies.getBody( "Mars" )->getGroundStation( "MSL" )->setTransmittingFrequencyCalculator( std::make_shared< ConstantFrequencyInterpolator >( 1.0 ) );

        // Set link ends for observation model
        LinkDefinition linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ reflector1 ] = groundStations[ 0 ];
        linkEnds[ receiver ] = groundStations[ 1 ];

        std::cout << "Case B " << std::endl;
        // Generate two-way doppler model
        std::shared_ptr< ObservationModel< 1 > > twoWayDopplerModel;
        std::vector< std::string > perturbingBodies;
        perturbingBodies.push_back( "Earth" );

        twoWayDopplerModel =
            observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                std::make_shared< TwoWayDopplerObservationSettings >
                    (  std::make_shared< OneWayDopplerObservationSettings >(
                           getUplinkFromTwoWayLinkEnds( linkEnds ),
                           std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                               perturbingBodies ),
                           std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" ),
                           std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                           nullptr, std::make_shared< LightTimeConvergenceCriteria>( ), true  ),
                       std::make_shared< OneWayDopplerObservationSettings >(
                           getDownlinkFromTwoWayLinkEnds( linkEnds ),
                           std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                               perturbingBodies ),
                           std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                           std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Mars" ) ,
                           nullptr, std::make_shared< LightTimeConvergenceCriteria>( ), true ), doppler_measured_frequency ), bodies );

        // Create parameter objects.
        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet;
        Eigen::VectorXd parameterPerturbationMultipliers = Eigen::Vector4d::Constant( 1.0 );
        fullEstimatableParameterSet = createEstimatableParameters( bodies, 1.1E7 );

        printEstimatableParameterEntries( fullEstimatableParameterSet );

        testObservationPartials< 1 >(
            twoWayDopplerModel, bodies, fullEstimatableParameterSet, linkEnds, doppler_measured_frequency, 1.0E-4, false, true,
            1.0, parameterPerturbationMultipliers, getDopplerMeasuredFrequencyAncilliarySettings(
                std::vector< FrequencyBands >{ x_band, x_band } ), 1.1E7, 100.0 );

    }
}

BOOST_AUTO_TEST_SUITE_END( )


}// namespace unit_tests

}// namespace tudat
