/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <tudat/simulation/estimation.h>
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/simulation/environment_setup/createRadiationPressureInterface.h"
#include "tudat/astro/electromagnetism/radiationPressureAcceleration.h"
#include "tudat/astro/electromagnetism/radiationPressureTargetModel.h"
#include <tudat/io/applicationOutput.h>

//! Execute propagation of orbits of Vehicle and Obelix around the Earth.
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
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
    using namespace tudat::electromagnetism;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );

    // Specify initial time
    double initialEphemerisTime = double( 1.0E7 );
    double finalEphemerisTime = initialEphemerisTime + 0.5 * 86400.0;

    // Create bodies needed in simulation
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - 3600.0,finalEphemerisTime + 3600.0 );
    setSimpleRotationSettingsFromSpice( bodySettings, "Earth", initialEphemerisTime );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create spacecraft object.
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );

    // Define constant rotational ephemeris
    //bodies.at( "Vehicle" )->setRotationalEphemeris(
    //            createRotationModel(
    //                orbitalStateBasedRotationSettings( "Sun", false, false, "ECLIPJ2000", "VehicleFixed" ),
    //                "Vehicle", bodies ));
    Eigen::Vector7d rotationalStateVehicle;
    rotationalStateVehicle.segment( 0, 4 ) = linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( Eigen::Matrix3d::Identity() ));
    rotationalStateVehicle.segment( 4, 3 ) = Eigen::Vector3d::Zero();
    bodies.at( "Vehicle" )->setRotationalEphemeris(
                std::make_shared< ConstantRotationalEphemeris >(
                    rotationalStateVehicle, "ECLIPJ2000", "VehicleFixed" ) );
    //bodySettings.at( "Vehicle" )->rotationModelSettings = std::make_shared< SynchronousRotationModelSettings >( "Jupiter", "ECLIPJ2000", "VehicleFixed" );

    std::vector < std::shared_ptr< system_models::VehicleExteriorPanel > > panels;
    panels = {
                    std::make_shared< system_models::VehicleExteriorPanel >(9.9, Eigen::Vector3d::UnitX(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.35, 0.20)),
                    std::make_shared< system_models::VehicleExteriorPanel >(9.9, Eigen::Vector3d::UnitY(),
                                reflectionLawFromSpecularAndDiffuseReflectivity(0.35, 0.25)),
            };
    const std::string panelTypeId = "SolarPanel";
    panels.at(0)->setPanelTypeId(panelTypeId);
    panels.at(1)->setPanelTypeId(panelTypeId);

    bodies.at( "Vehicle" )->setRadiationPressureTargetModels(
            { std::make_shared<PaneledRadiationPressureTargetModel>(panels) } );


    //const auto bodyShape = std::make_shared< FullPanelledBodySettings > bodyWingPanelledGeometry(2,2,2,2,0.1,0.2,0.1,0.2);


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE GROUND STATIONS               //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create ground stations from geodetic positions.
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    groundStationNames.push_back( "Station3" );

    createGroundStation( bodies.at( "Earth" ), "Station1",
                         ( Eigen::Vector3d( ) << 0.0, 1.25, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2",
                         ( Eigen::Vector3d( ) << 0.0, -1.55, 2.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station3",
                         ( Eigen::Vector3d( ) << 0.0, 0.8, 4.0 ).finished( ), geodetic_position );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationSettingsList;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;

    accelerationsOfVehicle[ "Earth" ] = {
            sphericalHarmonicAcceleration( 32, 32 )};
    accelerationsOfVehicle[ "Sun" ] = {
            pointMassGravityAcceleration( ),
            radiationPressureAcceleration( ) };

    accelerationSettingsList[ "Vehicle" ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate = { "Vehicle" };
    std::vector< std::string > centralBodies = { "Earth" };

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationSettingsList, bodiesToIntegrate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Vehicle.
    Eigen::Vector6d vehicleInitialStateInKeplerianElements;
    vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    vehicleInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    // Set initial state.
    double earthGravitationalParameter = getBodyGravitationalParameter( bodies, "Earth" );
    Eigen::Matrix< double, 6, 1 > systemInitialState = convertKeplerianToCartesianElements(
                vehicleInitialStateInKeplerianElements, earthGravitationalParameter );

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, double( finalEphemerisTime ) );

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                double( initialEphemerisTime ), 40.0,
                rungeKuttaFehlberg78,
                40.0, 40.0, 1.0, 1.0 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             DEFINE LINK ENDS FOR OBSERVATIONS            //////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Define link ends.
    std::vector< LinkDefinition > stationReceiverLinkEnds;
    std::vector< LinkDefinition > stationTransmitterLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkDefinition linkEnds;
        linkEnds[ transmitter ] = std::pair< std::string, std::string >( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Vehicle", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds[ receiver ] = std::pair< std::string, std::string >( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Vehicle", "" );
        stationReceiverLinkEnds.push_back( linkEnds );
    }

    // Define (arbitrarily) link ends to be used for 1-way range, 1-way doppler and angular position observables
    std::map< ObservableType, std::vector< LinkDefinition > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );

    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 2 ] );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////    DEFINE PARAMETERS THAT ARE TO BE ESTIMATED      ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define list of parameters to estimate.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames = getInitialStateParameterSettings< double >( propagatorSettings, bodies );
    const std::string dummy = "SolarPanel";
    parameterNames.push_back(std::make_shared<EstimatableParameterSettings>("Vehicle", specular_reflectivity, dummy));
    parameterNames.push_back(std::make_shared<EstimatableParameterSettings>("Vehicle", diffuse_reflectivity, dummy));

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodies );

    // Print identifiers and indices of parameters to terminal.
    printEstimatableParameterEntries( parametersToEstimate );

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    for( std::map< ObservableType, std::vector< LinkDefinition > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;


        std::vector< LinkDefinition > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            observationSettingsList.push_back(
                        std::make_shared< ObservationModelSettings >(
                            currentObservable, currentLinkEndsList.at( i ), std::shared_ptr< LightTimeCorrectionSettings >( ) ) );
        }
    }

    // Create orbit determination object.
    OrbitDeterminationManager< double, double > orbitDeterminationManager = OrbitDeterminationManager< double, double >(
                bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );

    // Compute list of observation times.
    std::vector< double > baseTimeList;
    double observationTime = initialEphemerisTime + 1000.0;
    double  observationInterval = 60.0;
    while(observationTime < finalEphemerisTime - 1000.0)
    {
        baseTimeList.push_back( observationTime);
        observationTime += observationInterval;
    }

    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkDefinition > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkDefinition > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput.push_back(
                        std::make_shared< TabulatedObservationSimulationSettings< > >(
                            currentObservable, currentLinkEndsList[ i ], baseTimeList, receiver ) );
        }
    }

    // Simulate observations.
    std::shared_ptr< ObservationCollection< > > observationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );


    // Set typedefs for POD input (observation types, observation link ends, observation values, associated times with
    // reference link ends.
    typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > >
            SingleObservablePodInputType;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////    PERTURB PARAMETER VECTOR AND ESTIMATE PARAMETERS     ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Perturb parameter estimate.
    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< double >( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    Eigen::Matrix< double, Eigen::Dynamic, 1 > parameterPerturbation =
            Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );

    // Perturbe initial state estimate.
    parameterPerturbation.segment( 0, 3 ) = Eigen::Vector3d::Constant( 1.0 );
    parameterPerturbation.segment( 3, 3 ) = Eigen::Vector3d::Constant( 1.E-3 );
    parameterPerturbation.segment( 6, 1 ) = Eigen::Vector1d::Constant( 0.05 );
    parameterPerturbation.segment( 7, 1 ) = Eigen::Vector1d::Constant( 0.05 );

    initialParameterEstimate += parameterPerturbation;
    parametersToEstimate->resetParameterValues( initialParameterEstimate );


    // Define estimation input
    std::shared_ptr< EstimationInput< double, double  > > estimationInput = std::make_shared< EstimationInput< double, double > >(
            observationsAndTimes );

    std::map< observation_models::ObservableType, double > weightPerObservable;
    weightPerObservable[ one_way_range ] = 1.0 / ( 1.0 * 1.0 );
    weightPerObservable[ one_way_doppler ] = 1.0 / ( 1.0E-11 * 1.0E-11 * physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT  );

    estimationInput->setConstantPerObservableWeightsMatrix( weightPerObservable );
    estimationInput->defineEstimationSettings( true, true, true, true, true );
    estimationInput->setConvergenceChecker(
            std::make_shared< EstimationConvergenceChecker >( 10 ) );

    // Perform estimation
    std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters(
            estimationInput );

    input_output::writeMatrixToFile( estimationOutput->getParameterHistoryMatrix( ),
                                     "MichaelearthOrbitParameterHistory.dat", 16,
                                     tudat_applications::getOutputPath( )  );

    Eigen::VectorXd parameterEstimate = estimationOutput->parameterEstimate_;
    std::cout <<"parameter estimate: "<< ( parameterEstimate ).transpose( ) << std::endl;

    Eigen::VectorXd estimationError = parameterEstimate - truthParameters;
    std::cout <<"estimation error: "<< ( estimationError ).transpose( ) << std::endl;


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}