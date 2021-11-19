/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef ORBITDETERMINATIONTESTCASES_H
#define ORBITDETERMINATIONTESTCASES_H

#include <boost/make_shared.hpp>

#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/simulation.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"

namespace tudat
{
namespace unit_tests
{

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


Eigen::VectorXd getDefaultInitialParameterPerturbation( );

template< typename TimeType = double, typename StateScalarType  = double >
std::pair< std::shared_ptr< PodOutput< StateScalarType, TimeType > >, Eigen::VectorXd > executePlanetaryParameterEstimation(
        const int observableType = 1,
        Eigen::VectorXd parameterPerturbation = getDefaultInitialParameterPerturbation( ),
        Eigen::MatrixXd inverseAPrioriCovariance  = Eigen::MatrixXd::Zero( 7, 7 ),
        const double weight = 1.0 )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    //Define setting for total number of bodies and those which need to be integrated numerically.
    //The first numberOfNumericalBodies from the bodyNames vector will be integrated numerically.

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Saturn" );

    // Specify initial time
    TimeType initialEphemerisTime = TimeType( 1.0E7 );
    TimeType finalEphemerisTime = TimeType( 3.0E7 );
    double maximumTimeStep = 3600.0;

    double buffer = 10.0 * maximumTimeStep;

    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames,initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Moon" )->ephemerisSettings->resetFrameOrigin( "Sun" );

    // Create bodies needed in simulation
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );


    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
    accelerationsOfEarth[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfEarth[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfEarth[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfEarth[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfEarth[ "Saturn" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

    accelerationMap[ "Earth" ] = accelerationsOfEarth;


    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToEstimate;
    bodiesToEstimate.push_back( "Earth" );
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Earth" );
    unsigned int numberOfNumericalBodies = bodiesToIntegrate.size( );

    // Define propagator settings.
    std::vector< std::string > centralBodies;
    std::map< std::string, std::string > centralBodyMap;

    centralBodies.resize( numberOfNumericalBodies );
    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        centralBodies[ i ] = "SSB";
        centralBodyMap[ bodiesToIntegrate[ i ] ] = centralBodies[ i ];
    }

    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, centralBodyMap );

    // Set parameters that are to be estimated.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                                  "Earth", propagators::getInitialStateOfBody< TimeType, StateScalarType >(
                                      "Earth", centralBodyMap.at( "Earth" ), bodies, initialEphemerisTime ),
                              centralBodyMap.at( "Earth" ) ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Moon", gravitational_parameter ) );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate< StateScalarType >( parameterNames, bodies );


    // Define integrator settings.
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            std::make_shared< IntegratorSettings< TimeType > >(
                rungeKutta4, TimeType( initialEphemerisTime - 4.0 * maximumTimeStep ), 900.0 );


    std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate,
              getInitialStateVectorOfBodiesToEstimate( parametersToEstimate ),
              TimeType( finalEphemerisTime + 4.0 * maximumTimeStep ),
              cowell, std::shared_ptr< DependentVariableSaveSettings >( ) );


    // Define link ends
    LinkEnds linkEnds;
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;

    if( observableType == 0 )
    {
        linkEnds[ observed_body ] = std::make_pair( "Earth", "" );
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                                           position_observable, linkEnds ) );
    }
    else if( observableType == 5 )
    {
        linkEnds[ observed_body ] = std::make_pair( "Earth", "" );
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                                          velocity_observable, linkEnds ) );
    }
    else
    {
        linkEnds[ transmitter ] = std::make_pair( "Earth", "" );
        linkEnds[ receiver ] = std::make_pair( "Mars", "" );

        if( observableType == 1 )
        {

            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                                               one_way_range, linkEnds ) );
        }
        else if( observableType == 2 )
        {
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                                               angular_position, linkEnds ) );
        }
        else if( observableType == 3 )
        {
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                                               one_way_doppler, linkEnds ) );
        }
        else if( observableType == 4 )
        {
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                                               one_way_range, linkEnds ) );
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                                               one_way_doppler, linkEnds ) );
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                                               angular_position, linkEnds ) );
        }
    }



    // Create orbit determination object.
    OrbitDeterminationManager< StateScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< StateScalarType, TimeType >(
                bodies, parametersToEstimate, observationSettingsList,
                integratorSettings, propagatorSettings );


    // Define observation times.
    double observationTimeStep = 1000.0;
    TimeType observationTime = Time( initialEphemerisTime + 10.0E4 );
    int numberOfObservations = 18000;

    std::vector< TimeType > initialObservationTimes;
    initialObservationTimes.resize( numberOfObservations );

    for( int i = 0; i < numberOfObservations; i++ )
    {
        initialObservationTimes[ i ] = observationTime;
        observationTime += observationTimeStep;
    }

    // Define observation simulation settings
    std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > measurementSimulationInput;
    initialObservationTimes = utilities::addScalarToVector( initialObservationTimes, 30.0 );
    if( observableType == 0 )
    {
        measurementSimulationInput.push_back(
                    std::make_shared< TabulatedObservationSimulationSettings< > >(
                        position_observable, linkEnds, initialObservationTimes, observed_body ) );
    }
    else if( observableType == 5 )
    {
        measurementSimulationInput.push_back(
                    std::make_shared< TabulatedObservationSimulationSettings< > >(
                        velocity_observable, linkEnds, initialObservationTimes, observed_body ) );
    }
    else
    {
        if( observableType == 1 )
        {
            measurementSimulationInput.push_back(
                        std::make_shared< TabulatedObservationSimulationSettings< > >(
                            one_way_range, linkEnds, initialObservationTimes, transmitter ) );
        }
        else if( observableType == 2 )
        {
            measurementSimulationInput.push_back(
                        std::make_shared< TabulatedObservationSimulationSettings< > >(
                            angular_position, linkEnds, initialObservationTimes, transmitter ) );
        }
        else if( observableType == 3 )
        {
            measurementSimulationInput.push_back(
                        std::make_shared< TabulatedObservationSimulationSettings< > >(
                            one_way_doppler, linkEnds, initialObservationTimes, transmitter ) );
        }
        else if( observableType == 4 )
        {
            measurementSimulationInput.push_back(
                        std::make_shared< TabulatedObservationSimulationSettings< > >(
                            one_way_range, linkEnds, initialObservationTimes, transmitter ) );
            measurementSimulationInput.push_back(
                        std::make_shared< TabulatedObservationSimulationSettings< > >(
                            angular_position, linkEnds, initialObservationTimes, transmitter ) );
            measurementSimulationInput.push_back(
                        std::make_shared< TabulatedObservationSimulationSettings< > >(
                            one_way_doppler, linkEnds, initialObservationTimes, transmitter ) );
        }
    }

    // Simulate observations
    std::shared_ptr< ObservationCollection< StateScalarType, TimeType > > simulatedObservations = simulateObservations< StateScalarType, TimeType >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Perturb parameter estimate
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    if( parameterPerturbation.rows( ) == 0 )
    {
        parameterPerturbation = Eigen::VectorXd::Zero( 7 );
    }
    for( unsigned int i = 0; i < 7; i++ )
    {
        initialParameterEstimate( i ) += parameterPerturbation( i );
    }

    // Define estimation input
    std::shared_ptr< PodInput< StateScalarType, TimeType > > podInput =
            std::make_shared< PodInput< StateScalarType, TimeType > >(
                simulatedObservations, initialParameterEstimate.rows( ), inverseAPrioriCovariance,
                initialParameterEstimate - truthParameters );
    if( observableType == 4 )
    {
        std::map< observation_models::ObservableType, double > weightPerObservable;
        weightPerObservable[ one_way_range ] = 1.0 / ( 1.0 * 1.0 );
        weightPerObservable[ angular_position ] = 1.0 / ( 1.0E-9 * 1.0E-9 );
        weightPerObservable[ one_way_doppler ] = 1.0 / ( 1.0E-12 * 1.0E-12 );

        podInput->setConstantPerObservableWeightsMatrix( weightPerObservable );
    }
    else
    {
        podInput->setConstantWeightsMatrix( weight );
    }
    podInput->defineEstimationSettings( true, true, false, false, false );

    // Perform estimation
    std::shared_ptr< PodOutput< StateScalarType, TimeType > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( ) );

    return std::make_pair( podOutput,
                           ( podOutput->parameterEstimate_.template cast< double >( ) -
                             truthParameters .template cast< double >( ) ) );
}


extern template std::pair< std::shared_ptr< PodOutput< double > >, Eigen::VectorXd > executePlanetaryParameterEstimation< double, double >(
        const int observableType,
        Eigen::VectorXd parameterPerturbation,
        Eigen::MatrixXd inverseAPrioriCovariance,
        const double weight );

#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template std::pair< std::shared_ptr< PodOutput< long double > >, Eigen::VectorXd > executePlanetaryParameterEstimation< double, long double >(
        const int observableType,
        Eigen::VectorXd parameterPerturbation,
        Eigen::MatrixXd inverseAPrioriCovariance,
        const double weight );
extern template std::pair< std::shared_ptr< PodOutput< double, Time > >, Eigen::VectorXd > executePlanetaryParameterEstimation< Time, double >(
        const int observableType ,
        Eigen::VectorXd parameterPerturbation,
        Eigen::MatrixXd inverseAPrioriCovariance,
        const double weight );
extern template std::pair< std::shared_ptr< PodOutput< long double, Time > >, Eigen::VectorXd > executePlanetaryParameterEstimation< Time, long double >(
        const int observableType,
        Eigen::VectorXd parameterPerturbation,
        Eigen::MatrixXd inverseAPrioriCovariance,
        const double weight );
#endif


template< typename TimeType = double, typename StateScalarType = double >
Eigen::VectorXd executeEarthOrbiterParameterEstimation(
        std::pair< std::shared_ptr< PodOutput< StateScalarType > >,
        std::shared_ptr< PodInput< StateScalarType, TimeType > > >& podData,
        const TimeType startTime = TimeType( 1.0E7 ),
        const int numberOfDaysOfData = 3,
        const int numberOfIterations = 5,
        const bool useFullParameterSet = true )
{

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );

    // Specify initial time
    TimeType initialEphemerisTime = startTime;
    TimeType finalEphemerisTime = initialEphemerisTime + numberOfDaysOfData * 86400.0;

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "ECLIPJ2000", "IAU_Earth",
                spice_interface::computeRotationQuaternionBetweenFrames(
                    "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
                initialEphemerisTime, 2.0 * mathematical_constants::PI /
                ( physical_constants::JULIAN_DAY ) );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 4.0;
    double aerodynamicCoefficient = 1.2;
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceArea, aerodynamicCoefficient * ( Eigen::Vector3d( ) << 1.2, -0.1, -0.4 ).finished( ), 1, 1 );

    // Create and set aerodynamic coefficients object
    bodies.at( "Vehicle" )->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodies.at( "Vehicle" )->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    asterixRadiationPressureSettings, "Vehicle", bodies ) );

    bodies.at( "Vehicle" )->setEphemeris( std::make_shared< TabulatedCartesianEphemeris< > >(
                                            std::shared_ptr< interpolators::OneDimensionalInterpolator
                                            < double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );


    // Creatre ground stations: same position, but different representation
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    groundStationNames.push_back( "Station3" );

    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station3", ( Eigen::Vector3d( ) << 0.0, 0.05, 4.0 ).finished( ), geodetic_position );

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 8, 8 ) );
    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                    basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                    basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::cannon_ball_radiation_pressure ) );
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::aerodynamic ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Set (perturbed) initial state.
    Eigen::Matrix< StateScalarType, 6, 1 > systemInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState,
              TimeType( finalEphemerisTime ), cowell );

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettings< TimeType > >
            ( TimeType( initialEphemerisTime ), 40.0,
              RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78,
              40.0, 40.0, 1.0, 1.0 );

    // Define parameters.
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ receiver ] = std::make_pair( "Vehicle", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = std::make_pair( "Vehicle", "" );
        stationReceiverLinkEnds.push_back( linkEnds );
    }

    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );

    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 2 ] );

    linkEndsPerObservable[ angular_position ].push_back( stationReceiverLinkEnds[ 2 ] );
    linkEndsPerObservable[ angular_position ].push_back( stationTransmitterLinkEnds[ 1 ] );

    std::cout << "Link ends " << getLinkEndsString( linkEndsPerObservable[ one_way_doppler ].at( 0 ) ) << std::endl;

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back(
                std::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                    "Vehicle", systemInitialState, "Earth" ) );

    if( useFullParameterSet )
    {
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", radiation_pressure_coefficient ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", constant_drag_coefficient ) );
        parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                      linkEndsPerObservable.at( one_way_range ).at( 0 ), one_way_range, true ) );
        parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                      linkEndsPerObservable.at( one_way_range ).at( 0 ), one_way_range, false ) );
        parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                      linkEndsPerObservable.at( one_way_range ).at( 1 ), one_way_range, false ) );

        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                      2, 0, 2, 2, "Earth", spherical_harmonics_cosine_coefficient_block ) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                                      2, 1, 2, 2, "Earth", spherical_harmonics_sine_coefficient_block ) );
        parameterNames.push_back(  std::make_shared< EstimatableParameterSettings >
                                   ( "Earth", rotation_pole_position ) );
        parameterNames.push_back(  std::make_shared< EstimatableParameterSettings >
                                   ( "Earth", ground_station_position, "Station1" ) );
    }

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodies );

    printEstimatableParameterEntries( parametersToEstimate );

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;

    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;


        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            std::shared_ptr< ObservationBiasSettings > biasSettings;
            if( ( currentObservable == one_way_range ) && ( i == 0 ) )
            {
                std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;

                biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >(
                                                Eigen::Vector1d::Zero( ), true ) );
                biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >(
                                                Eigen::Vector1d::Zero( ), false ) );
                biasSettings = std::make_shared< MultipleObservationBiasSettings >(
                            biasSettingsList );
            }
            else if( ( currentObservable == one_way_range ) && ( i == 1 ) )
            {
                biasSettings = std::make_shared< ConstantObservationBiasSettings >(
                            Eigen::Vector1d::Zero( ), false );
            }

            observationSettingsList.push_back(
                                        std::make_shared< ObservationModelSettings >(
                                            currentObservable, currentLinkEndsList.at( i ), std::shared_ptr< LightTimeCorrectionSettings >( ),
                                            biasSettings ) );
        }
    }

    // Create orbit determination object.
    OrbitDeterminationManager< StateScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< StateScalarType, TimeType >(
                bodies, parametersToEstimate, observationSettingsList,
                integratorSettings, propagatorSettings );

    std::vector< TimeType > baseTimeList;
    double observationTimeStart = initialEphemerisTime + 1000.0;
    double  observationInterval = 20.0;
    for( int i = 0; i < numberOfDaysOfData; i++ )
    {
        for( unsigned int j = 0; j < 500; j++ )
        {
            baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 86400.0 +
                                    static_cast< double >( j ) * observationInterval );
        }
    }

    std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > measurementSimulationInput =
            getObservationSimulationSettings< TimeType >(
                linkEndsPerObservable, baseTimeList, receiver );

    // Simulate observations
    std::shared_ptr< ObservationCollection< StateScalarType, TimeType > > simulatedObservations =
            simulateObservations< StateScalarType, TimeType >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Perturb parameter estimate
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > parameterPerturbation =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );

    if( numberOfIterations > 0 )
    {
        parameterPerturbation.segment( 0, 3 ) = Eigen::Vector3d::Constant( 1.0 );
        parameterPerturbation.segment( 3, 3 ) = Eigen::Vector3d::Constant( 1.E-3 );

        if( useFullParameterSet )
        {
            parameterPerturbation( 6 ) = 0.05;
            parameterPerturbation( 7 ) = 0.05;
        }
        initialParameterEstimate += parameterPerturbation;
    }


    // Define estimation input
    std::shared_ptr< PodInput< StateScalarType, TimeType  > > podInput =
            std::make_shared< PodInput< StateScalarType, TimeType > >(
                simulatedObservations, initialParameterEstimate.rows( ),
                Eigen::MatrixXd::Zero( truthParameters.rows( ), truthParameters.rows( ) ),
                initialParameterEstimate - truthParameters );

    std::map< observation_models::ObservableType, double > weightPerObservable;
    weightPerObservable[ one_way_range ] = 1.0 / ( 1.0 * 1.0 );
    weightPerObservable[ angular_position ] = 1.0 / ( 1.0E-5 * 1.0E-5 );
    weightPerObservable[ one_way_doppler ] = 1.0 / ( 1.0E-11 * 1.0E-11 );

    podInput->setConstantPerObservableWeightsMatrix( weightPerObservable );
    podInput->defineEstimationSettings( true, true, true, true, false );

    // Perform estimation
    std::shared_ptr< PodOutput< StateScalarType > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( numberOfIterations ) );

    Eigen::VectorXd estimationError = podOutput->parameterEstimate_ - truthParameters;
    std::cout <<"estimation error: "<< ( estimationError ).transpose( ) << std::endl;

    podData = std::make_pair( podOutput, podInput );

    return estimationError;
}

extern template Eigen::VectorXd executeEarthOrbiterParameterEstimation< double, double >(
        std::pair< std::shared_ptr< PodOutput< double > >, std::shared_ptr< PodInput< double, double > > >& podData,
        const double startTime,
        const int numberOfDaysOfData,
        const int numberOfIterations,
        const bool useFullParameterSet );




//! Test the estimation of observation biases.
template< typename TimeType = double, typename StateScalarType  = double >
std::pair< Eigen::VectorXd, bool > executeEarthOrbiterBiasEstimation(
        const bool estimateRangeBiases = true,
        const bool estimateTwoWayBiases = false,
        const bool useSingleBiasModel = true,
        const bool estimateAbsoluteBiases = true,
        const bool omitRangeData = false,
        const bool useMultiArcBiases = false )
{

    const int numberOfDaysOfData = 1;
    const int numberOfIterations = 3;

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );

    // Specify initial time
    TimeType initialEphemerisTime = 1.0E7;
    TimeType finalEphemerisTime = initialEphemerisTime + numberOfDaysOfData * 86400.0;

    std::vector< double > biasArcs;
    biasArcs.push_back( initialEphemerisTime );
    biasArcs.push_back( initialEphemerisTime + 4.0 * 3600.0 );
    biasArcs.push_back( initialEphemerisTime + 12.0 * 3600.0 );

    std::vector< Eigen::VectorXd > biasPerArc;
    biasPerArc.push_back( Eigen::Vector1d::Zero( ) );
    biasPerArc.push_back( Eigen::Vector1d::Zero( ) );
    biasPerArc.push_back( Eigen::Vector1d::Zero( ) );

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setEphemeris( std::make_shared< TabulatedCartesianEphemeris< > >(
                                            std::shared_ptr< interpolators::OneDimensionalInterpolator
                                            < double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );


    // Creatre ground stations: same position, but different representation
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );

    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ),
                         geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ),
                         geodetic_position );

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );
    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::Matrix< StateScalarType, 6, 1 > systemInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState,
              TimeType( finalEphemerisTime ), cowell );

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettings< TimeType > >
            ( TimeType( initialEphemerisTime ), 120.0,
              RungeKuttaCoefficients::CoefficientSets::rungeKuttaFehlberg78,
              120.0, 120.0, 1.0, 1.0 );


    // Define parameters.
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;
    std::vector< LinkEnds > stationTwoWayLinkEnds;
    std::vector< LinkEnds > stationTwoWayInverseLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ receiver ] = std::make_pair( "Vehicle", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = std::make_pair( "Vehicle", "" );
        stationReceiverLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ reflector1 ] = std::make_pair( "Vehicle", "" );
        linkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        stationTwoWayLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = std::make_pair( "Vehicle", "" );
        linkEnds[ reflector1 ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = std::make_pair( "Vehicle", "" );
        stationTwoWayInverseLinkEnds.push_back( linkEnds );
    }

    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    if( estimateRangeBiases )
    {
        linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
        linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
        linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );
        linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 1 ] );

        linkEndsPerObservable[ n_way_range ].push_back( stationTwoWayLinkEnds[ 0 ] );
        linkEndsPerObservable[ n_way_range ].push_back( stationTwoWayInverseLinkEnds[ 0 ] );
        linkEndsPerObservable[ n_way_range ].push_back( stationTwoWayLinkEnds[ 1 ] );
        linkEndsPerObservable[ n_way_range ].push_back( stationTwoWayInverseLinkEnds[ 1 ] );
    }

    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 1 ] );

    linkEndsPerObservable[ two_way_doppler ].push_back( stationTwoWayLinkEnds[ 0 ] );
    linkEndsPerObservable[ two_way_doppler ].push_back( stationTwoWayInverseLinkEnds[ 0 ] );
    linkEndsPerObservable[ two_way_doppler ].push_back( stationTwoWayLinkEnds[ 1 ] );
    linkEndsPerObservable[ two_way_doppler ].push_back( stationTwoWayInverseLinkEnds[ 1 ] );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back(
                std::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                    "Vehicle", systemInitialState, "Earth" ) );

    if( useMultiArcBiases )
    {
        if( estimateRangeBiases )
        {
            parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                                          linkEndsPerObservable.at( one_way_range ).at( 0 ), one_way_range, biasArcs, transmitter, estimateAbsoluteBiases ) );
            parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                                          linkEndsPerObservable.at( one_way_range ).at( 1 ), one_way_range, biasArcs, transmitter, estimateAbsoluteBiases ) );
            parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                                          linkEndsPerObservable.at( one_way_range ).at( 2 ), one_way_range, biasArcs, transmitter, estimateAbsoluteBiases ) );
            parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                                          linkEndsPerObservable.at( one_way_range ).at( 3 ), one_way_range, biasArcs, transmitter, estimateAbsoluteBiases ) );
            if( estimateTwoWayBiases )
            {
                parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                                              linkEndsPerObservable.at( n_way_range ).at( 0 ), n_way_range, biasArcs, transmitter, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                                              linkEndsPerObservable.at( n_way_range ).at( 1 ), n_way_range, biasArcs, transmitter, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                                              linkEndsPerObservable.at( n_way_range ).at( 2 ), n_way_range, biasArcs, transmitter, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                                              linkEndsPerObservable.at( n_way_range ).at( 3 ), n_way_range, biasArcs, transmitter, estimateAbsoluteBiases ) );
            }
        }
        else
        {
            parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                                          linkEndsPerObservable.at( one_way_doppler ).at( 0 ), one_way_doppler, biasArcs, transmitter, estimateAbsoluteBiases ) );
            parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                                          linkEndsPerObservable.at( one_way_doppler ).at( 1 ), one_way_doppler, biasArcs, transmitter, estimateAbsoluteBiases ) );
            parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                                          linkEndsPerObservable.at( one_way_doppler ).at( 2 ), one_way_doppler, biasArcs, transmitter, estimateAbsoluteBiases ) );
            parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                                          linkEndsPerObservable.at( one_way_doppler ).at( 3 ), one_way_doppler, biasArcs, transmitter, estimateAbsoluteBiases ) );
            if( estimateTwoWayBiases )
            {
                parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                                              linkEndsPerObservable.at( two_way_doppler ).at( 0 ), two_way_doppler, biasArcs, transmitter, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                                              linkEndsPerObservable.at( two_way_doppler ).at( 1 ), two_way_doppler, biasArcs, transmitter, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                                              linkEndsPerObservable.at( two_way_doppler ).at( 2 ), two_way_doppler, biasArcs, transmitter, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                                              linkEndsPerObservable.at( two_way_doppler ).at( 3 ), two_way_doppler, biasArcs, transmitter, estimateAbsoluteBiases ) );
            }
        }
    }
    else
    {
        if( estimateRangeBiases )
        {
            parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                          linkEndsPerObservable.at( one_way_range ).at( 0 ), one_way_range, estimateAbsoluteBiases ) );
            parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                          linkEndsPerObservable.at( one_way_range ).at( 1 ), one_way_range, estimateAbsoluteBiases ) );
            parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                          linkEndsPerObservable.at( one_way_range ).at( 2 ), one_way_range, estimateAbsoluteBiases ) );
            parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                          linkEndsPerObservable.at( one_way_range ).at( 3 ), one_way_range, estimateAbsoluteBiases ) );
            if( estimateTwoWayBiases )
            {
                parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                              linkEndsPerObservable.at( n_way_range ).at( 0 ), n_way_range, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                              linkEndsPerObservable.at( n_way_range ).at( 1 ), n_way_range, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                              linkEndsPerObservable.at( n_way_range ).at( 2 ), n_way_range, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                              linkEndsPerObservable.at( n_way_range ).at( 3 ), n_way_range, estimateAbsoluteBiases ) );
            }
        }
        else
        {
            parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                          linkEndsPerObservable.at( one_way_doppler ).at( 0 ), one_way_doppler, estimateAbsoluteBiases ) );
            parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                          linkEndsPerObservable.at( one_way_doppler ).at( 1 ), one_way_doppler, estimateAbsoluteBiases ) );
            parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                          linkEndsPerObservable.at( one_way_doppler ).at( 2 ), one_way_doppler, estimateAbsoluteBiases ) );
            parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                          linkEndsPerObservable.at( one_way_doppler ).at( 3 ), one_way_doppler, estimateAbsoluteBiases ) );
            if( estimateTwoWayBiases )
            {
                parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                              linkEndsPerObservable.at( two_way_doppler ).at( 0 ), two_way_doppler, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                              linkEndsPerObservable.at( two_way_doppler ).at( 1 ), two_way_doppler, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                              linkEndsPerObservable.at( two_way_doppler ).at( 2 ), two_way_doppler, estimateAbsoluteBiases ) );
                parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
                                              linkEndsPerObservable.at( two_way_doppler ).at( 3 ), two_way_doppler, estimateAbsoluteBiases ) );
            }
        }
    }



    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate( parameterNames, bodies );

    printEstimatableParameterEntries( parametersToEstimate );

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            std::shared_ptr< ObservationBiasSettings > biasSettings;
            if( useMultiArcBiases )
            {
                if( useSingleBiasModel )
                {
                    if( estimateAbsoluteBiases )
                    {
                        biasSettings = std::make_shared< ArcWiseConstantObservationBiasSettings >(
                                    biasArcs, biasPerArc, transmitter, true );
                    }
                    else
                    {
                        biasSettings = std::make_shared< ArcWiseConstantObservationBiasSettings >(
                                    biasArcs, biasPerArc, transmitter, false );
                    }
                }
                else
                {
                    std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;

                    biasSettingsList.push_back( std::make_shared< ArcWiseConstantObservationBiasSettings >(
                                                    biasArcs, biasPerArc, transmitter, true ) );
                    biasSettingsList.push_back( std::make_shared< ArcWiseConstantObservationBiasSettings >(
                                                    biasArcs, biasPerArc, transmitter, false ) );
                    biasSettings = std::make_shared< MultipleObservationBiasSettings >(
                                biasSettingsList );
                }
            }
            else
            {
                if( useSingleBiasModel )
                {
                    if( estimateAbsoluteBiases )
                    {
                        biasSettings = std::make_shared< ConstantObservationBiasSettings >(
                                    Eigen::Vector1d::Zero( ), true );
                    }
                    else
                    {
                        biasSettings = std::make_shared< ConstantObservationBiasSettings >(
                                    Eigen::Vector1d::Zero( ), false );
                    }
                }
                else
                {
                    std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;

                    biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >(
                                                    Eigen::Vector1d::Zero( ), true ) );
                    biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >(
                                                    Eigen::Vector1d::Zero( ), false ) );
                    biasSettings = std::make_shared< MultipleObservationBiasSettings >(
                                biasSettingsList );
                }
            }
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                            currentObservable, currentLinkEndsList.at( i ), std::shared_ptr< LightTimeCorrectionSettings >( ),
                                            biasSettings ) );
        }
    }

    // Create orbit determination object.
    OrbitDeterminationManager< StateScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< StateScalarType, TimeType >(
                bodies, parametersToEstimate, observationSettingsList,
                integratorSettings, propagatorSettings );

    std::vector< TimeType > baseTimeList;
    double observationTimeStart = initialEphemerisTime + 600.0;
    double  observationInterval = 600.0;
    for( int i = 0; i < numberOfDaysOfData; i++ )
    {
        for( unsigned int j = 0; j < 100; j++ )
        {
            baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 86400.0 +
                                    static_cast< double >( j ) * observationInterval );
        }
    }

    std::cout<<"Final time: "<<baseTimeList.at( baseTimeList.size( ) - 1 )<<std::endl;
    std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        if( !( omitRangeData && ( ( currentObservable == one_way_range ) ||
                                  ( currentObservable == n_way_range ) ) ) )
        {
            std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
            for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
            {
                measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< TimeType > >(
                            currentObservable, currentLinkEndsList.at( i ), baseTimeList, receiver ) );
            }
        }
    }


    // Simulate observations
    std::shared_ptr< ObservationCollection< StateScalarType, TimeType > > simulatedObservations = simulateObservations< StateScalarType, TimeType >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Perturb parameter estimate
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > parameterPerturbation =
            Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );

    if( numberOfIterations > 0 )
    {
        parameterPerturbation.segment( 0, 3 ) = Eigen::Vector3d::Constant( 1.0 );
        parameterPerturbation.segment( 3, 3 ) = Eigen::Vector3d::Constant( 1.E-3 );

        for( unsigned int i = 6; i < initialParameterEstimate.rows( ); i++ )
        {
            if( estimateAbsoluteBiases )
            {
                if( estimateRangeBiases )
                {
                    parameterPerturbation( i ) = 1.0;
                }
                else
                {
                    parameterPerturbation( i ) = 1.0E-8;
                }
            }
            else
            {
                parameterPerturbation( i ) = 1.0E-6;
            }
        }
        initialParameterEstimate += parameterPerturbation;
    }


    // Define estimation input
    std::shared_ptr< PodInput< StateScalarType, TimeType  > > podInput =
            std::make_shared< PodInput< StateScalarType, TimeType > >(
                simulatedObservations, initialParameterEstimate.rows( ),
                Eigen::MatrixXd::Zero( truthParameters.rows( ), truthParameters.rows( ) ),
                initialParameterEstimate - truthParameters );

    std::map< observation_models::ObservableType, double > weightPerObservable;
    weightPerObservable[ one_way_range ] = 1.0 / ( 1.0 * 1.0 );
    weightPerObservable[ n_way_range ] = 1.0 / ( 1.0 * 1.0 );

    weightPerObservable[ one_way_doppler ] = 1.0 / ( 1.0E-12 * 1.0E-12 );
    weightPerObservable[ two_way_doppler ] = 1.0 / ( 1.0E-12 * 1.0E-12 );

    podInput->setConstantPerObservableWeightsMatrix( weightPerObservable );
    podInput->defineEstimationSettings( true, false, false, true, false );

    // Perform estimation
    std::shared_ptr< PodOutput< StateScalarType > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( numberOfIterations ) );

    Eigen::VectorXd estimationError = podOutput->parameterEstimate_ - truthParameters;
    std::cout <<"estimation error: "<< ( estimationError ).transpose( ) << std::endl<< std::endl;

    return std::make_pair( estimationError,
                           ( podOutput->exceptionDuringInversion_ ||
                             !( podOutput->getUnnormalizedCovarianceMatrix( ) == podOutput->getUnnormalizedCovarianceMatrix( ) ) ) );
}

extern template std::pair< Eigen::VectorXd, bool > executeEarthOrbiterBiasEstimation< double, double >(
        const bool estimateRangeBiases,
        const bool estimateTwoWayBiases,
        const bool useSingleBiasModel,
        const bool estimateAbsoluteBiases,
        const bool omitRangeData,
        const bool useMultiArcBiases );

}

}
#endif // ORBITDETERMINATIONTESTCASES_H
