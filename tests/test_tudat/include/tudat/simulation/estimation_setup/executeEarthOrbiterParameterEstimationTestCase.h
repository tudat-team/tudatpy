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

#ifndef EXECUTEEARTHORBITERPARAMETERESTIMATIONTESTCASE_H
#define EXECUTEEARTHORBITERPARAMETERESTIMATIONTESTCASE_H

#include "tudat/simulation/estimation_setup/orbitDeterminationTestCaseUtilities.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"

namespace tudat
{
namespace unit_tests
{

using namespace ephemerides;
using namespace basic_astrodynamics;
using namespace simulation_setup;
using namespace orbital_element_conversions;
using namespace coordinate_conversions;
using namespace physical_constants;

template< typename TimeType = double, typename StateScalarType = double >
Eigen::VectorXd executeEarthOrbiterParameterEstimation(
        std::pair< std::shared_ptr< EstimationOutput< StateScalarType > >,
                   std::shared_ptr< EstimationInput< StateScalarType, TimeType > > >& podData,
        const TimeType startTime = TimeType( 1.0E7 ),
        const int numberOfDaysOfData = 3,
        const int numberOfIterations = 5,
        const bool useFullParameterSet = true,
        const bool saveDesignMatrix = true )
{
    // Load spice kernels.
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
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
            "ECLIPJ2000",
            "IAU_Earth",
            spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
            initialEphemerisTime,
            2.0 * mathematical_constants::PI / ( physical_constants::JULIAN_DAY ) );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 4.0;
    double aerodynamicCoefficient = 1.2;
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                    referenceArea,
                    aerodynamicCoefficient * ( Eigen::Vector3d( ) << 1.2, -0.1, -0.4 ).finished( ),
                    negative_aerodynamic_frame_coefficients );

    // Create and set aerodynamic coefficients object
    bodies.at( "Vehicle" )
            ->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle", bodies ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodies.at( "Vehicle" )
            ->setRadiationPressureInterface( "Sun",
                                             createRadiationPressureInterface( asterixRadiationPressureSettings, "Vehicle", bodies ) );

    bodies.at( "Vehicle" )
            ->setEphemeris( std::make_shared< TabulatedCartesianEphemeris<> >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );

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
    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Sun" ].push_back(
            std::make_shared< AccelerationSettings >( basic_astrodynamics::cannon_ball_radiation_pressure ) );
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( basic_astrodynamics::aerodynamic ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Set (perturbed) initial state.
    Eigen::Matrix< StateScalarType, 6, 1 > systemInitialState =
            convertKeplerianToCartesianElements( asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType, TimeType > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< StateScalarType, TimeType > >(
                    centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, TimeType( finalEphemerisTime ), cowell );

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            std::make_shared< RungeKuttaVariableStepSizeSettings< TimeType > >(
                    TimeType( initialEphemerisTime ), 40.0, CoefficientSets::rungeKuttaFehlberg78, 40.0, 40.0, 1.0, 1.0 );

    // Define parameters.
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = LinkEndId( "Earth", groundStationNames.at( i ) );
        linkEnds[ receiver ] = LinkEndId( "Vehicle", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = LinkEndId( "Vehicle", "" );
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
    parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
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
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", rotation_pole_position ) );
        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Earth", ground_station_position, "Station1" ) );
    }

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate< StateScalarType, TimeType >( parameterNames, bodies );

    printEstimatableParameterEntries( parametersToEstimate );

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;

    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( );
         linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            std::shared_ptr< ObservationBiasSettings > biasSettings;
            if( ( currentObservable == one_way_range ) && ( i == 0 ) )
            {
                std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;

                biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d::Zero( ), true ) );
                biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d::Zero( ), false ) );
                biasSettings = std::make_shared< MultipleObservationBiasSettings >( biasSettingsList );
            }
            else if( ( currentObservable == one_way_range ) && ( i == 1 ) )
            {
                biasSettings = std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d::Zero( ), false );
            }

            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                    currentObservable, currentLinkEndsList.at( i ), std::shared_ptr< LightTimeCorrectionSettings >( ), biasSettings ) );
        }
    }

    // Create orbit determination object.
    OrbitDeterminationManager< StateScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< StateScalarType, TimeType >(
                    bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );

    std::vector< TimeType > baseTimeList;
    double observationTimeStart = initialEphemerisTime + 1000.0;
    double observationInterval = 20.0;
    for( int i = 0; i < numberOfDaysOfData; i++ )
    {
        for( unsigned int j = 0; j < 500; j++ )
        {
            baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 86400.0 +
                                    static_cast< double >( j ) * observationInterval );
        }
    }

    std::vector< std::shared_ptr< ObservationSimulationSettings< TimeType > > > measurementSimulationInput =
            getObservationSimulationSettings< TimeType >( linkEndsPerObservable, baseTimeList, receiver );

    // Simulate observations
    std::shared_ptr< ObservationDataset< StateScalarType, TimeType > > simulatedObservations =
            simulateObservationDataset< StateScalarType, TimeType >(
                    measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    for( int setId = 0; setId < simulatedObservations->getNumberOfObservationSets( ); ++setId )
    {
        const ObservableType currentObservable = simulatedObservations->getObservationSetMetadata( setId ).observableType_;
        if( currentObservable == one_way_range )
        {
            simulatedObservations->setConstantWeightForSet( setId, 1.0 / ( 1.0 * 1.0 ) );
        }
        else if( currentObservable == angular_position )
        {
            simulatedObservations->setConstantWeightForSet( setId, 1.0 / ( 1.0E-5 * 1.0E-5 ) );
        }
        else if( currentObservable == one_way_doppler )
        {
            simulatedObservations->setConstantWeightForSet( setId, 1.0 / ( 1.0E-11 * 1.0E-11 * SPEED_OF_LIGHT * SPEED_OF_LIGHT ) );
        }
    }

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
            //            parameterPerturbation( 7 ) = 0.05;
        }
        initialParameterEstimate += parameterPerturbation;
    }
    parametersToEstimate->resetParameterValues( initialParameterEstimate );

    // Define estimation input
    std::shared_ptr< EstimationInput< StateScalarType, TimeType > > estimationInput =
            std::make_shared< EstimationInput< StateScalarType, TimeType > >( simulatedObservations );
    estimationInput->defineEstimationSettings( true, true, saveDesignMatrix, true, true, true );
    estimationInput->setConvergenceChecker( std::make_shared< EstimationConvergenceChecker >( numberOfIterations ) );

    // Perform estimation
    std::shared_ptr< EstimationOutput< StateScalarType > > estimationOutput =
            orbitDeterminationManager.estimateParameters( estimationInput );

    Eigen::VectorXd estimationError = estimationOutput->parameterEstimate_ - truthParameters;
    std::cout << "estimation error: " << ( estimationError ).transpose( ) << std::endl;

    podData = std::make_pair( estimationOutput, estimationInput );

    return estimationError;
}

}  // namespace unit_tests

}  // namespace tudat

#endif  // EXECUTEEARTHORBITERPARAMETERESTIMATIONTESTCASE_H
