/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
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
#include <boost/test/unit_test.hpp>
#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationTestCases.h"


namespace tudat
{
namespace unit_tests
{

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
using namespace tudat::physical_constants;
using namespace tudat::ground_stations;


BOOST_AUTO_TEST_SUITE( test_time_bias_partials )

//! Test partial derivatives of angular position observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testTimeBiasPartials )
//int main( )
{
    const int numberOfDaysOfData = 1;
    int numberOfIterations = 10;

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

//    for( int stepIndex = -2; stepIndex < 2; stepIndex++ )
    {
//        std::cout<<"STEP INDEX "<<stepIndex<<" ================================"<<std::endl;

        // Define bodies in simulation
        std::vector< std::string > bodyNames;
        bodyNames.push_back( "Earth" );

        // Specify initial time
        double initialTime = 1.0e7;
        double finalTime = initialTime + numberOfDaysOfData * 86400.0;

        std::vector< double > arcs;
        arcs.push_back( initialTime );
        arcs.push_back( initialTime + 2100.0 );
        arcs.push_back( initialTime + 4100.0 );

        std::vector< double > timeBiasesPerArc = { 0.0, 0.0, 0.0 };

        std::vector< Eigen::VectorXd > biasesPerArc;
        biasesPerArc.push_back( Eigen::Vector1d::Zero( ) );
        biasesPerArc.push_back( Eigen::Vector1d::Zero( ) );
        biasesPerArc.push_back( Eigen::Vector1d::Zero( ) );

        // Create bodies needed in simulation
        BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );
        bodies.createEmptyBody( "Vehicle" );
        bodies.at( "Vehicle" )->setEphemeris( std::make_shared< TabulatedCartesianEphemeris< > >(
            std::shared_ptr< interpolators::OneDimensionalInterpolator < double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );
        bodies.createEmptyBody( "Vehicle2" );
        bodies.at( "Vehicle2" )->setEphemeris( std::make_shared< TabulatedCartesianEphemeris< > >(
            std::shared_ptr< interpolators::OneDimensionalInterpolator < double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );

        std::shared_ptr< system_models::VehicleSystems > vehicleSystems = std::make_shared< system_models::VehicleSystems >( );
        vehicleSystems->setTransponderTurnaroundRatio( &getDsnDefaultTurnaroundRatios );
        bodies.at( "Vehicle" )->setVehicleSystems( vehicleSystems );

        // Create ground station
        createGroundStation( bodies.at( "Earth" ), "Station", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
        createGroundStation( bodies.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.5, 0.8 ).finished( ), geodetic_position );

        bodies.at( "Earth" )->getGroundStation( "Station" )->setTransmittingFrequencyCalculator( std::make_shared< ConstantFrequencyInterpolator >( 7E9 ) );
        bodies.at( "Earth" )->getGroundStation( "Station2" )->setTransmittingFrequencyCalculator( std::make_shared< ConstantFrequencyInterpolator >( 7E9 ) );

        // Set accelerations on Vehicle that are to be taken into account.
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
        accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
        accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
        accelerationMap[ "Vehicle2" ] = accelerationsOfVehicle;

        // Set bodies for which initial state is to be estimated and integrated.
        std::vector< std::string > bodiesToIntegrate;
        std::vector< std::string > centralBodies;
        bodiesToIntegrate.push_back( "Vehicle" );
        bodiesToIntegrate.push_back( "Vehicle2" );
        centralBodies.push_back( "Earth" );
        centralBodies.push_back( "Earth" );

        // Create acceleration models
        AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

        // Set Keplerian elements
        Eigen::Vector6d initialKeplerianState;
        initialKeplerianState( semiMajorAxisIndex ) = 7200.0E3;
        initialKeplerianState( eccentricityIndex ) = 0.05;
        initialKeplerianState( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
        initialKeplerianState( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
        initialKeplerianState( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
        initialKeplerianState( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );
        double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        Eigen::Vector6d initialState1 = convertKeplerianToCartesianElements( initialKeplerianState, earthGravitationalParameter );

        initialKeplerianState( trueAnomalyIndex ) -= unit_conversions::convertDegreesToRadians( 10.0 );
        Eigen::Vector6d initialState2 = convertKeplerianToCartesianElements( initialKeplerianState, earthGravitationalParameter );
        Eigen::VectorXd initialState = Eigen::VectorXd::Zero( 12 );
        initialState << initialState1, initialState2;

        // Create integrator settings
        std::shared_ptr< IntegratorSettings< double > > integratorSettings = rungeKuttaFixedStepSettings(
            120.0, CoefficientSets::rungeKuttaFehlberg78 );

        // Create propagator settings
        std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double, double > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate, initialState, initialTime, integratorSettings,
              propagationTimeTerminationSettings( finalTime ) );

        // Define link ends.
        LinkEnds oneWayLinkEnds;
        oneWayLinkEnds[ receiver ] = LinkEndId( "Earth", "Station" );
        oneWayLinkEnds[ transmitter ] = LinkEndId( "Vehicle", "" );

        LinkEnds twoWayLinkEnds;
        twoWayLinkEnds[ receiver ] = LinkEndId( "Earth", "Station" );
        twoWayLinkEnds[ retransmitter ] = LinkEndId( "Vehicle", "" );
        twoWayLinkEnds[ transmitter ] = LinkEndId( "Earth", "Station" );

        LinkEnds threeWayLinkEnds;
        threeWayLinkEnds[ receiver ] = LinkEndId( "Earth", "Station" );
        threeWayLinkEnds[ retransmitter ] = LinkEndId( "Vehicle", "" );
        threeWayLinkEnds[ transmitter ] = LinkEndId( "Earth", "Station2" );

        LinkEnds differencedTransmitterLinkEnds;
        differencedTransmitterLinkEnds[ receiver ] = LinkEndId( "Earth", "Station" );
        differencedTransmitterLinkEnds[ transmitter ] = LinkEndId( "Vehicle", "" );
        differencedTransmitterLinkEnds[ transmitter2 ] = LinkEndId( "Vehicle2", "" );


        LinkEnds differencedTransmitterLinkEndsInverse;
        differencedTransmitterLinkEndsInverse[ receiver ] = LinkEndId( "Earth", "Station" );
        differencedTransmitterLinkEndsInverse[ transmitter2 ] = LinkEndId( "Vehicle", "" );
        differencedTransmitterLinkEndsInverse[ transmitter ] = LinkEndId( "Vehicle2", "" );

        for( unsigned int observableCase = 0; observableCase < 11; observableCase++ )
        {
            std::cout<<"OBSERVABLE CASE "<<observableCase<<" ================================"<<std::endl;
            ObservableType currentObservable = undefined_observation_model;
            LinkEnds testLinkEnds = oneWayLinkEnds;
            LinkEndType referenceLinkEnd = receiver;

            switch( observableCase )
            {
            case 0:
                currentObservable = one_way_range;
                break;
            case 1:
                currentObservable = angular_position;
                break;
            case 2:
                currentObservable = one_way_differenced_range;
                break;
            case 3:
                currentObservable = relative_angular_position;
                testLinkEnds = differencedTransmitterLinkEnds;
                break;
            case 4:
                currentObservable = relative_angular_position;
                testLinkEnds = differencedTransmitterLinkEndsInverse;
                break;
            case 5:
                currentObservable = n_way_range;
                testLinkEnds = twoWayLinkEnds;
                break;
            case 6:
                currentObservable = n_way_range;
                testLinkEnds = threeWayLinkEnds;
                break;
            case 7:
                currentObservable = n_way_differenced_range;
                testLinkEnds = twoWayLinkEnds;
                break;
            case 8:
                currentObservable = n_way_differenced_range;
                testLinkEnds = threeWayLinkEnds;
                break;
            case 9:
                currentObservable = dsn_n_way_averaged_doppler;
                testLinkEnds = twoWayLinkEnds;
                break;
            case 10:
                currentObservable = dsn_n_way_averaged_doppler;
                testLinkEnds = threeWayLinkEnds;
                break;


    //        case 2:
    //            currentObservable = position_observable;
    //            break;
    //        case 3:
    //            currentObservable = velocity_observable;
    //            break;
            default:
                throw std::runtime_error( "Error when testing time bias partials; observable not implemented" );
            }
            for ( unsigned int testCase = 0 ; testCase < 2 ; testCase++ )
            {
                bool multiArcBiases = testCase;

                // Add bias parameters to estimation
                std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
                parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >( "Vehicle", initialState1, "Earth" ) );
                parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >( "Vehicle2", initialState2, "Earth" ) );

                if ( !multiArcBiases )
                {
                    parameterNames.push_back( std::make_shared< ConstantTimeBiasEstimatableParameterSettings >(
                        testLinkEnds, currentObservable, referenceLinkEnd ) );
                }
                else
                {
                    parameterNames.push_back( std::make_shared< ArcWiseTimeBiasEstimatableParameterSettings >(
                        testLinkEnds, currentObservable, arcs, referenceLinkEnd ) );
                }

                // Create parameters
                std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                    createParametersToEstimate< double, double >( parameterNames, bodies );
    //            printEstimatableParameterEntries( parametersToEstimate );

                // Define bias settings for observation model
                std::shared_ptr< ObservationBiasSettings > biasSettings;
                if ( !multiArcBiases )
                {
                    std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;
                    biasSettingsList.push_back( std::make_shared< ConstantTimeBiasSettings >( 0.0, referenceLinkEnd ) );
                    biasSettings = std::make_shared< MultipleObservationBiasSettings >( biasSettingsList );
                }
                else
                {
                    std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;
                    biasSettingsList.push_back( std::make_shared< ArcWiseTimeBiasSettings >( arcs, timeBiasesPerArc, referenceLinkEnd ) );
                    biasSettings = std::make_shared< MultipleObservationBiasSettings >( biasSettingsList );
                }
                // Define observation model settings
                std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
                if( currentObservable == n_way_differenced_range )
                {
                    observationSettingsList.push_back( std::make_shared<NWayDifferencedRangeObservationSettings>(
                        LinkDefinition( testLinkEnds ), std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ), biasSettings ) );
                }
                else if( currentObservable == dsn_n_way_averaged_doppler )
                {
                    observationSettingsList.push_back( std::make_shared<DsnNWayAveragedDopplerObservationSettings>(
                        LinkDefinition( testLinkEnds ), std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ), biasSettings ) );
                }
                else
                {
                    observationSettingsList.push_back( std::make_shared<ObservationModelSettings>(
                        currentObservable, testLinkEnds, std::shared_ptr<LightTimeCorrectionSettings>( ), biasSettings ));
                }

                // Create orbit determination object.
                OrbitDeterminationManager< double, double > orbitDeterminationManager = OrbitDeterminationManager< double, double >(
                    bodies, parametersToEstimate, observationSettingsList, propagatorSettings );

                // Define list of observation times
                std::vector< double > baseTimeList;
                double observationInterval = 600.0;
                int numberOfObservations = 10;
                for( int j = 0; j < numberOfObservations; j++ )
                {
                    baseTimeList.push_back( ( initialTime + 600.0 ) + static_cast< double >( j ) * observationInterval );
                }

                // Define observation simulation times

                std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
                measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
                    currentObservable, testLinkEnds, baseTimeList, referenceLinkEnd ) );

                std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings;
                if( currentObservable == one_way_differenced_range || currentObservable == n_way_differenced_range )
                {
                    ancilliarySettings = std::make_shared< ObservationAncilliarySimulationSettings >( );
                    ancilliarySettings->setAncilliaryDoubleData( doppler_integration_time, 60.0 );
                    measurementSimulationInput.at( 0 )->setAncilliarySettings( ancilliarySettings );
                }
                else if( currentObservable == dsn_n_way_averaged_doppler )
                {
                    ancilliarySettings = getDsnNWayAveragedDopplerAncillarySettings(
                        std::vector< FrequencyBands >{ x_band, x_band }, x_band, 8.5e9, 60.0 );
                    measurementSimulationInput.at( 0 )->setAncilliarySettings( ancilliarySettings );
                }

                // Compute observations and partials
                Eigen::VectorXd observations;
                Eigen::MatrixXd partials;
                switch( currentObservable )
                {
                case one_way_differenced_range:
                case one_way_range:
                case position_observable:
                case n_way_range:
                case n_way_differenced_range:
                case dsn_n_way_averaged_doppler:
                    orbitDeterminationManager.computePartialsAndObservations< 1 >(
                        testLinkEnds, ancilliarySettings, currentObservable, referenceLinkEnd, baseTimeList, observations, partials );
                    break;
                case angular_position:
                case relative_angular_position:
                    orbitDeterminationManager.computePartialsAndObservations< 2 >(
                        testLinkEnds, ancilliarySettings, currentObservable, referenceLinkEnd, baseTimeList, observations, partials );
                    break;
                default:
                    throw std::runtime_error( "Error when testing time bias partials; observable partial size not defined" );
                }

                // Get nominal parameter values
                Eigen::VectorXd originalParameters = parametersToEstimate->getFullParameterValues< double >( );

                for( int parameterIndex = 12; parameterIndex < parametersToEstimate->getParameterSetSize( ); parameterIndex++ )
                {
                    // Define perturbation for numerical partial
                    double timeBiasPerturbation = 1.0;//std::pow( 10, static_cast< double >( stepIndex ) / 2 );
                    if( observableCase > 6 )
                    {
                        timeBiasPerturbation *= 10.0;
                    }
//                    double timeBiasPerturbation = std::pow( 10, static_cast< double >( stepIndex ) );

                    // Perturb bias value up and recompute observations
                    Eigen::VectorXd perturbedParameters = originalParameters;
                    perturbedParameters( parameterIndex ) += timeBiasPerturbation;
                    parametersToEstimate->resetParameterValues( perturbedParameters );
                    std::shared_ptr< ObservationCollection< double, double > > upperturbedSimulatedObservations = simulateObservations< double, double >(
                        measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

                    // Perturb bias value down and recompute observations
                    perturbedParameters = originalParameters;
                    perturbedParameters( parameterIndex ) -= timeBiasPerturbation;
                    parametersToEstimate->resetParameterValues( perturbedParameters );
                    std::shared_ptr< ObservationCollection< double, double > > downperturbedSimulatedObservations = simulateObservations< double, double >(
                        measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

                    // Reset to nominal values
                    parametersToEstimate->resetParameterValues( perturbedParameters );

                    // Compute numerical partials
                    Eigen::VectorXd numericalPartials = ( upperturbedSimulatedObservations->getObservationVector( ) - downperturbedSimulatedObservations->getObservationVector( ) ) /
                        ( 2.0 * timeBiasPerturbation );
                    Eigen::VectorXd analyticalPartials =  partials.block( 0, parameterIndex, partials.rows( ), 1 );

                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( numericalPartials, analyticalPartials, 1.0E-3 );
//                    std::cout<<"Numerical "<<numericalPartials.transpose( )<<std::endl;
//                    std::cout<<"Analytical "<<analyticalPartials.transpose( )<<std::endl;
//                    std::cout<<"Relative error "<<( numericalPartials - analyticalPartials ).cwiseQuotient( analyticalPartials ).transpose( )<<std::endl;

                }
            }
        }
        std::cout<<std::endl;
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





