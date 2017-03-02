/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <string>
#include <thread>
#include <omp.h>
#include <limits>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/ObservationModels/simulateObservations.h"
#include "Tudat/Astrodynamics/OrbitDetermination/orbitDeterminationManager.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_estimation_from_positions )

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


template< typename TimeType = double, typename StateScalarType  = double >
Eigen::VectorXd  executeParameterEstimation(
        const int observableType )
{
    //Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

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

    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames,initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings[ "Moon" ]->ephemerisSettings->resetFrameOrigin( "Sun" );

    // Create bodies needed in simulation
    NamedBodyMap bodyMap = createBodies( bodySettings );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
    accelerationsOfEarth[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfEarth[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfEarth[ "Mars" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfEarth[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfEarth[ "Saturn" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );

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
                bodyMap, accelerationMap, centralBodyMap );

    // Set parameters that are to be estimated.
    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( boost::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                                  "Earth", propagators::getInitialStateOfBody< TimeType, StateScalarType >(
                                      "Earth", centralBodyMap[ "Earth" ], bodyMap, initialEphemerisTime ),
                              centralBodyMap[ "Earth" ] ) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Moon", gravitational_parameter ) );

    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate< StateScalarType >( parameterNames, bodyMap );


    // Define integrator settings.
    boost::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            boost::make_shared< IntegratorSettings< TimeType > >(
                rungeKutta4, TimeType( initialEphemerisTime - 4.0 * maximumTimeStep ), 900.0 );


    boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate,
              getInitialStateVectorOfBodiesToEstimate( parametersToEstimate ),
              TimeType( finalEphemerisTime + 4.0 * maximumTimeStep ),
              cowell, boost::shared_ptr< DependentVariableSaveSettings >( ), 600.0 );

    // Define light-time corrections.
    std::map< ObservableType, std::map< LinkEnds, std::vector< boost::shared_ptr< LightTimeCorrectionSettings > > > >
            observableCorrections;

    std::vector< LinkEnds > linkEnds;
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsMap;

    if(observableType == 0 )
    {
        linkEnds.resize( 1 );
        linkEnds[ 0 ][ observed_body ] = std::make_pair( "Earth", "" );
        linkEndsMap[ position_observable ] = linkEnds;
    }
    else
    {
        linkEnds.resize( 1 );
        linkEnds[ 0 ][ transmitter ] = std::make_pair( "Earth", "" );
        linkEnds[ 0 ][ receiver ] = std::make_pair( "Mars", "" );

        if( observableType == 1 )
        {
            linkEndsMap[ oneWayRange ] = linkEnds;
        }
        else if( observableType == 2 )
        {
            linkEndsMap[ angular_position ] = linkEnds;
        }
        else if( observableType == 3 )
        {
            linkEndsMap[ oneWayDoppler ] = linkEnds;
        }
    }



    // Create orbit determination object.
    OrbitDeterminationManager< StateScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< StateScalarType, TimeType >(
                bodyMap, parametersToEstimate,
                linkEndsMap, integratorSettings, propagatorSettings, observableCorrections );


    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );


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

    std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > > measurementSimulationInput;
    std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > singleObservableSimulationInput;


    initialObservationTimes = utilities::addScalarToVector( initialObservationTimes, 30.0 );
    if( observableType == 0 )
    {
        singleObservableSimulationInput[ linkEnds[ 0 ] ] = std::make_pair( initialObservationTimes, observed_body );
        measurementSimulationInput[ position_observable ] = singleObservableSimulationInput;
    }
    else
    {
        singleObservableSimulationInput[ linkEnds[ 0 ] ] = std::make_pair( initialObservationTimes, transmitter );

        if( observableType == 1 )
        {
            measurementSimulationInput[ oneWayRange ] = singleObservableSimulationInput;
        }
        else if( observableType == 2 )
        {
            measurementSimulationInput[ angular_position ] = singleObservableSimulationInput;
        }
        else if( observableType == 3 )
        {
            measurementSimulationInput[ oneWayDoppler ] = singleObservableSimulationInput;
        }
    }

    singleObservableSimulationInput.clear( );

    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< TimeType >, LinkEndType > > > SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

    PodInputDataType observationsAndTimes = simulateObservations< StateScalarType, TimeType >(
                measurementSimulationInput, orbitDeterminationManager.getObservationManagers( ) );


    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;

    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        initialParameterEstimate[ 0 + 6 * i ] += 1.0E3;
        initialParameterEstimate[ 1 + 6 * i ] += 1.0E3;
        initialParameterEstimate[ 2 + 6 * i ] += 1.0E3;
        initialParameterEstimate[ 3 + 6 * i ] += 1.0E-2;
        initialParameterEstimate[ 4 + 6 * i ] += 1.0E-2;
        initialParameterEstimate[ 5 + 6 * i ] += 1.0E-2;
    }
    initialParameterEstimate( 6 ) *= ( 1.0 + 1.0E-6 );


    boost::shared_ptr< PodInput< StateScalarType, TimeType > > podInput =
            boost::make_shared< PodInput< StateScalarType, TimeType > >(
                observationsAndTimes, initialParameterEstimate.rows( ), Eigen::MatrixXd::Zero( 0, 0 ),
                initialParameterEstimate - truthParameters );

    boost::shared_ptr< PodOutput< StateScalarType > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, boost::make_shared< EstimationConvergenceChecker >( ), true, true, false, false );

    return( podOutput->parameterEstimate_.template cast< double >( ) - truthParameters .template cast< double >( ) );
}

BOOST_AUTO_TEST_CASE( test_EstimationFromPosition )
{
    for( int simulationType = 0; simulationType < 4; simulationType++ )
    {
        for( unsigned int i = 0; i < 4; i++ )
        {
            std::cout<<"=============================================== Running Case: "<<i<<" "<<simulationType<<std::endl;

            Eigen::VectorXd totalError;
            if( i == 0 )
            {
                totalError = executeParameterEstimation< double, double >( simulationType );
            }
            else if( i == 1 )
            {
                totalError = executeParameterEstimation< double, long double >( simulationType );
            }
            else if( i == 2 )
            {
                totalError = executeParameterEstimation< Time, double >( simulationType );
            }
            else if( i == 3 )
            {
                totalError = executeParameterEstimation< Time, long double >( simulationType );
            }

            double toleranceMultiplier = 1.0;
            if( i % 2 == 1 )
            {
                toleranceMultiplier *= 1.0E-3;

                if( simulationType > 0 )
                {
                    toleranceMultiplier *= 100.0;
                }
            }
            else if( simulationType > 0 )
            {
                toleranceMultiplier *= 20.0;
            }

            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( totalError( j ), toleranceMultiplier * 5.0E-3 );
            }

            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( totalError( j + 3 ), toleranceMultiplier * 1.0E-7 );
            }

            BOOST_CHECK_SMALL( totalError( 6 ), toleranceMultiplier * 1.0E3 );
            std::cout<<totalError.transpose( )<<std::endl;
        }
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}


