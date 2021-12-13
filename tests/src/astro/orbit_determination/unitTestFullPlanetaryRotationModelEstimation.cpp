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

#include <string>
#include <thread>

#include <limits>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/periodicSpinVariation.h"


namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_full_planetary_rotational_parameters_estimation )

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
using namespace tudat;


//! Unit test to check if periodic spin variation, polar motion, and free-core factor/ampliture
//! (for a full planetary rotational model) are estimated correctly. Translatuonal state estimation is included for interface
//! consistency only
BOOST_AUTO_TEST_CASE( test_FullPlanetaryRotationalParameters )
{    
    //Load spice kernels.
    std::string kernelsPath = paths::getSpiceKernelPath( );
    spice_interface::loadStandardSpiceKernels( );
    
    //Define environment settings
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    
    // Specify total time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = 1000.0 * 86400.0;
    double maximumTimeStep = 86400.0;
    double buffer = 10.0 * maximumTimeStep;
    // Set-up different cases with various parameters to estimate.
    for ( int testCase = 0 ; testCase < 3 ; testCase++ )
    {
        // Create body objects; Mars with high-accuracy rotation model
        BodyListSettings bodySettings =
                getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
        bodySettings.at( "Mars" )->rotationModelSettings = getHighAccuracyMarsRotationModel( );
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );


        // Create ground stations
        std::pair< std::string, std::string > grazStation = std::pair< std::string, std::string >( "Earth", "" );
        std::pair< std::string, std::string > mslStation = std::pair< std::string, std::string >( "Mars", "MarsStation" );
        createGroundStation( bodies.at( "Mars" ), "MarsStation", ( Eigen::Vector3d( ) << 100.0, 0.5, 2.1 ).finished( ),
                             coordinate_conversions::geodetic_position );

        // Set accelerations between bodies that are to be taken into account.
        SelectedAccelerationMap accelerationMap;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
        accelerationsOfEarth[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
        accelerationsOfEarth[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
        accelerationMap[ "Earth" ] = accelerationsOfEarth;

        // Set bodies for which initial state is to be estimated and integrated.
        std::vector< std::string > bodiesToEstimate;
        bodiesToEstimate.push_back( "Earth" );
        std::vector< std::string > bodiesToIntegrate;
        bodiesToIntegrate.push_back( "Earth" );

        // Define propagator settings.
        std::vector< std::string > centralBodies; centralBodies.push_back( "SSB" );
        AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );
        Eigen::VectorXd initialState = getInitialStateOfBody< double, double>(
                    bodiesToIntegrate.at( 0 ), centralBodies.at( 0 ), bodies, initialEphemerisTime );
        std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToIntegrate, initialState,
                  finalEphemerisTime );

        // Define integrator settings.
        std::shared_ptr< IntegratorSettings< double > > integratorSettings = std::make_shared< IntegratorSettings< double > >
                ( rungeKutta4, initialEphemerisTime, maximumTimeStep );

        // Define links in simulation.
        std::vector< LinkEnds > linkEnds; linkEnds.resize( 1 );
        linkEnds[ 0 ][ transmitter ] = grazStation;
        linkEnds[ 0 ][ receiver ] = mslStation;

        // Define observation model settings
        std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( one_way_range, linkEnds[ 0 ] ) );

        // Define observation times.
        double observationTime;
        double timeBuffer = maximumTimeStep * 5.0;
        double observationCadence = 3600.0;
        std::vector< double > observationTimes;
        observationTime = initialEphemerisTime + timeBuffer;
        while( observationTime < finalEphemerisTime - timeBuffer )
        {
            observationTimes.push_back( observationTime );
            observationTime += observationCadence;
        }

        // Create observation simulation settings
        std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
        measurementSimulationInput.push_back(
                    std::make_shared< TabulatedObservationSimulationSettings< > >(
                        one_way_range, linkEnds[ 0 ], observationTimes, receiver ) );

        // Set parameters that are to be estimated.
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
                getInitialStateParameterSettings< double >( propagatorSettings, bodies );

        // Estimate core factor and free core nutation rate
        if ( testCase == 0 )
        {
            parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Mars", core_factor ) );
            parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Mars", free_core_nutation_rate ) );
        }

        // Estimate periodic spin variation
        else if ( testCase == 1 )
        {
            parameterNames.push_back(  std::make_shared< EstimatableParameterSettings >( "Mars", periodic_spin_variation ) );
        }

        // Estimate polar motion amplitude
        else if ( testCase == 2 )
        {

            parameterNames.push_back(  std::make_shared< EstimatableParameterSettings >( "Mars", polar_motion_amplitude ) );
        }

        // Create parameters
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double >( parameterNames, bodies );
        printEstimatableParameterEntries( parametersToEstimate );

        // Create orbit determination object.
        OrbitDeterminationManager< double, double > orbitDeterminationManager = OrbitDeterminationManager< double, double >(
                    bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );
        
        // Define initial parameter estimate.
        Eigen::VectorXd initialParameterEstimate =
                parametersToEstimate->template getFullParameterValues< double >( );
        
        // Simulate observations
        std::shared_ptr< ObservationCollection< > > observationsAndTimes = simulateObservations< double, double >(
                    measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );
        
        
        // Define perturbation of parameter estimate
        Eigen::VectorXd truthParameters = initialParameterEstimate;
        if ( testCase == 0 )
        {
            initialParameterEstimate[ 6 ] += 1.0E-4;
            initialParameterEstimate[ 7 ] += 1.0E-10;
        }
        else if ( testCase == 1 )
        {
            for( int i = 6 + 0 ; i < static_cast< int >( initialParameterEstimate.rows( ) ) ; i++ )
            {
                initialParameterEstimate[ i ] += 1.0E-6;
            }
        }
        else if ( testCase == 2 )
        {
            for( int i = 6 ; i < static_cast< int >( initialParameterEstimate.rows( ) ) ; i++ )
            {
                initialParameterEstimate[ i ] += 1.0E-6;
            }
        }
        parametersToEstimate->resetParameterValues( initialParameterEstimate );

        // Set strong a priori covariance for translational state (test only rotational variation estimation)
        int numberOfParameters = initialParameterEstimate.rows( );
        Eigen::MatrixXd inverseAprioriCovariance = Eigen::MatrixXd::Zero( numberOfParameters, numberOfParameters );
        for( int i = 0; i < 3; i++ )
        {
            inverseAprioriCovariance( i, i ) = 1.0 / ( 1.0E-3 * 1.0E-3 );
            inverseAprioriCovariance( i + 3, i + 3 ) = 1.0 / ( 1.0E-6 * 1.0E-6 );

        }

        // Create estimation input
        std::shared_ptr< PodInput< double, double > > podInput = std::make_shared< PodInput< double, double > >(
                    observationsAndTimes, numberOfParameters, inverseAprioriCovariance );
        podInput->defineEstimationSettings( false, false );

        // Perform state estimation
        std::shared_ptr< PodOutput< double, double > > podOutput = orbitDeterminationManager.estimateParameters(
                    podInput, std::make_shared< EstimationConvergenceChecker >( 3 ) );
        
        
        // Retrieve estimated parameter, and compare against true values
        Eigen::VectorXd parameterError = podOutput->parameterEstimate_ - truthParameters;
        if ( testCase == 0 ) {
            
            BOOST_CHECK_SMALL( std::fabs( parameterError( 6 ) ), 1.0E-7 );
            BOOST_CHECK_SMALL( std::fabs( parameterError( 6 + 1 ) ), 1.0E-12 );
            
        }
        else
        {
            for( int i = 6 ; i < static_cast< int >( initialParameterEstimate.rows( ) ); i++ )
            {
                if ( testCase == 1 )
                {
                    BOOST_CHECK_SMALL( std::fabs( parameterError( i ) ), 1.0E-12 );
                }
                else if ( testCase == 2 )
                {
                    BOOST_CHECK_SMALL( std::fabs( parameterError( i ) ), 1.0E-12 );
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
