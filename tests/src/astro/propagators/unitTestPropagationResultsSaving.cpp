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


#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/directTidalTimeLag.h"


namespace tudat
{
    namespace unit_tests
    {
        BOOST_AUTO_TEST_SUITE( test_arcwise_environment )

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


//! Unit test to check if tidal time lag parameters are estimated correctly
            BOOST_AUTO_TEST_CASE( test_ArcwiseEnvironmentParameters )
            {
                //Load spice kernels.
                spice_interface::loadStandardSpiceKernels( );

                // Define bodies in simulation
                std::vector< std::string > bodyNames;
                bodyNames.push_back( "Earth" );
                bodyNames.push_back( "Sun" );

                // Specify initial time
                double initialEphemerisTime = double( 1.0E7 );
                double finalEphemerisTime = initialEphemerisTime + 1.0 * 86400.0;

                // Create bodies needed in simulation
                BodyListSettings bodySettings =
                        getDefaultBodySettings( bodyNames );

                SystemOfBodies bodies = createSystemOfBodies( bodySettings );
                bodies.createEmptyBody( "Vehicle" );


                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                // Set accelerations on Vehicle that are to be taken into account.
                SelectedAccelerationMap accelerationMap;
                std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
                accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                        basic_astrodynamics::point_mass_gravity ) );
                accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

                // Set bodies for which initial state is to be estimated and integrated.
                std::vector< std::string > bodiesToIntegrate;
                std::vector< std::string > centralBodies;
                bodiesToIntegrate.push_back( "Vehicle" );
                centralBodies.push_back( "Earth" );

                // Create acceleration models
                AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                        bodies, accelerationMap, bodiesToIntegrate, centralBodies );

                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                // Set Keplerian elements for Asterix.
                Eigen::Vector6d asterixInitialStateInKeplerianElements;
                asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 15000.0E3;
                asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.9;
                asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
                asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
                asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
                asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

                double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

                // Set (perturbed) initial state.
                Eigen::Matrix< double, 6, 1 > systemInitialState = convertKeplerianToCartesianElements(
                        asterixInitialStateInKeplerianElements, earthGravitationalParameter );

                // Create propagator settings
                std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
                dependentVariables.push_back(
                        std::make_shared< SingleDependentVariableSaveSettings >(
                                altitude_dependent_variable, "Vehicle", "Earth" ) );


                // Create integrator settings
                std::shared_ptr< IntegratorSettings< double > > integratorSettings =
                        std::make_shared< RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                                double( initialEphemerisTime ), 70.0,
                                CoefficientSets::rungeKuttaFehlberg78,
                                0.01, 3600.0, 1.0E-12, 1.0E-12 );


                std::vector< std::map< double, Eigen::VectorXd > > numericalResultsVector;
                for( unsigned int i = 0; i < 4; i++ )
                {

                    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                                    centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, initialEphemerisTime, integratorSettings,
                                    std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime ),
                                    cowell, dependentVariables);
                    if( i == 1 )
                    {
                        propagatorSettings->getOutputSettings( )->setResultsSaveFrequencyInSteps( 2 );
                    }

                    if( i == 2 )
                    {
                        propagatorSettings->getOutputSettings( )->setResultsSaveFrequencyInSteps( 0 );
                        propagatorSettings->getOutputSettings( )->setResultsSaveFrequencyInSeconds( 60 );
                    }

                    if( i == 3 )
                    {
                        propagatorSettings->getOutputSettings( )->setResultsSaveFrequencyInSteps( 2 );
                        propagatorSettings->getOutputSettings( )->setResultsSaveFrequencyInSeconds( 60 );
                    }
                    SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, propagatorSettings );
                    std::map< double, Eigen::VectorXd > numericalResults = dynamicsSimulator.getSingleArcPropagationResults( )->getEquationsOfMotionNumericalSolution( );
                    numericalResultsVector.push_back( numericalResults );
                    std::map< double, double > timeStepHistory = getTimeStepHistory( numericalResults );
                }

                {
                    auto nominalIterator = numericalResultsVector.at( 0 ).begin( );
                    auto halfFrequencyIterator = numericalResultsVector.at( 1 ).begin( );

                    while( halfFrequencyIterator != numericalResultsVector.at( 1 ).end( ) )
                    {
                        BOOST_CHECK_EQUAL( nominalIterator->first, halfFrequencyIterator->first );
                        nominalIterator++;
                        nominalIterator++;
                        halfFrequencyIterator++;
                    }

                }

                {
                    auto nominalIterator = numericalResultsVector.at( 0 ).begin( );
                    double previousOutputTime = nominalIterator->first;
                    nominalIterator++;
                    auto perMinuteIterator = numericalResultsVector.at( 2 ).begin( );
                    perMinuteIterator++;

                    while( nominalIterator != numericalResultsVector.at( 0 ).end( ) && perMinuteIterator != numericalResultsVector.at( 2 ).end( ) )
                    {

                        if( nominalIterator->first - previousOutputTime < 60.0 )
                        {
                            nominalIterator++;
                        }
                        else
                        {
                            BOOST_CHECK_EQUAL( nominalIterator->first, perMinuteIterator->first );

                            previousOutputTime = nominalIterator->first;
                            nominalIterator++;
                            perMinuteIterator++;
                        }
                    }

                }

                {
                    auto nominalIterator = numericalResultsVector.at( 0 ).begin( );
                    double previousOutputTime = nominalIterator->first;
                    nominalIterator++;
                    int stepsSinceLastSave = 1;
                    auto perMinuteIterator = numericalResultsVector.at( 3 ).begin( );
                    perMinuteIterator++;

                    while( nominalIterator != numericalResultsVector.at( 0 ).end( ) && perMinuteIterator != numericalResultsVector.at( 3 ).end( ) )
                    {

                        if( nominalIterator->first - previousOutputTime < 60.0 && stepsSinceLastSave < 2 )
                        {
                            nominalIterator++;
                        }
                        else
                        {
                            BOOST_CHECK_EQUAL( nominalIterator->first, perMinuteIterator->first );

                            previousOutputTime = nominalIterator->first;
                            stepsSinceLastSave = 0;
                            nominalIterator++;
                            perMinuteIterator++;
                        }
                        stepsSinceLastSave++;
                    }

                }
            }



        BOOST_AUTO_TEST_SUITE_END( )

    }

}
