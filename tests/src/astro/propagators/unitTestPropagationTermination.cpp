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

#include "tudat/simulation/estimation.h"
#include <boost/test/unit_test.hpp>


using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::gravitation;
using namespace tudat::numerical_integrators;
using namespace tudat::estimatable_parameters;

//! Test suite for astro functions.
BOOST_AUTO_TEST_SUITE( test_propagation_termination )


bool customTerminationConditionOnStateElement( const double time, const Eigen::MatrixXd& state )
{
    if( state( 0 ) < 0 )
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool customTerminationConditionOnAltitude( const double time, const SystemOfBodies& bodies )
{
    if( bodies.at( "Asterix" )->getFlightConditions( )->getCurrentAltitude( ) < 1.0E6 )
    {
        return true;
    }
    else
    {
        return false;
    }
}


BOOST_AUTO_TEST_CASE( testCustomPropagationTermination )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Define body settings for simulation.
    std::vector<std::string> bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );

    // Create body objects.
    BodyListSettings bodySettings =
        getDefaultBodySettings( bodiesToCreate, "SSB", "J2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodies.createEmptyBody( "Asterix" );
    bodies.at( "Asterix" )->setConstantBodyMass( 400.0 );

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 4.0;
    double aerodynamicCoefficient = 1.2;
    std::shared_ptr<AerodynamicCoefficientSettings> aerodynamicCoefficientSettings =
        std::make_shared<ConstantAerodynamicCoefficientSettings>(
            referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ),
            negative_aerodynamic_frame_coefficients );

    // Create and set aerodynamic coefficients object
    bodies.at( "Asterix" )->setAerodynamicCoefficientInterface(
        createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Asterix", bodies ));

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector<std::string> bodiesToPropagate;
    std::vector<std::string> centralBodies;

    // Define propagation settings.
    std::map<std::string, std::vector<std::shared_ptr<AccelerationSettings> > > accelerationsOfAsterix;
    accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared<AccelerationSettings>(
        basic_astrodynamics::point_mass_gravity ));
    accelerationsOfAsterix[ "Sun" ].push_back( std::make_shared<AccelerationSettings>(
        basic_astrodynamics::point_mass_gravity ));
    accelerationsOfAsterix[ "Moon" ].push_back( std::make_shared<AccelerationSettings>(
        basic_astrodynamics::point_mass_gravity ));
    accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared<AccelerationSettings>(
        basic_astrodynamics::aerodynamic ));

    accelerationMap[ "Asterix" ] = accelerationsOfAsterix;
    bodiesToPropagate.push_back( "Asterix" );
    centralBodies.push_back( "Earth" );

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
        bodies, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
        = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
        = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d asterixInitialState = convertKeplerianToCartesianElements(
        asterixInitialStateInKeplerianElements, earthGravitationalParameter );
    std::cout << asterixInitialState << std::endl;

    const double fixedStepSize = 30.0;
    std::shared_ptr<IntegratorSettings<> > integratorSettings =
        rungeKuttaFixedStepSettings<double>( fixedStepSize, numerical_integrators::rungeKutta4Classic );

    std::vector<std::shared_ptr<SingleDependentVariableSaveSettings> > dependentVariablesList;
    dependentVariablesList.push_back(
        std::make_shared<SingleDependentVariableSaveSettings>(
            altitude_dependent_variable, "Asterix", "Earth" ));

    for ( int test = 0; test < 2; test++ )
    {

        std::shared_ptr< PropagationTerminationSettings > terminationSettings;
        if( test == 0 )
        {
            std::function< bool( const double, const Eigen::MatrixXd& ) > terminationFunction =
                std::bind( &customTerminationConditionOnStateElement, std::placeholders::_1, std::placeholders::_2 );
            terminationSettings = popagationCustomTerminationSettingsFromFullState( terminationFunction );
        }
        else if( test == 1 )
        {
            terminationSettings = popagationCustomTerminationSettings( std::bind( &customTerminationConditionOnAltitude, std::placeholders::_1, bodies ) );
        }

        std::shared_ptr<TranslationalStatePropagatorSettings<double> > propagatorSettings =
            std::make_shared<TranslationalStatePropagatorSettings<double> >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, asterixInitialState, simulationStartEpoch,
                  integratorSettings,
                  terminationSettings, cowell,
                  dependentVariablesList );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


        // Create simulation object (but do not propagate dynamics).
        SingleArcDynamicsSimulator<> dynamicsSimulator( bodies, propagatorSettings );

        BOOST_CHECK( ( dynamicsSimulator.getEquationsOfMotionNumericalSolution().size( ) > 50 ) );
        if( test == 0 )
        {
            int checkTermination = 0;
            bool terminationReachedAtEnd = false;
            std::map< double, Eigen::VectorXd > stateHistory = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
            for ( std::map< double, Eigen::VectorXd >::iterator it = stateHistory.begin( ); it != stateHistory.end( ); it++ )
            {
                if( it->second( 0 ) < 0.0 )
                {
                    if( it->first == stateHistory.rbegin( )->first )
                    {
                        terminationReachedAtEnd = true;
                    }
                    checkTermination++;
                }
            }
            BOOST_CHECK( checkTermination == 1 );
            BOOST_CHECK( terminationReachedAtEnd );
        }

        if( test == 1 )
        {
            int checkTermination = 0;
            bool terminationReachedAtEnd = false;
            std::map< double, Eigen::VectorXd > dependentVariableHistory = dynamicsSimulator.getDependentVariableHistory( );
            for ( std::map< double, Eigen::VectorXd >::iterator it = dependentVariableHistory.begin( ); it != dependentVariableHistory.end( ); it++ )
            {
                if( it->second( 0 ) < 1.0E6 )
                {
                    if( it->first == dependentVariableHistory.rbegin( )->first )
                    {
                        terminationReachedAtEnd = true;
                    }
                    checkTermination++;
                }
            }
            BOOST_CHECK( checkTermination == 1 );
            BOOST_CHECK( terminationReachedAtEnd );
        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

