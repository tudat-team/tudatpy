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

#include "tudat/simulation/simulation.h"
#include <boost/test/unit_test.hpp>

//! Test for the setting assessTerminationOnMinorSteps.
BOOST_AUTO_TEST_SUITE( test_assess_propagation_termination_condition_during_integration_substeps )

//! Unit test description:
//! - Use a "simple" termination condition, only limited by end epoch.
//! - Use constant step-size integrator (RK4).
//! - Check that the propagation stops after (before) reaching the termination condition when
//!   `assessTerminationOnMinorSteps` is off (on).
//!
//! For example: initial epoch = 0, step-size = 50 s, end epoch = 1020 s.
//! Expected outcome:
//! - Final epoch 1050 s when `assessTerminationOnMinorSteps` is off.
//! - Final epoch 1000 s when `assessTerminationOnMinorSteps` is on
//!   (because epoch will be 1025 > 1020 when computing k2).
BOOST_AUTO_TEST_CASE( testassessTerminationOnMinorStepsRKFixedStepSize )
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              ///////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::numerical_integrators;


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       ///////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for ( unsigned int assessDuringSubsteps = 0; assessDuringSubsteps <= 1; assessDuringSubsteps++ )
    {
        // Load Spice kernels.
        spice_interface::loadStandardSpiceKernels( );

        // Set simulation time settings.
        const double simulationStartEpoch =  0.0;
        const double simulationEndEpoch = 1020.0;

        // Define body settings for simulation.
        std::vector< std::string > bodiesToCreate;
        bodiesToCreate.push_back( "Sun" );
        bodiesToCreate.push_back( "Earth" );
        bodiesToCreate.push_back( "Moon" );

        // Create body objects.
        BodyListSettings bodySettings =
                getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0,
                                        "SSB", "J2000" );
        for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            bodySettings.at( bodiesToCreate.at( i ) )->ephemerisSettings->resetFrameOrientation( "J2000" );
            bodySettings.at( bodiesToCreate.at( i ) )->rotationModelSettings->resetOriginalFrame( "J2000" );
        }
        SystemOfBodies bodies = createSystemOfBodies( bodySettings );

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE VEHICLE            //////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create spacecraft object.
        bodies.createEmptyBody( "Asterix" );
        bodies.at( "Asterix" )->setConstantBodyMass( 400.0 );

        // Create aerodynamic coefficient interface settings.
        double referenceArea = 4.0;
        double aerodynamicCoefficient = 1.2;
        std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
                std::make_shared< ConstantAerodynamicCoefficientSettings >(
                    referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

        // Create and set aerodynamic coefficients object
        bodies.at( "Asterix" )->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Asterix" ) );

        // Create radiation pressure settings
        double referenceAreaRadiation = 4.0;
        double radiationPressureCoefficient = 1.2;
        std::vector< std::string > occultingBodies;
        occultingBodies.push_back( "Earth" );
        std::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
                std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

        // Create and set radiation pressure settings
        bodies.at( "Asterix" )->setRadiationPressureInterface(
                    "Sun", createRadiationPressureInterface(
                        asterixRadiationPressureSettings, "Asterix", bodies ) );

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ACCELERATIONS          ///////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define propagation settings.
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
        accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

        accelerationsOfAsterix[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::point_mass_gravity ) );
        accelerationsOfAsterix[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                        basic_astrodynamics::point_mass_gravity ) );
        accelerationsOfAsterix[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::cannon_ball_radiation_pressure ) );
        accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::aerodynamic ) );

        accelerationMap[ "Asterix" ] = accelerationsOfAsterix;
        bodiesToPropagate.push_back( "Asterix" );
        centralBodies.push_back( "Earth" );

        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToPropagate, centralBodies );

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            /////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////

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


        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, asterixInitialState, simulationEndEpoch );

        const double fixedStepSize = 50.0;
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >( rungeKutta4, 0.0, fixedStepSize, 1, assessDuringSubsteps );

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            /////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create simulation object (but do not propagate dynamics).
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodies, integratorSettings, propagatorSettings, true, false, false );

        double finalPropagatedEpoch = ( --dynamicsSimulator.getEquationsOfMotionNumericalSolution().end() )->first;
        BOOST_CHECK( finalPropagatedEpoch == ( assessDuringSubsteps ? 1000.0 : 1050.0 ) );
    }
}


//! Unit test description:
//! - Use a termination condition based on altitude (and time limit sufficiently large so that it's not reached).
//! - Use variable step-size integrator (RK78).
//! - Check that the propagation stops after (before) reaching the termination condition when
//!   `assessTerminationOnMinorSteps` is off (on).
//!
//! For example: altitude limit 100 km.
//! Expected outcome:
//! - Final altitude below 100 km when `assessTerminationOnMinorSteps` is off.
//! - Final altitude above 100 km when `assessTerminationOnMinorSteps` is on.
BOOST_AUTO_TEST_CASE( testassessTerminationOnMinorStepsRKVariableStepSize )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::numerical_integrators;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for ( unsigned int assessDuringSubsteps = 0; assessDuringSubsteps <= 1; assessDuringSubsteps++ )
    {
        // Load Spice kernels.
        spice_interface::loadStandardSpiceKernels( );

        // Set simulation time settings.
        const double simulationStartEpoch = 0.0;
        const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

        // Define body settings for simulation.
        std::vector< std::string > bodiesToCreate;
        bodiesToCreate.push_back( "Sun" );
        bodiesToCreate.push_back( "Earth" );
        bodiesToCreate.push_back( "Moon" );

        // Create body objects.
        BodyListSettings bodySettings =
                getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0,
                                         "SSB", "J2000" );
        for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
        {
            bodySettings.at( bodiesToCreate.at( i ) )->ephemerisSettings->resetFrameOrientation( "J2000" );
            bodySettings.at( bodiesToCreate.at( i ) )->rotationModelSettings->resetOriginalFrame( "J2000" );
        }
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
        std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
                std::make_shared< ConstantAerodynamicCoefficientSettings >(
                    referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

        // Create and set aerodynamic coefficients object
        bodies.at( "Asterix" )->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Asterix" ) );

        // Create radiation pressure settings
        double referenceAreaRadiation = 4.0;
        double radiationPressureCoefficient = 1.2;
        std::vector< std::string > occultingBodies;
        occultingBodies.push_back( "Earth" );
        std::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
                std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                    "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

        // Create and set radiation pressure settings
        bodies.at( "Asterix" )->setRadiationPressureInterface(
                    "Sun", createRadiationPressureInterface(
                        asterixRadiationPressureSettings, "Asterix", bodies ) );



        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Define propagator settings variables.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;

        // Define propagation settings.
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
        accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

        accelerationsOfAsterix[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::point_mass_gravity ) );
        accelerationsOfAsterix[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                        basic_astrodynamics::point_mass_gravity ) );
        accelerationsOfAsterix[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                       basic_astrodynamics::cannon_ball_radiation_pressure ) );
        accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                         basic_astrodynamics::aerodynamic ) );

        accelerationMap[ "Asterix" ] = accelerationsOfAsterix;
        bodiesToPropagate.push_back( "Asterix" );
        centralBodies.push_back( "Earth" );

        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToPropagate, centralBodies );

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Set Keplerian elements for Asterix.
        Eigen::Vector6d asterixInitialStateInKeplerianElements;      // initial altitude ~129 km
        asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6500.0E3;
        asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.0;
        asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
        asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                = unit_conversions::convertDegreesToRadians( 235.7 );
        asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                = unit_conversions::convertDegreesToRadians( 23.4 );
        asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

        double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
        const Eigen::Vector6d asterixInitialState = convertKeplerianToCartesianElements(
                    asterixInitialStateInKeplerianElements, earthGravitationalParameter );


        // Termination condition

        std::vector< std::shared_ptr< propagators::PropagationTerminationSettings > > constituentSettings;

        // Time limit
        constituentSettings.push_back(
                    std::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ) );

        // Altitude limit (100 km)
        std::shared_ptr< PropagationTerminationSettings > altitudeTerminationSettings =
                std::make_shared< propagators::PropagationDependentVariableTerminationSettings >(
                    std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                        propagators::altitude_dependent_variable, "Asterix", "Earth" ), 100.0E3, 1 );

        constituentSettings.push_back( altitudeTerminationSettings );

        // Stop if ANY of the two is met
        std::shared_ptr< PropagationTerminationSettings > terminationSettings = std::make_shared<
                propagators::PropagationHybridTerminationSettings >( constituentSettings, 1 );

        // Save altitude dependent variable
        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
           dependentVariables.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                              altitude_dependent_variable, "Asterix", "Earth" ) );
        std::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings =
                std::make_shared< DependentVariableSaveSettings >( dependentVariables, 0 );

        // Create propagator settings
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, asterixInitialState,
                  terminationSettings, cowell, dependentVariableSaveSettings );


        const double initialStepSize = 30.0;
        const double minStepSize = 30.0;
        const double maxStepSize = 30.0;
        const double tolerance = 1.0E-11;
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< RungeKuttaVariableStepSizeSettings< > >
                ( simulationStartEpoch, initialStepSize,
                  rungeKuttaFehlberg78, minStepSize, maxStepSize, tolerance, tolerance, 1,
                  assessDuringSubsteps );



        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Create simulation object (but do not propagate dynamics).
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodies, integratorSettings, propagatorSettings, true, false, false );

        double finalAltitude = ( --dynamicsSimulator.getDependentVariableHistory().end() )->second( 0 );
        BOOST_CHECK( assessDuringSubsteps ? finalAltitude > 100.0E+3 : finalAltitude < 100.0E+3 );
    }
}


BOOST_AUTO_TEST_SUITE_END( )

