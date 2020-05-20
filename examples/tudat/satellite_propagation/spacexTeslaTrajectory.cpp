/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <tudat/simulation/simulation.h>

#include "tudat/io/applicationOutput.h"

//! Execute propagation of orbit of TeslaRoadster around the Earth.
int main( )
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


    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( { input_output::getSpiceKernelPath( ) + "de430.bsp" } );

    // Set simulation time settings.
    const double simulationStartEpoch =
            ( 2458163.500000000 - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * 86400.0;
    const double simulationEndEpoch = 500.0 * physical_constants::JULIAN_YEAR;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Venus" );
    bodiesToCreate.push_back( "Jupiter" );
    bodiesToCreate.push_back( "Saturn" );
    bodiesToCreate.push_back( "Mercury" );

    // Create body objects.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate );
//    bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( "J2000" );
//    bodySettings[ "Moon" ]->ephemerisSettings->resetFrameOrientation( "J2000" );
//    bodySettings[ "Earth" ]->ephemerisSettings->resetFrameOrientation( "J2000" );

//    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
//    {
//        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
//    }

    bodySettings[ "Mars" ]->ephemerisSettings = std::make_shared< ApproximatePlanetPositionSettings >(
                ephemerides::ApproximatePlanetPositionsBase::mars, false );
    bodySettings[ "Jupiter" ]->ephemerisSettings = std::make_shared< ApproximatePlanetPositionSettings >(
                ephemerides::ApproximatePlanetPositionsBase::mars, false );
    bodySettings[ "Saturn" ]->ephemerisSettings = std::make_shared< ApproximatePlanetPositionSettings >(
                ephemerides::ApproximatePlanetPositionsBase::mars, false );
    bodySettings[ "Venus" ]->ephemerisSettings = std::make_shared< ApproximatePlanetPositionSettings >(
                ephemerides::ApproximatePlanetPositionsBase::mars, false );
    bodySettings[ "Mercury" ]->ephemerisSettings = std::make_shared< ApproximatePlanetPositionSettings >(
                ephemerides::ApproximatePlanetPositionsBase::mars, false );

    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "TeslaRoadster" ] = std::make_shared< simulation_setup::Body >( );

    bodyMap[ "TeslaRoadster" ]->setConstantBodyMass( 1500.0 );

    // Create radiation pressure settings
    double referenceAreaRadiation = 8.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    std::shared_ptr< RadiationPressureInterfaceSettings > teslaRoadsterRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "TeslaRoadster" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    teslaRoadsterRadiationPressureSettings, "TeslaRoadster", bodyMap ) );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;
    accelerationsOfAsterix[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfAsterix[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                   basic_astrodynamics::central_gravity ) );
    accelerationsOfAsterix[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfAsterix[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfAsterix[ "Venus" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfAsterix[ "Mercury" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfAsterix[ "Saturn" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfAsterix[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfAsterix[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::cannon_ball_radiation_pressure ) );


    accelerationMap[  "TeslaRoadster" ] = accelerationsOfAsterix;
    bodiesToPropagate.push_back( "TeslaRoadster" );
    centralBodies.push_back( "Sun" );

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for TeslaRoadster.
    Eigen::Vector6d teslaRoadsterInitialStateJ2000;
    teslaRoadsterInitialStateJ2000 << -1.484544094935328E+06, -1.530028575781120E+06, -4.100899567817330E+05,
    -2.358201939417945E+00, -2.459124989364402E+00, -6.370785284397021E-01;
    teslaRoadsterInitialStateJ2000 *= 1000.0;
    teslaRoadsterInitialStateJ2000 += spice_interface::getBodyCartesianStateAtEpoch(
                "Earth", "SSB", "J2000", "None", simulationStartEpoch );

    Eigen::Quaterniond toEclipJ2000 = spice_interface::computeRotationQuaternionBetweenFrames(
                                 "J2000", "ECLIPJ2000", 0.0 );

    Eigen::Vector6d teslaRoadsterInitialState;
    teslaRoadsterInitialState.segment( 0, 3 ) = toEclipJ2000 * teslaRoadsterInitialStateJ2000.segment( 0, 3 );
    teslaRoadsterInitialState.segment( 3, 3 ) = toEclipJ2000 * teslaRoadsterInitialStateJ2000.segment( 3, 3 );
\




    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    relative_position_dependent_variable, "Earth", "Sun" ) );
    dependentVariablesList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    relative_position_dependent_variable, "Mars", "Sun" ) );
    dependentVariablesList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    relative_position_dependent_variable, "Venus", "Sun" ) );
    dependentVariablesList.push_back(
                std::make_shared< SingleDependentVariableSaveSettings >(
                    relative_distance_dependent_variable, "TeslaRoadster", "Earth" ) );

    // Create object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, teslaRoadsterInitialState, simulationEndEpoch, cowell,
              dependentVariablesToSave );

    std::shared_ptr< IntegratorSettings< > > integratorSettings
            = std::make_shared< tudat::numerical_integrators::BulirschStoerIntegratorSettings< > >(
                simulationStartEpoch, 3600.0,
                numerical_integrators::bulirsch_stoer_sequence, 4,
                std::numeric_limits< double >::epsilon( ), std::numeric_limits< double >::infinity( ),
                1.0E-15, 1.0E-12 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableResult = dynamicsSimulator.getDependentVariableHistory( );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::string outputSubFolder = "SpacexTeslaExample/";

    // Write perturbed satellite propagation history to file.
    input_output::writeDataMapToTextFile( integrationResult,
                                          "spacexTeslaPropagationHistory.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    input_output::writeDataMapToTextFile( dependentVariableResult,
                                          "spacexTeslaDependentVariableHistory.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}

