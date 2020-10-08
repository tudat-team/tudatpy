/*    Copyright (c) 2010-2019, Delft University of Technology
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

//! Execute propagation of orbit of Asterix around the Earth.
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
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;
    const double environmentTimeBuffer = 300.0;

    // Create bodies in simulation
    std::vector< std::string > bodiesToCreate = { "Sun", "Earth", "Moon", "Mars", "Venus" };
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - environmentTimeBuffer,
                                    simulationEndEpoch + environmentTimeBuffer,  "Earth", "J2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodies.createEmptyBody( "Asterix" );
    bodies.at( "Asterix" )->setConstantBodyMass( 400.0 );

    // Create and add aerodynamic coefficient interface
    double referenceArea = 4.0;
    double aerodynamicCoefficient = 1.2;
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            constantAerodynamicCoefficientSettings(
                referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );
    addAerodynamicCoefficientInterface(
                bodies, "Asterix", aerodynamicCoefficientSettings );

    // Create and add radiation pressure interace
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies = { "Earth" };
    std::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
            cannonBallRadiationPressureSettings(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );
    addRadiationPressureInterface(
                bodies, "Asterix", asterixRadiationPressureSettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;

    using namespace tudat::basic_astrodynamics;
    accelerationsOfAsterix[ "Earth" ] = {
            sphericalHarmonicAcceleration( 5, 5 ),
            aerodynamicAcceleration() };

    accelerationsOfAsterix[ "Sun" ] = {
            pointMassGravityAcceleration( ),
            cannonBallRadiationPressureAcceleration( ) };

    accelerationsOfAsterix[ "Mars" ] = {
            pointMassGravityAcceleration( ) };

    accelerationsOfAsterix[ "Venus" ] = {
            pointMassGravityAcceleration( ) };

    accelerationsOfAsterix[ "Moon" ] = {
            pointMassGravityAcceleration( ) };

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationSettingsList;
    accelerationSettingsList[ "Asterix" ] = accelerationsOfAsterix;

    std::vector< std::string > bodiesToPropagate = { "Asterix" };
    std::vector< std::string > centralBodies = { "Earth" };

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationSettingsList, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = getBodyGravitationalParameter( bodies, "Earth" );
    const Eigen::Vector6d asterixInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, asterixInitialState, simulationEndEpoch );

    const double fixedStepSize = 10.0;
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, 0.0, fixedStepSize );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator( bodies, integratorSettings, propagatorSettings );
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::string outputSubFolder = "PerturbedSatelliteExample/";

    Eigen::VectorXd finalIntegratedState = (--integrationResult.end( ) )->second;
    // Print the position (in km) and the velocity (in km/s) at t = 0.
    std::cout << "Single Earth-Orbiting Satellite Example." << std::endl <<
                 "The initial position vector of Asterix is [km]:" << std::endl <<
                 asterixInitialState.segment( 0, 3 ) / 1E3 << std::endl <<
                 "The initial velocity vector of Asterix is [km/s]:" << std::endl <<
                 asterixInitialState.segment( 3, 3 ) / 1E3 << std::endl;

    // Print the position (in km) and the velocity (in km/s) at t = 86400.
    std::cout << "After " << simulationEndEpoch <<
                 " seconds, the position vector of Asterix is [km]:" << std::endl <<
                 finalIntegratedState.segment( 0, 3 ) / 1E3 << std::endl <<
                 "And the velocity vector of Asterix is [km/s]:" << std::endl <<
                 finalIntegratedState.segment( 3, 3 ) / 1E3 << std::endl;

    // Write perturbed satellite propagation history to file.
    input_output::writeDataMapToTextFile( integrationResult,
                                          "singlePerturbedSatellitePropagationHistory.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );


    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;
}
