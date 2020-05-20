/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <tudat/simulation/simulation.h>
#include <tudat/Optimization/missionLinker.h>

#include <pagmo/rng.hpp> //<-Only needed for setting the seed

#include "tudat/io/applicationOutput.h"

//! Execute optimization of burn time of a shuttle model chaser
//! in order to reach an ISS model target
int main( )
{

    //Set seed for reproducible results

    pagmo::random_device::set_seed(123);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::basic_mathematics;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::unit_conversions;
    using namespace tudat::input_output;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLES      //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             Mission segment 1               ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set simulation start epoch.
    const double simulationStartEpoch = 0.0;

    // Set preliminary simulation end epoch.
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Set numerical integration fixed step size for the first simulation.
    double fixedStepSize = 0.1;

    // Define simulation body settings.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( { "Earth" }, simulationStartEpoch - 300,
                                    simulationEndEpoch + 300 );

    bodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );

    bodySettings[ "Earth" ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    bodySettings[ "Earth" ]->atmosphereSettings = NULL;
    bodySettings[ "Earth" ]->shapeModelSettings = NULL;

    // Create Earth object for first simulation
    simulation_setup::NamedBodyMap bodyMap_1 = simulation_setup::createBodies( bodySettings );

    // Create vehicle objects for first simulation.
    bodyMap_1[ "Target" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap_1[ "Chaser" ] = std::make_shared< simulation_setup::Body >( );

    // Set initial body mass of chaser (Space shuttle weight without
    // payload and rendez vous + de-orbit fuel )
    double chaserMass = 90116.4; //kg
    bodyMap_1[ "Chaser" ]->setConstantBodyMass( chaserMass );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap_1, "SSB", "J2000" );

    // Define propagator settings variables for first simulation.
    SelectedAccelerationMap accelerationMap_1;
    std::vector< std::string > bodiesToPropagate_1;
    std::vector< std::string > centralBodies_1;

    // Define target acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfTarget_1;
    accelerationsOfTarget_1[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
    accelerationMap_1[  "Target" ] = accelerationsOfTarget_1;
    bodiesToPropagate_1.push_back( "Target" );
    centralBodies_1.push_back( "Earth" );

    // Define thrust settings of chaser

    // Space shuttle OMS thrust (AJ10-190 engine)
    double thrustMagnitude = 26700.0; //N
    double specificImpulse = 316.0; //s

    // Thrust opposite to direction of motion
    std::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings =
            std::make_shared< ThrustDirectionFromStateGuidanceSettings >(
                "Earth", true, true);
    std::shared_ptr< ThrustEngineSettings > thrustMagnitudeSettings =
            std::make_shared< ConstantThrustEngineSettings >(
                thrustMagnitude, specificImpulse );

    // Define chaser acceleration model settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfChaser_1;
    accelerationsOfChaser_1[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
    accelerationsOfChaser_1[ "Chaser" ].push_back(
                std::make_shared< ThrustAccelerationSettings >( thrustDirectionGuidanceSettings, thrustMagnitudeSettings) );
    accelerationMap_1[  "Chaser" ] = accelerationsOfChaser_1;
    bodiesToPropagate_1.push_back( "Chaser" );
    centralBodies_1.push_back( "Earth" );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap_1 = createAccelerationModelsMap(
                bodyMap_1, accelerationMap_1, bodiesToPropagate_1, centralBodies_1 );

    // Set initial Keplerian elements for target (ISS)
    Eigen::Vector6d targetInitialStateInKeplerianElements;
    targetInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6782650.0;
    targetInitialStateInKeplerianElements( eccentricityIndex ) = 0.0005;
    targetInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 51.640 );
    targetInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 235.7 );
    targetInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 23.4 );
    targetInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.8 );

    // Set Keplerian elements for Chaser close to target (calculated separation ~105 km)
    Eigen::Vector6d chaserInitialStateInKeplerianElements( 6 );
    chaserInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6782650.0;
    chaserInitialStateInKeplerianElements( eccentricityIndex ) = 0.00025;
    chaserInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 51.635 );
    chaserInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 84.1 );
    chaserInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 23.4 );
    chaserInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 290.5 );


    // Convert initial states from Keplerian to Cartesian elements.
    double earthGravitationalParameter = bodyMap_1.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Convert target state from Keplerian elements to Cartesian elements.
    const Eigen::Vector6d targetInitialState = convertKeplerianToCartesianElements(
                targetInitialStateInKeplerianElements,
                earthGravitationalParameter );

    // Convert chaser state from Keplerian elements to Cartesian elements.
    const Eigen::Vector6d chaserInitialState = convertKeplerianToCartesianElements(
                chaserInitialStateInKeplerianElements,
                earthGravitationalParameter );


    // Calculate preliminary information
    double spacecraftSeparation = (targetInitialState.segment(0,3) - chaserInitialState.segment(0,3)).norm();
    std::cout << "Initial spacecraft seaparation = " << spacecraftSeparation/1000 << "km\n\n";
    double meanOrbitalMotion = pow( earthGravitationalParameter/pow(
            targetInitialStateInKeplerianElements( semiMajorAxisIndex ), 3 ), 0.5 );

    // Clohessy Wiltshire solution for estimated impulse (with x_0 = 0 and in y-direction)
    double deltaV = fabs(spacecraftSeparation*meanOrbitalMotion/(6*mathematical_constants::PI));
    std::cout << "Required impulse = " << deltaV << "m/s\n\n";

    // Tsiolkowsky solution for estimated burn time
    double massRatio = exp( deltaV/( specificImpulse*physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION ) );
    double estimatedBurnTime = specificImpulse*physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION*
            chaserMass*(1 - 1/massRatio )/thrustMagnitude;
    std::cout << "Estimated burn time = " << estimatedBurnTime << "s\n\n";

    // Clohessy Wiltshire solution for estimated time to rendez-vous
    double estimatedRendezVousTime = 2.0 *mathematical_constants::PI/meanOrbitalMotion;

    // Set initial state for first simulation
    Eigen::VectorXd systemInitialState = Eigen::VectorXd( 12 );
    systemInitialState.segment( 0, 6 ) = targetInitialState;
    systemInitialState.segment( 6, 6 ) = chaserInitialState;

    // Preliminary time termination settings for first simulation
    std::shared_ptr< PropagationTimeTerminationSettings > timeTerminationSettings_1 =
            std::make_shared< propagators::PropagationTimeTerminationSettings >( estimatedBurnTime );

    // State propagator for first simulation
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies_1, accelerationModelMap_1, bodiesToPropagate_1, systemInitialState, timeTerminationSettings_1 );

    // Create mass rate model for the chaser
    std::shared_ptr< MassRateModelSettings > massRateModelSettings =
            std::make_shared< FromThrustMassModelSettings >( 1 );
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Chaser" ] = createMassRateModel(
            "Chaser", massRateModelSettings, bodyMap_1, accelerationModelMap_1 );

    // Create mass propagator settings for the chaser
    std::vector< std::string > bodiesWithMassToPropagate;
    bodiesWithMassToPropagate.push_back( "Chaser" );

    Eigen::VectorXd initialBodyMasses = Eigen::VectorXd( 1 );
    initialBodyMasses( 0 ) = chaserMass;

    std::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
            std::make_shared< MassPropagatorSettings< double > >(
                    bodiesWithMassToPropagate, massRateModels, initialBodyMasses, timeTerminationSettings_1 );

    // Create multi-type propagator settings (state + chaser mass propagator)
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
                    propagatorSettingsVector.push_back( translationalPropagatorSettings);
                    propagatorSettingsVector.push_back( massPropagatorSettings );

    std::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings_1 =
            std::make_shared< MultiTypePropagatorSettings< double > >(
                                    propagatorSettingsVector, timeTerminationSettings_1 );

    // Create integrator settings for first simulation
    std::shared_ptr< IntegratorSettings< > > integratorSettings_1 =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    // Create dynamics simulator for first simulation
    std::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator_1
            = std::make_shared< SingleArcDynamicsSimulator< > >(
                bodyMap_1, integratorSettings_1, propagatorSettings_1, false, false, false );

    // Set the simulation time as decision variable, with boundaries +-5 seconds from
    // the estimated burn time.
    std::vector< std::shared_ptr< optimization::SingleDecisionVariableSettings > > decisionVariableList;
    decisionVariableList.push_back( std::make_shared< optimization::SingleDecisionVariableSettings >(
                    optimization::simulation_time_decision_variable,
                            estimatedBurnTime - 10.0, estimatedBurnTime + 10.0 ) );
    decisionVariableList.push_back(
                std::make_shared< optimization::SingleCartesianComponentDecisionVariableSettings >(
                    orbital_elements::xCartesianVelocity,
                            chaserInitialState( 3 ) - 10.0, chaserInitialState( 3 ) + 10.0, "Chaser" ) );
    decisionVariableList.push_back(
                std::make_shared< optimization::SingleCartesianComponentDecisionVariableSettings >(
                    orbital_elements::yCartesianVelocity,
                    chaserInitialState( 4 ) - 10.0, chaserInitialState( 4 ) + 10.0, "Chaser" ) );
    decisionVariableList.push_back(
                std::make_shared< optimization::SingleCartesianComponentDecisionVariableSettings >(
                    orbital_elements::zCartesianVelocity,
                    chaserInitialState( 5 ) - 10.0, chaserInitialState( 5 ) + 10.0, "Chaser" ) );

    std::cout<<"Boundaries 0 "<<estimatedBurnTime - 10.0<<" "<<estimatedBurnTime + 10.0<<std::endl;
    std::cout<<"Boundaries 1 "<<chaserInitialState( 3 ) - 10.<<" "<<chaserInitialState( 3 ) + 10.<<std::endl;
    std::cout<<"Boundaries 2 "<<chaserInitialState( 4 ) - 10.<<" "<<chaserInitialState( 4 ) + 10.<<std::endl;
    std::cout<<"Boundaries 3 "<<chaserInitialState( 5 ) - 10.<<" "<<chaserInitialState( 5 ) + 10.<<std::endl;

    std::shared_ptr< optimization::DecisionVariableSettings > decisionVariableSettings = std::make_shared<
            optimization::DecisionVariableSettings >( decisionVariableList );

    // Create first mission segment with the dynamics simulator and the decision variable
    std::shared_ptr< optimization::MissionSegmentSettings > missionSegment_1 =
            std::make_shared< optimization::MissionSegmentSettings >( dynamicsSimulator_1, decisionVariableSettings );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             Mission segment 2               ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create Earth object for second simulation
    simulation_setup::NamedBodyMap bodyMap_2 = simulation_setup::createBodies( bodySettings );

    // Create target and chaser bodies for second simulation
    bodyMap_2[ "Target" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap_2[ "Chaser" ] = std::make_shared< simulation_setup::Body >( );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap_2, "SSB", "J2000" );

    // Define propagator settings variables for second simulation
    SelectedAccelerationMap accelerationMap_2;
    std::vector< std::string > bodiesToPropagate_2;
    std::vector< std::string > centralBodies_2;

    fixedStepSize = 5.0; //s

    // Define acceleration model settings for target and chaser fro second simulation (unpowered flight for both).
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfTarget_2;
    accelerationsOfTarget_2[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
    accelerationMap_2[ "Target" ] = accelerationsOfTarget_2;
    bodiesToPropagate_2.push_back( "Target" );
    centralBodies_2.push_back( "Earth" );
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfChaser_2;
    accelerationsOfChaser_2[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 4, 0 ) );
    accelerationMap_2[  "Chaser" ] = accelerationsOfChaser_2;
    bodiesToPropagate_2.push_back( "Chaser" );
    centralBodies_2.push_back( "Earth" );


    // Create acceleration model for second simulation.
    basic_astrodynamics::AccelerationMap accelerationModelMap_2 = createAccelerationModelsMap(
                bodyMap_2, accelerationMap_2, bodiesToPropagate_2, centralBodies_2 );

    // Create dependent variable save settings to retrieve
    // the distance between the two spacecraft
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariableSaveSettingsVector;
    dependentVariableSaveSettingsVector.push_back( std::make_shared< SingleDependentVariableSaveSettings > (
            relative_distance_dependent_variable, "Chaser", "Target" ) );
    std::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings =
            std::make_shared< DependentVariableSaveSettings >( dependentVariableSaveSettingsVector, false );

    // Create hybrid termination settings to simulate until estimated
    // rendez vous + margin or until minimum separation
    const double minimumSeparation = 1000; //m

    std::vector< std::shared_ptr< propagators::PropagationTerminationSettings > > multiTerminationSettings;

    multiTerminationSettings.push_back( std::make_shared< propagators::PropagationTimeTerminationSettings >(
                                            estimatedBurnTime + estimatedRendezVousTime + 10.0 * fixedStepSize ) );

    multiTerminationSettings.push_back( std::make_shared<
            propagators::PropagationDependentVariableTerminationSettings >(
                    dependentVariableSaveSettingsVector[0], minimumSeparation, true ) );

    std::shared_ptr< propagators::PropagationHybridTerminationSettings > hybridTerminationSettings =
        std::make_shared< propagators::PropagationHybridTerminationSettings >( multiTerminationSettings, true );

    // Create propagator settings for second simulation
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings_2 =
            std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies_2, accelerationModelMap_2, bodiesToPropagate_2, Eigen::VectorXd::Zero(12),
                            hybridTerminationSettings, cowell, dependentVariableSaveSettings );

    // Create integrator settings for second simulation
    std::shared_ptr< IntegratorSettings< > > integratorSettings_2 =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, simulationStartEpoch, fixedStepSize );

    // Create dynamics simulator for second simulation
    std::shared_ptr< SingleArcDynamicsSimulator< > > dynamicsSimulator_2 =
            std::make_shared< SingleArcDynamicsSimulator< > >( bodyMap_2, integratorSettings_2,
                    propagatorSettings_2, false );

    // Use the minimum separation between spacecraft as objective function

    const int maxNumberOfEvolutions = 100;
    const double objectiveValue = 0.0;
    std::shared_ptr< optimization::ObjectiveFunctionFromMinOrMaxDependentVariableSettings > objectiveFunction =
            std::make_shared< optimization::ObjectiveFunctionFromMinOrMaxDependentVariableSettings >(
                    dependentVariableSaveSettingsVector[0], objectiveValue, 0, true,
                            minimumSeparation, maxNumberOfEvolutions );

    // Create mission segment for second simulation containing the second
    // dynamics simulator and the objective function
    std::shared_ptr< optimization::MissionSegmentSettings > missionSegment_2 =
            std::make_shared< optimization::MissionSegmentSettings >( dynamicsSimulator_2,
                    objectiveFunction );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////                        MISSION LINKING AND OPTIMIZATION                ////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create object to link the two mission segments and optimize the mission constraint

    std::vector< std::shared_ptr< optimization::MissionSegmentSettings > > missionSegments;

    missionSegments.push_back( missionSegment_1 );
    missionSegments.push_back( missionSegment_2 );

    optimization::MissionLinker missionLinker( missionSegments,
            std::make_shared< optimization::OptimizationSettings >( optimization::global_optimization,
                    16, true, 1 ), false  );

    // Start optimization
    missionLinker.optimize();

    // Show the optimum value. Since the decision variable is the simulation time of the
    // first simulation, it is stored inside the termination settings.
    std::cout << "\n Calculated burn time = " << timeTerminationSettings_1->terminationTime_ << "s\n";



    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT OF FIRST SIMULATION TO FILE         /////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::map< double, Eigen::VectorXd > integrationResult_1 = dynamicsSimulator_1->getEquationsOfMotionNumericalSolution();

    // Retrieve numerically integrated states of vehicles.
    std::map< double, Eigen::VectorXd > chaserPropagationHistory_1;
    std::map< double, Eigen::VectorXd > targetPropagationHistory_1;
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult_1.begin( );
         stateIterator != integrationResult_1.end( ); stateIterator++ )
    {
        targetPropagationHistory_1[ stateIterator->first ] = stateIterator->second.segment( 0, 6 );
        chaserPropagationHistory_1[ stateIterator->first ] = stateIterator->second.segment( 6, 7 );
    }


    // Write Asterix propagation history to file.
    writeDataMapToTextFile( chaserPropagationHistory_1,
                            "chaserPropagationHistory_1.dat",
                            tudat_applications::getOutputPath( ),
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

    // Write obelix propagation history to file.
    writeDataMapToTextFile( targetPropagationHistory_1,
                            "targetPropagationHistory_1.dat",
                            tudat_applications::getOutputPath( ),
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////        PROVIDE OUTPUT OF SECOND SIMULATION TO FILE                        //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    std::map< double, Eigen::VectorXd > integrationResult_2 = dynamicsSimulator_2->getEquationsOfMotionNumericalSolution();

    // Retrieve numerically integrated states of vehicles.
    std::map< double, Eigen::VectorXd > chaserPropagationHistory_2;
    std::map< double, Eigen::VectorXd > targetPropagationHistory_2;
    for( std::map< double, Eigen::VectorXd >::const_iterator stateIterator = integrationResult_2.begin( );
         stateIterator != integrationResult_2.end( ); stateIterator++ )
    {
        targetPropagationHistory_2[ stateIterator->first ] = stateIterator->second.segment( 0, 6 );
        chaserPropagationHistory_2[ stateIterator->first ] = stateIterator->second.segment( 6, 6 );
    }

    // Write Asterix propagation history to file.
    writeDataMapToTextFile( chaserPropagationHistory_2,
                            "chaserPropagationHistory_2.dat",
                            tudat_applications::getOutputPath( ),
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

    // Write obelix propagation history to file.
    writeDataMapToTextFile( targetPropagationHistory_2,
                            "targetPropagationHistory_2.dat",
                            tudat_applications::getOutputPath( ),
                            "",
                            std::numeric_limits< double >::digits10,
                            std::numeric_limits< double >::digits10,
                            "," );

    // Final statement.
    // The exit code EXIT_SUCCESS indicates that the program was successfully executed.
    return EXIT_SUCCESS;}
