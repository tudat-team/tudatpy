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
#include <tudat/io/applicationOutput.h>

//! Execute propagation of orbit of vehicle around the Earth. The vehicle is subject to a thrustforce, which is specified in
//! the nonconstantThrust.txt file. In that file, the first column is time in seconds, the last three columns give the x, y
//! and z components of the thrust force in the J2000 (?) frame.
int main( )
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::numerical_integrators;
    using namespace tudat::interpolators;
    using namespace tudat::unit_conversions;
    using namespace tudat::propagators;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define body settings for simulation.
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Earth" ] = std::make_shared< BodySettings >( );
    bodySettings[ "Earth" ]->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodySettings[ "Earth" ]->gravityFieldSettings = std::make_shared< GravityFieldSettings >( central_spice );

    // Create Earth object
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Create spacecraft object.
    double vehicleMass = 2000.0;
    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 4.0;
    double aerodynamicCoefficient = 1.2;
    std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            std::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

    // Create and set aerodynamic coefficients object
    bodyMap[ "Vehicle" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle" ) );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define data to be used for thrust as a function of time.
    const std::string cppFilePath( __FILE__ );
    const std::string cppFolder = cppFilePath.substr( 0 , cppFilePath.find_last_of("/\\")+1 );
    std::shared_ptr< FromFileDataMapSettings< Eigen::Vector3d > > thrustDataSettings =
            std::make_shared< FromFileDataMapSettings< Eigen::Vector3d > >( cppFolder + "testThrustValues.txt" );

    // Define interpolator settings.
    std::shared_ptr< InterpolatorSettings > thrustInterpolatorSettings =
            std::make_shared< InterpolatorSettings >( linear_interpolator );

    // Create data interpolation settings
    std::shared_ptr< DataInterpolationSettings< double, Eigen::Vector3d > > thrustDataInterpolatorSettings =
            std::make_shared< DataInterpolationSettings< double, Eigen::Vector3d > >(
                thrustDataSettings, thrustInterpolatorSettings );

    // Define specific impulse
    double constantSpecificImpulse = 3000.0;

    // Define propagation settings.
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
                                                     basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared< ThrustAccelerationSettings >(
                                                       thrustDataInterpolatorSettings, constantSpecificImpulse,
                                                       lvlh_thrust_frame, "Earth" ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Set initial conditions for the vehicle satellite that will be propagated in this simulation.
    // The initial conditions are given in Keplerian elements and later on converted to Cartesian
    // elements.
    // Set Keplerian elements for vehicle.
    Eigen::Vector6d vehicleInitialStateInKeplerianElements;
    vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 72130.0e3;
    vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.6;
    vehicleInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 169.0 );
    vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 45.0 );
    vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = convertDegreesToRadians( 80.0 );
    vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 15.0 );

    // Convert vehicle state from Keplerian elements to Cartesian elements.
    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::VectorXd systemInitialState = convertKeplerianToCartesianElements( vehicleInitialStateInKeplerianElements,
                                                                              earthGravitationalParameter );

    // Define propagation termination conditions (stop after 2 weeks).
    std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            std::make_shared< propagators::PropagationTimeTerminationSettings >( 4.0E6 );

    // Define settings for propagation of translational dynamics.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings );

    // Crete mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Vehicle" ] = createMassRateModel( "Vehicle", std::make_shared< FromThrustMassModelSettings >( 1 ),
                                                       bodyMap, accelerationModelMap );

    // Create settings for propagating the mass of the vehicle
    std::shared_ptr< MassPropagatorSettings< double > > massPropagatorSettings =
            std::make_shared< MassPropagatorSettings< double > >( std::vector< std::string >{ "Vehicle" }, massRateModels,
                                                                  ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ),
                                                                  terminationSettings );

    // Create list of propagation settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsVector;
    propagatorSettingsVector.push_back( translationalPropagatorSettings );
    propagatorSettingsVector.push_back( massPropagatorSettings );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                                          basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );
    dependentVariablesList.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
                                          lvlh_to_inertial_frame_rotation_dependent_variable, "Vehicle", "Earth" ) );

    // Create object with list of dependent variables
    std::shared_ptr< DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< DependentVariableSaveSettings >( dependentVariablesList );

    // Create propagation settings for mass and translational dynamics concurrently
    std::shared_ptr< PropagatorSettings< > > propagatorSettings =
            std::make_shared< MultiTypePropagatorSettings< double > >(
                propagatorSettingsVector, terminationSettings, dependentVariablesToSave );

    // Define integrator settings
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >( rungeKutta4, 0.0, 30.0 );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT AND PRINT OUTPUT TO FILE         //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings );
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableResult = dynamicsSimulator.getDependentVariableHistory( );

    // Manually add thrust force in LVLH frame to output
    for( std::map< double, Eigen::VectorXd >::iterator outputIterator = dependentVariableResult.begin( );
         outputIterator != dependentVariableResult.end( ); outputIterator++ )
    {
        Eigen::Matrix3d currentRotationMatrix =
                getMatrixFromVectorRotationRepresentation( outputIterator->second.segment( 3, 9 ) );
        Eigen::Vector3d currentThrust = outputIterator->second.segment( 0, 3 );
        Eigen::VectorXd newOutput = Eigen::VectorXd( 15 );
        newOutput.segment( 0, 12 ) = outputIterator->second;
        newOutput.segment( 12, 3 ) =
                integrationResult.at( outputIterator->first )( 6 ) *
                ( currentRotationMatrix.transpose( ) * currentThrust );
        dependentVariableResult[ outputIterator->first ] = newOutput;
    }

    std::string outputSubFolder = "ThrustFromFileExample/";

    // Write propagation history to file.
    input_output::writeDataMapToTextFile( integrationResult,
                                          "thrustExampleFromFilePropagationHistory.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );

    // Write dependent variable history to file.
    input_output::writeDataMapToTextFile( dependentVariableResult,
                                          "thrustExampleFromFileDependentVariableHistory.dat",
                                          tudat_applications::getOutputPath( ) + outputSubFolder,
                                          "",
                                          std::numeric_limits< double >::digits10,
                                          std::numeric_limits< double >::digits10,
                                          "," );
}
