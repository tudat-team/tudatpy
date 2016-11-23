/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>

#include "applicationOutput.h"
#include "Tudat/InputOutput/basicInputOutput.h"


namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_aerodynamic_acceleration_force_moment_models )

BOOST_AUTO_TEST_CASE( testTabulatedDragCoefficient )
{
    double tolerance = 1.0E-12;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace std;
    using namespace tudat;
    using namespace simulation_setup;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace orbital_element_conversions;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;
    using namespace gravitation;
    using namespace numerical_integrators;
    using namespace interpolators;
    using namespace input_output;


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Define body settings for simulation.
    vector< string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );

    // Create body objects.
    map< string, boost::shared_ptr< BodySettings > > bodySettings = getDefaultBodySettings( bodiesToCreate, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0 );

    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ ) {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }

    // EARTH
    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice );
    bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< AtmosphereSettings >( nrlmsise00 );

    // MOON
    bodySettings[ "Moon" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice );

    NamedBodyMap bodyMap = createBodies( bodySettings );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Vehicle" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( 400 );
    double referenceArea = 10;

    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings;

    // Read CD as a function of altitude
    Eigen::MatrixXd CDh = readMatrixFromFile( tudat_applications::getOutputPath( ) + "tabulatedDragCoefficient.txt" );
    Eigen::VectorXd altitudes_km = CDh.col(0);
    vector< double > altitudes;
    for ( unsigned int i = 0; i < altitudes_km.size( ); i++ ) {
        altitudes.push_back( altitudes_km(i) * 1.0E3 );
    }
    Eigen::VectorXd CDs = CDh.col(1);
    vector< double > dragCoefficients( CDs.data(), CDs.data() + CDs.size() );
    vector< Eigen::Vector3d > aerodynamicCoefficients;
    for ( unsigned int i = 0; i < dragCoefficients.size(); i++ ) {
        aerodynamicCoefficients.push_back( dragCoefficients[i] * Eigen::Vector3d::UnitX() );
    }

    // Create interpolator
    boost::shared_ptr< InterpolatorSettings > interpolatorSettings = boost::make_shared< InterpolatorSettings >( OneDimensionalInterpolatorTypes::linear_interpolator );

    // Tabulated aerodynamic settings
    aerodynamicCoefficientSettings = boost::make_shared< TabulatedAerodynamicCoefficientSettings< 1 > >( altitudes, aerodynamicCoefficients, referenceArea, aerodynamics::AerodynamicCoefficientsIndependentVariables::altitude_dependent, interpolatorSettings, 1, 1 );

    // Aerodynamics interface
    bodyMap[ "Vehicle" ]->setAerodynamicCoefficientInterface( createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle" ) );

    // SRP interface
    boost::shared_ptr< RadiationPressureInterfaceSettings > vehicleRadiationPressureSettings = boost::make_shared< CannonBallRadiationPressureInterfaceSettings >( "Sun", referenceArea, 2.2, boost::assign::list_of( "Earth" ) );
    bodyMap[ "Vehicle" ]->setRadiationPressureInterface( "Sun", createRadiationPressureInterface( vehicleRadiationPressureSettings, "Vehicle", bodyMap ) );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    vector< string > bodiesToPropagate;
    vector< string > centralBodies;

    // Define propagation settings.
    map< string, vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;

    accelerationsOfVehicle[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );
    accelerationsOfVehicle[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );

    accelerationsOfVehicle[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( basic_astrodynamics::cannon_ball_radiation_pressure ) );

    accelerationsOfVehicle[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( basic_astrodynamics::aerodynamic ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodyMap, accelerationMap, bodiesToPropagate, centralBodies );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double meanEarthRadius = 6371E3;
    double h_p = 150*1E3;
    double h_a = 35780*1E3;
    double a = (h_a + h_p)/2 + meanEarthRadius;
    double e = (h_a - h_p)/(h_a + h_p + 2*meanEarthRadius);

    // Set Keplerian elements for Vehicle.
    Vector6d vehicleInitialStateInKeplerianElements;
    vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = a;
    vehicleInitialStateInKeplerianElements( eccentricityIndex ) = e;
    vehicleInitialStateInKeplerianElements( inclinationIndex ) = 23.4;

    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Vector6d vehicleInitialState = convertKeplerianToCartesianElements( vehicleInitialStateInKeplerianElements, earthGravitationalParameter );

    // Hybrid termination conditions
    vector< boost::shared_ptr< propagators::PropagationTerminationSettings > > constituentSettings;

    // Time limit
    constituentSettings.push_back( boost::make_shared< propagators::PropagationTimeTerminationSettings >( simulationEndEpoch ) );

    // Altitude limit
    boost::shared_ptr< PropagationTerminationSettings > altitudeTerminationSettings = boost::make_shared< propagators::PropagationDependentVariableTerminationSettings >( boost::make_shared< propagators::SingleDependentVariableSaveSettings >( propagators::altitude_dependent_variable, "Vehicle" ), 100E3, 1 );
    constituentSettings.push_back( altitudeTerminationSettings );

    // Stop if ANY of the two is met
    boost::shared_ptr< PropagationTerminationSettings > terminationSettings = boost::make_shared< propagators::PropagationHybridTerminationSettings >( constituentSettings, 1 );

    // Save dependent variables
    vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesToSave;
    dependentVariablesToSave.push_back( boost::make_shared< SingleDependentVariableSaveSettings >( altitude_dependent_variable, "Vehicle" ) );
    dependentVariablesToSave.push_back( boost::make_shared< SingleDependentVariableSaveSettings >( aerodynamic_moment_coefficients_dependent_variable, "Vehicle" ) );

    boost::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings = boost::make_shared< DependentVariableSaveSettings >( dependentVariablesToSave, 0 ) ;

    // Translational propagator settings
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings = boost::make_shared< TranslationalStatePropagatorSettings< double > > ( centralBodies, accelerationModelMap, bodiesToPropagate, vehicleInitialState, terminationSettings );

    // Create multiple propagator settings
    vector< boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsList;
    propagatorSettingsList.push_back( translationalPropagatorSettings );

    boost::shared_ptr< MultiTypePropagatorSettings< double > > propagatorSettings = boost::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsList, terminationSettings, dependentVariableSaveSettings );

    boost::shared_ptr< IntegratorSettings< > > integratorSettings = boost::make_shared< IntegratorSettings< > > ( rungeKutta4, simulationStartEpoch, 60 );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator( bodyMap, integratorSettings, propagatorSettings, true, false, false );
    map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > dependentVariableOutput = dynamicsSimulator.getDependentVariableHistory( );


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////       CHECK DRAG COEFFICIENTS         ////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Retrieve altitudes and drag coefficients from propagation results

    vector< double > altitudes_propagation_unsorted;
    vector< double > dragCoefficients_propagation_unsorted;
    for ( map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::iterator iter = dependentVariableOutput.begin(); iter != dependentVariableOutput.end(); ++iter ) {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > value = iter->second;
        altitudes_propagation_unsorted.push_back( value[0] );
        dragCoefficients_propagation_unsorted.push_back( value[1] );
    }


    // Sort altitude and drag coefficients by altitude for proper interpolation

    vector< pair< double, double > > altitudeDragPairs;
    for ( int i = 0; i < altitudes_propagation_unsorted.size(); i++ ) {
        altitudeDragPairs.push_back( { altitudes_propagation_unsorted[i] , dragCoefficients_propagation_unsorted[i] } );
    }

    auto sort_by_scores = [](const pair<double,double>& _lhs, const pair<double,double>& _rhs) { return _lhs.first < _rhs.first; };
    sort( altitudeDragPairs.begin(), altitudeDragPairs.end(), sort_by_scores );

    vector< double > altitudes_propagation;
    vector< double > dragCoefficients_propagation;
    for ( int i = 0; i < altitudeDragPairs.size(); i++ ) {
        altitudes_propagation.push_back( altitudeDragPairs[i].first );
        dragCoefficients_propagation.push_back( altitudeDragPairs[i].second );
    }


    // Compare to original input from txt file using interpolators

    LinearInterpolatorDoublePointer interpolator = boost::make_shared<LinearInterpolatorDouble>( LinearInterpolatorDouble( altitudes, dragCoefficients ) );

    LinearInterpolatorDoublePointer interpolator_propagation = boost::make_shared<LinearInterpolatorDouble>( altitudes_propagation, dragCoefficients_propagation );

    BOOST_CHECK_CLOSE_FRACTION( interpolator->interpolate( h_p ), interpolator_propagation->interpolate( h_p ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( interpolator->interpolate( 2*h_p ), interpolator_propagation->interpolate( 2*h_p ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( interpolator->interpolate( h_a ), interpolator_propagation->interpolate( h_a ), tolerance );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
