/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
//
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>




#include <memory>

#include "tudat/astro/aerodynamics/testApolloCapsuleCoefficients.h"
#include "tudat/astro/basic_astro/sphericalStateConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ephemerides/directionBasedRotationalEphemeris.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/basics/testMacros.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/multiDimensionalArrayReader.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/propagation_setup/createMassRateModels.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createSystemModel.h"
#include <limits>
#include <string>

#include <Eigen/Core>

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_body_mass_properties )



double dummyMassFunction( const double time )
{
    return 5.0E3 - 0.5 * ( time - 1.0E7 ) - 0.5 * std::sin( ( time - 1.0E7 ) / 300.0 );
}

Eigen::Vector3d dummyCenterOfMassFunction( const double time )
{
    return ( Eigen::Vector3d( ) << 1.0, 2.0, 3.0 ).finished( ) - ( Eigen::Vector3d( ) << 0.1, 0.3, 0.5 ).finished( ) / 1000.0 * ( time - 1.0E7 );
}

Eigen::Matrix3d dummyInertiaTensorFunction( const double time )
{
    return 3.0 * Eigen::Matrix3d::Identity( ) - ( Eigen::Vector3d::Ones( ) *  Eigen::Vector3d::Ones( ).transpose( ) ) / 1000.0 * ( time - 1.0E7 ) ;
}

Eigen::Vector3d dummyCenterOfMassFunction2( const double mass )
{
    return ( Eigen::Vector3d( ) << 1.0, 2.0, 3.0 ).finished( ) - ( Eigen::Vector3d( ) << 0.1, 0.3, 0.5 ).finished( ) * ( 5.0E3 - mass ) / 1000.0 ;
}

Eigen::Matrix3d dummyInertiaTensorFunction2( const double mass )
{
    return 3.0 * Eigen::Matrix3d::Identity( ) - ( Eigen::Vector3d::Ones( ) *  Eigen::Vector3d::Ones( ).transpose( ) ) * ( 5.0E3 - mass ) / 1000.0 ;
}


BOOST_AUTO_TEST_CASE( testDirectBodyMassProperties )
{
    using namespace tudat;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;
    using namespace orbital_element_conversions;


    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );



    for( unsigned int propagateMass = 0; propagateMass < 2; propagateMass++ )
    {
        unsigned int maximumNumberOfRuns = 2;
        if( propagateMass == 1 )
        {
            maximumNumberOfRuns = 3;
        }
        for( unsigned int i = 0; i < maximumNumberOfRuns; i++ )
        {
            double initialTime = 1.0E7;
            // Create vehicle objects.
            double vehicleMass = 5.0E3;
            double dryVehicleMass = 2.0E3;

            // Define simulation body settings.
            BodyListSettings bodySettings =
                getDefaultBodySettings( { "Earth", "Moon", "Sun" }, "Earth", "ECLIPJ2000" );

            // Create Earth object
            simulation_setup::SystemOfBodies bodies = simulation_setup::createSystemOfBodies( bodySettings );


            bodies.createEmptyBody( "Vehicle" );

            std::shared_ptr< BodyMassPropertiesSettings > bodyMassProperties;
            if ( i == 0 )
            {
                bodyMassProperties = std::make_shared< ConstantBodyMassPropertiesSettings >(
                    vehicleMass, ( Eigen::Vector3d( ) << 1.0, 2.0, 3.0 ).finished( ),
                    3.0 * Eigen::Matrix3d::Identity( ) );
            }
            else if( i == 1 )
            {
                bodyMassProperties = std::make_shared< FromFunctionBodyMassPropertiesSettings >(
                    &dummyMassFunction,
                    &dummyCenterOfMassFunction,
                    &dummyInertiaTensorFunction );
            }
            else if( i == 2 )
            {
                bodyMassProperties = std::make_shared< MassDependentMassDistributionSettings >(
                    vehicleMass,
                    &dummyCenterOfMassFunction2,
                    &dummyInertiaTensorFunction2 );
            }
            addBodyMassProperties(
                bodies, "Vehicle", bodyMassProperties );

            double thrustMagnitude1 = 1.0E3;
            double specificImpulse1 = 250.0;
            double massFlow1 = propulsion::computePropellantMassRateFromSpecificImpulse(
                thrustMagnitude1, specificImpulse1 );

            addEngineModel(
                "Vehicle", "Engine",
                std::make_shared<ConstantThrustMagnitudeSettings>(
                    thrustMagnitude1, specificImpulse1 ), bodies );

            bodies.at( "Vehicle" )->setRotationalEphemeris(
                createRotationModel(
                    orbitalStateBasedRotationSettings( "Earth", true, false, "ECLIPJ2000", "VehicleFixed" ),
                    "Vehicle", bodies ));

            // Define propagator settings variables.
            SelectedAccelerationMap accelerationMap;
            std::vector<std::string> bodiesToPropagate;
            std::vector<std::string> centralBodies;

            std::map<std::string, std::vector<std::shared_ptr<AccelerationSettings> > > accelerationsOfVehicle;
            accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared<AccelerationSettings>(
                basic_astrodynamics::point_mass_gravity ));
            accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared<AccelerationSettings>(
                basic_astrodynamics::point_mass_gravity ));
            accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared<AccelerationSettings>(
                basic_astrodynamics::point_mass_gravity ));
            accelerationsOfVehicle[ "Vehicle" ].push_back( std::make_shared<ThrustAccelerationSettings>(
                "Engine" ));

            accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

            bodiesToPropagate.push_back( "Vehicle" );
            centralBodies.push_back( "Earth" );

            // Create acceleration models and propagation settings.
            basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToPropagate, centralBodies );

            // Set Keplerian elements for Vehicle.
            Eigen::Vector6d vehicleInitialStateInKeplerianElements;
            vehicleInitialStateInKeplerianElements( semiMajorAxisIndex ) = 8000.0E3;
            vehicleInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
            vehicleInitialStateInKeplerianElements( inclinationIndex ) =
                unit_conversions::convertDegreesToRadians( 85.3 );
            vehicleInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                = unit_conversions::convertDegreesToRadians( 235.7 );
            vehicleInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                = unit_conversions::convertDegreesToRadians( 23.4 );
            vehicleInitialStateInKeplerianElements( trueAnomalyIndex ) =
                unit_conversions::convertDegreesToRadians( 139.87 );

            double earthGravitationalParameter =
                bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
            const Eigen::Vector6d vehicleInitialState = convertKeplerianToCartesianElements(
                vehicleInitialStateInKeplerianElements, earthGravitationalParameter );

            std::vector<std::shared_ptr<SingleDependentVariableSaveSettings> > dependentVariablesToSave;
            dependentVariablesToSave.push_back(
                std::make_shared<SingleDependentVariableSaveSettings>(
                    current_body_mass_dependent_variable, "Vehicle" ));
            dependentVariablesToSave.push_back(
                std::make_shared<SingleDependentVariableSaveSettings>(
                    body_center_of_mass, "Vehicle" ));
            dependentVariablesToSave.push_back(
                std::make_shared<SingleDependentVariableSaveSettings>(
                    body_inertia_tensor, "Vehicle" ));

            std::shared_ptr<IntegratorSettings<> > integratorSettings =
                std::make_shared<IntegratorSettings<> >
                    ( rungeKutta4, 0.0, 1.0 );
            std::shared_ptr<PropagationTimeTerminationSettings> terminationSettings =
                std::make_shared<propagators::PropagationTimeTerminationSettings>( initialTime + 1000.0 );
            std::shared_ptr<TranslationalStatePropagatorSettings<double> > translationalPropagatorSettings =
                std::make_shared<TranslationalStatePropagatorSettings<double> >
                    ( centralBodies, accelerationModelMap, bodiesToPropagate, vehicleInitialState, initialTime,
                      integratorSettings, terminationSettings,
                      cowell, dependentVariablesToSave );


            std::map<std::string, std::vector<std::shared_ptr<basic_astrodynamics::MassRateModel> > > massRateModels;
            massRateModels[ "Vehicle" ].push_back(
                createMassRateModel( "Vehicle", std::make_shared<FromThrustMassRateSettings>( 1 ),
                                     bodies, accelerationModelMap ));

            std::shared_ptr<SingleArcPropagatorSettings<double> > massPropagatorSettings =
                std::make_shared<MassPropagatorSettings<double> >(
                    std::vector<std::string>{ "Vehicle" }, massRateModels,
                    ( Eigen::Matrix<double, 1, 1>( ) << vehicleMass ).finished( ), initialTime, integratorSettings,
                    terminationSettings,
                    dependentVariablesToSave );

            std::vector<std::shared_ptr<SingleArcPropagatorSettings<double> > > propagatorSettingsVector;
            propagatorSettingsVector.push_back( translationalPropagatorSettings );
            if( propagateMass == 1 )
            {
                propagatorSettingsVector.push_back( massPropagatorSettings );
            }

            std::shared_ptr<SingleArcPropagatorSettings<double> > propagatorSettings =
                std::make_shared<MultiTypePropagatorSettings<double> >( propagatorSettingsVector, integratorSettings,
                                                                        initialTime, terminationSettings,
                                                                        dependentVariablesToSave );

            // Create simulation object and propagate dynamics.
            auto dynamicsSimulator = std::dynamic_pointer_cast<SingleArcDynamicsSimulator<> >(
                createDynamicsSimulator<double, double>(
                    bodies, propagatorSettings ));


            std::map<double, Eigen::Matrix<double, Eigen::Dynamic, 1> > dependentVariableSolution =
                dynamicsSimulator->getDependentVariableHistory( );
            for ( auto it : dependentVariableSolution )
            {
                double currentMass = it.second( 0 );
                Eigen::Vector3d currentCenterOfMass = it.second.segment( 1, 3 );
                Eigen::Matrix3d currentInertiaTensor = getMatrixFromVectorRotationRepresentation( it.second.segment( 4, 9 ) );
                if( propagateMass == 0 )
                {
                    if( i == 0 )
                    {
                        BOOST_CHECK_CLOSE_FRACTION( currentMass, vehicleMass, std::numeric_limits<double>::epsilon( ));
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentCenterOfMass,
                                                           (( Eigen::Vector3d( ) << 1.0, 2.0, 3.0 ).finished( )),
                                                           std::numeric_limits<double>::epsilon( ));
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentInertiaTensor, ( 3.0 * Eigen::Matrix3d::Identity( )),
                                                           std::numeric_limits<double>::epsilon( ));
                    }
                    else if( i == 1 )
                    {
                        BOOST_CHECK_CLOSE_FRACTION( currentMass, dummyMassFunction( it.first ), std::numeric_limits<double>::epsilon( ));
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentCenterOfMass,
                                                           ( dummyCenterOfMassFunction( it.first ) ),
                                                           std::numeric_limits<double>::epsilon( ));
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentInertiaTensor,
                                                           ( dummyInertiaTensorFunction( it.first ) ),
                                                           std::numeric_limits<double>::epsilon( ));
                    }
                }
                else if( propagateMass == 1 )
                {
                    BOOST_CHECK_CLOSE_FRACTION( currentMass, vehicleMass - ( it.first - initialTime ) * massFlow1,
                                                1.0E3 * std::numeric_limits<double>::epsilon( ));
                    if( i == 0 )
                    {
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentCenterOfMass,
                                                           (( Eigen::Vector3d( ) << 1.0, 2.0, 3.0 ).finished( )),
                                                           std::numeric_limits<double>::epsilon( ));
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentInertiaTensor, ( 3.0 * Eigen::Matrix3d::Identity( )),
                                                           std::numeric_limits<double>::epsilon( ));
                    }
                    else if( i == 1 )
                    {
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentCenterOfMass,
                                                           ( dummyCenterOfMassFunction( it.first ) ),
                                                           std::numeric_limits<double>::epsilon( ));
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentInertiaTensor,
                                                           ( dummyInertiaTensorFunction( it.first ) ),
                                                           std::numeric_limits<double>::epsilon( ));
                    }
                    else if( i == 2 )
                    {
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentCenterOfMass,
                                                           ( dummyCenterOfMassFunction2( currentMass ) ),
                                                           std::numeric_limits<double>::epsilon( ));
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( currentInertiaTensor,
                                                           ( dummyInertiaTensorFunction2( currentMass ) ),
                                                           std::numeric_limits<double>::epsilon( ));
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
