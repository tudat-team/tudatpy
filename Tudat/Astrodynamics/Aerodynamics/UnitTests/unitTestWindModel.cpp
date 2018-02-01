/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#define BOOST_TEST_MAIN

#include <limits>

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Astrodynamics/Aerodynamics/customAerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicAcceleration.h"
#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"
#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#include "Tudat/External/SpiceInterface/spiceRotationalEphemeris.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/PropagationSetup/createNumericalSimulator.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_wind_models )

//! Function to compute a (non-physical) wind vector
Eigen::Vector3d getCustomWindVector(
        const double altitude, const double longitude, const double latitude, const double time )
{
   return ( Eigen::Vector3d( ) << 200.0, -120.0, 75.0 ).finished( ) *
            ( longitude / ( 2.0 * mathematical_constants::PI ) ) *
              ( latitude / mathematical_constants::PI ) * ( altitude - 250.0E3 ) / 10E3 * ( time - 500.0 ) / 500.0;
}

BOOST_AUTO_TEST_CASE( testWindModelInPropagation )
{
    using namespace tudat;
    using namespace aerodynamics;
    using namespace simulation_setup;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create Earth object
    std::map< std::string, boost::shared_ptr< BodySettings > > defaultBodySettings =
            getDefaultBodySettings( boost::assign::list_of( "Earth" ), -1.0E6, 1.0E6 );
    defaultBodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< ConstantEphemerisSettings >(
                Eigen::Vector6d::Zero( ) );
    defaultBodySettings[ "Earth" ]->atmosphereSettings->setWindSettings(
                boost::make_shared< CustomWindModelSettings >(
                    boost::bind( &getCustomWindVector, _1, _2, _3, _4 ) ) );
    NamedBodyMap bodyMap = createBodies( defaultBodySettings );

    // Create vehicle object.
    double vehicleMass = 5.0E3;
    bodyMap[ "Vehicle" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );

    // Set aerodynamic coefficients.
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                2.0, 4.0, 1.5, Eigen::Vector3d::Zero( ), Eigen::Vector3d::UnitX( ), Eigen::Vector3d::Zero( ), 1, 1 );
    bodyMap[ "Vehicle" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle" ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Define acceleration model settings.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfVehicle[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Set initial state
    Eigen::Vector6d systemInitialState = Eigen::Vector6d::Zero( );
    systemInitialState( 0 ) = 6.8E6;
    systemInitialState( 4 ) = 7.5E3;

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    // Set variables to save
    boost::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings;
    std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    altitude_dependent_variable, "Vehicle", "Earth" ) );
    dependentVariables.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "Vehicle", reference_frames::longitude_angle ) );
    dependentVariables.push_back(
                boost::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                    "Vehicle", reference_frames::latitude_angle ) );
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    local_density_dependent_variable, "Vehicle", "Earth" ) );
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    airspeed_dependent_variable, "Vehicle", "Earth" ) );
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_airspeed_based_velocity_variable, "Vehicle" ) );
    dependentVariables.push_back(
                boost::make_shared< SingleDependentVariableSaveSettings >(
                    body_fixed_groundspeed_based_velocity_variable, "Vehicle" ) );
    dependentVariables.push_back(
                boost::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    aerodynamic, "Vehicle", "Earth", 0 ) );
    dependentVariableSaveSettings = boost::make_shared< DependentVariableSaveSettings >( dependentVariables );

    // Set propagation/integration settings
    boost::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            boost::make_shared< propagators::PropagationTimeTerminationSettings >( 1000.0 );
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings,
              cowell, dependentVariableSaveSettings );
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, 0.0, 5.0 );

    // Create simulation object and propagate dynamics.
    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, translationalPropagatorSettings, true, false, false );

    // Retrieve numerical solutions for state and dependent variables
    std::map< double, Eigen::VectorXd > numericalSolution =
            dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableOutput =
            dynamicsSimulator.getDependentVariableHistory( );


    // Get function to transform aerodynamic coefficients to inertial frame
    boost::function< Eigen::Vector3d( const Eigen::Vector3d& ) > toPropagationFrameTransformation;
    toPropagationFrameTransformation =
            reference_frames::getAerodynamicForceTransformationFunction(
                bodyMap.at( "Vehicle" )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                reference_frames::aerodynamic_frame,
                boost::bind( &Body::getCurrentRotationToGlobalFrame, bodyMap.at( "Earth" ) ),
                reference_frames::inertial_frame );

    // Iterate over all time steps and compute test quantities
    for( std::map< double, Eigen::VectorXd >::const_iterator dataIterator = dependentVariableOutput.begin( );
         dataIterator != dependentVariableOutput.end( ); dataIterator++ )
    {
        // Update environment
        bodyMap.at( "Vehicle" )->setState( numericalSolution.at( dataIterator->first ) );
        bodyMap.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( dataIterator->first );
        bodyMap.at( "Vehicle" )->getFlightConditions( )->updateConditions( dataIterator->first );

        // Retrieve dependent variables
        double altitude = dataIterator->second( 0 );
        double longitude = dataIterator->second( 1 );
        double latitude = dataIterator->second( 2 );
        double localDensity = dataIterator->second( 3 );
        double airspeed = dataIterator->second( 4 );
        Eigen::Vector3d airspeedBasedVelocity = dataIterator->second.segment( 5, 3 );
        Eigen::Vector3d groundspeedBasedVelocity = dataIterator->second.segment( 8, 3 );
        Eigen::Vector3d aerodynamicAcceleration = dataIterator->second.segment( 11, 3 );

        // Compute wind vector from dependent variables and directly from function.
        Eigen::Vector3d expectedWindVelocity = getCustomWindVector(
                    altitude, longitude, latitude, dataIterator->first );
        Eigen::Vector3d computedWindVector = airspeedBasedVelocity - groundspeedBasedVelocity;

        // Manually compute aerodynamic acceleration vector
        Eigen::Vector3d airSpeedVelocityUnitVectorInInertialFrame =
                reference_frames::transformVectorFunctionFromVectorFunctions(
                    boost::lambda::constant( Eigen::Vector3d::UnitX( ) ), toPropagationFrameTransformation );
        Eigen::Vector3d expectedAerodynamicAcceleration = -0.5 * localDensity * airspeed * airspeed *
                airSpeedVelocityUnitVectorInInertialFrame * 4.0 / vehicleMass;

        for( unsigned int i = 0; i < 3; i++ )
        {
            // Test wind velocity vector
            BOOST_CHECK_SMALL( std::fabs( expectedWindVelocity( i ) - computedWindVector( i ) ),
                               systemInitialState( 4 ) * 5.0 * std::numeric_limits< double >::epsilon( ) );

            // Test aerodynamic acceleration unit vector
            BOOST_CHECK_SMALL( std::fabs( airSpeedVelocityUnitVectorInInertialFrame.normalized( )( i ) +
                                          aerodynamicAcceleration.normalized( )( i ) ),
                               5.0 * std::numeric_limits< double >::epsilon( ) );

            // Test aerodynamic acceleration
            BOOST_CHECK_SMALL( std::fabs( aerodynamicAcceleration( i ) - expectedAerodynamicAcceleration( i ) ),
                               aerodynamicAcceleration.norm( ) * 5.0 * std::numeric_limits< double >::epsilon( ) );
        }

        // Test airspeed
        BOOST_CHECK_SMALL( std::fabs( airspeedBasedVelocity.norm( ) - airspeed ),
                           systemInitialState( 4 ) * 5.0 * std::numeric_limits< double >::epsilon( ) );


    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
