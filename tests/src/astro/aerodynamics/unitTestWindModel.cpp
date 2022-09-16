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

#include <limits>

#include <boost/make_shared.hpp>
#include <memory>
#include <boost/lambda/lambda.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/aerodynamics/customAerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/aerodynamicAcceleration.h"
#include "tudat/astro/reference_frames/aerodynamicAngleCalculator.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

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

    for( int test = 0; test < 2; test++ )
    {
        // Create Earth object
        BodyListSettings defaultBodySettings =
                getDefaultBodySettings( { "Earth" }, -1.0E6, 1.0E6 );
        defaultBodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >(
                    Eigen::Vector6d::Zero( ) );
        if( test == 0 )
        {
            defaultBodySettings.at( "Earth" )->atmosphereSettings->setWindSettings(
                        std::make_shared< CustomWindModelSettings >(
                            std::bind( &getCustomWindVector,
                                       std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4 ),
                            reference_frames::corotating_frame ) );
        }
        else
        {
            defaultBodySettings.at( "Earth" )->atmosphereSettings->setWindSettings(
                        std::make_shared< CustomWindModelSettings >(
                            std::bind( &getCustomWindVector,
                                       std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4 ),
                            reference_frames::vertical_frame ) );
        }
        SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

        // Create vehicle object.
        double vehicleMass = 5.0E3;
        bodies.createEmptyBody( "Vehicle" );
        bodies.at( "Vehicle" )->setConstantBodyMass( vehicleMass );

        // Set aerodynamic coefficients.
        std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
                std::make_shared< ConstantAerodynamicCoefficientSettings >(
                    2.0, 4.0, 1.5, Eigen::Vector3d::Zero( ), Eigen::Vector3d::UnitX( ), Eigen::Vector3d::Zero( ), 1, 1 );
        bodies.at( "Vehicle" )->setAerodynamicCoefficientInterface(
                    createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Vehicle" ) );

        // Define acceleration model settings.
        SelectedAccelerationMap accelerationMap;
        std::vector< std::string > bodiesToPropagate;
        std::vector< std::string > centralBodies;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
        accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
        accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
        accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;
        bodiesToPropagate.push_back( "Vehicle" );
        centralBodies.push_back( "Earth" );

        // Set initial state
        Eigen::Vector6d systemInitialState = Eigen::Vector6d::Zero( );
        systemInitialState( 0 ) = 6.8E6;
        systemInitialState( 4 ) = 7.5E3;

        // Create acceleration models and propagation settings.
        basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodies, accelerationMap, bodiesToPropagate, centralBodies );

        // Set variables to save
        std::shared_ptr< DependentVariableSaveSettings > dependentVariableSaveSettings;
        std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        altitude_dependent_variable, "Vehicle", "Earth" ) );
        dependentVariables.push_back(
                    std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                        "Vehicle", reference_frames::longitude_angle ) );
        dependentVariables.push_back(
                    std::make_shared< BodyAerodynamicAngleVariableSaveSettings >(
                        "Vehicle", reference_frames::latitude_angle ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        local_density_dependent_variable, "Vehicle", "Earth" ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        airspeed_dependent_variable, "Vehicle", "Earth" ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        body_fixed_airspeed_based_velocity_variable, "Vehicle" ) );
        dependentVariables.push_back(
                    std::make_shared< SingleDependentVariableSaveSettings >(
                        body_fixed_groundspeed_based_velocity_variable, "Vehicle" ) );
        dependentVariables.push_back(
                    std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                        aerodynamic, "Vehicle", "Earth", 0 ) );
        dependentVariableSaveSettings = std::make_shared< DependentVariableSaveSettings >( dependentVariables );

        // Set propagation/integration settings
        std::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
                std::make_shared< propagators::PropagationTimeTerminationSettings >( 1000.0 );
        std::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
                std::make_shared< TranslationalStatePropagatorSettings< double > >
                ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings,
                  cowell, dependentVariableSaveSettings );
        std::shared_ptr< IntegratorSettings< > > integratorSettings =
                std::make_shared< IntegratorSettings< > >
                ( rungeKutta4, 0.0, 5.0 );

        // Create simulation object and propagate dynamics.
        SingleArcDynamicsSimulator< > dynamicsSimulator(
                    bodies, integratorSettings, translationalPropagatorSettings, true, false, false );

        // Retrieve numerical solutions for state and dependent variables
        std::map< double, Eigen::VectorXd > numericalSolution =
                dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
        std::map< double, Eigen::VectorXd > dependentVariableOutput =
                dynamicsSimulator.getDependentVariableHistory( );


        // Get function to transform aerodynamic coefficients to inertial frame
        std::function< Eigen::Vector3d( const Eigen::Vector3d& ) > toPropagationFrameTransformation;
        toPropagationFrameTransformation =
                reference_frames::getAerodynamicForceTransformationFunction(
                    bodies.at( "Vehicle" )->getFlightConditions( )->getAerodynamicAngleCalculator( ),
                    reference_frames::aerodynamic_frame,
                    std::bind( &Body::getCurrentRotationToGlobalFrame, bodies.at( "Earth" ) ),
                    reference_frames::inertial_frame );

        // Iterate over all time steps and compute test quantities
        for( std::map< double, Eigen::VectorXd >::const_iterator dataIterator = dependentVariableOutput.begin( );
             dataIterator != dependentVariableOutput.end( ); dataIterator++ )
        {
            // Update environment
            bodies.at( "Vehicle" )->setState( numericalSolution.at( dataIterator->first ) );
            bodies.at( "Earth" )->setState( Eigen::Vector6d::Zero( ) );
            bodies.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( dataIterator->first );
            bodies.at( "Vehicle" )->getFlightConditions( )->updateConditions( dataIterator->first );

            // Retrieve dependent variables
            double altitude = dataIterator->second( 0 );
            double longitude = dataIterator->second( 1 );
            double latitude = dataIterator->second( 2 );
            double localDensity = dataIterator->second( 3 );
            double airspeed = dataIterator->second( 4 );
            Eigen::Vector3d airspeedBasedVelocity = dataIterator->second.segment( 5, 3 );
            Eigen::Vector3d groundspeedBasedVelocity = dataIterator->second.segment( 8, 3 );
            Eigen::Vector3d aerodynamicAcceleration = dataIterator->second.segment( 11, 3 );
            Eigen::Quaterniond rotationToVerticalFrame =
                    reference_frames::getLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion(
                        longitude, latitude ).inverse( );

            // Compute wind vector from dependent variables and directly from function.
            Eigen::Vector3d expectedWindVelocity = getCustomWindVector(
                        altitude, longitude, latitude, dataIterator->first );
            Eigen::Vector3d computedWindVector = groundspeedBasedVelocity - airspeedBasedVelocity;
            if( test == 1 )
            {
                computedWindVector = rotationToVerticalFrame * computedWindVector;
            }

            // Manually compute aerodynamic acceleration vector
            Eigen::Vector3d airSpeedVelocityUnitVectorInInertialFrame =
                    reference_frames::transformVectorFunctionFromVectorFunctions(
                        [ & ]( ){ return Eigen::Vector3d::UnitX( ); }, toPropagationFrameTransformation );
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
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
