
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
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include <Tudat/Basics/testMacros.h>
#include <Tudat/Astrodynamics/Propagators/dynamicsSimulator.h>
#include <Tudat/SimulationSetup/body.h>
#include <Tudat/SimulationSetup/createAccelerationModels.h>
#include <Tudat/SimulationSetup/createMassRateModels.h>
#include <Tudat/SimulationSetup/defaultBodies.h>

#include <iostream>
#include <limits>
#include <string>

#include <Eigen/Core>

namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_thrust_acceleration )

BOOST_AUTO_TEST_CASE( testConstantThrustAcceleration )
{
    using namespace tudat;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace basic_astrodynamics;
    using namespace propagators;
    using namespace basic_mathematics;
    using namespace basic_astrodynamics;

    // Create Earth object
    simulation_setup::NamedBodyMap bodyMap;

    // Create vehicle objects.
    double vehicleMass = 5.0E3;
    bodyMap[ "Vehicle" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Vehicle" ]->setConstantBodyMass( vehicleMass );
    bodyMap[ "Vehicle" ]->setEphemeris(
                boost::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                    boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, basic_mathematics::Vector6d  > >( ),
                    "SSB" ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    Eigen::Vector3d thrustDirection;
    thrustDirection << -1.4, 2.4, 5,6;

    double thrustMagnitude = 1.0E3;
    double specificImpulse = 250.0;

    double massRate = thrustMagnitude / ( specificImpulse * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION );
    // Define acceleration model settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Vehicle" ].push_back( boost::make_shared< ThrustAccelerationSettings >(
                                                    boost::make_shared< CustomThrustDirectionSettings >(
                                                        boost::lambda::constant( thrustDirection ) ),
                                                boost::make_shared< ConstantThrustMagnitudeSettings >(
                                                    thrustMagnitude, specificImpulse ) ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    bodiesToPropagate.push_back( "Vehicle" );
    centralBodies.push_back( "SSB" );

    // Set initial state
    basic_mathematics::Vector6d systemInitialState = basic_mathematics::Vector6d::Zero( );

    // Create acceleration models and propagation settings.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    boost::shared_ptr< PropagationTimeTerminationSettings > terminationSettings =
            boost::make_shared< propagators::PropagationTimeTerminationSettings >( 1000.0 );
    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > translationalPropagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, systemInitialState, terminationSettings );
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, 0.0, 0.1 );

    for( unsigned int i = 0; i < 2; i++ )
    {
        if( i == 0 )
        {
            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodyMap, integratorSettings, translationalPropagatorSettings, true, false, false );

            // Retrieve numerical solutions for state and dependent variables
            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
                    dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            Eigen::Vector3d constantAcceleration = thrustDirection.normalized( ) * thrustMagnitude / vehicleMass;
            for( std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::const_iterator outputIterator =
                 numericalSolution.begin( ); outputIterator != numericalSolution.end( ); outputIterator++ )
            {
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            ( outputIterator->second.segment( 0, 3 ) ),
                            ( 0.5 * constantAcceleration * std::pow( outputIterator->first, 2.0 ) ), 1.0E-12 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            ( outputIterator->second.segment( 3, 3 ) ),
                            ( constantAcceleration * outputIterator->first ), 1.0E-12 );
            }
        }
        else if( i == 1 )
        {
            std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
            massRateModels[ "Vehicle" ] = (
                        createMassRateModel( "Vehicle", boost::make_shared< FromThrustMassModelSettings >( 1 ),
                                                            bodyMap, accelerationModelMap ) );

            boost::shared_ptr< PropagatorSettings< double > > massPropagatorSettings =
                    boost::make_shared< MassPropagatorSettings< double > >(
                        boost::assign::list_of( "Vehicle" ), massRateModels,
                        ( Eigen::Matrix< double, 1, 1 >( ) << vehicleMass ).finished( ), terminationSettings );

            std::vector< boost::shared_ptr< PropagatorSettings< double > > > propagatorSettingsVector;
            propagatorSettingsVector.push_back( translationalPropagatorSettings );
            propagatorSettingsVector.push_back( massPropagatorSettings );

            boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
                    boost::make_shared< MultiTypePropagatorSettings< double > >( propagatorSettingsVector, terminationSettings );

            // Create simulation object and propagate dynamics.
            SingleArcDynamicsSimulator< > dynamicsSimulator(
                        bodyMap, integratorSettings, propagatorSettings, true, false, false );

            std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > > numericalSolution =
                    dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

            for( std::map< double, Eigen::Matrix< double, Eigen::Dynamic, 1 > >::const_iterator outputIterator =
                 numericalSolution.begin( ); outputIterator != numericalSolution.end( ); outputIterator++ )
            {
                double currentMass = vehicleMass - outputIterator->first * massRate;

                Eigen::Vector3d currentVelocity = thrustDirection.normalized( ) *
                        specificImpulse * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION *
                        std::log( vehicleMass / currentMass );
                BOOST_CHECK_CLOSE_FRACTION( outputIterator->second( 6 ), currentMass, 1.0E-12 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            ( outputIterator->second.segment( 3, 3 ) ), currentVelocity, 1.0E-12 );

            }
        }
    }

}

BOOST_AUTO_TEST_SUITE_END( )


}

}


