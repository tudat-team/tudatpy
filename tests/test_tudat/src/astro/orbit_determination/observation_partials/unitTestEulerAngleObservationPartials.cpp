/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/support/numericalObservationPartial.h"
#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

#include "tudat/support/observationPartialTestFunctions.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;

BOOST_AUTO_TEST_SUITE( test_euler_angle_partials )

//! Test partial derivatives of Euler angle observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testEulerAnglePartials )
{
    std::cout << " **************************************************************************************** " << std::endl;

    // Test partials with constant rotational ephemerides (with test of rotation state partial)
    {
        // Create environment
        SystemOfBodies bodies =
                setupEnvironment( std::vector< std::pair< std::string, std::string > >( ), 1.0E7, 1.2E7, 1.1E7, false, 1.0, true );

        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ observed_body ] = std::make_pair< std::string, std::string >( "Mars", "" );

        // Generate one-way range model
        std::shared_ptr< ObservationModel< 3 > > eulerAngleModel =
                observation_models::ObservationModelCreator< 3, double, double >::createObservationModel(
                        std::make_shared< observation_models::ObservationModelSettings >( euler_angle_313_observable, linkEnds ), bodies );

        // Create parameter objects.
        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                createEstimatableParameters( bodies, 1.1E7, false, true );

        testObservationPartials< 3 >(
                eulerAngleModel, bodies, fullEstimatableParameterSet, linkEnds, euler_angle_313_observable, 1.0E-6, false, false );
    }
}

//! Test all components of the partial w.r.t. the rotational state, created for the observation model, against central differences of the
//! observation model, for orientations with all combinations of the Euler angles below, and both signs of the quaternion.
BOOST_AUTO_TEST_CASE( testEulerAnglePartialWrtRotationalState )
{
    const double observationTime = 1.1E7;
    SystemOfBodies bodies =
            setupEnvironment( std::vector< std::pair< std::string, std::string > >( ), 1.0E7, 1.2E7, observationTime, false, 1.0, true );
    std::shared_ptr< ConstantRotationalEphemeris > rotationModel =
            std::dynamic_pointer_cast< ConstantRotationalEphemeris >( bodies.at( "Mars" )->getRotationalEphemeris( ) );

    LinkEnds linkEnds;
    linkEnds[ observed_body ] = LinkEndId( "Mars" );
    std::shared_ptr< ObservationModel< 3 > > eulerAngleModel = ObservationModelCreator< 3, double, double >::createObservationModel(
            std::make_shared< ObservationModelSettings >( euler_angle_313_observable, linkEnds ), bodies );

    // Translational states of Earth and Mars, and rotational state of Mars
    std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 3 > > > partials =
            ObservationPartialCreator< 3, double, double >::createObservationPartials(
                    eulerAngleModel, bodies, createEstimatableParameters( bodies, observationTime, false, true ) )
                    .first;

    // Only the rotational state of Mars
    BOOST_CHECK_EQUAL( partials.size( ), 1 );
    BOOST_CHECK_EQUAL( partials.begin( )->first.second, 7 );
    std::shared_ptr< ObservationPartial< 3 > > partialObject = partials.begin( )->second;
    BOOST_CHECK( partialObject->getParameterIdentifier( ).first == initial_rotational_body_state );

    // Quaternion perturbation [-]
    const double perturbation = 1.0E-6;
    for( double firstAngle : { 0.7, 2.6, -2.2 } )
    {
        for( double secondAngle : { 0.3, 1.4, 2.8 } )
        {
            for( double thirdAngle : { -1.9, 0.5 } )
            {
                for( double sign : { 1.0, -1.0 } )
                {
                    Eigen::Vector7d state = Eigen::Vector7d::Zero( );
                    state.segment( 0, 4 ) = sign *
                            linear_algebra::convertQuaternionToVectorFormat(
                                                    Eigen::Quaterniond( Eigen::AngleAxisd( firstAngle, Eigen::Vector3d::UnitZ( ) ) *
                                                                        Eigen::AngleAxisd( secondAngle, Eigen::Vector3d::UnitX( ) ) *
                                                                        Eigen::AngleAxisd( thirdAngle, Eigen::Vector3d::UnitZ( ) ) ) );
                    rotationModel->updateConstantState( state );

                    std::vector< double > linkEndTimes;
                    std::vector< Eigen::Vector6d > linkEndStates;
                    const Eigen::Vector3d observation = eulerAngleModel->computeObservationsWithLinkEndData(
                            observationTime, observed_body, linkEndTimes, linkEndStates );
                    const Eigen::Matrix< double, 3, Eigen::Dynamic > partial =
                            partialObject->calculatePartial( linkEndStates, linkEndTimes, observed_body, nullptr, observation )
                                    .at( 0 )
                                    .first;

                    // The state is a unit quaternion: compare along each state entry, projected on the tangent space
                    for( int i = 0; i < 4; i++ )
                    {
                        Eigen::Vector7d direction = Eigen::Vector7d::Zero( );
                        direction( i ) = 1.0;
                        direction.segment( 0, 4 ) -= state( i ) * state.segment( 0, 4 );

                        rotationModel->updateConstantState( state + perturbation * direction );
                        const Eigen::Vector3d upperAngles = eulerAngleModel->computeObservations( observationTime, observed_body );
                        rotationModel->updateConstantState( state - perturbation * direction );
                        const Eigen::Vector3d lowerAngles = eulerAngleModel->computeObservations( observationTime, observed_body );
                        rotationModel->updateConstantState( state );

                        const Eigen::Vector3d numericalChange = ( upperAngles - lowerAngles ) / ( 2.0 * perturbation );
                        const Eigen::Vector3d analyticalChange = partial.block( 0, 0, 3, 4 ) * direction.segment( 0, 4 );
                        BOOST_CHECK_SMALL( ( numericalChange - analyticalChange ).norm( ), 1.0E-8 );
                    }

                    // The observable does not depend on the angular velocity
                    BOOST_CHECK_EQUAL( partial.block( 0, 4, 3, 3 ).norm( ), 0.0 );
                }
            }
        }
    }
}

//! Test all components of the partials w.r.t. rotation model parameters (rotation rate and pole of Mars) against central differences of
//! the observation model, and check that no partials are created for the rotation model parameters of another body (Earth).
BOOST_AUTO_TEST_CASE( testEulerAnglePartialWrtRotationModelParameters )
{
    const double observationTime = 1.1E7;
    SystemOfBodies bodies =
            setupEnvironment( std::vector< std::pair< std::string, std::string > >( ), 1.0E7, 1.2E7, observationTime, false );

    LinkEnds linkEnds;
    linkEnds[ observed_body ] = LinkEndId( "Mars" );
    std::shared_ptr< ObservationModel< 3 > > eulerAngleModel = ObservationModelCreator< 3, double, double >::createObservationModel(
            std::make_shared< ObservationModelSettings >( euler_angle_313_observable, linkEnds ), bodies );

    // Translational states, rotation rates and poles of Earth and Mars
    std::shared_ptr< EstimatableParameterSet< double > > parameters = createEstimatableParameters( bodies, observationTime );
    std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 3 > > > partials =
            ObservationPartialCreator< 3, double, double >::createObservationPartials( eulerAngleModel, bodies, parameters ).first;

    // Only the rotation rate and pole of Mars
    BOOST_CHECK_EQUAL( partials.size( ), 2 );
    const Eigen::VectorXd nominalParameters = parameters->getFullParameterValues< double >( );
    for( auto partialIterator : partials )
    {
        BOOST_CHECK_EQUAL( partialIterator.second->getParameterIdentifier( ).second.first, "Mars" );
        BOOST_CHECK( isParameterRotationMatrixProperty( partialIterator.second->getParameterIdentifier( ).first ) );

        const Eigen::Matrix< double, 3, Eigen::Dynamic > analyticalPartial =
                partialIterator.second
                        ->calculatePartial( std::vector< Eigen::Vector6d >( { Eigen::Vector6d::Zero( ) } ),
                                            std::vector< double >( { observationTime } ),
                                            observed_body )
                        .at( 0 )
                        .first;
        BOOST_CHECK_EQUAL( analyticalPartial.cols( ), partialIterator.first.second );

        // Rotation rate [rad/s] and pole right ascension and declination [rad]
        const double perturbation =
                ( partialIterator.second->getParameterIdentifier( ).first == constant_rotation_rate ) ? 1.0E-12 : 1.0E-6;
        for( int i = 0; i < partialIterator.first.second; i++ )
        {
            Eigen::VectorXd perturbedParameters = nominalParameters;
            perturbedParameters( partialIterator.first.first + i ) += perturbation;
            parameters->resetParameterValues( perturbedParameters );
            const Eigen::Vector3d upperAngles = eulerAngleModel->computeObservations( observationTime, observed_body );

            perturbedParameters( partialIterator.first.first + i ) -= 2.0 * perturbation;
            parameters->resetParameterValues( perturbedParameters );
            const Eigen::Vector3d lowerAngles = eulerAngleModel->computeObservations( observationTime, observed_body );
            parameters->resetParameterValues( nominalParameters );

            const Eigen::Vector3d numericalPartial = ( upperAngles - lowerAngles ) / ( 2.0 * perturbation );
            BOOST_CHECK_SMALL( ( analyticalPartial.col( i ) - numericalPartial ).norm( ), 1.0E-6 * numericalPartial.norm( ) );
        }
    }
}

//! Test the partials w.r.t. the translational states of a synchronously rotating body and of its central body, created for the observation
//! model, against central differences of the observation model, and test their indices for arc-wise translational states.
BOOST_AUTO_TEST_CASE( testEulerAnglePartialWrtStateForSynchronousRotation )
{
    spice_interface::loadStandardSpiceKernels( );

    // Constant states of the Earth, at the origin, and the Moon, with approximately their relative state
    Eigen::Vector6d moonState;
    moonState << 3.0E8, -2.0E8, 0.5E8, 600.0, 750.0, -100.0;
    const std::vector< std::string > estimatedBodies = { "Earth", "Moon" };
    BodyListSettings bodySettings = getDefaultBodySettings( estimatedBodies, "SSB", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->ephemerisSettings = constantEphemerisSettings( Eigen::Vector6d::Zero( ), "SSB", "ECLIPJ2000" );
    bodySettings.at( "Moon" )->ephemerisSettings = constantEphemerisSettings( moonState, "SSB", "ECLIPJ2000" );
    bodySettings.at( "Moon" )->rotationModelSettings = synchronousRotationModelSettings( "Earth", "ECLIPJ2000", "IAU_Moon" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    LinkEnds linkEnds;
    linkEnds[ observed_body ] = LinkEndId( "Moon" );
    std::shared_ptr< ObservationModel< 3 > > eulerAngleModel = ObservationModelCreator< 3, double, double >::createObservationModel(
            std::make_shared< ObservationModelSettings >( euler_angle_313_observable, linkEnds ), bodies );

    std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > singleArcStates, arcWiseStates;
    const std::vector< double > arcStartTimes = { 0.0, 1.0E6 };
    const Eigen::VectorXd arcStates = Eigen::VectorXd::Zero( 12 );
    const std::vector< std::string > centralBodies = { "SSB", "SSB" };
    for( const std::string& bodyName : estimatedBodies )
    {
        singleArcStates.push_back( std::make_shared< InitialTranslationalStateParameter< double > >(
                bodyName, bodies.at( bodyName )->getEphemeris( )->getCartesianState( 0.0 ) ) );
        arcWiseStates.push_back( std::make_shared< ArcWiseInitialTranslationalStateParameter< double > >(
                bodyName, arcStartTimes, arcStates, centralBodies ) );
    }
    const std::vector< std::shared_ptr< EstimatableParameter< double > > > noDoubleParameters;
    const std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > noVectorParameters;

    // The index in the state transition matrix does not depend on the number of arcs
    std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 3 > > > arcWisePartials =
            ObservationPartialCreator< 3, double, double >::createObservationPartials(
                    eulerAngleModel,
                    bodies,
                    std::make_shared< EstimatableParameterSet< double > >( noDoubleParameters, noVectorParameters, arcWiseStates ) )
                    .first;
    BOOST_CHECK_EQUAL( arcWisePartials.size( ), 2 );
    BOOST_CHECK_EQUAL( arcWisePartials.count( std::make_pair( 0, 6 ) ), 1 );
    BOOST_CHECK_EQUAL( arcWisePartials.count( std::make_pair( 6, 6 ) ), 1 );

    std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 3 > > > partials =
            ObservationPartialCreator< 3, double, double >::createObservationPartials(
                    eulerAngleModel,
                    bodies,
                    std::make_shared< EstimatableParameterSet< double > >( noDoubleParameters, noVectorParameters, singleArcStates ) )
                    .first;
    BOOST_CHECK_EQUAL( partials.size( ), 2 );

    const double observationTime = 1.0E7;
    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;
    const Eigen::Vector3d observation =
            eulerAngleModel->computeObservationsWithLinkEndData( observationTime, observed_body, linkEndTimes, linkEndStates );

    for( unsigned int bodyIndex = 0; bodyIndex < estimatedBodies.size( ); bodyIndex++ )
    {
        const Eigen::Matrix< double, 3, Eigen::Dynamic > analyticalPartial =
                partials.at( std::make_pair( static_cast< int >( 6 * bodyIndex ), 6 ) )
                        ->calculatePartial( linkEndStates, linkEndTimes, observed_body, nullptr, observation )
                        .at( 0 )
                        .first;

        std::shared_ptr< Body > body = bodies.at( estimatedBodies.at( bodyIndex ) );
        std::shared_ptr< ConstantEphemeris > ephemeris = std::dynamic_pointer_cast< ConstantEphemeris >( body->getEphemeris( ) );
        const Eigen::Vector6d nominalState = ephemeris->getCartesianState( observationTime );

        // Position [m] and velocity [m/s] perturbations
        for( int i = 0; i < 6; i++ )
        {
            const double perturbation = ( i < 3 ) ? 10.0 : 1.0E-3;
            Eigen::Vector6d perturbedState = nominalState;

            perturbedState( i ) = nominalState( i ) + perturbation;
            ephemeris->updateConstantState( perturbedState );
            body->recomputeStateOnNextCall( );
            const Eigen::Vector3d upperAngles = eulerAngleModel->computeObservations( observationTime, observed_body );

            perturbedState( i ) = nominalState( i ) - perturbation;
            ephemeris->updateConstantState( perturbedState );
            body->recomputeStateOnNextCall( );
            const Eigen::Vector3d lowerAngles = eulerAngleModel->computeObservations( observationTime, observed_body );

            ephemeris->updateConstantState( nominalState );
            body->recomputeStateOnNextCall( );

            const Eigen::Vector3d numericalPartial = ( upperAngles - lowerAngles ) / ( 2.0 * perturbation );
            BOOST_CHECK_SMALL( ( analyticalPartial.col( i ) - numericalPartial ).norm( ), 1.0E-6 * numericalPartial.norm( ) );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
