/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <iostream>
#include <string>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"

namespace tudat
{

namespace unit_tests
{

using namespace tudat::estimatable_parameters;
using namespace tudat::simulation_setup;
using namespace tudat::propagators;
using namespace tudat::numerical_integrators;
using namespace tudat::ephemerides;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_astrodynamics;
using namespace tudat::observation_models;

BOOST_AUTO_TEST_SUITE( test_hybrid_arc_parameter_coupling )

//! Function to compare the analytical design matrix of position observations of a multi-arc vehicle with a numerical (central
//! difference) design matrix, for a hybrid-arc estimation with the Moon as single-arc body. The estimated parameters are the
//! initial states, the Earth gravitational parameter and a constant empirical acceleration of the Moon. The latter acts only on
//! the Moon, so that its influence on the vehicle is entirely through the Moon's trajectory.
//! \param vehicleCentralBody Central body of the vehicle propagation (Earth: Moon is a perturbing body; Moon: Moon is the central
//! body)
//! \return Relative difference ||H_analytical - H_numerical|| / ||H_numerical|| for each parameter
Eigen::VectorXd compareAnalyticalAndNumericalDesignMatrix( const std::string& vehicleCentralBody )
{
    const double initialTime = 1.0E7;
    const double arcDuration = 86400.0;
    const int numberOfArcs = 3;
    const double finalTime = initialTime + numberOfArcs * arcDuration;

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies(
            getDefaultBodySettings( { "Earth", "Moon", "Sun" }, initialTime - 3600.0, finalTime + 3600.0, "Earth", "ECLIPJ2000" ) );
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 1.0E3 );
    bodies.at( "Vehicle" )
            ->setEphemeris( std::make_shared< MultiArcEphemeris >(
                    std::map< double, std::shared_ptr< Ephemeris > >( ), vehicleCentralBody, "ECLIPJ2000" ) );
    bodies.processBodyFrameDefinitions( );

    // Single-arc propagation of the Moon, with an (estimated) empirical acceleration
    SelectedAccelerationMap moonAccelerationSettings;
    moonAccelerationSettings[ "Moon" ][ "Earth" ] = { pointMassGravityAcceleration( ), empiricalAcceleration( ) };
    moonAccelerationSettings[ "Moon" ][ "Sun" ] = { pointMassGravityAcceleration( ) };
    std::shared_ptr< SingleArcPropagatorSettings< double > > moonPropagatorSettings = translationalStatePropagatorSettings< double >(
            { "Earth" },
            createAccelerationModelsMap( bodies, moonAccelerationSettings, { "Moon" }, { "Earth" } ),
            { "Moon" },
            getInitialStatesOfBodies( { "Moon" }, { "Earth" }, bodies, initialTime ),
            initialTime,
            rungeKutta4Settings( 60.0 ),
            propagationTimeTerminationSettings( finalTime ) );

    // Multi-arc propagation of the vehicle
    SelectedAccelerationMap vehicleAccelerationSettings;
    vehicleAccelerationSettings[ "Vehicle" ][ "Earth" ] = { pointMassGravityAcceleration( ) };
    vehicleAccelerationSettings[ "Vehicle" ][ "Moon" ] = { pointMassGravityAcceleration( ) };
    AccelerationMap vehicleAccelerations =
            createAccelerationModelsMap( bodies, vehicleAccelerationSettings, { "Vehicle" }, { vehicleCentralBody } );

    double centralBodyGravitationalParameter = bodies.at( vehicleCentralBody )->getGravityFieldModel( )->getGravitationalParameter( );
    double semiMajorAxis = ( vehicleCentralBody == "Earth" ) ? 1.5E8 : 2.5E6;
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > vehicleArcSettings;
    std::vector< double > arcStartTimes;
    for( int i = 0; i < numberOfArcs; i++ )
    {
        arcStartTimes.push_back( initialTime + i * arcDuration );
        Eigen::Vector6d keplerianState = ( Eigen::Vector6d( ) << semiMajorAxis, 0.1, 0.5, 1.0 + i, 2.0 - i, 0.5 * i ).finished( );
        vehicleArcSettings.push_back( translationalStatePropagatorSettings< double >(
                { vehicleCentralBody },
                vehicleAccelerations,
                { "Vehicle" },
                convertKeplerianToCartesianElements( keplerianState, centralBodyGravitationalParameter ),
                arcStartTimes.back( ),
                rungeKutta4Settings( 60.0 ),
                propagationTimeTerminationSettings( arcStartTimes.back( ) + arcDuration ) ) );
    }
    std::shared_ptr< HybridArcPropagatorSettings< double > > propagatorSettings = std::make_shared< HybridArcPropagatorSettings< double > >(
            moonPropagatorSettings, std::make_shared< MultiArcPropagatorSettings< double > >( vehicleArcSettings ) );

    // Create parameters
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterSettings =
            getInitialStateParameterSettings< double >( propagatorSettings, bodies, arcStartTimes );
    parameterSettings.push_back( gravitationalParameter( "Earth" ) );
    parameterSettings.push_back( constantEmpiricalAccelerationMagnitudes( "Moon", "Earth" ) );
    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double >( parameterSettings, bodies, propagatorSettings );

    // Create orbit determination manager, with vehicle position observations
    LinkEnds linkEnds;
    linkEnds[ observed_body ] = LinkEndId( "Vehicle", "" );
    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies,
            parametersToEstimate,
            { std::make_shared< ObservationModelSettings >( position_observable, linkEnds ) },
            propagatorSettings );

    // Define observation times inside each arc
    std::vector< double > observationTimes;
    for( double time = initialTime + 1800.0; time < finalTime; time += 3600.0 )
    {
        observationTimes.push_back( time );
    }

    // Compute analytical design matrix
    Eigen::VectorXd nominalObservations;
    Eigen::MatrixXd analyticalDesignMatrix;
    orbitDeterminationManager.computePartialsAndObservations< 3 >(
            linkEnds, nullptr, position_observable, observed_body, observationTimes, nominalObservations, analyticalDesignMatrix );

    // Compute numerical design matrix, by central differences
    Eigen::VectorXd nominalParameters = parametersToEstimate->getFullParameterValues< double >( );
    int numberOfInitialStateParameters = parametersToEstimate->getInitialDynamicalStateParameterSize( );
    Eigen::MatrixXd numericalDesignMatrix = Eigen::MatrixXd::Zero( analyticalDesignMatrix.rows( ), analyticalDesignMatrix.cols( ) );
    for( int parameterIndex = 0; parameterIndex < nominalParameters.rows( ); parameterIndex++ )
    {
        double perturbation;
        if( parameterIndex < numberOfInitialStateParameters )
        {
            // Use larger perturbation for (weakly coupled) single-arc Moon state, to limit round-off error in numerical partials
            perturbation = ( ( parameterIndex % 6 < 3 ) ? 10.0 : 1.0E-3 ) * ( ( parameterIndex < 6 ) ? 100.0 : 1.0 );
        }
        else if( parameterIndex == numberOfInitialStateParameters )
        {
            perturbation = 1.0E-7 * nominalParameters( parameterIndex );
        }
        else
        {
            perturbation = 1.0E-7;
        }

        Eigen::VectorXd perturbedObservations[ 2 ];
        Eigen::MatrixXd unusedPartials;
        for( int sign = 0; sign < 2; sign++ )
        {
            Eigen::VectorXd perturbedParameters = nominalParameters;
            perturbedParameters( parameterIndex ) += ( sign == 0 ? 1.0 : -1.0 ) * perturbation;
            orbitDeterminationManager.resetParameterEstimate( perturbedParameters, false );
            orbitDeterminationManager.computePartialsAndObservations< 3 >( linkEnds,
                                                                           nullptr,
                                                                           position_observable,
                                                                           observed_body,
                                                                           observationTimes,
                                                                           perturbedObservations[ sign ],
                                                                           unusedPartials );
        }
        numericalDesignMatrix.col( parameterIndex ) = ( perturbedObservations[ 0 ] - perturbedObservations[ 1 ] ) / ( 2.0 * perturbation );
    }
    orbitDeterminationManager.resetParameterEstimate( nominalParameters, true );

    // Compute relative difference per parameter
    Eigen::VectorXd relativeDifferences = Eigen::VectorXd::Zero( nominalParameters.rows( ) );
    for( int parameterIndex = 0; parameterIndex < nominalParameters.rows( ); parameterIndex++ )
    {
        relativeDifferences( parameterIndex ) =
                ( analyticalDesignMatrix.col( parameterIndex ) - numericalDesignMatrix.col( parameterIndex ) ).norm( ) /
                numericalDesignMatrix.col( parameterIndex ).norm( );
        std::cout << "Vehicle central body " << vehicleCentralBody << ", parameter " << parameterIndex << ": relative difference "
                  << relativeDifferences( parameterIndex ) << " (numerical column norm "
                  << numericalDesignMatrix.col( parameterIndex ).norm( ) << ")" << std::endl;
    }
    return relativeDifferences;
}

BOOST_AUTO_TEST_CASE( testHybridArcParameterCoupling )
{
    spice_interface::loadStandardSpiceKernels( );

    for( std::string vehicleCentralBody : { "Earth", "Moon" } )
    {
        Eigen::VectorXd relativeDifferences = compareAnalyticalAndNumericalDesignMatrix( vehicleCentralBody );
        for( int parameterIndex = 0; parameterIndex < relativeDifferences.rows( ); parameterIndex++ )
        {
            BOOST_CHECK_SMALL( relativeDifferences( parameterIndex ), 5.0E-6 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
