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
//! difference) design matrix, for a hybrid-arc estimation with one or more natural satellites as single-arc bodies. The estimated
//! parameters are the initial states, the gravitational parameters of the central body and the first single-arc body, and a
//! constant empirical acceleration of the last single-arc body. The latter acts only on the single-arc bodies, so that its
//! influence on the vehicle is entirely through their trajectories.
//! \param centralBody Body w.r.t. which the single-arc bodies are propagated
//! \param singleArcBodies Bodies propagated in the single arc (mutually attracting)
//! \param vehicleCentralBody Central body of the vehicle propagation
//! \param semiMajorAxis Semi-major axis of the vehicle orbit
//! \param firstArcOffset Time between start of single-arc propagation and start of first vehicle arc
//! \param arcGap Time between end of one vehicle arc and start of the next
//! \return Relative difference ||H_analytical - H_numerical|| / ||H_numerical|| for each parameter
Eigen::VectorXd compareAnalyticalAndNumericalDesignMatrix( const std::string& centralBody,
                                                           const std::vector< std::string >& singleArcBodies,
                                                           const std::string& vehicleCentralBody,
                                                           const double semiMajorAxis,
                                                           const double firstArcOffset,
                                                           const double arcGap )
{
    const double initialTime = 1.0E7;
    const double arcDuration = 86400.0;
    const int numberOfArcs = 3;
    const double finalTime = initialTime + firstArcOffset + numberOfArcs * arcDuration + ( numberOfArcs - 1 ) * arcGap;

    // Create bodies
    std::vector< std::string > bodiesToCreate = singleArcBodies;
    bodiesToCreate.push_back( centralBody );
    bodiesToCreate.push_back( "Sun" );
    SystemOfBodies bodies = createSystemOfBodies(
            getDefaultBodySettings( bodiesToCreate, initialTime - 3600.0, finalTime + 3600.0, centralBody, "ECLIPJ2000" ) );
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 1.0E3 );
    bodies.at( "Vehicle" )
            ->setEphemeris( std::make_shared< MultiArcEphemeris >(
                    std::map< double, std::shared_ptr< Ephemeris > >( ), vehicleCentralBody, "ECLIPJ2000" ) );
    bodies.processBodyFrameDefinitions( );

    // Single-arc propagation of the natural satellites, with an (estimated) empirical acceleration on the last one
    SelectedAccelerationMap singleArcAccelerationSettings;
    for( const std::string& bodyName : singleArcBodies )
    {
        singleArcAccelerationSettings[ bodyName ][ centralBody ] = { pointMassGravityAcceleration( ) };
        singleArcAccelerationSettings[ bodyName ][ "Sun" ] = { pointMassGravityAcceleration( ) };
        for( const std::string& otherBodyName : singleArcBodies )
        {
            if( otherBodyName != bodyName )
            {
                singleArcAccelerationSettings[ bodyName ][ otherBodyName ] = { pointMassGravityAcceleration( ) };
            }
        }
    }
    singleArcAccelerationSettings[ singleArcBodies.back( ) ][ centralBody ].push_back( empiricalAcceleration( ) );
    std::vector< std::string > singleArcCentralBodies( singleArcBodies.size( ), centralBody );
    std::shared_ptr< SingleArcPropagatorSettings< double > > singleArcPropagatorSettings = translationalStatePropagatorSettings< double >(
            singleArcCentralBodies,
            createAccelerationModelsMap( bodies, singleArcAccelerationSettings, singleArcBodies, singleArcCentralBodies ),
            singleArcBodies,
            getInitialStatesOfBodies( singleArcBodies, singleArcCentralBodies, bodies, initialTime ),
            initialTime,
            rungeKutta4Settings( 60.0 ),
            propagationTimeTerminationSettings( finalTime ) );

    // Multi-arc propagation of the vehicle
    SelectedAccelerationMap vehicleAccelerationSettings;
    vehicleAccelerationSettings[ "Vehicle" ][ centralBody ] = { pointMassGravityAcceleration( ) };
    for( const std::string& bodyName : singleArcBodies )
    {
        vehicleAccelerationSettings[ "Vehicle" ][ bodyName ] = { pointMassGravityAcceleration( ) };
    }
    AccelerationMap vehicleAccelerations =
            createAccelerationModelsMap( bodies, vehicleAccelerationSettings, { "Vehicle" }, { vehicleCentralBody } );

    double centralBodyGravitationalParameter = bodies.at( vehicleCentralBody )->getGravityFieldModel( )->getGravitationalParameter( );
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > vehicleArcSettings;
    std::vector< double > arcStartTimes;
    std::vector< double > observationTimes;
    for( int i = 0; i < numberOfArcs; i++ )
    {
        arcStartTimes.push_back( initialTime + firstArcOffset + i * ( arcDuration + arcGap ) );
        Eigen::Vector6d keplerianState = ( Eigen::Vector6d( ) << semiMajorAxis, 0.1, 0.5, 1.0 + i, 2.0 - i, 0.5 * i ).finished( );
        vehicleArcSettings.push_back( translationalStatePropagatorSettings< double >(
                { vehicleCentralBody },
                vehicleAccelerations,
                { "Vehicle" },
                convertKeplerianToCartesianElements( keplerianState, centralBodyGravitationalParameter ),
                arcStartTimes.back( ),
                rungeKutta4Settings( 60.0 ),
                propagationTimeTerminationSettings( arcStartTimes.back( ) + arcDuration ) ) );

        // Define observation times inside arc
        for( double time = arcStartTimes.back( ) + 1800.0; time < arcStartTimes.back( ) + arcDuration; time += 3600.0 )
        {
            observationTimes.push_back( time );
        }
    }
    std::shared_ptr< HybridArcPropagatorSettings< double > > propagatorSettings = std::make_shared< HybridArcPropagatorSettings< double > >(
            singleArcPropagatorSettings, std::make_shared< MultiArcPropagatorSettings< double > >( vehicleArcSettings ) );

    // Create parameters
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterSettings =
            getInitialStateParameterSettings< double >( propagatorSettings, bodies, arcStartTimes );
    parameterSettings.push_back( gravitationalParameter( centralBody ) );
    parameterSettings.push_back( gravitationalParameter( singleArcBodies.front( ) ) );
    parameterSettings.push_back( constantEmpiricalAccelerationMagnitudes( singleArcBodies.back( ), centralBody ) );
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

    // Compute analytical design matrix
    Eigen::VectorXd nominalObservations;
    Eigen::MatrixXd analyticalDesignMatrix;
    orbitDeterminationManager.computePartialsAndObservations< 3 >(
            linkEnds, nullptr, position_observable, observed_body, observationTimes, nominalObservations, analyticalDesignMatrix );

    // Compute numerical design matrix, by central differences
    Eigen::VectorXd nominalParameters = parametersToEstimate->getFullParameterValues< double >( );
    int numberOfInitialStateParameters = parametersToEstimate->getInitialDynamicalStateParameterSize( );
    int numberOfSingleArcStates = 6 * singleArcBodies.size( );
    Eigen::MatrixXd numericalDesignMatrix = Eigen::MatrixXd::Zero( analyticalDesignMatrix.rows( ), analyticalDesignMatrix.cols( ) );
    for( int parameterIndex = 0; parameterIndex < nominalParameters.rows( ); parameterIndex++ )
    {
        // Use larger perturbations for parameters that are weakly coupled to the vehicle (single-arc states and gravitational
        // parameter, empirical acceleration), to limit round-off error in numerical partials
        double perturbation;
        if( parameterIndex < numberOfInitialStateParameters )
        {
            perturbation = ( ( parameterIndex % 6 < 3 ) ? 10.0 : 1.0E-3 ) * ( ( parameterIndex < numberOfSingleArcStates ) ? 1000.0 : 1.0 );
        }
        else if( parameterIndex == numberOfInitialStateParameters )
        {
            perturbation = 1.0E-7 * nominalParameters( parameterIndex );
        }
        else if( parameterIndex == numberOfInitialStateParameters + 1 )
        {
            perturbation = 1.0E-5 * nominalParameters( parameterIndex );
        }
        else
        {
            perturbation = 1.0E-5;
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

    // Test cases: single-arc Moon, with vehicle orbiting Earth (Moon is perturbing body) or Moon (Moon is central body); single-arc
    // Io and Europa (two mutually attracting single-arc bodies), with vehicle orbiting Jupiter or Europa (second single-arc body is
    // central body), with vehicle arcs that start after the single-arc initial time and have gaps between them.
    for( int testCase = 0; testCase < 4; testCase++ )
    {
        Eigen::VectorXd relativeDifferences;
        switch( testCase )
        {
            case 0:
                relativeDifferences = compareAnalyticalAndNumericalDesignMatrix( "Earth", { "Moon" }, "Earth", 1.5E8, 0.0, 0.0 );
                break;
            case 1:
                relativeDifferences = compareAnalyticalAndNumericalDesignMatrix( "Earth", { "Moon" }, "Moon", 2.5E6, 0.0, 0.0 );
                break;
            case 2:
                relativeDifferences =
                        compareAnalyticalAndNumericalDesignMatrix( "Jupiter", { "Io", "Europa" }, "Jupiter", 1.0E9, 43200.0, 21600.0 );
                break;
            default:
                relativeDifferences =
                        compareAnalyticalAndNumericalDesignMatrix( "Jupiter", { "Io", "Europa" }, "Europa", 2.5E6, 43200.0, 21600.0 );
                break;
        }
        for( int parameterIndex = 0; parameterIndex < relativeDifferences.rows( ); parameterIndex++ )
        {
            BOOST_CHECK_SMALL( relativeDifferences( parameterIndex ), 2.0E-5 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
