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

//! Function to compute observations and partials for a list of position observables
//! \param orbitDeterminationManager Object used to compute observations and partials
//! \param observationLinks List of observable types and link ends
//! \param observationTimes Observation times (identical for each observable)
//! \param observations Observations for each entry of observationLinks (returned by reference)
//! \param partials Partials of observations w.r.t. estimated parameters for each entry of observationLinks (returned by reference)
void computePositionObservationsAndPartials( OrbitDeterminationManager< double, double >& orbitDeterminationManager,
                                             const std::vector< std::pair< ObservableType, LinkEnds > >& observationLinks,
                                             const std::vector< double >& observationTimes,
                                             std::vector< Eigen::VectorXd >& observations,
                                             std::vector< Eigen::MatrixXd >& partials )
{
    observations.resize( observationLinks.size( ) );
    partials.resize( observationLinks.size( ) );
    for( unsigned int i = 0; i < observationLinks.size( ); i++ )
    {
        orbitDeterminationManager.computePartialsAndObservations< 3 >( observationLinks.at( i ).second,
                                                                       nullptr,
                                                                       observationLinks.at( i ).first,
                                                                       observed_body,
                                                                       observationTimes,
                                                                       observations.at( i ),
                                                                       partials.at( i ) );
    }
}

//! Function to compare the analytical design matrix of position observations with a numerical (central difference) design
//! matrix, for a hybrid-arc estimation of a multi-arc vehicle and one or more natural satellites as single-arc bodies, with the
//! single-arc central body or the SSB as global frame origin. The observations are the absolute positions of the vehicle and the single-arc bodies, and the positions
//! of the vehicle w.r.t. its central body and w.r.t. a third body. The estimated parameters are the initial states, the
//! gravitational parameters of the central body and the first single-arc body, and a constant empirical acceleration of the last
//! single-arc body. The latter acts only on the single-arc bodies, so that its influence on the vehicle is entirely through their
//! trajectories. The relative difference ||H_analytical - H_numerical|| / ||H_numerical|| is checked for each observable and
//! parameter (||H_analytical|| where H_numerical is zero).
//! \param centralBody Body w.r.t. which the single-arc bodies are propagated
//! \param singleArcBodies Bodies propagated in the single arc (mutually attracting)
//! \param vehicleCentralBody Central body of the vehicle propagation
//! \param vehicleThirdBody Body (other than vehicleCentralBody) w.r.t. which the vehicle position is also observed
//! \param semiMajorAxis Semi-major axis of the vehicle orbit
//! \param firstArcOffset Time between start of single-arc propagation and start of first vehicle arc
//! \param arcGap Time between end of one vehicle arc and start of the next
//! \param useBarycentricOrigin Boolean denoting whether the SSB (if true) or centralBody (if false) is the global frame origin
void compareAnalyticalAndNumericalDesignMatrix( const std::string& centralBody,
                                                const std::vector< std::string >& singleArcBodies,
                                                const std::string& vehicleCentralBody,
                                                const std::string& vehicleThirdBody,
                                                const double semiMajorAxis,
                                                const double firstArcOffset,
                                                const double arcGap,
                                                const bool useBarycentricOrigin )
{
    const double initialTime = 1.0E7;
    const double arcDuration = 86400.0;
    const int numberOfArcs = 3;
    const double finalTime = initialTime + firstArcOffset + numberOfArcs * arcDuration + ( numberOfArcs - 1 ) * arcGap;

    // Create bodies, with ephemerides of single-arc bodies w.r.t. their propagation central body
    std::vector< std::string > bodiesToCreate = singleArcBodies;
    bodiesToCreate.push_back( centralBody );
    bodiesToCreate.push_back( "Sun" );
    BodyListSettings bodySettings = getDefaultBodySettings(
            bodiesToCreate, initialTime - 3600.0, finalTime + 3600.0, useBarycentricOrigin ? "SSB" : centralBody, "ECLIPJ2000" );
    for( const std::string& bodyName : singleArcBodies )
    {
        bodySettings.at( bodyName )->ephemerisSettings->resetFrameOrigin( centralBody );
    }
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
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

    // Define observations: absolute positions of vehicle and single-arc bodies, and vehicle position w.r.t. its central body and
    // w.r.t. a third body
    std::vector< std::pair< ObservableType, LinkEnds > > observationLinks;
    LinkEnds vehicleLinkEnds;
    vehicleLinkEnds[ observed_body ] = LinkEndId( "Vehicle", "" );
    observationLinks.push_back( std::make_pair( position_observable, vehicleLinkEnds ) );
    for( const std::string& observerName : { vehicleCentralBody, vehicleThirdBody } )
    {
        vehicleLinkEnds[ observer ] = LinkEndId( observerName, "" );
        observationLinks.push_back( std::make_pair( relative_position_observable, vehicleLinkEnds ) );
    }
    for( const std::string& bodyName : singleArcBodies )
    {
        LinkEnds singleArcBodyLinkEnds;
        singleArcBodyLinkEnds[ observed_body ] = LinkEndId( bodyName, "" );
        observationLinks.push_back( std::make_pair( position_observable, singleArcBodyLinkEnds ) );
    }

    // Create orbit determination manager
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettings;
    for( unsigned int i = 0; i < observationLinks.size( ); i++ )
    {
        observationSettings.push_back(
                std::make_shared< ObservationModelSettings >( observationLinks.at( i ).first, observationLinks.at( i ).second ) );
    }
    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, observationSettings, propagatorSettings );

    // Compute analytical design matrices
    std::vector< Eigen::VectorXd > nominalObservations;
    std::vector< Eigen::MatrixXd > analyticalDesignMatrices;
    computePositionObservationsAndPartials(
            orbitDeterminationManager, observationLinks, observationTimes, nominalObservations, analyticalDesignMatrices );

    // Compute numerical design matrices by central differences with nominal and doubled perturbation, which are combined into a
    // fourth-order central difference. Their difference is used as an estimate of the error of the numerical partials, which is
    // large for weakly coupled parameters (due to round-off error)
    Eigen::VectorXd nominalParameters = parametersToEstimate->getFullParameterValues< double >( );
    int numberOfInitialStateParameters = parametersToEstimate->getInitialDynamicalStateParameterSize( );
    int numberOfSingleArcStates = 6 * singleArcBodies.size( );
    std::vector< Eigen::MatrixXd > numericalDesignMatrices[ 2 ];
    for( int stepCase = 0; stepCase < 2; stepCase++ )
    {
        for( unsigned int i = 0; i < observationLinks.size( ); i++ )
        {
            numericalDesignMatrices[ stepCase ].push_back(
                    Eigen::MatrixXd::Zero( analyticalDesignMatrices.at( i ).rows( ), analyticalDesignMatrices.at( i ).cols( ) ) );
        }
    }
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

        for( int stepCase = 0; stepCase < 2; stepCase++ )
        {
            double currentPerturbation = ( stepCase + 1 ) * perturbation;
            std::vector< Eigen::VectorXd > perturbedObservations[ 2 ];
            std::vector< Eigen::MatrixXd > unusedPartials;
            for( int sign = 0; sign < 2; sign++ )
            {
                Eigen::VectorXd perturbedParameters = nominalParameters;
                perturbedParameters( parameterIndex ) += ( sign == 0 ? 1.0 : -1.0 ) * currentPerturbation;
                orbitDeterminationManager.resetParameterEstimate( perturbedParameters, false );
                computePositionObservationsAndPartials(
                        orbitDeterminationManager, observationLinks, observationTimes, perturbedObservations[ sign ], unusedPartials );
            }
            for( unsigned int i = 0; i < observationLinks.size( ); i++ )
            {
                numericalDesignMatrices[ stepCase ].at( i ).col( parameterIndex ) =
                        ( perturbedObservations[ 0 ].at( i ) - perturbedObservations[ 1 ].at( i ) ) / ( 2.0 * currentPerturbation );
            }
        }
    }
    orbitDeterminationManager.resetParameterEstimate( nominalParameters, true );

    // Compare analytical and numerical partials, per observable and parameter, allowing for the estimated error of the numerical
    // partials (with a margin of 3, since the fourth-order difference is noisier than the difference of the two central
    // differences). Where the numerical partials are zero (single-arc body positions w.r.t. vehicle initial states), the analytical
    // partials must be zero as well
    for( unsigned int i = 0; i < observationLinks.size( ); i++ )
    {
        for( int parameterIndex = 0; parameterIndex < nominalParameters.rows( ); parameterIndex++ )
        {
            Eigen::VectorXd numericalPartial = ( 4.0 * numericalDesignMatrices[ 0 ].at( i ).col( parameterIndex ) -
                                                 numericalDesignMatrices[ 1 ].at( i ).col( parameterIndex ) ) /
                    3.0;
            double numericalError = ( numericalDesignMatrices[ 0 ].at( i ).col( parameterIndex ) -
                                      numericalDesignMatrices[ 1 ].at( i ).col( parameterIndex ) )
                                            .norm( );
            double difference = ( analyticalDesignMatrices.at( i ).col( parameterIndex ) - numericalPartial ).norm( );
            std::cout << ( useBarycentricOrigin ? "SSB" : centralBody ) << " origin, vehicle central body " << vehicleCentralBody << ", "
                      << getObservableName( observationLinks.at( i ).first, observationLinks.at( i ).second.size( ) ) << " "
                      << getLinkEndsString( observationLinks.at( i ).second ) << ", parameter " << parameterIndex
                      << ": relative difference " << difference / numericalPartial.norm( ) << ", relative numerical error "
                      << numericalError / numericalPartial.norm( ) << std::endl;
            BOOST_CHECK_LE( difference, 2.0E-5 * numericalPartial.norm( ) + 3.0 * numericalError );
        }
    }
}

BOOST_AUTO_TEST_CASE( testHybridArcParameterCoupling )
{
    spice_interface::loadStandardSpiceKernels( );

    // Test cases: single-arc Moon, with vehicle orbiting Earth (Moon is perturbing body) or Moon (Moon is central body); single-arc
    // Io and Europa (two mutually attracting single-arc bodies), with vehicle orbiting Jupiter or Europa (second single-arc body is
    // central body), with vehicle arcs that start after the single-arc initial time and have gaps between them. The vehicle position
    // is also observed w.r.t. a third body, which is a single-arc body (Moon, Io) or a body that is not propagated (Earth).
    // Each test case is run with the single-arc central body (Earth, Jupiter) and with the SSB as global frame origin.
    for( bool useBarycentricOrigin : { false, true } )
    {
        compareAnalyticalAndNumericalDesignMatrix( "Earth", { "Moon" }, "Earth", "Moon", 1.5E8, 0.0, 0.0, useBarycentricOrigin );
        compareAnalyticalAndNumericalDesignMatrix( "Earth", { "Moon" }, "Moon", "Earth", 2.5E6, 0.0, 0.0, useBarycentricOrigin );
        compareAnalyticalAndNumericalDesignMatrix(
                "Jupiter", { "Io", "Europa" }, "Jupiter", "Io", 1.0E9, 43200.0, 21600.0, useBarycentricOrigin );
        compareAnalyticalAndNumericalDesignMatrix(
                "Jupiter", { "Io", "Europa" }, "Europa", "Io", 2.5E6, 43200.0, 21600.0, useBarycentricOrigin );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
