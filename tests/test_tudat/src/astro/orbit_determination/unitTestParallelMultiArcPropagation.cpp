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

#include <string>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/simulatePseudoObservations.h"
#include "tudat/simulation/propagation_setup/multiArcDynamicsSimulator.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::ephemerides;
using namespace tudat::estimatable_parameters;
using namespace tudat::numerical_integrators;
using namespace tudat::observation_models;
using namespace tudat::propagators;
using namespace tudat::simulation_setup;

BOOST_AUTO_TEST_SUITE( test_parallel_multi_arc_propagation )

// Estimate arc-wise atmosphere parameters as in unitTestEstimationArcWiseAtmosphereParameters.cpp, using four arcs, and with
// environment models that require no SPICE calls during propagation. For numberOfThreads equal to 1, the acceleration models are
// created on the bodies, otherwise from the acceleration settings for each thread. Returns the estimated parameters and residuals.
std::pair< Eigen::VectorXd, Eigen::VectorXd > estimateArcWiseAtmosphereParameters( const unsigned int numberOfThreads )
{
    // Cassini T022 and T068 states (c/a-5min) from unitTestEstimationArcWiseAtmosphereParameters.cpp; the two additional arcs use
    // the same Titan-centered states one day later
    std::vector< double > arcStartTimes = {
        2.205723865899772644e+08, 2.205723865899772644e+08 + 86400.0, 3.275976261739959717e+08, 3.275976261739959717e+08 + 86400.0
    };
    std::vector< Eigen::Vector6d > arcInitialStates;
    arcInitialStates.push_back( Eigen::Vector6d{ 2.166092477145108860e+06,
                                                 -7.925218812290439382e+05,
                                                 3.553645666049534921e+06,
                                                 3.304561684998871897e+03,
                                                 2.420278071381293103e+03,
                                                 -4.244915762044645817e+03 } );
    arcInitialStates.push_back( arcInitialStates.at( 0 ) );
    arcInitialStates.push_back( Eigen::Vector6d{ -3.091501373934821691e+06,
                                                 9.134845736499809427e+05,
                                                 -2.880028589145255741e+06,
                                                 4.617009140001175183e+03,
                                                 3.555047142774454642e+03,
                                                 -4.711669772164187862e+02 } );
    arcInitialStates.push_back( arcInitialStates.at( 2 ) );

    double initialEnvironmentTime = arcStartTimes.front( ) - 10.0 * 3600.0;
    double finalEnvironmentTime = arcStartTimes.back( ) + 10.0 * 3600.0;
    BodyListSettings bodySettings =
            getDefaultBodySettings( { "Titan", "Saturn" }, initialEnvironmentTime, finalEnvironmentTime, "Titan", "ECLIPJ2000", 3600.0 );
    for( const std::string& bodyName : std::vector< std::string >{ "Titan", "Saturn" } )
    {
        // SPICE cannot be used during multi-threaded propagation
        bodySettings.at( bodyName )->rotationModelSettings =
                simpleRotationModelFromSpiceSettings( "ECLIPJ2000", "IAU_" + bodyName, "IAU_" + bodyName, initialEnvironmentTime );
    }
    bodySettings.at( "Titan" )->atmosphereSettings = std::make_shared< ExponentialAtmosphereSettings >( 7.58e04, 175, 3.0553447e-04 );
    bodySettings.addSettings( "Cassini" );
    bodySettings.at( "Cassini" )->ephemerisSettings =
            std::make_shared< TabulatedEphemerisSettings >( std::map< double, Eigen::Vector6d >( ), "Titan", "ECLIPJ2000" );
    bodySettings.at( "Cassini" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Cassini" )->constantMass = 2.0E3;
    bodySettings.at( "Cassini" )->aerodynamicCoefficientSettings = std::make_shared< ConstantAerodynamicCoefficientSettings >(
            20.0, Eigen::Vector3d( 2.5, -0.1, 0.05 ), aerodynamics::negative_aerodynamic_frame_coefficients );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Cassini" ][ "Titan" ].push_back( std::make_shared< AccelerationSettings >( aerodynamic ) );
    accelerationMap[ "Cassini" ][ "Titan" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Cassini" ][ "Saturn" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    std::vector< std::string > bodiesToEstimate = { "Cassini" };
    std::vector< std::string > centralBodies = { "Titan" };

    AccelerationMap accelerationModelMap;
    if( numberOfThreads == 1 )
    {
        accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToEstimate, centralBodies );
    }

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesList = { singleAccelerationDependentVariable(
            aerodynamic, "Cassini", "Titan" ) };

    std::vector< double > arcEndTimes;
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
    for( unsigned int j = 0; j < arcStartTimes.size( ); j++ )
    {
        arcEndTimes.push_back( arcStartTimes.at( j ) + 10 * 60 );
        if( numberOfThreads == 1 )
        {
            propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies,
                    accelerationModelMap,
                    bodiesToEstimate,
                    Eigen::VectorXd( arcInitialStates.at( j ) ),
                    arcStartTimes.at( j ),
                    rungeKuttaFixedStepSettings( 1.0, rungeKuttaFehlberg78 ),
                    propagationTimeTerminationSettings( arcEndTimes.at( j ) ),
                    cowell,
                    dependentVariablesList ) );
        }
        else
        {
            propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                    centralBodies,
                    accelerationMap,
                    bodiesToEstimate,
                    Eigen::VectorXd( arcInitialStates.at( j ) ),
                    arcStartTimes.at( j ),
                    rungeKuttaFixedStepSettings( 1.0, rungeKuttaFehlberg78 ),
                    propagationTimeTerminationSettings( arcEndTimes.at( j ) ),
                    cowell,
                    dependentVariablesList ) );
        }
    }
    std::shared_ptr< MultiArcPropagatorSettings< double > > propagatorSettings = std::make_shared< MultiArcPropagatorSettings< double > >(
            propagatorSettingsList, false, std::make_shared< MultiArcPropagatorProcessingSettings >( false, true ) );
    propagatorSettings->setNumberOfThreads( numberOfThreads );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialMultiArcParameterSettings< double >( propagatorSettings, bodies, arcStartTimes );
    parameterNames.push_back( arcwiseExponentialAtmosphereBaseDensity( "Titan", arcStartTimes ) );
    parameterNames.push_back( arcwiseExponentialAtmosphereScaleHeight( "Titan", arcStartTimes ) );

    // Propagate the true dynamics, and simulate observations of the resulting Cassini states
    MultiArcDynamicsSimulator<> dynamicsSimulator( bodies, propagatorSettings, true );
    std::vector< double > observationTimes;
    for( unsigned int k = 0; k < arcStartTimes.size( ); k++ )
    {
        double currentTime = arcStartTimes.at( k ) + 5.0;
        while( currentTime < arcEndTimes.at( k ) - 5.0 )
        {
            observationTimes.push_back( currentTime );
            currentTime += 1.0;
        }
    }
    std::pair< std::vector< std::shared_ptr< ObservationModelSettings > >, std::shared_ptr< ObservationCollection< double > > >
            observationCollectionAndModelSettings = simulatePseudoObservations( bodies, bodiesToEstimate, centralBodies, observationTimes );

    // Create the estimated parameters (from their settings), perturb the atmosphere parameters, and estimate
    OrbitDeterminationManager< double > orbitDeterminationManager(
            bodies, parameterNames, observationCollectionAndModelSettings.first, propagatorSettings );
    Eigen::VectorXd perturbation( 8 );
    perturbation << -3e-05, 6e-05, -3e-05, 6e-05, 1e03, -1.5e03, 1e03, -1.5e03;
    Eigen::VectorXd parameterVector = orbitDeterminationManager.getParametersToEstimate( )->getFullParameterValues< double >( );
    parameterVector.tail( perturbation.size( ) ) += perturbation;
    orbitDeterminationManager.resetParameterEstimate( parameterVector );
    std::shared_ptr< EstimationInput< double > > estimationInput =
            std::make_shared< EstimationInput< double > >( observationCollectionAndModelSettings.second );
    estimationInput->defineEstimationSettings( 0, 1, 0, 1, 1, 1 );
    std::shared_ptr< EstimationOutput<> > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );

    return std::make_pair( estimationOutput->parameterEstimate_, estimationOutput->residuals_ );
}

//! Test whether estimation with multi-threaded multi-arc propagation gives results equal (to numerical noise) to single-threaded propagation
BOOST_AUTO_TEST_CASE( testParallelMultiArcEstimation )
{
    spice_interface::loadStandardSpiceKernels( );

    std::pair< Eigen::VectorXd, Eigen::VectorXd > singleThreadResults = estimateArcWiseAtmosphereParameters( 1 );
    BOOST_CHECK( !singleThreadResults.first.hasNaN( ) );

    // Use fewer threads than arcs (several arcs per thread), and as many threads as arcs
    for( unsigned int numberOfThreads : { 2, 3, 4 } )
    {
        std::pair< Eigen::VectorXd, Eigen::VectorXd > multiThreadResults = estimateArcWiseAtmosphereParameters( numberOfThreads );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( multiThreadResults.first, singleThreadResults.first, 1.0E-8 );
        BOOST_CHECK_SMALL( ( multiThreadResults.second - singleThreadResults.second ).cwiseAbs( ).maxCoeff( ), 1.0E-6 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
