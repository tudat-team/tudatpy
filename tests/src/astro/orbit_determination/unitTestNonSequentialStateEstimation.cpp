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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <string>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/math/integrators/rungeKuttaCoefficients.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/simulation/simulation.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/simulation/estimation.h"

namespace tudat {

namespace unit_tests {

//! Using declarations.
using namespace interpolators;
using namespace numerical_integrators;
using namespace spice_interface;
using namespace ephemerides;
using namespace simulation_setup;
using namespace basic_astrodynamics;
using namespace orbital_element_conversions;
using namespace propagators;
using namespace tudat::estimatable_parameters;

BOOST_AUTO_TEST_SUITE( test_non_sequential_state_estimation )

BOOST_AUTO_TEST_CASE( testNonSequentialSingleArcStateEstimation )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Ganymede" );
    bodyNames.push_back( "Callisto" );

    // Specify initial time
    double initialEpoch = 0.0;
    double finalEpoch = 86400.0;
    double midArcEpoch = ( initialEpoch + finalEpoch ) / 2.0;

    std::string globalFrameOrientation = "ECLIPJ2000";
    std::string globalFrameOrigin = "Jupiter";

    // Create bodies.
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEpoch - 86400.0, finalEpoch + 86400.0 );
    bodySettings.at( "Io" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch, finalEpoch, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Europa" )->ephemerisSettings = std::make_shared<InterpolatedSpiceEphemerisSettings >(
            initialEpoch, finalEpoch, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Ganymede" )->ephemerisSettings = std::make_shared<InterpolatedSpiceEphemerisSettings >(
            initialEpoch, finalEpoch, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Callisto" )->ephemerisSettings = std::make_shared<InterpolatedSpiceEphemerisSettings >(
            initialEpoch, finalEpoch, 3600.0, globalFrameOrigin, globalFrameOrientation );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    //! Define bodies to propagate.
    std::vector< std::string > bodiesToPropagate = { "Io", "Europa", "Ganymede", "Callisto" };
    std::vector< std::string > centralBodies = { "Jupiter", "Jupiter", "Jupiter", "Jupiter" };

    // Set accelerations.
    SelectedAccelerationMap accelerationSettings;
    for ( unsigned int i = 0 ; i < bodiesToPropagate.size( ) ; i++ )
    {
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
        accelerationsOfSatellite[ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 8, 0, 2, 2 ) );
        for ( unsigned int j = 0 ; j < bodiesToPropagate.size( ) ; j++ )
        {
            if ( i != j )
            {
                accelerationsOfSatellite[ bodiesToPropagate[ j ] ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2, 8, 0 ) );
            }
        }
        accelerationSettings[ bodiesToPropagate[ i ] ] = accelerationsOfSatellite;
    }

    basic_astrodynamics::AccelerationMap accelerationsMap = createAccelerationModelsMap(
            bodies, accelerationSettings, bodiesToPropagate, centralBodies );

    // Define integrator settings
    double timeStep = 3600.0;
    std::shared_ptr< numerical_integrators::IntegratorSettings< > > forwardIntegratorSettings =
            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    midArcEpoch, timeStep, CoefficientSets::rungeKuttaFehlberg78, timeStep, timeStep, 1.0e3, 1.0e3 );

    std::shared_ptr< numerical_integrators::IntegratorSettings< > > backwardIntegratorSettings =
            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    midArcEpoch, - timeStep, CoefficientSets::rungeKuttaFehlberg78, - timeStep, - timeStep, 1.0e3, 1.0e3 );

    // Define initial states
    Eigen::VectorXd midArcStatesMoons = propagators::getInitialStatesOfBodies( bodiesToPropagate, centralBodies, bodies, midArcEpoch );

    // Define dependent variables to save
    std::vector<std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
            basic_astrodynamics::mutual_spherical_harmonic_gravity, "Io", "Jupiter" ) );
    for ( unsigned int k = 0 ; k < bodiesToPropagate.size( ) ; k++ )
    {
        if ( bodiesToPropagate.at( k ) != "Io" )
        {
            dependentVariables.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    basic_astrodynamics::mutual_spherical_harmonic_gravity, "Io", bodiesToPropagate.at( k ) ) );
        }
    }


    // Define propagator settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< > > forwardPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
            centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons, finalEpoch, cowell, dependentVariables );

    std::shared_ptr< TranslationalStatePropagatorSettings< > > backwardPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
            centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons, initialEpoch, cowell, dependentVariables );

    //! Create settings for non-sequential propagation
    std::shared_ptr< NonSequentialPropagationTerminationSettings > terminationSettings = std::make_shared< NonSequentialPropagationTerminationSettings >(
            std::make_shared< PropagationTimeTerminationSettings >( finalEpoch ), std::make_shared< PropagationTimeTerminationSettings >( initialEpoch ) );
    std::shared_ptr< TranslationalStatePropagatorSettings< > > nonsequentialPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
            centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons, terminationSettings, cowell, dependentVariables );

    // Define parameters to estimate for non-sequentiql propagation / estimation
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    for ( unsigned int i = 0 ; i < bodiesToPropagate.size( ) ; i++ )
    {
        parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
                bodiesToPropagate.at( i ), midArcStatesMoons.segment( i * 6, 6 ), centralBodies.at( i ) ) );

        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( bodiesToPropagate.at( i ), gravitational_parameter) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 0, 2, 2, bodiesToPropagate.at( i ), spherical_harmonics_cosine_coefficient_block ) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 1, 2, 2, bodiesToPropagate.at( i ), spherical_harmonics_sine_coefficient_block ) );
    }
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > nonSequentialPropagationParametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies, nonsequentialPropagatorSettings );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > nonSequentialParameterEstimate =
            nonSequentialPropagationParametersToEstimate->template getFullParameterValues< double >( );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > forwardPropagationParametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies, forwardPropagatorSettings );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > backwardPropagationParametersToEstimate =
            createParametersToEstimate< double >( parameterNames, bodies, backwardPropagatorSettings );


    // Define links and observations.
    std::vector< observation_models::LinkEnds > linkEndsList;
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationSettingsList;
    linkEndsList.resize( bodiesToPropagate.size( ) );
    for ( int i = 0 ; i < bodiesToPropagate.size( ) ; i++ )
    {
        linkEndsList[ i ][ observation_models::observed_body ] = observation_models::LinkEndId( bodiesToPropagate.at( i ), "" );
        observationSettingsList.push_back( std::make_shared< observation_models::ObservationModelSettings >(
                observation_models::position_observable, linkEndsList[ i ] ) );
    }

    // Define observation times
    std::vector< double > observationTimesForward;
    for ( double time = midArcEpoch + 3600.0 ; time < finalEpoch - 3600.0 ; time += 3.0 * 3600.0 )
    {
        observationTimesForward.push_back( time );
    }

    std::vector< double > observationTimesBackward;
    for ( double time = initialEpoch + 3600.0 ; time < midArcEpoch - 3600.0 ; time += 3.0 * 3600.0 )
    {
        observationTimesBackward.push_back( time );
    }

    std::vector< double > allObservationTimes = observationTimesBackward;
    for ( unsigned int i = 0 ; i < observationTimesForward.size( ) ; i++ )
    {
        allObservationTimes.push_back( observationTimesForward.at( i ) );
    }

    // Define observation settings
    std::vector< std::shared_ptr<ObservationSimulationSettings< double > > > measurementSimulationInputForward,
            measurementSimulationInputBackward, measurementSimulationInputAll;
    for ( unsigned int i = 0 ; i < bodiesToPropagate.size( ) ; i++ )
    {
        measurementSimulationInputForward.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
                observation_models::position_observable, linkEndsList[ i ], observationTimesForward, observation_models::observed_body ) );
        measurementSimulationInputBackward.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
                observation_models::position_observable, linkEndsList[ i ], observationTimesBackward, observation_models::observed_body ) );
        measurementSimulationInputAll.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
                observation_models::position_observable, linkEndsList[ i ], allObservationTimes, observation_models::observed_body ) );
    }


    // Create orbit determination object for forward propagation / estimation.
    OrbitDeterminationManager< > orbitDeterminationManagerForward = OrbitDeterminationManager< >(
            bodies, forwardPropagationParametersToEstimate, observationSettingsList, forwardIntegratorSettings, forwardPropagatorSettings );

    // Simulate observations for forward propagation / estimation
    std::shared_ptr< observation_models::ObservationCollection< > > observationsAndTimesForward = simulateObservations< >(
            measurementSimulationInputForward, orbitDeterminationManagerForward.getObservationSimulators( ), bodies );

    // Define estimation input for forward propagation / estimation
    std::shared_ptr< EstimationInput< double, double  > > estimationInputForward =
            std::make_shared< EstimationInput< double, double > >( observationsAndTimesForward );

    // Perform forward estimation
    std::shared_ptr< EstimationOutput< double, double > > estimationOutputForward = orbitDeterminationManagerForward.estimateParameters( estimationInputForward );


    // Create orbit determination object for backward propagation / estimation.
    OrbitDeterminationManager< > orbitDeterminationManagerBackward = OrbitDeterminationManager< >(
            bodies, backwardPropagationParametersToEstimate, observationSettingsList, backwardIntegratorSettings, backwardPropagatorSettings );

    // Simulate observations for backward propagation / estimation
    std::shared_ptr< observation_models::ObservationCollection< > > observationsAndTimesBackward = simulateObservations< >(
            measurementSimulationInputBackward, orbitDeterminationManagerBackward.getObservationSimulators( ), bodies );

    // Define POD input for backward propagation / estimation
    std::shared_ptr< EstimationInput< double, double  > > estimationInputBackward =
            std::make_shared< EstimationInput< double, double > >( observationsAndTimesBackward );

    // Perform backward estimation
    std::shared_ptr< EstimationOutput< double, double > > estimationOutputBackward = orbitDeterminationManagerBackward.estimateParameters( estimationInputBackward );


    // Create orbit determination object for non-sequential propagation / estimation.
    OrbitDeterminationManager< > orbitDeterminationManagerNonSequential = OrbitDeterminationManager< >(
            bodies, nonSequentialPropagationParametersToEstimate, observationSettingsList, forwardIntegratorSettings, nonsequentialPropagatorSettings );

    // Simulate observations for non-sequential propagation / estimation
    std::shared_ptr< observation_models::ObservationCollection< > > observationsAndTimesNonSequential = simulateObservations< >(
            measurementSimulationInputAll, orbitDeterminationManagerNonSequential.getObservationSimulators( ), bodies );

    // Define POD input for non-sequential propgation / estimation
    std::shared_ptr< EstimationInput< double, double  > > estimationInputNonSequential =
            std::make_shared< EstimationInput< double, double > >( observationsAndTimesNonSequential );

    // Perform non-sequential estimation
    std::shared_ptr< EstimationOutput< double, double > > estimationOutputNonSequential = orbitDeterminationManagerNonSequential.estimateParameters( estimationInputNonSequential );


    // Retrieve partials from each estimation.
    Eigen::MatrixXd partialsForwardEstimation = estimationOutputForward->getUnnormalizedDesignMatrix( );
    Eigen::MatrixXd partialsBackwardEstimation = estimationOutputBackward->getUnnormalizedDesignMatrix( );
    Eigen::MatrixXd partialsNonSequentialEstimation = estimationOutputNonSequential->getUnnormalizedDesignMatrix( );

    Eigen::MatrixXd combinedPartials = Eigen::MatrixXd::Zero( partialsForwardEstimation.rows( ) + partialsBackwardEstimation.rows( ), partialsForwardEstimation.cols( ) );

    // Combine forward and backward partials in proper order.
    unsigned int nbParameter = nonSequentialParameterEstimate.size( );
    for ( unsigned int j = 0 ; j < bodiesToPropagate.size( ) ; j++ )
    {
        std::vector< std::pair< int, int > > backwardPartialsIndices =
                estimationInputBackward->getObservationCollection( )->getObservationSetStartAndSize( ).at( observation_models::position_observable ).at( linkEndsList[ j ] );
        std::vector< std::pair< int, int > > forwardPartialsIndices =
                estimationInputForward->getObservationCollection( )->getObservationSetStartAndSize( ).at( observation_models::position_observable ).at( linkEndsList[ j ] );
        std::vector< std::pair< int, int > > nonSequentialPartialsIndices =
                estimationInputNonSequential->getObservationCollection( )->getObservationSetStartAndSize( ).at( observation_models::position_observable ).at(linkEndsList[ j ] );

        //! Add partials from backward estimation
        combinedPartials.block( nonSequentialPartialsIndices.at( 0 ).first, 0, backwardPartialsIndices.at( 0 ).second, nbParameter ) =
                partialsBackwardEstimation.block( backwardPartialsIndices.at( 0 ).first, 0, backwardPartialsIndices.at( 0 ).second, nbParameter );

        //! Add partials from forward estimation
        combinedPartials.block( nonSequentialPartialsIndices.at( 0 ).first + backwardPartialsIndices.at( 0 ).second, 0,
                                forwardPartialsIndices.at( 0 ).second, nbParameter ) =
                partialsForwardEstimation.block( forwardPartialsIndices.at( 0 ).first, 0, forwardPartialsIndices.at( 0 ).second, nbParameter );
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( combinedPartials, partialsNonSequentialEstimation, 1.0e-12 );

}

BOOST_AUTO_TEST_CASE( testNonSequentialMultiArcStateEstimation )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Ganymede" );
    bodyNames.push_back( "Callisto" );

    // Specify arc-wise start and end times
    double initialEpoch = 0.0;
    double finalEpoch = 86400.0;
    std::vector< double > arcStartTimes = { 0.0, 12.0 * 3600.0, 24.0 * 3600.0 };
    std::vector< double > arcEndTimes = { 12.0 * 3600.0, 24.0 * 3600.0, 36.0 * 3600.0 };
    std::vector< double > midArcTimes;
    for ( unsigned int i = 0 ; i < arcStartTimes.size( ) ; i++ )
    {
        midArcTimes.push_back( ( arcStartTimes.at( i ) + arcEndTimes.at( i ) ) / 2.0 );
    }

    unsigned int nbArcs = arcStartTimes.size( );

    std::string globalFrameOrientation = "ECLIPJ2000";
    std::string globalFrameOrigin = "Jupiter";

    // Create bodies.
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEpoch - 86400.0, finalEpoch + 86400.0 );
    bodySettings.at( "Io" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch, finalEpoch, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Europa" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch, finalEpoch, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Ganymede" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch, finalEpoch, 3600.0, globalFrameOrigin, globalFrameOrientation );
    bodySettings.at( "Callisto" )->ephemerisSettings = std::make_shared< InterpolatedSpiceEphemerisSettings >(
            initialEpoch, finalEpoch, 3600.0, globalFrameOrigin, globalFrameOrientation );

    bodySettings.at( "Io" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Europa" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Ganymede" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Callisto" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    //! Define bodies to propagate.
    std::vector< std::string > bodiesToPropagate = { "Io", "Europa", "Ganymede", "Callisto" };
    std::vector< std::string > centralBodies = { "Jupiter", "Jupiter", "Jupiter", "Jupiter" };

    // Set accelerations.
    SelectedAccelerationMap accelerationSettings;
    for ( unsigned int i = 0 ; i < bodiesToPropagate.size( ) ; i++ )
    {
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfSatellite;
        accelerationsOfSatellite[ "Jupiter" ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 8, 0, 2, 2  )  );
        for ( unsigned int j = 0 ; j < bodiesToPropagate.size( ) ; j++ )
        {
            if ( i != j )
            {
                accelerationsOfSatellite[ bodiesToPropagate[ j ] ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2, 8, 0 ) );
            }
        }
        accelerationSettings[ bodiesToPropagate[ i ] ] = accelerationsOfSatellite;
    }

    basic_astrodynamics::AccelerationMap accelerationsMap = createAccelerationModelsMap(
            bodies, accelerationSettings, bodiesToPropagate, centralBodies );

    // Define integrator settings
    double timeStep = 60.0;
    std::shared_ptr< numerical_integrators::IntegratorSettings< > > forwardIntegratorSettings =
            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    TUDAT_NAN, timeStep, CoefficientSets::rungeKuttaFehlberg78, timeStep, timeStep, 1.0e3, 1.0e3 );

    std::shared_ptr< numerical_integrators::IntegratorSettings< > > backwardIntegratorSettings =
            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    TUDAT_NAN, - timeStep, CoefficientSets::rungeKuttaFehlberg78, - timeStep, - timeStep, 1.0e3, 1.0e3 );

    // Define arc-wise initial states
    std::vector< Eigen::VectorXd > midArcStatesMoons;
    Eigen::VectorXd arcWiseMidStatesIo, arcWiseMidStatesEuropa, arcWiseMidStatesGanymede, arcWiseMidStatesCallisto;
    arcWiseMidStatesIo.resize( 6 * nbArcs );
    arcWiseMidStatesEuropa.resize( 6 * nbArcs );
    arcWiseMidStatesGanymede.resize( 6 * nbArcs );
    arcWiseMidStatesCallisto.resize( 6 * nbArcs );
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        midArcStatesMoons.push_back( propagators::getInitialStatesOfBodies( bodiesToPropagate, centralBodies, bodies, midArcTimes.at( i ) ) );
        arcWiseMidStatesIo.segment( i * 6, 6 ) = midArcStatesMoons[ i ].segment( 0, 6 );
        arcWiseMidStatesEuropa.segment( i * 6, 6 ) = midArcStatesMoons[ i ].segment( 6, 6 );
        arcWiseMidStatesGanymede.segment( i * 6, 6 ) = midArcStatesMoons[ i ].segment( 12, 6 );
        arcWiseMidStatesCallisto.segment( i * 6, 6 ) = midArcStatesMoons[ i ].segment( 18, 6 );
    }
    std::map< std::string, Eigen::VectorXd > concatenatedArcWiseStatesPerBody;
    concatenatedArcWiseStatesPerBody[ "Io" ] =  arcWiseMidStatesIo;
    concatenatedArcWiseStatesPerBody[ "Europa" ] =  arcWiseMidStatesEuropa;
    concatenatedArcWiseStatesPerBody[ "Ganymede" ] =  arcWiseMidStatesGanymede;
    concatenatedArcWiseStatesPerBody[ "Callisto" ] =  arcWiseMidStatesCallisto;

    // Define dependent variables to save
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
            basic_astrodynamics::mutual_spherical_harmonic_gravity, "Io", "Jupiter" ) );
    for ( unsigned int k = 0 ; k < bodiesToPropagate.size( ) ; k++ )
    {
        if ( bodiesToPropagate.at( k ) != "Io" )
        {
            dependentVariables.push_back( std::make_shared< SingleAccelerationDependentVariableSaveSettings >(
                    basic_astrodynamics::mutual_spherical_harmonic_gravity, "Io", bodiesToPropagate.at( k ) ) );
        }
    }

    // Define propagator settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > forwardPropagationSettingsList, backwardPropagationSettingsList, nonSequentialPropagationSettingsList;
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        forwardPropagationSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons.at( i ),  midArcTimes.at( i ), forwardIntegratorSettings,
                std::make_shared< PropagationTimeTerminationSettings >( arcEndTimes.at( i ) ), cowell, dependentVariables ) );
        backwardPropagationSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >(
                centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons.at( i ), midArcTimes.at( i ), backwardIntegratorSettings,
                std::make_shared< PropagationTimeTerminationSettings >( arcStartTimes.at( i ) ), cowell, dependentVariables ) );

        std::shared_ptr< NonSequentialPropagationTerminationSettings > terminationSettings = std::make_shared< NonSequentialPropagationTerminationSettings >(
                std::make_shared< PropagationTimeTerminationSettings >( arcEndTimes.at( i ) ),
                std::make_shared< PropagationTimeTerminationSettings >( arcStartTimes.at( i ) ) );
        std::shared_ptr< TranslationalStatePropagatorSettings< > > nonsequentialPropagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
                centralBodies, accelerationsMap, bodiesToPropagate, midArcStatesMoons.at( i ), midArcTimes.at( i ), forwardIntegratorSettings,
                terminationSettings, cowell, dependentVariables );
        nonSequentialPropagationSettingsList.push_back( nonsequentialPropagatorSettings );
    }

    std::shared_ptr< MultiArcPropagatorSettings< > > forwardPropagatorSettings = std::make_shared< MultiArcPropagatorSettings< > >( forwardPropagationSettingsList );

    std::shared_ptr< MultiArcPropagatorSettings< > > backwardPropagatorSettings = std::make_shared< MultiArcPropagatorSettings< > >( backwardPropagationSettingsList );

    std::shared_ptr< MultiArcPropagatorSettings< > > nonSequentialPropagatorSettings = std::make_shared< MultiArcPropagatorSettings< > >( nonSequentialPropagationSettingsList );


    // Define parameters to estimate for mid-arc propagation / estimation
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
    {
        parameterNames.push_back(
                std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                        bodiesToPropagate.at( i ), concatenatedArcWiseStatesPerBody.at( bodiesToPropagate.at( i ) ), midArcTimes, centralBodies.at( i ) ) );

        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >(
                bodiesToPropagate.at( i ), gravitational_parameter ) );

        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 0,  2, 2, bodiesToPropagate.at( i ), spherical_harmonics_cosine_coefficient_block ) );
        parameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 1,  2, 2, bodiesToPropagate.at( i ), spherical_harmonics_sine_coefficient_block ) );
    }
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > nonSequentialParameters =
            createParametersToEstimate< double >( parameterNames, bodies, nonSequentialPropagatorSettings );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > nominalParameters =
            nonSequentialParameters->template getFullParameterValues< double >( );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > forwardParameters =
            createParametersToEstimate< double >( parameterNames, bodies, forwardPropagatorSettings );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > backwardParameters =
            createParametersToEstimate< double >( parameterNames, bodies, backwardPropagatorSettings );


    // Define links and observations.
    std::vector< observation_models::LinkEnds > linkEndsList;
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationSettingsList;
    linkEndsList.resize( bodiesToPropagate.size( ) );
    for( int i = 0; i < bodiesToPropagate.size( ) ; i++ )
    {
        linkEndsList[ i ][ observation_models::observed_body ] = observation_models::LinkEndId( bodiesToPropagate.at( i ), "" );
        observationSettingsList.push_back( std::make_shared< observation_models::ObservationModelSettings >(
                observation_models::position_observable, linkEndsList[ i ] ) );
    }

    // Define observation times
    std::vector< double > observationTimesForward;
    std::vector< int > nbObservationsPerArcForward;
    int counterObservations = 0;
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        counterObservations = 0;
        for ( double time = midArcTimes.at( i ) + 3600.0 ; time < arcEndTimes.at( i ) - 3600.0 ; time += 3.0 * 3600.0 )
        {
            observationTimesForward.push_back( time );
            counterObservations++;
        }
        nbObservationsPerArcForward.push_back( counterObservations * 3 );
    }

    std::vector< double > observationTimesBackward;
    std::vector< int > nbObservationsPerArcBackward;
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        counterObservations = 0;
        for ( double time = arcStartTimes.at( i ) + 3600.0 ; time < midArcTimes.at( i ) - 3600.0 ; time += 3.0 * 3600.0 )
        {
            observationTimesBackward.push_back( time );
            counterObservations++;
        }
        nbObservationsPerArcBackward.push_back( counterObservations * 3 );
    }

    std::vector< double > allObservationTimes = observationTimesBackward;
    for ( unsigned int i = 0 ; i < observationTimesForward.size( ) ; i++ )
    {
        allObservationTimes.push_back( observationTimesForward.at( i ) );
    }

    // Define observation settings
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementInputForward, measurementInputBackward, measurementInputAll;
    for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
    {
        measurementInputForward.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
                observation_models::position_observable, linkEndsList[ i ], observationTimesForward, observation_models::observed_body ) );

        measurementInputBackward.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
                observation_models::position_observable, linkEndsList[ i ], observationTimesBackward, observation_models::observed_body ) );

        measurementInputAll.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
                observation_models::position_observable, linkEndsList[ i ], allObservationTimes, observation_models::observed_body ) );
    }


    // Create orbit determination object for forward propagation / estimation.
    OrbitDeterminationManager< > orbitDeterminationForward = OrbitDeterminationManager< >(
            bodies, forwardParameters, observationSettingsList, forwardPropagatorSettings );

    // Create orbit determination object for backward propagation / estimation.
    OrbitDeterminationManager< > orbitDeterminationBackward = OrbitDeterminationManager< >(
            bodies, backwardParameters, observationSettingsList, backwardPropagatorSettings );

    // Create orbit determination object for non-sequential propagation / estimation.
    nonSequentialParameters->resetParameterValues( nominalParameters );
    OrbitDeterminationManager< > orbitDeterminationNonSequential = OrbitDeterminationManager< >(
            bodies, nonSequentialParameters, observationSettingsList, nonSequentialPropagatorSettings );


    // Simulate observations for forward propagation / estimation
    std::shared_ptr< observation_models::ObservationCollection< > > observationsAndTimesForward = simulateObservations< >(
            measurementInputForward, orbitDeterminationForward.getObservationSimulators( ), bodies  );

    // Simulate observations for backward propagation / estimation
    std::shared_ptr< observation_models::ObservationCollection< > > observationsAndTimesBackward = simulateObservations< >(
            measurementInputBackward, orbitDeterminationBackward.getObservationSimulators( ), bodies  );

    // Simulate observations for non-sequential propagation / estimation
    std::shared_ptr< observation_models::ObservationCollection< > > observationsAndTimesNonSequential = simulateObservations< >(
            measurementInputAll, orbitDeterminationNonSequential.getObservationSimulators( ), bodies  );


    // Define estimation input for forward propagation / estimation
    std::shared_ptr< EstimationInput< double, double  > > estimationInputForward =
            std::make_shared< EstimationInput< double, double > >( observationsAndTimesForward, Eigen::MatrixXd::Zero( 0, 0 ), std::make_shared< EstimationConvergenceChecker >( 1 ) );

    // Define estimation input for backward propagation / estimation
    std::shared_ptr< EstimationInput< double, double  > > estimationInputBackward =
            std::make_shared< EstimationInput< double, double > >( observationsAndTimesBackward, Eigen::MatrixXd::Zero( 0, 0 ), std::make_shared< EstimationConvergenceChecker >( 1 ) );

    // Define POD input for non-sequential propagation / estimation
    std::shared_ptr< EstimationInput< double, double > > estimationInputNonSequential =
            std::make_shared< EstimationInput< double, double > >( observationsAndTimesNonSequential, Eigen::MatrixXd::Zero( 0, 0 ), std::make_shared< EstimationConvergenceChecker >( 1 ) );

    // Perform forward estimation
    std::shared_ptr< EstimationOutput< double, double > > estimationOutputForward = orbitDeterminationForward.estimateParameters( estimationInputForward );

    // Perform backward estimation
    std::shared_ptr< EstimationOutput< double, double > > estimationOutputBackward = orbitDeterminationBackward.estimateParameters( estimationInputBackward );

    // Perform non-sequential estimation
    std::shared_ptr< EstimationOutput< double, double > > estimationOutputNonSequential = orbitDeterminationNonSequential.estimateParameters( estimationInputNonSequential );

    // Retrieve partials from each estimation.
    Eigen::MatrixXd partialsForwardEstimation = estimationOutputForward->getUnnormalizedDesignMatrix( );
    Eigen::MatrixXd partialsBackwardEstimation = estimationOutputBackward->getUnnormalizedDesignMatrix( );
    Eigen::MatrixXd partialsNonSequentialEstimation = estimationOutputNonSequential->getUnnormalizedDesignMatrix( );

    Eigen::MatrixXd combinedPartials = Eigen::MatrixXd::Zero( partialsForwardEstimation.rows( ) + partialsBackwardEstimation.rows( ),
                                                              partialsForwardEstimation.cols( ) );

    // Combine forward and backward partials in proper order.
    unsigned int nbParameter = nominalParameters.size( );
    for ( unsigned int j = 0 ; j < bodiesToPropagate.size( ) ; j++ )
    {
        std::vector< std::pair< int, int > > backwardPartialsIndices =
                estimationInputBackward->getObservationCollection( )->getObservationSetStartAndSize( ).at( observation_models::position_observable ).at( linkEndsList[ j ] );
        std::vector< std::pair< int, int > > forwardPartialsIndices =
                estimationInputForward->getObservationCollection( )->getObservationSetStartAndSize( ).at( observation_models::position_observable ).at( linkEndsList[ j ] );
        std::vector< std::pair< int, int > > nonSequentialPartialsIndices =
                estimationInputNonSequential->getObservationCollection( )->getObservationSetStartAndSize( ).at( observation_models::position_observable ).at( linkEndsList[ j ] );

        int counterIndices = nonSequentialPartialsIndices.at( 0 ).first;
        int counterIndicesBackward = backwardPartialsIndices.at( 0 ).first;
        int counterIndicesForward = forwardPartialsIndices.at( 0 ).first;
        for ( unsigned int k = 0 ; k < nbObservationsPerArcBackward.size( ) ; k++ )
        {
            //! Add partials from backward estimation
            combinedPartials.block( counterIndices, 0, nbObservationsPerArcBackward.at( k ), nbParameter ) =
                    partialsBackwardEstimation.block( counterIndicesBackward, 0, nbObservationsPerArcBackward.at( k ), nbParameter );

            counterIndices += nbObservationsPerArcBackward.at( k );
            counterIndicesBackward += nbObservationsPerArcBackward.at( k );

            //! Add partials from forward estimation
            combinedPartials.block( counterIndices, 0, nbObservationsPerArcForward.at( k ), nbParameter ) =
                    partialsForwardEstimation.block( counterIndicesForward, 0, nbObservationsPerArcForward.at( k ), nbParameter );

            counterIndices += nbObservationsPerArcForward.at( k );
            counterIndicesForward += nbObservationsPerArcForward.at( k );
        }
    }

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( combinedPartials, partialsNonSequentialEstimation, 1.0e-12 );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}