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

BOOST_AUTO_TEST_SUITE( test_consider_parameters )

BOOST_AUTO_TEST_CASE( testConsiderParameters )
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
    double finalEpoch = 2.0 * 86400.0;
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
    std::shared_ptr< numerical_integrators::IntegratorSettings< > > integratorSettings =
            std::make_shared< numerical_integrators::RungeKuttaVariableStepSizeSettingsScalarTolerances< double > >(
                    TUDAT_NAN, timeStep, CoefficientSets::rungeKuttaFehlberg78, timeStep, timeStep, 1.0e3, 1.0e3 );

    // Define arc-wise initial states
    std::vector< Eigen::VectorXd > initialStatesMoons;
    Eigen::VectorXd initialStatesIo, initialStatesEuropa, initialStatesGanymede, initialStatesCallisto;
    initialStatesIo.resize( 6 * nbArcs );
    initialStatesEuropa.resize( 6 * nbArcs );
    initialStatesGanymede.resize( 6 * nbArcs );
    initialStatesCallisto.resize( 6 * nbArcs );
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        initialStatesMoons.push_back( propagators::getInitialStatesOfBodies( bodiesToPropagate, centralBodies, bodies, midArcTimes.at( i ) ) );
        initialStatesIo.segment( i * 6, 6 ) = initialStatesMoons[ i ].segment( 0, 6 );
        initialStatesEuropa.segment( i * 6, 6 ) = initialStatesMoons[ i ].segment( 6, 6 );
        initialStatesGanymede.segment( i * 6, 6 ) = initialStatesMoons[ i ].segment( 12, 6 );
        initialStatesCallisto.segment( i * 6, 6 ) = initialStatesMoons[ i ].segment( 18, 6 );
    }
    std::map< std::string, Eigen::VectorXd > initialStatesPerBody;
    initialStatesPerBody[ "Io" ] = initialStatesIo;
    initialStatesPerBody[ "Europa" ] = initialStatesEuropa;
    initialStatesPerBody[ "Ganymede" ] = initialStatesGanymede;
    initialStatesPerBody[ "Callisto" ] = initialStatesCallisto;

    // Define propagator settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagationSettingsList;
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        std::shared_ptr< NonSequentialPropagationTerminationSettings > terminationSettings = std::make_shared< NonSequentialPropagationTerminationSettings >(
                std::make_shared< PropagationTimeTerminationSettings >( arcEndTimes.at( i ) ),
                std::make_shared< PropagationTimeTerminationSettings >( arcStartTimes.at( i ) ) );
        std::shared_ptr< TranslationalStatePropagatorSettings< > > propagatorSettings = std::make_shared< TranslationalStatePropagatorSettings< > >(
                centralBodies, accelerationsMap, bodiesToPropagate, initialStatesMoons.at( i ), midArcTimes.at( i ), integratorSettings, terminationSettings );
        propagationSettingsList.push_back( propagatorSettings );
    }
    std::shared_ptr< MultiArcPropagatorSettings< > > propagatorSettings = std::make_shared< MultiArcPropagatorSettings< > >( propagationSettingsList );


    // Define parameters to estimate for non-sequential propagation / estimation
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames, considerParameterNames, parameterNamesAll;
    for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
    {
        parameterNames.push_back( std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                bodiesToPropagate.at( i ), initialStatesPerBody.at( bodiesToPropagate.at( i ) ), midArcTimes, centralBodies.at( i ) ) );

        parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( bodiesToPropagate.at( i ), gravitational_parameter ) );
    }

    for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
    {
        considerParameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 0,  2, 2, bodiesToPropagate.at( i ), spherical_harmonics_cosine_coefficient_block ) );
        considerParameterNames.push_back( std::make_shared< SphericalHarmonicEstimatableParameterSettings >(
                2, 1,  2, 2, bodiesToPropagate.at( i ), spherical_harmonics_sine_coefficient_block ) );
    }
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parameters =
            createParametersToEstimate< double >( parameterNames, bodies, propagatorSettings, considerParameterNames );
    printEstimatableParameterEntries( parameters );

    // Get nominal parameter values
    Eigen::VectorXd nominalParameters = parameters->getFullParameterValues< double >( );
    int nbEstimatedParameters = nominalParameters.size( );

    // Create parameters object with all parameters estimated
    parameterNamesAll = parameterNames;
    for ( unsigned int i = 0 ; i < considerParameterNames.size( ) ; i++ )
    {
        parameterNamesAll.push_back( considerParameterNames[ i ] );
    }
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersAll =
            createParametersToEstimate< double >( parameterNamesAll, bodies, propagatorSettings );
    printEstimatableParameterEntries( parametersAll );


    // Define links and observations.
    std::vector< observation_models::LinkEnds > linkEndsList;
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationSettingsList;
    linkEndsList.resize( bodiesToPropagate.size( ) );
    for( unsigned int i = 0; i < bodiesToPropagate.size( ) ; i++ )
    {
        linkEndsList[ i ][ observation_models::observed_body ] = observation_models::LinkEndId( bodiesToPropagate.at( i ), "" );
        observationSettingsList.push_back( std::make_shared< observation_models::ObservationModelSettings >(
                observation_models::position_observable, linkEndsList[ i ] ) );
    }

    // Define observation times
    std::vector< double > observationTimes;
    for ( unsigned int i = 0 ; i < nbArcs ; i++ )
    {
        for ( double time = arcStartTimes.at( i ) + 3600.0 ; time < arcEndTimes.at( i ) - 3600.0 ; time += 3.0 * 3600.0 )
        {
            observationTimes.push_back( time );
        }
    }

    // Define observation settings
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementInput;
    for( unsigned int i = 0; i < bodiesToPropagate.size( ); i++ )
    {
        measurementInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
                observation_models::position_observable, linkEndsList[ i ], observationTimes, observation_models::observed_body ) );
    }

    // Create orbit determination managers
    OrbitDeterminationManager< > orbitDeterminationManager = OrbitDeterminationManager< >(
            bodies, parameters, observationSettingsList, propagatorSettings );
    OrbitDeterminationManager< > orbitDeterminationManagerAll = OrbitDeterminationManager< >(
            bodies, parametersAll, observationSettingsList, propagatorSettings );

    // Simulate observations
    std::shared_ptr< observation_models::ObservationCollection< > > observationsAndTimes = simulateObservations< >(
            measurementInput, orbitDeterminationManager.getObservationSimulators( ), bodies  );
    std::shared_ptr< observation_models::ObservationCollection< > > observationsAndTimesAll = simulateObservations< >(
            measurementInput, orbitDeterminationManagerAll.getObservationSimulators( ), bodies  );

    // Define consider covariance
    int nbConsiderParameters = parameters->getConsiderParameters( )->getEstimatedParameterSetSize( );
    Eigen::VectorXd considerParametersValues = parameters->getConsiderParameters( )->getFullParameterValues< double >( );
    Eigen::MatrixXd considerCovariance = Eigen::MatrixXd::Zero( nbConsiderParameters, nbConsiderParameters );
    for ( int i = 0 ; i < nbConsiderParameters ; i++ )
    {
        considerCovariance( i, i ) = ( 0.1 * considerParametersValues[ i ] ) * ( 0.1 * considerParametersValues[ i ] );
    }

    // Define estimation input when estimating all parameters
    std::shared_ptr< EstimationInput< double, double  > > estimationInputAll = std::make_shared< EstimationInput< double, double > >(
            observationsAndTimesAll, Eigen::MatrixXd::Zero( 0, 0 ), std::make_shared< EstimationConvergenceChecker >( 1 ) );
    std::shared_ptr< CovarianceAnalysisInput< double, double  > > covarianceInputAll = std::make_shared< CovarianceAnalysisInput< double, double > >(
            observationsAndTimesAll, Eigen::MatrixXd::Zero( 0, 0 ) );

    // Perform estimation for all parameters
    std::shared_ptr< CovarianceAnalysisOutput< double, double > > covarianceOutputAll = orbitDeterminationManagerAll.computeCovariance( covarianceInputAll );
    std::shared_ptr< EstimationOutput< double, double > > estimationOutputAll = orbitDeterminationManagerAll.estimateParameters( estimationInputAll );


    // Retrieve covariance matrix when estimating all parameters
    Eigen::MatrixXd covarianceAll = covarianceOutputAll->getUnnormalizedCovarianceMatrix( );

    // Define estimation input with consider parameters
    Eigen::VectorXd considerParametersDeviations = 0.1 * considerParametersValues;
    std::shared_ptr< EstimationInput< double, double  > > estimationInput = std::make_shared< EstimationInput< double, double > >(
            observationsAndTimes, Eigen::MatrixXd::Zero( 0, 0 ), std::make_shared< EstimationConvergenceChecker >( 1 ), considerCovariance, considerParametersDeviations );
    std::shared_ptr< CovarianceAnalysisInput< double, double  > > covarianceInput = std::make_shared< CovarianceAnalysisInput< double, double > >(
            observationsAndTimes, Eigen::MatrixXd::Zero( 0, 0 ), considerCovariance );

    // Perform estimation with consider parameters
    std::shared_ptr< CovarianceAnalysisOutput< double, double> > covarianceOutput = orbitDeterminationManager.computeCovariance( covarianceInput );
    std::shared_ptr< EstimationOutput< double, double > > estimationOutput = orbitDeterminationManager.estimateParameters( estimationInput );
    Eigen::VectorXd updatedParameters = estimationOutput->parameterEstimate_;

    // Retrieve covariance matrix
    Eigen::MatrixXd covariance = covarianceOutput->getUnnormalizedCovarianceMatrix( );

    // Retrieve unnormalised partials estimated and consider parameters
    Eigen::MatrixXd fullPartials = estimationOutputAll->getUnnormalizedDesignMatrix( );
    unsigned int nbObservations = fullPartials.rows( );
    Eigen::MatrixXd considerPartials = fullPartials.block( 0, nbEstimatedParameters, nbObservations, nbConsiderParameters );
    Eigen::MatrixXd estimatedPartials = fullPartials.block( 0, 0, nbObservations, nbEstimatedParameters );
    Eigen::MatrixXd weightedEstimatedPartials = linear_algebra::multiplyDesignMatrixByDiagonalWeightMatrix(
            estimatedPartials, estimationInput->getWeightsMatrixDiagonals( ) );

    // Manually compute normalised inverse covariance (w/o consider parameters)
    Eigen::MatrixXd normalisedEstimatedPartials = estimationOutputAll->getNormalizedDesignMatrix( ).block( 0, 0, nbObservations, nbEstimatedParameters );
    Eigen::MatrixXd normalisedConsiderPartials = estimationOutputAll->getNormalizedDesignMatrix( ).block( 0, nbEstimatedParameters, nbObservations, nbConsiderParameters );
    Eigen::MatrixXd weightedNormalisedEstimatedPartials = linear_algebra::multiplyDesignMatrixByDiagonalWeightMatrix(
            normalisedEstimatedPartials, estimationInput->getWeightsMatrixDiagonals( ) );
    Eigen::MatrixXd computedNormalisedInvCovariance = normalisedEstimatedPartials.transpose( ) * weightedNormalisedEstimatedPartials;

    // Retrieve normalisation terms
    Eigen::VectorXd normalisationTerms = estimationOutputAll->designMatrixTransformationDiagonal_.segment( 0, nbEstimatedParameters );
    Eigen::VectorXd normalisationTermsConsider = estimationOutputAll->designMatrixTransformationDiagonal_.segment( nbEstimatedParameters, nbConsiderParameters );

    // Invert covariance matrix
    Eigen::MatrixXd computedNormalizedCovariance = computedNormalisedInvCovariance.inverse( );

    //  Check consistency normalized covariance matrix without consider contributions
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( covarianceOutput->normalizedCovarianceMatrix_, computedNormalizedCovariance, 1.0e-12 );

    // Compute unnormalized covariance matrix
    Eigen::MatrixXd computedUnnormalizedCovariance = computedNormalizedCovariance;
    for ( int i = 0 ; i < normalisationTerms.rows( ) ; i++ )
    {
        for ( int j = 0 ; j < normalisationTerms.rows( ); j++ )
        {
            computedUnnormalizedCovariance( i, j ) /= ( normalisationTerms( i ) * normalisationTerms( j ) );
        }
    }

    // Compute normalized consider covariance
    Eigen::MatrixXd normalizedConsiderCovariance = orbitDeterminationManagerAll.normalizeCovariance( considerCovariance, normalisationTermsConsider );

    // Compute (unnormalised) contribution of consider parameters
    Eigen::MatrixXd computedContributionConsiderCovariance = ( computedUnnormalizedCovariance * weightedEstimatedPartials.transpose( ) )
    * ( considerPartials * considerCovariance * considerPartials.transpose( ) )
    * ( computedUnnormalizedCovariance * weightedEstimatedPartials.transpose( ) ).transpose( );

    // Add contribution consider parameters
    Eigen::MatrixXd computedNormalizedCovarianceConsiderParameters = computedNormalizedCovariance
            + ( computedNormalizedCovariance * weightedNormalisedEstimatedPartials.transpose( ) )
            * ( normalisedConsiderPartials * normalizedConsiderCovariance * normalisedConsiderPartials.transpose( ) )
            * ( computedNormalizedCovariance * weightedNormalisedEstimatedPartials.transpose( ) ).transpose( );

    // Unnormalise covariance matrix (with consider parameters)
    Eigen::MatrixXd computedCovarianceConsiderParameters = computedNormalizedCovarianceConsiderParameters;
    for ( int i = 0 ; i < normalisationTerms.rows( ) ; i++ )
    {
        for ( int j = 0 ; j < normalisationTerms.rows( ); j++ )
        {
            computedCovarianceConsiderParameters( i, j ) /= ( normalisationTerms( i ) * normalisationTerms( j ) );
        }
    }

    // Compare with manually computed covariance
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( covarianceOutput->unnormalizedCovarianceWithConsiderParameters_, computedCovarianceConsiderParameters, 1.0e-12 );

    // Manually perform full estimation step
    Eigen::VectorXd rightHandSide = normalisedEstimatedPartials.transpose( ) * estimationInput->getWeightsMatrixDiagonals( ).cwiseProduct( normalisedConsiderPartials *
            considerParametersDeviations.cwiseProduct( normalisationTermsConsider ) );
    Eigen::MatrixXd leftHandSide = computedNormalisedInvCovariance;

    Eigen::VectorXd leastSquaresOutput = linear_algebra::solveSystemOfEquationsWithSvd( leftHandSide, rightHandSide, 1, 1.0e8 );
    Eigen::VectorXd computedUpdatedParameters = leastSquaresOutput.cwiseQuotient( normalisationTerms ) + nominalParameters;

    // Check consistency
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( updatedParameters, computedUpdatedParameters, 1.0e-12 );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}