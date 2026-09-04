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

#include <boost/test/included/unit_test.hpp>
#include <Eigen/LU>

#include <limits>

#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/simulation/estimation_setup/executePlanetaryParameterEstimationTestCase.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManagerHelpers.h"

namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_estimation_input_output )

//! Initial translational state with a user-defined hard linear constraint on each estimated correction.
class LinearlyConstrainedInitialTranslationalStateParameter : public estimatable_parameters::InitialTranslationalStateParameter< double >
{
public:
    LinearlyConstrainedInitialTranslationalStateParameter( const std::string& associatedBody,
                                                           const Eigen::Vector6d& initialState,
                                                           const std::string& centralBody ):
        InitialTranslationalStateParameter< double >( associatedBody, initialState, centralBody )
    {}

    int getConstraintSize( ) override
    {
        return 1;
    }

    Eigen::MatrixXd getConstraintStateMultipler( ) override
    {
        Eigen::MatrixXd constraintMultiplier = Eigen::MatrixXd::Zero( 1, 6 );
        constraintMultiplier( 0, 0 ) = 1.0;
        constraintMultiplier( 0, 4 ) = 1000.0;
        return constraintMultiplier;
    }

    Eigen::VectorXd getConstraintRightHandSide( ) override
    {
        return Eigen::VectorXd::Zero( 1 );
    }
};

//! Test that a physical-space constraint remains valid when design-matrix columns have mixed scales.
BOOST_AUTO_TEST_CASE( test_NormalizedLinearConstraints )
{
    Eigen::MatrixXd normalizedDesignMatrix( 3, 2 );
    normalizedDesignMatrix << 1.0, 0.0, 0.0, 1.0, 1.0, 1.0;
    const Eigen::VectorXd designMatrixNormalization = ( Eigen::Vector2d( ) << 1.0E6, 1.0E-3 ).finished( );
    const Eigen::VectorXd residuals = ( Eigen::Vector3d( ) << 0.25, -0.5, 0.1 ).finished( );
    const Eigen::VectorXd weights = Eigen::Vector3d::Ones( );

    Eigen::MatrixXd physicalConstraint( 1, 2 );
    physicalConstraint << 1.0, 1.0;
    const Eigen::VectorXd physicalConstraintRightHandSide = Eigen::VectorXd::Constant( 1, 1.0 );

    Eigen::MatrixXd normalizedConstraint = physicalConstraint;
    Eigen::VectorXd normalizedConstraintRightHandSide = physicalConstraintRightHandSide;
    simulation_setup::normalizeLinearConstraints( normalizedConstraint, normalizedConstraintRightHandSide, designMatrixNormalization );

    const auto leastSquaresOutput =
            linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix( normalizedDesignMatrix,
                                                                           residuals,
                                                                           weights,
                                                                           Eigen::Matrix2d::Zero( ),
                                                                           std::numeric_limits< double >::quiet_NaN( ),
                                                                           normalizedConstraint,
                                                                           normalizedConstraintRightHandSide );
    const Eigen::Vector2d physicalCorrection = leastSquaresOutput.first.head( 2 ).cwiseQuotient( designMatrixNormalization );

    BOOST_CHECK_SMALL( std::fabs( ( physicalConstraint * physicalCorrection - physicalConstraintRightHandSide )( 0 ) ), 1.0E-12 );
    BOOST_CHECK_CLOSE_FRACTION( normalizedConstraint.cwiseAbs( ).maxCoeff( ), 1.0, 1.0E-15 );
}

//! Test that a priori information constrains the total deviation from the initial parameter vector.
BOOST_AUTO_TEST_CASE( test_APrioriParameterDeviation )
{
    Eigen::MatrixXd designMatrix( 3, 2 );
    designMatrix << 1.0, 2.0, -1.0, 0.5, 0.0, 1.0;
    const Eigen::Vector3d residuals = ( Eigen::Vector3d( ) << 4.0, -2.0, 1.0 ).finished( );
    const Eigen::Vector3d weights = ( Eigen::Vector3d( ) << 2.0, 3.0, 5.0 ).finished( );
    Eigen::Matrix2d inverseAprioriCovariance;
    inverseAprioriCovariance << 4.0, 1.0, 1.0, 3.0;
    const Eigen::Vector2d aprioriParameterDeviation = ( Eigen::Vector2d( ) << 0.5, -0.25 ).finished( );

    const auto leastSquaresOutput =
            linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix( designMatrix,
                                                                           residuals,
                                                                           weights,
                                                                           inverseAprioriCovariance,
                                                                           std::numeric_limits< double >::quiet_NaN( ),
                                                                           Eigen::MatrixXd( 0, 0 ),
                                                                           Eigen::VectorXd( 0 ),
                                                                           Eigen::MatrixXd( 0, 0 ),
                                                                           Eigen::VectorXd( 0 ),
                                                                           Eigen::MatrixXd( 0, 0 ),
                                                                           Eigen::VectorXd( 0 ),
                                                                           aprioriParameterDeviation );

    const Eigen::Matrix2d expectedNormalMatrix =
            designMatrix.transpose( ) * weights.asDiagonal( ) * designMatrix + inverseAprioriCovariance;
    const Eigen::Vector2d expectedRightHandSide =
            designMatrix.transpose( ) * weights.asDiagonal( ) * residuals - inverseAprioriCovariance * aprioriParameterDeviation;
    const Eigen::Vector2d expectedParameterCorrection = expectedNormalMatrix.inverse( ) * expectedRightHandSide;

    BOOST_CHECK_SMALL( ( leastSquaresOutput.second - expectedNormalMatrix ).norm( ), 1.0E-14 );
    BOOST_CHECK_SMALL( ( leastSquaresOutput.first - expectedParameterCorrection ).norm( ), 1.0E-14 );

    // Exercise the same constraint through a complete iterative estimation. With zero acceleration and a position
    // observation at the initial epoch, the estimation problem is linear. The first update therefore reaches the MAP
    // solution. Although an observation residual remains, the second update must be zero because the growing deviation
    // from the fixed initial parameter vector supplies an equal and opposite a priori contribution.
    using namespace observation_models;

    const double initialTime = 1.0E7;
    const double finalTime = initialTime + 60.0;
    SystemOfBodies bodies( "SSB", "ECLIPJ2000" );
    bodies.createEmptyBody( "Satellite" );
    bodies.at( "Satellite" )->setConstantBodyMass( 100.0 );

    const std::vector< std::string > bodiesToPropagate = { "Satellite" };
    const std::vector< std::string > centralBodies = { "SSB" };
    AccelerationMap accelerationModels;
    accelerationModels[ "Satellite" ] = basic_astrodynamics::SingleBodyAccelerationMap( );

    Eigen::Vector6d trueInitialState;
    trueInitialState << 7.0E6, 2.0E6, -1.0E6, 100.0, 200.0, -50.0;
    const std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            rungeKuttaFixedStepSettings( 10.0, CoefficientSets::rungeKuttaFehlberg78 );
    const std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double, double > >( centralBodies,
                                                                                        accelerationModels,
                                                                                        bodiesToPropagate,
                                                                                        trueInitialState,
                                                                                        initialTime,
                                                                                        integratorSettings,
                                                                                        propagationTimeTerminationSettings( finalTime ),
                                                                                        cowell );

    const std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterSettings =
            getInitialStateParameterSettings< double, double >( propagatorSettings, bodies );
    const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterSettings, bodies );

    LinkEnds linkEnds;
    linkEnds[ observed_body ] = LinkEndId( "Satellite", "" );
    const std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettings = {
        std::make_shared< ObservationModelSettings >( position_observable, linkEnds )
    };
    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, observationModelSettings, propagatorSettings );

    const std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > observationSimulationSettings = {
        std::make_shared< TabulatedObservationSimulationSettings< double > >(
                position_observable, linkEnds, std::vector< double >{ finalTime }, observed_body )
    };
    const std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateObservations< double, double >(
            observationSimulationSettings, orbitDeterminationManager.getObservationSimulators( ), bodies );
    simulatedObservations->setConstantWeight( 1.0 );

    Eigen::Vector6d aprioriInitialState = trueInitialState;
    aprioriInitialState( 0 ) += 10.0;
    orbitDeterminationManager.resetParameterEstimate( aprioriInitialState, true );

    const Eigen::MatrixXd fullEstimationInverseAprioriCovariance = 1.0E4 * Eigen::MatrixXd::Identity( 6, 6 );
    const std::shared_ptr< EstimationInput< double, double > > estimationInput = std::make_shared< EstimationInput< double, double > >(
            simulatedObservations, fullEstimationInverseAprioriCovariance, estimationConvergenceChecker( 2, -1.0, -1.0, 100 ) );
    estimationInput->defineEstimationSettings( true, true, true, false, true, false );

    const std::shared_ptr< EstimationOutput< double, double > > estimationOutput =
            orbitDeterminationManager.estimateParameters( estimationInput );
    BOOST_REQUIRE( !estimationOutput->exceptionDuringInversion_ );
    BOOST_REQUIRE( !estimationOutput->exceptionDuringPropagation_ );
    BOOST_REQUIRE_EQUAL( estimationOutput->residualHistory_.size( ), 2 );

    const Eigen::MatrixXd parameterHistory = estimationOutput->getParameterHistoryMatrix( );
    BOOST_REQUIRE_EQUAL( parameterHistory.rows( ), 6 );
    BOOST_REQUIRE_EQUAL( parameterHistory.cols( ), 3 );
    const Eigen::Vector6d deviationAfterFirstIteration = parameterHistory.col( 1 ) - parameterHistory.col( 0 );
    const Eigen::Vector6d secondParameterCorrection = parameterHistory.col( 2 ) - parameterHistory.col( 1 );

    BOOST_CHECK_GT( deviationAfterFirstIteration.norm( ), 1.0E-3 );
    BOOST_CHECK_GT( estimationOutput->residualHistory_.at( 1 ).norm( ), 1.0 );
    BOOST_CHECK_LT( estimationOutput->residualHistory_.at( 1 ).norm( ), estimationOutput->residualHistory_.at( 0 ).norm( ) );
    BOOST_CHECK_SMALL( secondParameterCorrection.norm( ), 1.0E-6 * deviationAfterFirstIteration.norm( ) );
}

//! Test a constrained least-squares orbit estimation with position and velocity corrections at different scales.
BOOST_AUTO_TEST_CASE( test_LinearConstraintInEarthSatelliteEstimation )
{
    using namespace observation_models;
    using namespace orbital_element_conversions;

    spice_interface::loadStandardSpiceKernels( );

    const double initialTime = 1.0E7;
    const double finalTime = initialTime + 3.0 * 3600.0;

    BodyListSettings bodySettings = getDefaultBodySettings( { "Earth", "Sun", "Moon" }, "Earth", "ECLIPJ2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Satellite" );
    bodies.at( "Satellite" )->setConstantBodyMass( 400.0 );

    // Use a modestly perturbed Earth-orbit model: Earth degree/order 2 gravity and point-mass Sun/Moon gravity.
    SelectedAccelerationMap accelerationSettings;
    accelerationSettings[ "Satellite" ][ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
    accelerationSettings[ "Satellite" ][ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationSettings[ "Satellite" ][ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );

    const std::vector< std::string > bodiesToPropagate = { "Satellite" };
    const std::vector< std::string > centralBodies = { "Earth" };
    const AccelerationMap accelerationModels =
            createAccelerationModelsMap( bodies, accelerationSettings, bodiesToPropagate, centralBodies );

    Eigen::Vector6d initialKeplerianState;
    initialKeplerianState << 7200.0E3, 0.01, unit_conversions::convertDegreesToRadians( 55.0 ),
            unit_conversions::convertDegreesToRadians( 40.0 ), unit_conversions::convertDegreesToRadians( 20.0 ),
            unit_conversions::convertDegreesToRadians( 10.0 );
    const Eigen::Vector6d trueInitialState = convertKeplerianToCartesianElements(
            initialKeplerianState, bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( ) );

    const std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            rungeKuttaFixedStepSettings( 60.0, CoefficientSets::rungeKuttaFehlberg78 );
    const std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double, double > >( centralBodies,
                                                                                        accelerationModels,
                                                                                        bodiesToPropagate,
                                                                                        trueInitialState,
                                                                                        initialTime,
                                                                                        integratorSettings,
                                                                                        propagationTimeTerminationSettings( finalTime ),
                                                                                        cowell );

    const std::shared_ptr< LinearlyConstrainedInitialTranslationalStateParameter > constrainedInitialState =
            std::make_shared< LinearlyConstrainedInitialTranslationalStateParameter >( "Satellite", trueInitialState, "Earth" );
    constrainedInitialState->addStateClosureFunctions(
            [ propagatorSettings ]( ) { return propagatorSettings->getStateOfBody( 0 ); },
            [ propagatorSettings ]( const Eigen::VectorXd& state ) { propagatorSettings->setStateOfBody( 0, state ); } );

    const std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > scalarParameters;
    const std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > vectorParameters;
    const std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > initialStateParameters = {
        constrainedInitialState
    };
    const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            std::make_shared< estimatable_parameters::EstimatableParameterSet< double > >(
                    scalarParameters, vectorParameters, initialStateParameters );

    LinkEnds linkEnds;
    linkEnds[ observed_body ] = LinkEndId( "Satellite", "" );
    const std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettings = {
        std::make_shared< ObservationModelSettings >( position_observable, linkEnds )
    };

    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, observationModelSettings, propagatorSettings );

    std::vector< double > observationTimes;
    for( double observationTime = initialTime; observationTime <= finalTime; observationTime += 300.0 )
    {
        observationTimes.push_back( observationTime );
    }
    const std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > observationSimulationSettings = {
        std::make_shared< TabulatedObservationSimulationSettings< double > >(
                position_observable, linkEnds, observationTimes, observed_body )
    };
    const std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateObservations< double, double >(
            observationSimulationSettings, orbitDeterminationManager.getObservationSimulators( ), bodies );

    Eigen::Vector6d initialStatePerturbation;
    initialStatePerturbation << 1.0E3, 1.0E3, 1.0E3, 1.0, 1.0, 1.0;
    const Eigen::Vector6d perturbedInitialState = trueInitialState + initialStatePerturbation;
    orbitDeterminationManager.resetParameterEstimate( perturbedInitialState, true );

    const std::shared_ptr< EstimationInput< double, double > > estimationInput = std::make_shared< EstimationInput< double, double > >(
            simulatedObservations, Eigen::MatrixXd::Zero( 0, 0 ), estimationConvergenceChecker( 4, 0.0, 0.0, 100 ) );
    estimationInput->defineEstimationSettings( true, true, true, true, true );

    const std::shared_ptr< EstimationOutput< double, double > > estimationOutput =
            orbitDeterminationManager.estimateParameters( estimationInput );
    BOOST_REQUIRE( !estimationOutput->exceptionDuringInversion_ );
    BOOST_REQUIRE( !estimationOutput->exceptionDuringPropagation_ );

    const Eigen::MatrixXd initialStateHistory = estimationOutput->getParameterHistoryMatrix( );
    BOOST_REQUIRE_EQUAL( initialStateHistory.rows( ), 6 );
    BOOST_REQUIRE_GE( initialStateHistory.cols( ), 2 );
    for( int stateIndex = 0; stateIndex < 6; stateIndex++ )
    {
        BOOST_CHECK_SMALL( initialStateHistory( stateIndex, 0 ) - perturbedInitialState( stateIndex ), 1.0E-12 );
    }

    for( int iteration = 1; iteration < initialStateHistory.cols( ); iteration++ )
    {
        const Eigen::Vector6d stateCorrection = initialStateHistory.col( iteration ) - initialStateHistory.col( iteration - 1 );
        const double constraintResidual = stateCorrection( 0 ) + 1000.0 * stateCorrection( 4 );
        const double constraintScale =
                std::max( 1.0, std::max( std::fabs( stateCorrection( 0 ) ), std::fabs( 1000.0 * stateCorrection( 4 ) ) ) );
        BOOST_CHECK_SMALL( constraintResidual, 1.0E-9 * constraintScale );
    }

    const Eigen::Vector6d totalStateCorrection = initialStateHistory.col( initialStateHistory.cols( ) - 1 ) - initialStateHistory.col( 0 );
    BOOST_CHECK_SMALL( totalStateCorrection( 0 ) + 1000.0 * totalStateCorrection( 4 ), 1.0E-7 );
}

//! This test checks whether the input/output of the estimation (weights, a priori covariance, unscaled covariance) are
//! correctly handed
BOOST_AUTO_TEST_CASE( test_EstimationInputAndOutput )
{
    int simulationType = 0;

    Eigen::VectorXd parameterPerturbation = getDefaultInitialParameterPerturbation( );

    // Define stringent a priori covariance
    Eigen::MatrixXd inverseAPrioriCovariance = 1.0E32 * Eigen::MatrixXd::Identity( 7, 7 );

    // Define moderate a priori covariance
    Eigen::MatrixXd moderateInverseAPriopriCovariance = Eigen::MatrixXd::Zero( 7, 7 );
    for( unsigned int i = 0; i < 7; i++ )
    {
        moderateInverseAPriopriCovariance( i, i ) = 1.0 / ( 1.0E-6 * parameterPerturbation( i ) * parameterPerturbation( i ) );
    }

    // Run estimation with strong a priori covariance
    std::pair< std::shared_ptr< EstimationOutput< double > >, Eigen::VectorXd > estimationOutputWithAprioriCovariance =
            executePlanetaryParameterEstimation< double, double >( simulationType, parameterPerturbation, inverseAPrioriCovariance );

    int numberOfSavedParameterVectors = estimationOutputWithAprioriCovariance.first->parameterHistory_.size( );
    int numberOfSavedResidualVectors = estimationOutputWithAprioriCovariance.first->residualHistory_.size( );

    BOOST_CHECK_EQUAL( numberOfSavedParameterVectors, numberOfSavedResidualVectors );

    // Run estimation with effectively zero covariance
    std::pair< std::shared_ptr< EstimationOutput< double > >, Eigen::VectorXd > estimationOutputWithSmallAprioriCovariance =
            executePlanetaryParameterEstimation< double, double >(
                    simulationType, parameterPerturbation, 1.0E-64 * inverseAPrioriCovariance );

    // Run estimation with moderate a priori covariance
    std::pair< std::shared_ptr< EstimationOutput< double > >, Eigen::VectorXd > estimationOutputWithModerateAprioriCovariance =
            executePlanetaryParameterEstimation< double, double >(
                    simulationType, parameterPerturbation, moderateInverseAPriopriCovariance );

    // Run estimation without a priori covariance
    std::pair< std::shared_ptr< EstimationOutput< double > >, Eigen::VectorXd > estimationOutputWithoutAprioriCovariance =
            executePlanetaryParameterEstimation< double, double >( simulationType, parameterPerturbation );

    // Run estimation without a priori covariance and increased weights
    double constantWeight = 100.0;
    std::pair< std::shared_ptr< EstimationOutput< double > >, Eigen::VectorXd > estimationOutputWithoutAprioriCovarianceAndWeakWeight =
            executePlanetaryParameterEstimation< double, double >(
                    simulationType, parameterPerturbation, Eigen::MatrixXd::Zero( 7, 7 ), constantWeight );

    // Retrieve estimation errors and a priori covariances
    Eigen::MatrixXd tightConstraintInverseCovariance =
            estimationOutputWithAprioriCovariance.first->getUnnormalizedInverseCovarianceMatrix( );
    Eigen::MatrixXd weakConstraintInverseCovariance =
            estimationOutputWithSmallAprioriCovariance.first->getUnnormalizedInverseCovarianceMatrix( );
    Eigen::MatrixXd moderateConstraintInverseCovariance =
            estimationOutputWithModerateAprioriCovariance.first->getUnnormalizedInverseCovarianceMatrix( );
    Eigen::MatrixXd noConstraintInverseCovariance =
            estimationOutputWithoutAprioriCovariance.first->getUnnormalizedInverseCovarianceMatrix( );
    Eigen::MatrixXd noConstraintInverseCovarianceWithWeakWeight =
            estimationOutputWithoutAprioriCovarianceAndWeakWeight.first->getUnnormalizedInverseCovarianceMatrix( );

    Eigen::VectorXd tightConstraintError = estimationOutputWithAprioriCovariance.second;
    Eigen::VectorXd weakConstraintError = estimationOutputWithSmallAprioriCovariance.second;
    Eigen::VectorXd moderateConstraintError = estimationOutputWithModerateAprioriCovariance.second;
    Eigen::VectorXd noConstraintError = estimationOutputWithoutAprioriCovariance.second;
    Eigen::VectorXd noConstraintWeakWeightError = estimationOutputWithoutAprioriCovarianceAndWeakWeight.second;

    // Check if (effectively) unconstrained solutions converge at expected level
    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( weakConstraintError( i ) ), 1.0E-2 );
        BOOST_CHECK_SMALL( std::fabs( weakConstraintError( i + 3 ) ), 1.0E-7 );

        BOOST_CHECK_SMALL( std::fabs( noConstraintError( i ) ), 1.0E-2 );
        BOOST_CHECK_SMALL( std::fabs( noConstraintError( i + 3 ) ), 1.0E-7 );

        BOOST_CHECK_SMALL( std::fabs( noConstraintWeakWeightError( i ) ), 1.0E-2 );
        BOOST_CHECK_SMALL( std::fabs( noConstraintWeakWeightError( i + 3 ) ), 1.0E-7 );
    }

    BOOST_CHECK_SMALL( std::fabs( weakConstraintError( 6 ) ), 500.0 );
    BOOST_CHECK_SMALL( std::fabs( noConstraintError( 6 ) ), 500.0 );
    BOOST_CHECK_SMALL( std::fabs( noConstraintWeakWeightError( 6 ) ), 500.0 );

    for( unsigned int i = 0; i < 7; i++ )
    {
        // Check if moderately constrained solution has intermediate accuracy
        BOOST_CHECK_EQUAL( std::fabs( moderateConstraintError( i ) ) > std::fabs( noConstraintError( i ) ), true );
        BOOST_CHECK_EQUAL( std::fabs( moderateConstraintError( i ) ) < std::fabs( tightConstraintError( i ) ), true );

        // Check if very tightly constrained solution has not differed from a priori error
        BOOST_CHECK_CLOSE_FRACTION( tightConstraintError( i ), parameterPerturbation( i ), 1.0E-8 );

        for( unsigned int j = 0; j < 7; j++ )
        {
            // Check if weights are correctly processed into covarince
            BOOST_CHECK_CLOSE_FRACTION(
                    constantWeight * noConstraintInverseCovariance( i, j ), noConstraintInverseCovarianceWithWeakWeight( i, j ), 1.0E-8 );

            // Check if tight a priori constraints are processed correctly to a posteriori covariance
            if( i == j )
            {
                BOOST_CHECK_CLOSE_FRACTION( tightConstraintInverseCovariance( i, j ), 1.0E32, 1.0E-10 );
            }
            else
            {
                BOOST_CHECK_SMALL( tightConstraintInverseCovariance( i, j ) / tightConstraintInverseCovariance( i, i ), 1.0E-10 );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
