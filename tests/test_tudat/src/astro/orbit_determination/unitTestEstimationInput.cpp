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

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/executePlanetaryParameterEstimationTestCase.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"

namespace tudat
{
namespace unit_tests
{

Eigen::MatrixXd createStrictlyDominantWeightMatrix( const int size, const double baseDiagonal, const double couplingScale )
{
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero( size, size );
    for( int row = 0; row < size; ++row )
    {
        for( int column = 0; column < size; ++column )
        {
            if( row != column )
            {
                matrix( row, column ) = couplingScale * ( 1.0 + 0.05 * static_cast< double >( row + column ) );
            }
        }
    }

    for( int row = 0; row < size; ++row )
    {
        matrix( row, row ) = baseDiagonal + matrix.row( row ).cwiseAbs( ).sum( );
    }
    return matrix;
}

BOOST_AUTO_TEST_SUITE( test_estimation_input_output )

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

/*!
 * Verifies off-diagonal weights through estimation and covariance analysis.
 *
 * Test outline: simulates range and angular observations inserted out of
 * observable/link-end order, assigns per-observation blocks, set-level blocks
 * and cross-set blocks through ObservationDataset, and checks the assembled
 * sparse weight matrix. It then verifies that these weights are used in the
 * parameter update, final estimation inverse covariance and covariance-analysis
 * inverse covariance.
 */
BOOST_AUTO_TEST_CASE( test_OffDiagonalWeightsInEstimationAndCovariance )
{
    using namespace observation_models;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace estimatable_parameters;
    using namespace basic_astrodynamics;

    spice_interface::loadStandardSpiceKernels( );

    const double initialEphemerisTime = 1.0E7;
    const double finalEphemerisTime = initialEphemerisTime + 7200.0;
    const int numberOfObservationsPerSet = 3;

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
            "ECLIPJ2000",
            "IAU_Earth",
            spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
            initialEphemerisTime,
            2.0 * mathematical_constants::PI / physical_constants::JULIAN_DAY );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );
    bodies.at( "Vehicle" )
            ->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris<> >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );

    createGroundStation( bodies.at( "Earth" ),
                         "Station1",
                         ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ),
                         coordinate_conversions::geodetic_position );
    createGroundStation( bodies.at( "Earth" ),
                         "Station2",
                         ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ),
                         coordinate_conversions::geodetic_position );
    createGroundStation( bodies.at( "Earth" ),
                         "Station3",
                         ( Eigen::Vector3d( ) << 0.0, 0.05, 4.0 ).finished( ),
                         coordinate_conversions::geodetic_position );

    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Vehicle" );
    std::vector< std::string > centralBodies;
    centralBodies.push_back( "Earth" );
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    Eigen::Vector6d initialStateInKeplerianElements;
    initialStateInKeplerianElements( orbital_element_conversions::semiMajorAxisIndex ) = 7200.0E3;
    initialStateInKeplerianElements( orbital_element_conversions::eccentricityIndex ) = 0.03;
    initialStateInKeplerianElements( orbital_element_conversions::inclinationIndex ) = unit_conversions::convertDegreesToRadians( 78.0 );
    initialStateInKeplerianElements( orbital_element_conversions::argumentOfPeriapsisIndex ) =
            unit_conversions::convertDegreesToRadians( 35.0 );
    initialStateInKeplerianElements( orbital_element_conversions::longitudeOfAscendingNodeIndex ) =
            unit_conversions::convertDegreesToRadians( 15.0 );
    initialStateInKeplerianElements( orbital_element_conversions::trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 120.0 );
    const double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d truthState = orbital_element_conversions::convertKeplerianToCartesianElements( initialStateInKeplerianElements,
                                                                                                         earthGravitationalParameter );

    std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double, double > >(
                    centralBodies, accelerationModelMap, bodiesToIntegrate, truthState, finalEphemerisTime, cowell );
    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            std::make_shared< IntegratorSettings< double > >( rungeKutta4, initialEphemerisTime, 20.0 );

    LinkEnds rangeLinkEndsStation1, rangeLinkEndsStation2, angularLinkEndsStation1, angularLinkEndsStation2, angularLinkEndsStation3;
    rangeLinkEndsStation1[ transmitter ] = LinkEndId( "Vehicle", "" );
    rangeLinkEndsStation1[ receiver ] = LinkEndId( "Earth", "Station1" );
    rangeLinkEndsStation2[ transmitter ] = LinkEndId( "Vehicle", "" );
    rangeLinkEndsStation2[ receiver ] = LinkEndId( "Earth", "Station2" );
    angularLinkEndsStation1 = rangeLinkEndsStation1;
    angularLinkEndsStation2 = rangeLinkEndsStation2;
    angularLinkEndsStation3[ transmitter ] = LinkEndId( "Vehicle", "" );
    angularLinkEndsStation3[ receiver ] = LinkEndId( "Earth", "Station3" );

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( one_way_range, rangeLinkEndsStation1 ) );
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( one_way_range, rangeLinkEndsStation2 ) );
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( angular_position, angularLinkEndsStation1 ) );
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( angular_position, angularLinkEndsStation2 ) );
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( angular_position, angularLinkEndsStation3 ) );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back(
            std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >( "Vehicle", truthState, "Earth" ) );
    std::shared_ptr< EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies );

    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );

    std::vector< double > observationTimes;
    for( int i = 0; i < numberOfObservationsPerSet; ++i )
    {
        observationTimes.push_back( initialEphemerisTime + 600.0 + i * 120.0 );
    }

    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
            angular_position, angularLinkEndsStation3, observationTimes, receiver ) );
    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
            one_way_range, rangeLinkEndsStation2, observationTimes, receiver ) );
    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
            angular_position, angularLinkEndsStation1, observationTimes, receiver ) );
    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
            one_way_range, rangeLinkEndsStation1, observationTimes, receiver ) );
    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
            angular_position, angularLinkEndsStation2, observationTimes, receiver ) );

    std::shared_ptr< ObservationDataset< double, double > > simulatedObservations = simulateObservationDataset< double, double >(
            measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    const std::vector< unsigned int > orderedSetIds = simulatedObservations->getSetIdsInOrderedFlattenedDataOrder( );

    // The simulated dataset is intentionally inserted out of observable/link-end order.
    BOOST_REQUIRE_EQUAL( orderedSetIds.size( ), simulatedObservations->getNumberOfObservationSets( ) );
    BOOST_CHECK( orderedSetIds.front( ) != 0 );

    std::map< unsigned int, Eigen::MatrixXd > expectedSetWeightMatrices;
    int numberOfObservationBlockSets = 0;
    int numberOfSetBlockSets = 0;
    for( unsigned int setId = 0; setId < simulatedObservations->getNumberOfObservationSets( ); ++setId )
    {
        const ObservationSetMetadata< double, double >& metadata = simulatedObservations->getObservationSetMetadata( setId );
        const LinkEnds& currentLinkEnds = simulatedObservations->getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_;
        const std::string stationName = currentLinkEnds.at( receiver ).stationName_;
        const int setSize = static_cast< int >( simulatedObservations->getTotalScalarSizeForSet( setId ) );

        if( metadata.observableType_ == angular_position && ( stationName == "Station1" || stationName == "Station2" ) )
        {
            Eigen::MatrixXd expectedSetMatrix = Eigen::MatrixXd::Zero( setSize, setSize );
            const double stationOffset = ( stationName == "Station1" ? 0.0 : 1.0 );
            const std::vector< unsigned int >& observationIds = simulatedObservations->getObservationIdsForSet( setId );
            for( unsigned int i = 0; i < observationIds.size( ); ++i )
            {
                Eigen::Matrix2d currentBlock;
                currentBlock << 8.0 + stationOffset + 0.2 * static_cast< double >( i ), 0.2 + 0.02 * static_cast< double >( i ),
                        0.2 + 0.02 * static_cast< double >( i ), 6.0 + stationOffset + 0.15 * static_cast< double >( i );
                simulatedObservations->setWeightMatrixForObservation( observationIds.at( i ), currentBlock );
                expectedSetMatrix.block( 2 * static_cast< int >( i ), 2 * static_cast< int >( i ), 2, 2 ) = currentBlock;
            }
            expectedSetWeightMatrices[ setId ] = expectedSetMatrix;
            ++numberOfObservationBlockSets;
        }
        else
        {
            const double baseDiagonal = metadata.observableType_ == one_way_range ? 4.0 : 7.0;
            const double couplingScale = metadata.observableType_ == one_way_range ? 0.02 : 0.01;
            const Eigen::MatrixXd setWeightMatrix = createStrictlyDominantWeightMatrix( setSize, baseDiagonal, couplingScale );
            simulatedObservations->setWeightMatrixForSet( setId, setWeightMatrix );
            expectedSetWeightMatrices[ setId ] = setWeightMatrix;
            ++numberOfSetBlockSets;
        }
    }

    // The simulated dataset must exercise both per-observation blocks and full set-level blocks.
    BOOST_CHECK_GE( numberOfObservationBlockSets, 2 );
    BOOST_CHECK_GE( numberOfSetBlockSets, 2 );

    const int totalObservationSize = static_cast< int >( simulatedObservations->getTotalScalarSize( ) );
    Eigen::MatrixXd expectedFullWeightsMatrix = Eigen::MatrixXd::Zero( totalObservationSize, totalObservationSize );
    int currentStartIndex = 0;
    std::map< unsigned int, int > orderedSetStartIndex;
    for( const unsigned int setId : orderedSetIds )
    {
        const int currentSetSize = static_cast< int >( simulatedObservations->getTotalScalarSizeForSet( setId ) );
        orderedSetStartIndex[ setId ] = currentStartIndex;
        expectedFullWeightsMatrix.block( currentStartIndex, currentStartIndex, currentSetSize, currentSetSize ) =
                expectedSetWeightMatrices.at( setId );
        currentStartIndex += currentSetSize;
    }

    // The reference dense matrix must span the complete ordered flattened data vector.
    BOOST_CHECK_EQUAL( currentStartIndex, totalObservationSize );

    std::vector< unsigned int > rangeSetIds;
    for( unsigned int setId = 0; setId < simulatedObservations->getNumberOfObservationSets( ); ++setId )
    {
        if( simulatedObservations->getObservationSetMetadata( setId ).observableType_ == one_way_range )
        {
            rangeSetIds.push_back( setId );
        }
    }

    // The test requires at least two range sets so the arbitrary block crosses unrelated observation sets.
    BOOST_REQUIRE_GE( rangeSetIds.size( ), 2 );
    const std::vector< unsigned int > rowBlockObservationIds = {
        simulatedObservations->getObservationIdsForSet( rangeSetIds.at( 0 ) ).at( 0 ),
        simulatedObservations->getObservationIdsForSet( rangeSetIds.at( 0 ) ).at( 2 )
    };
    const std::vector< unsigned int > columnBlockObservationIds = {
        simulatedObservations->getObservationIdsForSet( rangeSetIds.at( 1 ) ).at( 1 ),
        simulatedObservations->getObservationIdsForSet( rangeSetIds.at( 1 ) ).at( 2 )
    };
    Eigen::Matrix2d crossSetWeightBlock;
    crossSetWeightBlock << 0.35, 0.04, 0.06, 0.31;
    simulatedObservations->setWeightBlock( rowBlockObservationIds, columnBlockObservationIds, crossSetWeightBlock, {}, {}, true );

    auto getOrderedFlattenedDataIndex = [ &simulatedObservations, &orderedSetStartIndex ]( const unsigned int observationId,
                                                                                           const unsigned int componentIndex ) {
        const ObservationDatasetRow< double >& row = simulatedObservations->getObservationRow( observationId );
        return orderedSetStartIndex.at( row.setId_ ) + static_cast< int >( row.indexInSet_ * row.scalarSize_ + componentIndex );
    };
    for( unsigned int i = 0; i < rowBlockObservationIds.size( ); ++i )
    {
        for( unsigned int j = 0; j < columnBlockObservationIds.size( ); ++j )
        {
            const int rowIndex = getOrderedFlattenedDataIndex( rowBlockObservationIds.at( i ), 0 );
            const int columnIndex = getOrderedFlattenedDataIndex( columnBlockObservationIds.at( j ), 0 );
            expectedFullWeightsMatrix( rowIndex, columnIndex ) = crossSetWeightBlock( i, j );
            expectedFullWeightsMatrix( columnIndex, rowIndex ) = crossSetWeightBlock( i, j );
        }
    }

    const FlattenedObservationData< double, double > weightData = simulatedObservations->createOrderedFlattenedObservationData( );

    // Dataset flattened data must contain the exact sparse off-diagonal matrix and expose its diagonal as the compact vector.
    BOOST_CHECK( weightData.hasOffDiagonalWeights( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( weightData.getWeightMatrix( ).toDense( ), expectedFullWeightsMatrix, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( weightData.getWeightVector( ), expectedFullWeightsMatrix.diagonal( ), 1.0E-15 );

    const Eigen::VectorXd truthStateVector = parametersToEstimate->template getFullParameterValues< double >( );
    Eigen::VectorXd perturbedState = truthStateVector;
    Eigen::VectorXd statePerturbation( 6 );
    statePerturbation << 50.0, -40.0, 30.0, 0.05, -0.04, 0.03;
    perturbedState += statePerturbation;

    parametersToEstimate->resetParameterValues( perturbedState );
    std::shared_ptr< EstimationInput< double, double > > singleStepEstimationInput =
            std::make_shared< EstimationInput< double, double > >( simulatedObservations );
    singleStepEstimationInput->defineEstimationSettings( true, true, true, false, true, false );
    singleStepEstimationInput->setConvergenceChecker( estimationConvergenceChecker( 1, 0.0, 1.0E-20, 4 ) );
    std::shared_ptr< EstimationOutput< double, double > > singleStepEstimationOutput =
            orbitDeterminationManager.estimateParameters( singleStepEstimationInput );

    // A one-iteration run must expose the initial and updated parameters plus the residuals used for the step.
    BOOST_REQUIRE_EQUAL( singleStepEstimationOutput->parameterHistory_.size( ), 2 );
    BOOST_REQUIRE_EQUAL( singleStepEstimationOutput->residualHistory_.size( ), 1 );
    const Eigen::MatrixXd singleStepDesignMatrix = singleStepEstimationOutput->getNormalizedDesignMatrix( );
    const Eigen::VectorXd singleStepResiduals = singleStepEstimationOutput->residualHistory_.at( 0 );
    const Eigen::VectorXd expectedSingleStepParameterUpdate =
            linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix(
                    singleStepDesignMatrix, singleStepResiduals, weightData.getWeightMatrix( ) )
                    .first.cwiseQuotient( singleStepEstimationOutput->getNormalizationTerms( ) );
    const Eigen::VectorXd actualSingleStepParameterUpdate =
            singleStepEstimationOutput->parameterHistory_.at( 1 ) - singleStepEstimationOutput->parameterHistory_.at( 0 );

    // The first estimation step must use the sparse weight matrix in the normal-equation solve.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( actualSingleStepParameterUpdate, expectedSingleStepParameterUpdate, 1.0E-10 );

    parametersToEstimate->resetParameterValues( perturbedState );
    std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( simulatedObservations );
    estimationInput->defineEstimationSettings( true, true, true, false, true, false );
    estimationInput->setConvergenceChecker( estimationConvergenceChecker( 8, 0.0, 1.0E-20, 4 ) );

    std::shared_ptr< EstimationOutput< double, double > > estimationOutput =
            orbitDeterminationManager.estimateParameters( estimationInput );
    const Eigen::VectorXd estimationError = estimationOutput->parameterEstimate_ - truthStateVector;

    // A multi-iteration estimation with off-diagonal weights must still converge close to the truth state.
    BOOST_CHECK_SMALL( estimationError.segment( 0, 3 ).norm( ), 2.5 );
    BOOST_CHECK_SMALL( estimationError.segment( 3, 3 ).norm( ), 5.0E-3 );

    const Eigen::MatrixXd designMatrixFromEstimation = estimationOutput->getUnnormalizedDesignMatrix( );
    const Eigen::MatrixXd expectedInverseCovarianceFromEstimation =
            designMatrixFromEstimation.transpose( ) * expectedFullWeightsMatrix * designMatrixFromEstimation;

    // The final estimation inverse covariance must be H^T W H using the full off-diagonal weight matrix.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            estimationOutput->getUnnormalizedInverseCovarianceMatrix( ), expectedInverseCovarianceFromEstimation, 1.0E-13 );

    parametersToEstimate->template resetParameterValues< double >( estimationOutput->parameterEstimate_ );
    std::shared_ptr< CovarianceAnalysisInput< double, double > > covarianceInput =
            std::make_shared< CovarianceAnalysisInput< double, double > >( simulatedObservations );
    covarianceInput->defineCovarianceSettings( true, true, true, false );

    std::shared_ptr< CovarianceAnalysisOutput< double, double > > covarianceOutput =
            orbitDeterminationManager.computeCovariance( covarianceInput );
    const Eigen::MatrixXd designMatrixFromCovariance = covarianceOutput->getUnnormalizedDesignMatrix( );
    const Eigen::MatrixXd expectedInverseCovarianceFromCovariance =
            designMatrixFromCovariance.transpose( ) * expectedFullWeightsMatrix * designMatrixFromCovariance;

    // A standalone covariance analysis must use the same full weight matrix and match the estimation covariance.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            covarianceOutput->getUnnormalizedInverseCovarianceMatrix( ), expectedInverseCovarianceFromCovariance, 1.0E-13 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( covarianceOutput->getUnnormalizedInverseCovarianceMatrix( ),
                                       estimationOutput->getUnnormalizedInverseCovarianceMatrix( ),
                                       1.0E-12 );

    const unsigned int rejectedObservationId = simulatedObservations->getObservationIdsForSet( rangeSetIds.at( 0 ) ).at( 1 );
    const unsigned int rejectedObservationSize = simulatedObservations->getObservationRow( rejectedObservationId ).scalarSize_;
    const ObservationSelectionCondition< double, double > rejectedObservationSelectionCondition(
            [ rejectedObservationId ]( const ObservationDataset< double, double >&, const unsigned int observationId ) {
                return observationId == rejectedObservationId;
            } );
    simulatedObservations->rejectObservations( rejectedObservationSelectionCondition, "excluded from estimation system" );
    simulatedObservations->setResidualVector(
            Eigen::VectorXd::Constant( static_cast< int >( simulatedObservations->getTotalScalarSize( ) ), -12345.0 ) );
    const FlattenedObservationData< double, double > activeData = simulatedObservations->createOrderedFlattenedObservationData( false );

    // Rejecting one observation must remove only its scalar rows from the estimator-facing flattened data.
    BOOST_CHECK_EQUAL( activeData.getObservationVector( ).size( ), totalObservationSize - rejectedObservationSize );

    parametersToEstimate->resetParameterValues( perturbedState );
    std::shared_ptr< EstimationInput< double, double > > rejectedEstimationInput =
            std::make_shared< EstimationInput< double, double > >( simulatedObservations );
    rejectedEstimationInput->defineEstimationSettings( true, true, true, false, true, false );
    rejectedEstimationInput->setConvergenceChecker( estimationConvergenceChecker( 1, 0.0, 1.0E-20, 4 ) );
    std::shared_ptr< EstimationOutput< double, double > > rejectedEstimationOutput =
            orbitDeterminationManager.estimateParameters( rejectedEstimationInput );

    // Differential correction must use only active rows in its residual history and saved design matrix.
    BOOST_REQUIRE_EQUAL( rejectedEstimationOutput->residualHistory_.size( ), 1 );
    BOOST_CHECK_EQUAL( rejectedEstimationOutput->residualHistory_.at( 0 ).size( ), activeData.getObservationVector( ).size( ) );
    BOOST_CHECK_EQUAL( rejectedEstimationOutput->getUnnormalizedDesignMatrix( ).rows( ), activeData.getObservationVector( ).size( ) );

    // The rejected observation remains stored and must still receive an updated residual during estimation.
    BOOST_CHECK( !simulatedObservations->getObservationRow( rejectedObservationId ).isActive_ );
    BOOST_CHECK_GT( std::fabs( simulatedObservations->getResidualValue( rejectedObservationId )( 0 ) + 12345.0 ), 1.0 );

    const Eigen::MatrixXd rejectedSingleStepDesignMatrix = rejectedEstimationOutput->getNormalizedDesignMatrix( );
    const Eigen::VectorXd rejectedSingleStepResiduals = rejectedEstimationOutput->residualHistory_.at( 0 );
    const Eigen::VectorXd expectedRejectedStepParameterUpdate =
            linear_algebra::performLeastSquaresAdjustmentFromDesignMatrix(
                    rejectedSingleStepDesignMatrix, rejectedSingleStepResiduals, activeData.getSparseWeightMatrix( ) )
                    .first.cwiseQuotient( rejectedEstimationOutput->getNormalizationTerms( ) );
    const Eigen::VectorXd actualRejectedStepParameterUpdate =
            rejectedEstimationOutput->parameterHistory_.at( 1 ) - rejectedEstimationOutput->parameterHistory_.at( 0 );

    // The rejected row must not contribute to the differential-correction step.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( actualRejectedStepParameterUpdate, expectedRejectedStepParameterUpdate, 1.0E-10 );

    parametersToEstimate->template resetParameterValues< double >( rejectedEstimationOutput->parameterEstimate_ );
    std::shared_ptr< CovarianceAnalysisInput< double, double > > rejectedCovarianceInput =
            std::make_shared< CovarianceAnalysisInput< double, double > >( simulatedObservations );
    rejectedCovarianceInput->defineCovarianceSettings( true, true, true, false );
    std::shared_ptr< CovarianceAnalysisOutput< double, double > > rejectedCovarianceOutput =
            orbitDeterminationManager.computeCovariance( rejectedCovarianceInput );

    // Covariance analysis must use the same active-only flattened data as differential correction.
    BOOST_CHECK_EQUAL( rejectedCovarianceOutput->getUnnormalizedDesignMatrix( ).rows( ), activeData.getObservationVector( ).size( ) );
    const Eigen::MatrixXd expectedRejectedInverseCovariance = rejectedCovarianceOutput->getUnnormalizedDesignMatrix( ).transpose( ) *
            activeData.getSparseWeightMatrix( ).toDense( ) * rejectedCovarianceOutput->getUnnormalizedDesignMatrix( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            rejectedCovarianceOutput->getUnnormalizedInverseCovarianceMatrix( ), expectedRejectedInverseCovariance, 1.0E-13 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
