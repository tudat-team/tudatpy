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

#include <limits>
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/simulation/estimation_setup/executePlanetaryParameterEstimationTestCase.h"
#include "tudat/simulation/estimation_setup/executeEarthOrbiterParameterEstimationTestCase.h"

#include "tudat/simulation/estimation_setup/podProcessing.h"

namespace tudat
{
namespace unit_tests
{

Eigen::MatrixXd buildBlockDiagonalMatrix( const std::vector< Eigen::MatrixXd >& blockWeights )
{
    const int setSize = static_cast< int >( blockWeights.size( ) ) * blockWeights.at( 0 ).rows( );
    Eigen::MatrixXd blockDiagonalMatrix = Eigen::MatrixXd::Zero( setSize, setSize );
    const int blockSize = blockWeights.at( 0 ).rows( );
    for( unsigned int i = 0; i < blockWeights.size( ); i++ )
    {
        blockDiagonalMatrix.block( static_cast< int >( i ) * blockSize,
                                   static_cast< int >( i ) * blockSize,
                                   blockSize,
                                   blockSize ) = blockWeights.at( i );
    }
    return blockDiagonalMatrix;
}

Eigen::MatrixXd createFullDenseWeightMatrix( const int size, const double baseDiagonal, const double couplingScale )
{
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero( size, size );
    for( int i = 0; i < size; i++ )
    {
        for( int j = 0; j < size; j++ )
        {
            if( i == j )
            {
                continue;
            }
            matrix( i, j ) = couplingScale * ( 1.0 + 0.05 * static_cast< double >( i + j ) );
        }
    }

    // Enforce strict diagonal dominance to guarantee positive definiteness.
    for( int i = 0; i < size; i++ )
    {
        matrix( i, i ) = baseDiagonal + matrix.row( i ).cwiseAbs( ).sum( );
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

//! Test whether the covariance is correctly computed as a function of time
BOOST_AUTO_TEST_CASE( test_CovarianceAsFunctionOfTime )
{
    std::pair< std::shared_ptr< EstimationOutput< double > >, std::shared_ptr< EstimationInput< double, double > > > podData;

    // Simulate covariances directly by propagating to different final tomes
    std::map< int, Eigen::MatrixXd > manualCovarianes;
    for( unsigned int i = 1; i < 5; i++ )
    {
        executeEarthOrbiterParameterEstimation< double, double >( podData, 1.0E7, i, 0, false );
        manualCovarianes[ i ] = podData.first->getUnnormalizedCovarianceMatrix( );
    }

    // Use final calculations to compute covariance as a function of time
    std::map< double, Eigen::MatrixXd > automaticCovariances =
            simulation_setup::calculateCovarianceUsingDataUpToEpoch( podData.second, podData.first, 86400.0 - 1.0 );

    // Check consistency
    int counter = 1;
    for( std::map< double, Eigen::MatrixXd >::const_iterator covarianceIterator = automaticCovariances.begin( );
         covarianceIterator != automaticCovariances.end( );
         covarianceIterator++ )
    {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( covarianceIterator->second, manualCovarianes.at( counter ), 2.0E-8 );
        counter++;
    }
}

BOOST_AUTO_TEST_CASE( test_DesignMatrixSaving )
{
    // Simulate covariances directly by propagating to different final tomes
    for( unsigned int i = 0; i < 2; i++ )
    {
        std::vector< Eigen::MatrixXd > designMatrices;
        std::pair< std::shared_ptr< EstimationOutput< double > >, std::shared_ptr< EstimationInput< double, double > > > podData;
        executeEarthOrbiterParameterEstimation< double, double >( podData, 1.0E7, 1, 0, false, static_cast< bool >( i ) );
        designMatrices.push_back( podData.first->getNormalizedDesignMatrix( ) );
        designMatrices.push_back( podData.first->getUnnormalizedDesignMatrix( ) );
        designMatrices.push_back( podData.first->getNormalizedWeightedDesignMatrix( ) );
        designMatrices.push_back( podData.first->getUnnormalizedWeightedDesignMatrix( ) );

        if( i == 1 )
        {
            int numberOfParameters = podData.first->parameterEstimate_.rows( );
            int numberOfObservations = podData.second->getObservationCollection( )->getTotalObservableSize( );

            for( unsigned int j = 0; j < designMatrices.size( ); j++ )
            {
                BOOST_CHECK_EQUAL( designMatrices.at( j ).rows( ), numberOfObservations );
                BOOST_CHECK_EQUAL( designMatrices.at( j ).cols( ), numberOfParameters );
            }
        }

        if( i == 0 )
        {
            for( unsigned int j = 0; j < designMatrices.size( ); j++ )
            {
                BOOST_CHECK_EQUAL( designMatrices.at( j ).rows( ), 0 );
                BOOST_CHECK_EQUAL( designMatrices.at( j ).cols( ), 0 );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( test_WeightDefinitions )

{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    const double startTime = double( 1.0E7 );
    const int numberOfDaysOfData = 3;

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );

    // Specify initial time
    double initialEphemerisTime = startTime;
    double finalEphemerisTime = initialEphemerisTime + numberOfDaysOfData * 86400.0;

    // Create bodies needed in simulation
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
            "ECLIPJ2000",
            "IAU_Earth",
            spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
            initialEphemerisTime,
            2.0 * mathematical_constants::PI / ( physical_constants::JULIAN_DAY ) );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );

    bodies.at( "Vehicle" )
            ->setEphemeris( std::make_shared< TabulatedCartesianEphemeris<> >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );

    // Creatre ground stations: same position, but different representation
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    groundStationNames.push_back( "Station3" );

    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station3", ( Eigen::Vector3d( ) << 0.0, 0.05, 4.0 ).finished( ), geodetic_position );

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 8, 8 ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Set (perturbed) initial state.
    Eigen::Matrix< double, 6, 1 > systemInitialState =
            convertKeplerianToCartesianElements( asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double, double > >(
                    centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, double( finalEphemerisTime ), cowell );

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< double > >(
            double( initialEphemerisTime ), 40.0, CoefficientSets::rungeKuttaFehlberg78, 40.0, 40.0, 1.0, 1.0 );

    // Define parameters.
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = LinkEndId( "Earth", groundStationNames.at( i ) );
        linkEnds[ receiver ] = LinkEndId( "Vehicle", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = LinkEndId( "Vehicle", "" );
        stationReceiverLinkEnds.push_back( linkEnds );
    }

    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );

    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 2 ] );

    linkEndsPerObservable[ angular_position ].push_back( stationReceiverLinkEnds[ 2 ] );
    linkEndsPerObservable[ angular_position ].push_back( stationTransmitterLinkEnds[ 1 ] );

    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back(
            std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >( "Vehicle", systemInitialState, "Earth" ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies );

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;

    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( );
         linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;

        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            observationSettingsList.push_back(
                    std::make_shared< ObservationModelSettings >( currentObservable, currentLinkEndsList.at( i ) ) );
        }
    }

    // Create orbit determination object.
    OrbitDeterminationManager< double, double > orbitDeterminationManager = OrbitDeterminationManager< double, double >(
            bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );

    std::vector< double > baseTimeList;
    double observationTimeStart = initialEphemerisTime + 1000.0;
    double observationInterval = 20.0;
    unsigned int nbObsPerDay = 50;
    for( int i = 0; i < numberOfDaysOfData; i++ )
    {
        for( unsigned int j = 0; j < nbObsPerDay; j++ )
        {
            baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 86400.0 +
                                    static_cast< double >( j ) * observationInterval );
        }
    }

    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput =
            getObservationSimulationSettings< double >( linkEndsPerObservable, baseTimeList, receiver );

    // Simulate observations
    std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateObservations< double, double >(
            measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Define estimation input
    std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( simulatedObservations );

    std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize =
            simulatedObservations->getObservationTypeStartAndSize( );

    {
        simulatedObservations->setConstantWeight( 0.1 );

        // Define estimation input
        std::shared_ptr< EstimationInput< double, double > > estimationInput =
                std::make_shared< EstimationInput< double, double > >( simulatedObservations );
        std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize =
                simulatedObservations->getObservationTypeStartAndSize( );

        Eigen::VectorXd totalWeights = estimationInput->getWeightsMatrixDiagonals( );

        for( unsigned int i = 0; i < totalWeights.rows( ); i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( totalWeights( i ), 0.1, std::numeric_limits< double >::epsilon( ) );
        }
    }

    {
        std::map< std::shared_ptr< observation_models::ObservationCollectionParser >, double > weightPerObservationParser;
        weightPerObservationParser[ observationParser( one_way_range ) ] = 1.0 / ( 3.0 * 3.0 );
        weightPerObservationParser[ observationParser( angular_position ) ] = 1.0 / ( 1.0E-5 * 1.0E-5 );
        weightPerObservationParser[ observationParser( one_way_doppler ) ] = 1.0 / ( 1.0E-11 * 1.0E-11 * SPEED_OF_LIGHT * SPEED_OF_LIGHT );
        simulatedObservations->setConstantWeightPerObservable( weightPerObservationParser );

        // Define estimation input
        std::shared_ptr< EstimationInput< double, double > > estimationInput =
                std::make_shared< EstimationInput< double, double > >( simulatedObservations );
        std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize =
                simulatedObservations->getObservationTypeStartAndSize( );

        Eigen::VectorXd totalWeights = estimationInput->getWeightsMatrixDiagonals( );

        for( auto it : weightPerObservationParser )
        {
            ObservableType observableType =
                    std::dynamic_pointer_cast< ObservationCollectionObservableTypeParser >( it.first )->getObservableTypes( ).at( 0 );
            for( int i = 0; i < observationTypeStartAndSize.at( observableType ).second; i++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( totalWeights( observationTypeStartAndSize.at( observableType ).first + i ),
                                            it.second,
                                            std::numeric_limits< double >::epsilon( ) );
            }
        }
    }

    {
        Eigen::Vector2d angularPositionWeight;
        angularPositionWeight << 0.1, 0.2;
        simulatedObservations->setConstantWeight( 2.0 );
        simulatedObservations->setConstantWeight( angularPositionWeight, observationParser( angular_position ) );

        // Define estimation input
        std::shared_ptr< EstimationInput< double, double > > estimationInput =
                std::make_shared< EstimationInput< double, double > >( simulatedObservations );
        std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize =
                simulatedObservations->getObservationTypeStartAndSize( );

        Eigen::VectorXd totalWeights = estimationInput->getWeightsMatrixDiagonals( );

        std::pair< int, int > startEndIndex = observationTypeStartAndSize.at( angular_position );

        for( int i = 0; i < startEndIndex.first; i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( totalWeights( i ), 2.0, std::numeric_limits< double >::epsilon( ) );
        }

        for( int i = 0; i < startEndIndex.second; i++ )
        {
            if( i % 2 == 0 )
            {
                BOOST_CHECK_CLOSE_FRACTION(
                        totalWeights( startEndIndex.first + i ), angularPositionWeight( 0 ), std::numeric_limits< double >::epsilon( ) );
            }
            else
            {
                BOOST_CHECK_CLOSE_FRACTION(
                        totalWeights( startEndIndex.first + i ), angularPositionWeight( 1 ), std::numeric_limits< double >::epsilon( ) );
            }
        }

        for( unsigned int i = startEndIndex.first + startEndIndex.second; i < totalWeights.rows( ); i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( totalWeights( i ), 2.0, std::numeric_limits< double >::epsilon( ) );
        }
    }

    // Test tabulated weights
    {
        // Set same tabulated weights to each range observation set
        int sizeRangeObsPerObsSet = nbObsPerDay * numberOfDaysOfData;
        Eigen::VectorXd singleSetRangeWeights =
                Eigen::VectorXd::LinSpaced( sizeRangeObsPerObsSet, 1.0 / ( 3.0 * 3.0 ), 1.0 / ( 4.0 * 4.0 ) );

        // Compute full range weight vector
        unsigned int nbRangeObsSets = simulatedObservations->getSingleObservationSets( observationParser( one_way_range ) ).size( );
        Eigen::VectorXd rangeWeights = Eigen::VectorXd::Zero( nbRangeObsSets * sizeRangeObsPerObsSet );
        for( unsigned int k = 0; k < nbRangeObsSets; k++ )
        {
            rangeWeights.segment( k * sizeRangeObsPerObsSet, sizeRangeObsPerObsSet ) = singleSetRangeWeights;
        }

        // Set total tabulated weights for all Doppler observation sets
        int totalSizeDopplerObs = simulatedObservations->getSingleObservationSets( observationParser( one_way_doppler ) ).size( ) *
                nbObsPerDay * numberOfDaysOfData;
        Eigen::VectorXd dopplerWeights = Eigen::VectorXd::LinSpaced( totalSizeDopplerObs,
                                                                     1.0 / ( 1.0e-11 * SPEED_OF_LIGHT * 1.0e-11 * SPEED_OF_LIGHT ),
                                                                     1.0 / ( 1.5e-11 * SPEED_OF_LIGHT * 1.5e-11 * SPEED_OF_LIGHT ) );

        // Default angular position weights set to 1
        int totalSizeAngularPositionObs = 2.0 *
                simulatedObservations->getSingleObservationSets( observationParser( angular_position ) ).size( ) * nbObsPerDay *
                numberOfDaysOfData;
        Eigen::VectorXd angularPositionWeights = Eigen::VectorXd::Ones( totalSizeAngularPositionObs );

        // Concatenate tabulated weights per observable type (default weihts for angular_position observables)
        std::map< std::shared_ptr< observation_models::ObservationCollectionParser >, Eigen::VectorXd > weightPerObservationParser;
        weightPerObservationParser[ observationParser( one_way_range ) ] = singleSetRangeWeights;
        weightPerObservationParser[ observationParser( one_way_doppler ) ] = dopplerWeights;
        weightPerObservationParser[ observationParser( angular_position ) ] = angularPositionWeights;
        simulatedObservations->setTabulatedWeights( weightPerObservationParser );

        // Define estimation input
        std::shared_ptr< EstimationInput< double, double > > estimationInput =
                std::make_shared< EstimationInput< double, double > >( simulatedObservations );
        std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize =
                simulatedObservations->getObservationTypeStartAndSize( );
        Eigen::VectorXd totalWeights = estimationInput->getWeightsMatrixDiagonals( );

        // Define expected weights per observable
        std::map< ObservableType, Eigen::VectorXd > expectedWeights;
        expectedWeights[ one_way_range ] = rangeWeights;
        expectedWeights[ one_way_doppler ] = dopplerWeights;
        expectedWeights[ angular_position ] = angularPositionWeights;

        for( auto it : expectedWeights )
        {
            for( int i = 0; i < observationTypeStartAndSize.at( it.first ).second; i++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( totalWeights( observationTypeStartAndSize.at( it.first ).first + i ),
                                            it.second( i ),
                                            std::numeric_limits< double >::epsilon( ) );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( test_OffDiagonalWeightsHandling )
{
    using namespace observation_models;
    using namespace propagators;
    using namespace numerical_integrators;
    using namespace simulation_setup;
    using namespace estimatable_parameters;
    using namespace basic_astrodynamics;

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    const double initialEphemerisTime = 1.0E7;
    const double finalEphemerisTime = initialEphemerisTime + 7200.0;
    const int numberOfObservationsPerSet = 3;

    // Define and create bodies.
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
            2.0 * mathematical_constants::PI / ( physical_constants::JULIAN_DAY ) );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );
    bodies.at( "Vehicle" )
            ->setEphemeris( std::make_shared< TabulatedCartesianEphemeris<> >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );

    // Create ground stations.
    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station3", ( Eigen::Vector3d( ) << 0.0, 0.05, 4.0 ).finished( ), geodetic_position );

    // Define dynamics.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    Eigen::Vector6d initialStateInKeplerianElements;
    initialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    initialStateInKeplerianElements( eccentricityIndex ) = 0.03;
    initialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 78.0 );
    initialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 35.0 );
    initialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 15.0 );
    initialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 120.0 );
    const double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::Matrix< double, 6, 1 > truthState =
            convertKeplerianToCartesianElements( initialStateInKeplerianElements, earthGravitationalParameter );

    std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< double, double > >(
                    centralBodies, accelerationModelMap, bodiesToIntegrate, truthState, finalEphemerisTime, cowell );
    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
            std::make_shared< IntegratorSettings< double > >( rungeKutta4, initialEphemerisTime, 20.0 );

    // Define link ends (multiple per observable type).
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
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies );

    OrbitDeterminationManager< double, double > orbitDeterminationManager(
            bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );

    // Simulate observations.
    std::vector< double > observationTimes;
    for( int i = 0; i < numberOfObservationsPerSet; i++ )
    {
        observationTimes.push_back( initialEphemerisTime + 600.0 + i * 120.0 );
    }

    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
            one_way_range, rangeLinkEndsStation1, observationTimes, receiver ) );
    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
            one_way_range, rangeLinkEndsStation2, observationTimes, receiver ) );
    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
            angular_position, angularLinkEndsStation1, observationTimes, receiver ) );
    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
            angular_position, angularLinkEndsStation2, observationTimes, receiver ) );
    measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
            angular_position, angularLinkEndsStation3, observationTimes, receiver ) );

    std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateObservations< double, double >(
            measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    simulatedObservations->setConstantWeight( 1.0 );
    std::vector< std::shared_ptr< SingleObservationSet< double, double > > > singleObservationSets =
            simulatedObservations->getSingleObservationSets( );
    BOOST_CHECK_EQUAL( singleObservationSets.size( ), 5 );

    // We store the exact expected matrix per single set, then assemble the full matrix from these blocks.
    std::map< SingleObservationSet< double, double >*, Eigen::MatrixXd > expectedSetWeightMatrices;
    int numberOfBlockDiagonalAngularSets = 0;
    int numberOfFullAngularSets = 0;
    int numberOfFullRangeSets = 0;

    for( auto singleObservationSet : singleObservationSets )
    {
        const ObservableType observableType = singleObservationSet->getObservableType( );
        const LinkEnds currentLinkEnds = singleObservationSet->getLinkEnds( ).linkEnds_;
        const std::string stationName = currentLinkEnds.at( receiver ).stationName_;
        const int setSize = static_cast< int >( singleObservationSet->getTotalObservationSetSize( ) );
        const int singleObservationSize = static_cast< int >( singleObservationSet->getSingleObservableSize( ) );
        Eigen::MatrixXd baseWeightMatrix = Eigen::MatrixXd::Zero( setSize, setSize );
        Eigen::MatrixXd fullWeightMatrixContribution = Eigen::MatrixXd::Zero( 0, 0 );
        bool fullWeightMatrixContributionIsSet = false;

        if( observableType == angular_position )
        {
            // Two angular sets get block-diagonal base weights plus a full contribution.
            if( stationName == "Station1" || stationName == "Station2" )
            {
                std::vector< Eigen::MatrixXd > blockWeights;
                blockWeights.reserve( static_cast< unsigned int >( numberOfObservationsPerSet ) );
                const double offset = ( stationName == "Station1" ? 0.0 : 1.0 );
                for( int i = 0; i < numberOfObservationsPerSet; i++ )
                {
                    Eigen::MatrixXd currentBlock = Eigen::MatrixXd::Zero( singleObservationSize, singleObservationSize );
                    currentBlock( 0, 0 ) = 8.0 + offset + 0.2 * static_cast< double >( i );
                    currentBlock( 1, 1 ) = 6.0 + offset + 0.15 * static_cast< double >( i );
                    currentBlock( 0, 1 ) = 0.2 + 0.02 * static_cast< double >( i );
                    currentBlock( 1, 0 ) = currentBlock( 0, 1 );
                    blockWeights.push_back( currentBlock );
                }
                singleObservationSet->setBlockDiagonalWeights( blockWeights );
                baseWeightMatrix = buildBlockDiagonalMatrix( blockWeights );
                numberOfBlockDiagonalAngularSets++;

                fullWeightMatrixContribution = createFullDenseWeightMatrix(
                        setSize, ( stationName == "Station1" ? 0.45 : 0.55 ), ( stationName == "Station1" ? 0.03 : 0.025 ) );
                singleObservationSet->setFullWeightMatrix( fullWeightMatrixContribution );
                fullWeightMatrixContributionIsSet = true;
                numberOfFullAngularSets++;
            }
            else if( stationName == "Station3" )
            {
                // One angular set remains pure diagonal to ensure not all angular sets are block-diagonal.
                Eigen::VectorXd diagonalWeights = Eigen::VectorXd::Zero( setSize );
                diagonalWeights << 12.0, 13.0, 14.0, 15.0, 16.0, 17.0;
                singleObservationSet->setTabulatedWeights( diagonalWeights );
                baseWeightMatrix.diagonal( ) = diagonalWeights;
            }
            else
            {
                BOOST_FAIL( "Unexpected station for angular position observable." );
            }
        }
        else if( observableType == one_way_range )
        {
            // Both range sets use diagonal base weights plus a full contribution.
            Eigen::VectorXd diagonalWeights = Eigen::VectorXd::Zero( setSize );
            if( stationName == "Station1" )
            {
                diagonalWeights << 5.0, 6.0, 7.0;
                fullWeightMatrixContribution = createFullDenseWeightMatrix( setSize, 0.35, 0.02 );
            }
            else if( stationName == "Station2" )
            {
                diagonalWeights << 6.5, 7.5, 8.5;
                fullWeightMatrixContribution = createFullDenseWeightMatrix( setSize, 0.4, 0.015 );
            }
            else
            {
                BOOST_FAIL( "Unexpected station for one-way range observable." );
            }
            singleObservationSet->setTabulatedWeights( diagonalWeights );
            baseWeightMatrix.diagonal( ) = diagonalWeights;
            singleObservationSet->setFullWeightMatrix( fullWeightMatrixContribution );
            fullWeightMatrixContributionIsSet = true;
            numberOfFullRangeSets++;
        }
        else
        {
            BOOST_FAIL( "Unexpected observable type in off-diagonal weights test." );
        }

        // Single-set expected matrix: base contribution + optional full contribution.
        Eigen::MatrixXd expectedCombinedWeightMatrix = baseWeightMatrix;
        if( fullWeightMatrixContributionIsSet )
        {
            expectedCombinedWeightMatrix += fullWeightMatrixContribution;
        }
        expectedSetWeightMatrices[ singleObservationSet.get( ) ] = expectedCombinedWeightMatrix;

        // Validate that each set returns exactly the weights we just assigned.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( singleObservationSet->getWeightMatrix( ), expectedCombinedWeightMatrix, 1.0E-15 );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( singleObservationSet->getWeightsDiagonalVector( ),
                                           expectedCombinedWeightMatrix.diagonal( ),
                                           1.0E-15 );
        if( fullWeightMatrixContributionIsSet )
        {
            BOOST_CHECK( singleObservationSet->hasFullWeightMatrixContribution( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    singleObservationSet->getFullWeightMatrix( ), fullWeightMatrixContribution, 1.0E-15 );
        }
    }

    const int numberOfAngularSets = static_cast< int >(
            simulatedObservations->getSingleObservationSets( observationParser( angular_position ) ).size( ) );
    BOOST_CHECK_GE( numberOfBlockDiagonalAngularSets, 2 );
    BOOST_CHECK_LT( numberOfBlockDiagonalAngularSets, numberOfAngularSets );
    BOOST_CHECK_GE( numberOfFullAngularSets, 2 );
    BOOST_CHECK_GE( numberOfFullRangeSets, 2 );
    BOOST_CHECK( simulatedObservations->hasOffDiagonalWeights( ) );

    // Build the expected dense full weights matrix by placing each set matrix on the global diagonal.
    const int totalObservationSize = simulatedObservations->getTotalObservableSize( );
    Eigen::MatrixXd expectedFullWeightsMatrix = Eigen::MatrixXd::Zero( totalObservationSize, totalObservationSize );
    int currentStartIndex = 0;
    for( auto singleObservationSet : simulatedObservations->getSingleObservationSets( ) )
    {
        const int currentSetSize = static_cast< int >( singleObservationSet->getTotalObservationSetSize( ) );
        expectedFullWeightsMatrix.block( currentStartIndex, currentStartIndex, currentSetSize, currentSetSize ) =
                expectedSetWeightMatrices.at( singleObservationSet.get( ) );
        currentStartIndex += currentSetSize;
    }
    BOOST_CHECK_EQUAL( currentStartIndex, totalObservationSize );

    // Global matrix from the ObservationCollection must match the exact matrix we assembled manually.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            Eigen::MatrixXd( simulatedObservations->getConcatenatedWeightMatrix( ) ), expectedFullWeightsMatrix, 1.0E-15 );

    // EstimationInput should pass through the exact same matrix and diagonal.
    std::shared_ptr< EstimationInput< double, double > > estimationInput =
            std::make_shared< EstimationInput< double, double > >( simulatedObservations );
    estimationInput->defineEstimationSettings( true, true, true, false, true, false );
    estimationInput->setConvergenceChecker( estimationConvergenceChecker( 8, 0.0, 1.0E-20, 4 ) );
    BOOST_CHECK( estimationInput->hasOffDiagonalWeights( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            Eigen::MatrixXd( estimationInput->getWeightsMatrix( ) ), expectedFullWeightsMatrix, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            estimationInput->getWeightsMatrixDiagonals( ), expectedFullWeightsMatrix.diagonal( ), 1.0E-15 );

    // Perturb the initial state estimate.
    Eigen::VectorXd truthStateVector = parametersToEstimate->template getFullParameterValues< double >( );
    Eigen::VectorXd perturbedState = truthStateVector;
    Eigen::VectorXd statePerturbation = Eigen::VectorXd::Zero( 6 );
    statePerturbation << 50.0, -40.0, 30.0, 0.05, -0.04, 0.03;
    perturbedState += statePerturbation;
    parametersToEstimate->resetParameterValues( perturbedState );

    // Perform estimation and test convergence.
    std::shared_ptr< EstimationOutput< double, double > > estimationOutput =
            orbitDeterminationManager.estimateParameters( estimationInput );
    Eigen::VectorXd estimatedState = estimationOutput->parameterEstimate_;
    Eigen::VectorXd estimationError = estimatedState - truthStateVector;
    BOOST_CHECK_SMALL( estimationError.segment( 0, 3 ).norm( ), 2.5 );
    BOOST_CHECK_SMALL( estimationError.segment( 3, 3 ).norm( ), 5.0E-3 );

    // For estimation output: inverse covariance should equal H^T W H with our exact dense W.
    Eigen::MatrixXd designMatrixFromEstimation = estimationOutput->getUnnormalizedDesignMatrix( );
    Eigen::MatrixXd expectedInverseCovarianceFromEstimation =
            designMatrixFromEstimation.transpose( ) * expectedFullWeightsMatrix * designMatrixFromEstimation;
    Eigen::MatrixXd expectedCovarianceFromEstimation = expectedInverseCovarianceFromEstimation.inverse( );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            estimationOutput->getUnnormalizedInverseCovarianceMatrix( ), expectedInverseCovarianceFromEstimation, 1.0E-14 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            estimationOutput->getUnnormalizedCovarianceMatrix( ), expectedCovarianceFromEstimation, 5.0E-8 );

    // Repeat the same check through the standalone covariance path.
    parametersToEstimate->template resetParameterValues< double >( estimationOutput->parameterEstimate_ );
    std::shared_ptr< CovarianceAnalysisInput< double, double > > covarianceInput =
            std::make_shared< CovarianceAnalysisInput< double, double > >( simulatedObservations );
    covarianceInput->defineCovarianceSettings( true, true, true, false );
    BOOST_CHECK( covarianceInput->hasOffDiagonalWeights( ) );

    std::shared_ptr< CovarianceAnalysisOutput< double, double > > covarianceOutput =
            orbitDeterminationManager.computeCovariance( covarianceInput );
    Eigen::MatrixXd designMatrixFromCovariance = covarianceOutput->getUnnormalizedDesignMatrix( );
    Eigen::MatrixXd expectedInverseCovarianceFromCovariance =
            designMatrixFromCovariance.transpose( ) * expectedFullWeightsMatrix * designMatrixFromCovariance;
    Eigen::MatrixXd expectedCovarianceFromCovariance = expectedInverseCovarianceFromCovariance.inverse( );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            covarianceOutput->getUnnormalizedInverseCovarianceMatrix( ), expectedInverseCovarianceFromCovariance, 1.0E-14 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            covarianceOutput->getUnnormalizedCovarianceMatrix( ), expectedCovarianceFromCovariance, 5.0E-8 );

    // Ensure estimate and standalone covariance are consistent.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( covarianceOutput->getUnnormalizedInverseCovarianceMatrix( ),
                                       estimationOutput->getUnnormalizedInverseCovarianceMatrix( ),
                                       1.0E-14 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
