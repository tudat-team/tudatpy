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

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/simulation/estimation_setup/orbitDeterminationTestCases.h"
#include "tudat/simulation/estimation_setup/podProcessing.h"


namespace tudat
{
namespace unit_tests
{
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
            executePlanetaryParameterEstimation< double, double >(
                simulationType, parameterPerturbation, inverseAPrioriCovariance );

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
                simulationType, parameterPerturbation,  moderateInverseAPriopriCovariance );

    // Run estimation without a priori covariance
    std::pair< std::shared_ptr< EstimationOutput< double > >, Eigen::VectorXd > estimationOutputWithoutAprioriCovariance =
            executePlanetaryParameterEstimation< double, double >(
                simulationType, parameterPerturbation );

    // Run estimation without a priori covariance and increased weights
    double constantWeight = 100.0;
    std::pair< std::shared_ptr< EstimationOutput< double > >, Eigen::VectorXd > estimationOutputWithoutAprioriCovarianceAndWeakWeight =
            executePlanetaryParameterEstimation< double, double >(
                simulationType, parameterPerturbation, Eigen::MatrixXd::Zero( 7, 7 ), constantWeight);

    // Retrieve estimation errors and a priori covariances
    Eigen::MatrixXd tightConstraintInverseCovariance  =
            estimationOutputWithAprioriCovariance.first->getUnnormalizedInverseCovarianceMatrix( );
    Eigen::MatrixXd weakConstraintInverseCovariance  =
            estimationOutputWithSmallAprioriCovariance.first->getUnnormalizedInverseCovarianceMatrix( );
    Eigen::MatrixXd moderateConstraintInverseCovariance  =
            estimationOutputWithModerateAprioriCovariance.first->getUnnormalizedInverseCovarianceMatrix( );
    Eigen::MatrixXd noConstraintInverseCovariance  =
            estimationOutputWithoutAprioriCovariance.first->getUnnormalizedInverseCovarianceMatrix( );
    Eigen::MatrixXd noConstraintInverseCovarianceWithWeakWeight  =
            estimationOutputWithoutAprioriCovarianceAndWeakWeight.first->getUnnormalizedInverseCovarianceMatrix( );

    Eigen::VectorXd tightConstraintError  =
            estimationOutputWithAprioriCovariance.second;
    Eigen::VectorXd weakConstraintError  =
            estimationOutputWithSmallAprioriCovariance.second;
    Eigen::VectorXd moderateConstraintError  =
            estimationOutputWithModerateAprioriCovariance.second;
    Eigen::VectorXd noConstraintError  =
            estimationOutputWithoutAprioriCovariance.second;
    Eigen::VectorXd noConstraintWeakWeightError  =
            estimationOutputWithoutAprioriCovarianceAndWeakWeight.second;

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
            BOOST_CHECK_CLOSE_FRACTION( constantWeight * noConstraintInverseCovariance( i, j ),
                                        noConstraintInverseCovarianceWithWeakWeight( i, j ), 1.0E-8 );

            // Check if tight a priori constraints are processed correctly to a posteriori covariance
            if( i == j )
            {
                BOOST_CHECK_CLOSE_FRACTION(
                            tightConstraintInverseCovariance( i, j ), 1.0E32, 1.0E-10 );
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
        executeEarthOrbiterParameterEstimation< double, double >(
                    podData, 1.0E7, i, 0, false );
        manualCovarianes[ i ] = podData.first->getUnnormalizedCovarianceMatrix( );
    }

    // Use final calculations to compute covariance as a function of time
    std::map< double, Eigen::MatrixXd > automaticCovariances = simulation_setup::calculateCovarianceUsingDataUpToEpoch(
                podData.second, podData.first, 86400.0 - 1.0 );

    // Check consistency
    int counter = 1;
    for( std::map< double, Eigen::MatrixXd >::const_iterator covarianceIterator = automaticCovariances.begin( );
         covarianceIterator != automaticCovariances.end( ); covarianceIterator++ )
    {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( covarianceIterator->second, manualCovarianes.at( counter ), 1.0E-8 );
        counter++;
    }
}

BOOST_AUTO_TEST_CASE( test_WeightDefinitions )

{

    //Load spice kernels.
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
    BodyListSettings bodySettings =
        getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
    bodySettings.at( "Earth" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
        "ECLIPJ2000", "IAU_Earth",
        spice_interface::computeRotationQuaternionBetweenFrames(
            "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
        initialEphemerisTime, 2.0 * mathematical_constants::PI /
                              ( physical_constants::JULIAN_DAY ) );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setConstantBodyMass( 400.0 );


    bodies.at( "Vehicle" )->setEphemeris( std::make_shared< TabulatedCartesianEphemeris< > >(
        std::shared_ptr< interpolators::OneDimensionalInterpolator
            < double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );


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
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
        bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
        = unit_conversions::convertDegreesToRadians( 235.7 );
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
        = unit_conversions::convertDegreesToRadians( 23.4 );
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Set (perturbed) initial state.
    Eigen::Matrix< double, 6, 1 > systemInitialState = convertKeplerianToCartesianElements(
        asterixInitialStateInKeplerianElements, earthGravitationalParameter );

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
        std::make_shared< TranslationalStatePropagatorSettings< double, double > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState,
              double( finalEphemerisTime ), cowell );

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings =
        std::make_shared< RungeKuttaVariableStepSizeSettings< double > >
            ( double( initialEphemerisTime ), 40.0,
              CoefficientSets::rungeKuttaFehlberg78,
              40.0, 40.0, 1.0, 1.0 );

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
        std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >(
            "Vehicle", systemInitialState, "Earth" ) );

    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
        createParametersToEstimate< double, double >( parameterNames, bodies );

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;

    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;


        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            observationSettingsList.push_back(
                std::make_shared< ObservationModelSettings >(
                    currentObservable, currentLinkEndsList.at( i ) ) );
        }
    }

    // Create orbit determination object.
    OrbitDeterminationManager< double, double > orbitDeterminationManager =
        OrbitDeterminationManager< double, double >(
            bodies, parametersToEstimate, observationSettingsList,
            integratorSettings, propagatorSettings );

    std::vector< double > baseTimeList;
    double observationTimeStart = initialEphemerisTime + 1000.0;
    double  observationInterval = 20.0;
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
        getObservationSimulationSettings< double >(
            linkEndsPerObservable, baseTimeList, receiver );

    // Simulate observations
    std::shared_ptr< ObservationCollection< double, double > > simulatedObservations =
        simulateObservations< double, double >(
            measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );


    // Define estimation input
    std::shared_ptr< EstimationInput< double, double  > > estimationInput =
        std::make_shared< EstimationInput< double, double > >(
            simulatedObservations );

    std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize =
        simulatedObservations->getObservationTypeStartAndSize( );

    {
        simulatedObservations->setConstantWeight( 0.1 );

        // Define estimation input
        std::shared_ptr< EstimationInput< double, double  > > estimationInput = std::make_shared< EstimationInput< double, double > >( simulatedObservations );
        std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize = simulatedObservations->getObservationTypeStartAndSize( );

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
        std::shared_ptr< EstimationInput< double, double  > > estimationInput = std::make_shared< EstimationInput< double, double > >( simulatedObservations );
        std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize = simulatedObservations->getObservationTypeStartAndSize( );

        Eigen::VectorXd totalWeights = estimationInput->getWeightsMatrixDiagonals( );

        for( auto it : weightPerObservationParser )
        {
            ObservableType observableType = std::dynamic_pointer_cast< ObservationCollectionObservableTypeParser >( it.first )->getObservableTypes( ).at( 0 );
            for( int i = 0; i < observationTypeStartAndSize.at( observableType ).second; i++ )
            {
                BOOST_CHECK_CLOSE_FRACTION( totalWeights( observationTypeStartAndSize.at( observableType ).first + i ), it.second, std::numeric_limits< double >::epsilon( ) );
            }
        }
    }

    {
        Eigen::Vector2d angularPositionWeight;
        angularPositionWeight << 0.1, 0.2;
        simulatedObservations->setConstantWeight( 2.0 );
        simulatedObservations->setConstantWeight( angularPositionWeight, observationParser( angular_position ) );

        // Define estimation input
        std::shared_ptr< EstimationInput< double, double  > > estimationInput = std::make_shared< EstimationInput< double, double > >( simulatedObservations );
        std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize = simulatedObservations->getObservationTypeStartAndSize( );

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
                BOOST_CHECK_CLOSE_FRACTION( totalWeights( startEndIndex.first + i ), angularPositionWeight( 0 ), std::numeric_limits< double >::epsilon( ) );
            }
            else
            {
                BOOST_CHECK_CLOSE_FRACTION( totalWeights( startEndIndex.first + i ), angularPositionWeight( 1 ), std::numeric_limits< double >::epsilon( ) );
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
        Eigen::VectorXd singleSetRangeWeights = Eigen::VectorXd::LinSpaced( sizeRangeObsPerObsSet, 1.0 / ( 3.0 * 3.0 ), 1.0 / ( 4.0 * 4.0 ) );

        // Compute full range weight vector
        unsigned int nbRangeObsSets = simulatedObservations->getSingleObservationSets( observationParser( one_way_range ) ).size( );
        Eigen::VectorXd rangeWeights = Eigen::VectorXd::Zero( nbRangeObsSets * sizeRangeObsPerObsSet );
        for ( unsigned int k = 0 ; k < nbRangeObsSets ; k++ )
        {
            rangeWeights.segment( k * sizeRangeObsPerObsSet, sizeRangeObsPerObsSet ) = singleSetRangeWeights;
        }

        // Set total tabulated weights for all Doppler observation sets
        int totalSizeDopplerObs = simulatedObservations->getSingleObservationSets( observationParser( one_way_doppler ) ).size( ) * nbObsPerDay * numberOfDaysOfData;
        Eigen::VectorXd dopplerWeights = Eigen::VectorXd::LinSpaced(
                totalSizeDopplerObs, 1.0 / ( 1.0e-11 * SPEED_OF_LIGHT * 1.0e-11 * SPEED_OF_LIGHT ), 1.0 / ( 1.5e-11 * SPEED_OF_LIGHT * 1.5e-11 * SPEED_OF_LIGHT ) );

        // Default angular position weights set to 1
        int totalSizeAngularPositionObs = 2.0 * simulatedObservations->getSingleObservationSets( observationParser( angular_position ) ).size( ) * nbObsPerDay * numberOfDaysOfData;
        Eigen::VectorXd angularPositionWeights = Eigen::VectorXd::Ones( totalSizeAngularPositionObs );

        // Concatenate tabulated weights per observable type (default weihts for angular_position observables)
        std::map< std::shared_ptr< observation_models::ObservationCollectionParser >, Eigen::VectorXd > weightPerObservationParser;
        weightPerObservationParser[ observationParser( one_way_range ) ] = singleSetRangeWeights;
        weightPerObservationParser[ observationParser( one_way_doppler ) ] = dopplerWeights;
        weightPerObservationParser[ observationParser( angular_position ) ] = angularPositionWeights;
        simulatedObservations->setTabulatedWeights( weightPerObservationParser );

        // Define estimation input
        std::shared_ptr< EstimationInput< double, double  > > estimationInput = std::make_shared< EstimationInput< double, double > >( simulatedObservations );
        std::map< ObservableType, std::pair< int, int > > observationTypeStartAndSize = simulatedObservations->getObservationTypeStartAndSize( );
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
                BOOST_CHECK_CLOSE_FRACTION( totalWeights( observationTypeStartAndSize.at( it.first ).first + i ), it.second( i ), std::numeric_limits< double >::epsilon( ) );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}



