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

#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/simulation/estimation_setup/observationDataset.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManagerHelpers.h"
#include "tudat/simulation/estimation_setup/outlierRejection.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;
using namespace tudat::simulation_setup;

BOOST_AUTO_TEST_SUITE( test_outlier_rejection )

LinkDefinition createTestLinkDefinition( const std::string& stationName )
{
    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( "Earth", stationName );
    linkEnds[ receiver ] = LinkEndId( "Vehicle", "" );
    return LinkDefinition( linkEnds );
}

//! Create a dataset with three one-way range observations and two angular position observations.
/*!
 * The dataset contains both a single-component and a two-component observable, so that the bookkeeping between
 * observations and the rows they occupy is exercised. Observation ids 0 to 2 are the range observations, ids 3 and 4
 * are the angular position observations, which occupy two rows each.
 */
std::shared_ptr< ObservationDataset< double, double > > createTestDataset( )
{
    std::shared_ptr< ObservationDataset< double, double > > dataset = std::make_shared< ObservationDataset< double, double > >( );

    dataset->addObservationSet( one_way_range,
                                createTestLinkDefinition( "StationA" ),
                                { ( Eigen::VectorXd( 1 ) << 10.0 ).finished( ),
                                  ( Eigen::VectorXd( 1 ) << 11.0 ).finished( ),
                                  ( Eigen::VectorXd( 1 ) << 12.0 ).finished( ) },
                                { 1.0, 2.0, 3.0 },
                                receiver,
                                std::vector< Eigen::VectorXd >( ),
                                nullptr,
                                nullptr,
                                { ( Eigen::VectorXd( 1 ) << 4.0 ).finished( ),
                                  ( Eigen::VectorXd( 1 ) << 5.0 ).finished( ),
                                  ( Eigen::VectorXd( 1 ) << 8.0 ).finished( ) } );

    dataset->addObservationSet( angular_position,
                                createTestLinkDefinition( "StationB" ),
                                { ( Eigen::Vector2d( ) << 20.0, 21.0 ).finished( ), ( Eigen::Vector2d( ) << 22.0, 23.0 ).finished( ) },
                                { 4.0, 5.0 },
                                receiver,
                                std::vector< Eigen::VectorXd >( ),
                                nullptr,
                                nullptr,
                                { ( Eigen::Vector2d( ) << 2.0, 10.0 ).finished( ), ( Eigen::Vector2d( ) << 20.0, 40.0 ).finished( ) } );

    return dataset;
}

//! Outlier rejection algorithm that rejects a fixed list of observations, used to test the base class machinery.
class FixedListOutlierRejection : public OutlierRejection< double, double >
{
public:
    FixedListOutlierRejection( const std::shared_ptr< ObservationDataset< double, double > >& observationDataset,
                               const std::vector< unsigned int >& observationsToReject ):
        OutlierRejection< double, double >( OutlierRejectionType::carpino_outlier_rejection, observationDataset ),
        observationsToReject_( observationsToReject )
    { }

    //! Number of rows that the algorithm was given in the last call, used to check that rejected observations are
    //! provided to the algorithm as well.
    int numberOfRowsInInput_ = -1;

protected:
    void computeRejectionStatus( const OutlierRejectionInput< double, double >& outlierRejectionInput ) override
    {
        numberOfRowsInInput_ = static_cast< int >( outlierRejectionInput.flattenedObservationData_.getObservationVector( ).size( ) );

        std::fill( isRejected_.begin( ), isRejected_.end( ), false );
        for( const unsigned int observationId : observationsToReject_ )
        {
            isRejected_.at( observationId ) = true;
        }
    }

private:
    std::vector< unsigned int > observationsToReject_;
};

//! Check that the settings are validated on creation, and that they are stored on the estimation input.
BOOST_AUTO_TEST_CASE( test_OutlierRejectionSettings )
{
    const std::shared_ptr< OutlierRejectionSettings > settings = carpinoOutlierRejectionSettings( 9.0, 8.0, 0.3, 2 );

    const std::shared_ptr< CarpinoOutlierRejectionSettings > carpinoSettings =
            std::dynamic_pointer_cast< CarpinoOutlierRejectionSettings >( settings );
    BOOST_REQUIRE( carpinoSettings != nullptr );
    BOOST_CHECK( ( carpinoSettings->getOutlierRejectionType( ) == OutlierRejectionType::carpino_outlier_rejection ) );
    BOOST_CHECK_CLOSE_FRACTION( carpinoSettings->getChi2RejectionThreshold( ), 9.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( carpinoSettings->getChi2RecoveryThreshold( ), 8.0, 1.0E-15 );
    BOOST_CHECK_CLOSE_FRACTION( carpinoSettings->getMaximumRejectedFraction( ), 0.3, 1.0E-15 );
    BOOST_CHECK_EQUAL( carpinoSettings->getFirstIterationWithRejection( ), 2 );

    // The recovery threshold must be smaller than the rejection threshold, and the rejected fraction must be a fraction
    BOOST_CHECK_THROW( carpinoOutlierRejectionSettings( 8.0, 9.0 ), std::runtime_error );
    BOOST_CHECK_THROW( carpinoOutlierRejectionSettings( 9.0, 8.0, 1.5 ), std::runtime_error );

    // No outlier rejection is performed by default; settings can be provided on construction or afterwards
    const std::shared_ptr< ObservationDataset< double, double > > dataset = createTestDataset( );
    BOOST_CHECK( ( EstimationInput< double, double >( dataset ).getOutlierRejectionSettings( ) == nullptr ) );

    EstimationInput< double, double > estimationInput( dataset,
                                                      Eigen::MatrixXd::Zero( 0, 0 ),
                                                      std::make_shared< EstimationConvergenceChecker >( ),
                                                      Eigen::MatrixXd::Zero( 0, 0 ),
                                                      Eigen::VectorXd::Zero( 0 ),
                                                      true,
                                                      settings );
    BOOST_CHECK( estimationInput.getOutlierRejectionSettings( ) == settings );

    estimationInput.setOutlierRejectionSettings( nullptr );
    BOOST_CHECK( estimationInput.getOutlierRejectionSettings( ) == nullptr );
}

//! Check that the outlier rejection object created from the settings is of the corresponding type.
BOOST_AUTO_TEST_CASE( test_OutlierRejectionCreation )
{
    const std::shared_ptr< ObservationDataset< double, double > > dataset = createTestDataset( );

    // Without settings, no outlier rejection object is created
    BOOST_CHECK( ( createOutlierRejection< double, double >( nullptr, dataset ) == nullptr ) );

    const std::shared_ptr< OutlierRejection< double, double > > outlierRejection =
            createOutlierRejection< double, double >( carpinoOutlierRejectionSettings( ), dataset );
    BOOST_REQUIRE( outlierRejection != nullptr );
    BOOST_CHECK( ( std::dynamic_pointer_cast< CarpinoOutlierRejection< double, double > >( outlierRejection ) != nullptr ) );
    BOOST_CHECK( ( outlierRejection->getOutlierRejectionType( ) == OutlierRejectionType::carpino_outlier_rejection ) );

    // All observations start out as accepted, since none of them is rejected in the dataset
    BOOST_CHECK_EQUAL( outlierRejection->getRejectionStatus( ).size( ), dataset->getNumberOfObservations( ) );
    BOOST_CHECK_EQUAL( outlierRejection->getNumberOfRejectedObservations( ), 0 );

    // Observations that are already rejected in the dataset start out as rejected
    dataset->rejectObservations( ObservationSelectionCondition< double, double >::observableType( angular_position ) );
    const std::shared_ptr< OutlierRejection< double, double > > rejectionFromRejectedDataset =
            createOutlierRejection< double, double >( carpinoOutlierRejectionSettings( ), dataset );
    BOOST_CHECK_EQUAL( rejectionFromRejectedDataset->getNumberOfRejectedObservations( ), 2 );

    // The algorithm itself is not implemented yet
    const Eigen::MatrixXd emptyMatrix = Eigen::MatrixXd::Zero( 0, 0 );
    const Eigen::VectorXd emptyVector = Eigen::VectorXd::Zero( 0 );
    const FlattenedObservationData< double, double > flattenedData = dataset->createOrderedFlattenedObservationData( true );
    const ObservationCovarianceInterface< double, double > observationCovarianceInterface( flattenedData );
    const OutlierRejectionInput< double, double > outlierRejectionInput(
            0, flattenedData, observationCovarianceInterface, emptyVector, emptyMatrix, emptyMatrix, emptyVector );
    BOOST_CHECK_THROW( outlierRejection->updateRejectionStatus( outlierRejectionInput ), std::runtime_error );
}

//! Check that the rejection status computed by an algorithm is applied to the observation dataset.
BOOST_AUTO_TEST_CASE( test_OutlierRejectionStatusIsAppliedToDataset )
{
    const std::shared_ptr< ObservationDataset< double, double > > dataset = createTestDataset( );

    const FlattenedObservationData< double, double > flattenedData = dataset->createOrderedFlattenedObservationData( true );
    const ObservationCovarianceInterface< double, double > observationCovarianceInterface( flattenedData );
    const Eigen::VectorXd residuals = flattenedData.getResidualVector( );
    const Eigen::MatrixXd designMatrix = Eigen::MatrixXd::Ones( residuals.rows( ), 2 );
    const Eigen::MatrixXd parameterCovariance = Eigen::MatrixXd::Identity( 2, 2 );
    const Eigen::VectorXd parameterCorrection = Eigen::VectorXd::Zero( 2 );

    const OutlierRejectionInput< double, double > outlierRejectionInput(
            0, flattenedData, observationCovarianceInterface, residuals, designMatrix, parameterCovariance, parameterCorrection );

    // Reject one range observation and one angular position observation
    FixedListOutlierRejection outlierRejection( dataset, { 1, 4 } );
    outlierRejection.updateRejectionStatus( outlierRejectionInput );

    BOOST_CHECK_EQUAL( outlierRejection.getNumberOfRejectedObservations( ), 2 );
    BOOST_CHECK_EQUAL( dataset->getObservationRow( 0 ).isActive_, true );
    BOOST_CHECK_EQUAL( dataset->getObservationRow( 1 ).isActive_, false );
    BOOST_CHECK_EQUAL( dataset->getObservationRow( 4 ).isActive_, false );

    // The rejected observations must be excluded from the data used in the estimation. The rejected angular position
    // observation removes two rows, since both of its components are rejected.
    BOOST_CHECK_EQUAL( dataset->createOrderedFlattenedObservationData( false ).getObservationVector( ).size( ), 4 );
    BOOST_CHECK_EQUAL( dataset->createOrderedFlattenedObservationData( true ).getObservationVector( ).size( ), 7 );

    // Rejected observations must be provided to the algorithm, so that they can be recovered
    outlierRejection.updateRejectionStatus( outlierRejectionInput );
    BOOST_CHECK_EQUAL( outlierRejection.numberOfRowsInInput_, 7 );

    // An empty list of outliers recovers all observations
    FixedListOutlierRejection recoveringOutlierRejection( dataset, { } );
    recoveringOutlierRejection.updateRejectionStatus( outlierRejectionInput );

    BOOST_CHECK_EQUAL( recoveringOutlierRejection.getNumberOfRejectedObservations( ), 0 );
    BOOST_CHECK_EQUAL( dataset->createOrderedFlattenedObservationData( false ).getObservationVector( ).size( ), 7 );
}

//! Check that the observation covariance is the inverse of the weight matrix, also with off-diagonal weights.
BOOST_AUTO_TEST_CASE( test_ObservationCovarianceInterface )
{
    const std::shared_ptr< ObservationDataset< double, double > > dataset = createTestDataset( );
    dataset->setConstantSingleObservationScalarWeight( ObservationSelectionCondition< double, double >::observableType( one_way_range ),
                                                       4.0 );
    dataset->setConstantSingleObservationDiagonalWeight(
            ObservationSelectionCondition< double, double >::observableType( angular_position ),
            ( Eigen::Vector2d( ) << 2.0, 8.0 ).finished( ) );

    {
        const FlattenedObservationData< double, double > flattenedData = dataset->createOrderedFlattenedObservationData( true );
        const ObservationCovarianceInterface< double, double > observationCovarianceInterface( flattenedData );

        // For diagonal weights, the covariance of an observation is the inverse of its own weights
        BOOST_CHECK_CLOSE_FRACTION( observationCovarianceInterface.getObservationCovariance( 0 )( 0, 0 ), 1.0 / 4.0, 1.0E-12 );

        const Eigen::MatrixXd angularCovariance = observationCovarianceInterface.getObservationCovariance( 3 );
        BOOST_CHECK_EQUAL( angularCovariance.rows( ), 2 );
        BOOST_CHECK_EQUAL( angularCovariance.cols( ), 2 );
        BOOST_CHECK_CLOSE_FRACTION( angularCovariance( 0, 0 ), 1.0 / 2.0, 1.0E-12 );
        BOOST_CHECK_CLOSE_FRACTION( angularCovariance( 1, 1 ), 1.0 / 8.0, 1.0E-12 );
        BOOST_CHECK_SMALL( angularCovariance( 0, 1 ), 1.0E-12 );

        // Uncorrelated observations have a zero cross-covariance
        BOOST_CHECK_SMALL( observationCovarianceInterface.getObservationCovariance( 0, 1 )( 0, 0 ), 1.0E-12 );
    }

    {
        // Correlate the first two range observations. The covariance of a single observation is then a block of the
        // inverse of the complete weight matrix, and not the inverse of that observation's own weight block.
        Eigen::MatrixXd weightBlock = Eigen::MatrixXd::Zero( 2, 2 );
        weightBlock << 4.0, 1.0, 1.0, 4.0;
        dataset->setWeightBlock( { 0, 1 }, { 0, 1 }, weightBlock );

        const FlattenedObservationData< double, double > flattenedData = dataset->createOrderedFlattenedObservationData( true );
        const ObservationCovarianceInterface< double, double > observationCovarianceInterface( flattenedData );

        const Eigen::MatrixXd expectedCovariance = weightBlock.inverse( );
        BOOST_CHECK_CLOSE_FRACTION( observationCovarianceInterface.getObservationCovariance( 0 )( 0, 0 ), expectedCovariance( 0, 0 ), 1.0E-12 );
        BOOST_CHECK_CLOSE_FRACTION( observationCovarianceInterface.getObservationCovariance( 0, 1 )( 0, 0 ), expectedCovariance( 0, 1 ), 1.0E-12 );

        // The uncorrelated observations are unaffected
        BOOST_CHECK_CLOSE_FRACTION( observationCovarianceInterface.getObservationCovariance( 2 )( 0, 0 ), 1.0 / 4.0, 1.0E-12 );
    }
}

//! Check that the rows of the observations used in the estimation are correctly extracted from the complete data.
BOOST_AUTO_TEST_CASE( test_EstimationObservationRowExtraction )
{
    const std::shared_ptr< ObservationDataset< double, double > > dataset = createTestDataset( );
    dataset->rejectObservations( ObservationSelectionCondition< double, double >::timeGreaterThan( 1.0 ) &&
                                 ObservationSelectionCondition< double, double >::timeLessThan( 4.5 ) );

    const FlattenedObservationData< double, double > computationData = dataset->createOrderedFlattenedObservationData( true );
    const FlattenedObservationData< double, double > estimationData = dataset->createOrderedFlattenedObservationData( false );

    // Observations 1, 2 and 3 are rejected, so 3 of the 7 rows remain
    BOOST_REQUIRE_EQUAL( computationData.getObservationVector( ).size( ), 7 );
    BOOST_REQUIRE_EQUAL( estimationData.getObservationVector( ).size( ), 3 );

    // Fill each row of the computed data with the number of the row it belongs to, so that the origin of every
    // extracted row can be verified
    Eigen::MatrixXd computedRows = Eigen::MatrixXd::Zero( 7, 2 );
    for( int row = 0; row < computedRows.rows( ); row++ )
    {
        computedRows( row, 0 ) = static_cast< double >( row );
        computedRows( row, 1 ) = 100.0 + static_cast< double >( row );
    }

    const Eigen::MatrixXd estimationRows = extractEstimationObservationRows( computationData, estimationData, computedRows );
    BOOST_REQUIRE_EQUAL( estimationRows.rows( ), 3 );

    for( int row = 0; row < estimationRows.rows( ); row++ )
    {
        const unsigned int observationId = estimationData.getObservationIds( ).at( row );
        const int componentIndex = row - estimationData.getFirstFlattenedRowForObservation( observationId );
        const int expectedComputationRow = computationData.getFirstFlattenedRowForObservation( observationId ) + componentIndex;

        BOOST_CHECK_CLOSE_FRACTION( estimationRows( row, 0 ), static_cast< double >( expectedComputationRow ), 1.0E-15 );
        BOOST_CHECK_CLOSE_FRACTION( estimationRows( row, 1 ), 100.0 + static_cast< double >( expectedComputationRow ), 1.0E-15 );
    }

    // The same extraction must work for a vector of residuals
    const Eigen::VectorXd estimationResiduals =
            extractEstimationObservationRows( computationData, estimationData, Eigen::VectorXd( computationData.getResidualVector( ) ) );
    BOOST_REQUIRE_EQUAL( estimationResiduals.rows( ), 3 );
    for( int row = 0; row < estimationResiduals.rows( ); row++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( estimationResiduals( row ), estimationData.getResidualVector( )( row ), 1.0E-15 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
