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

#include <cmath>
#include <limits>
#include <vector>

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////  Helper functions and objects, shared by all test cases in this file          ////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

LinkDefinition createTestLinkDefinition( const std::string& stationName )
{
    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( "Earth", stationName );
    linkEnds[ receiver ] = LinkEndId( "Vehicle", "" );
    return LinkDefinition( linkEnds );
}

//! Weights of the test datasets, which are the inverse of the square of an uncertainty of the observations.
/*!
 * The algorithm of Carpino et al. (2003) derives the covariance of a residual from the weight of the observation, and
 * therefore requires weights that represent the uncertainty of the observations. It refuses to run on a dataset whose
 * weights are all equal to one, since those are the weights that a dataset has when no weights were set at all. The
 * test datasets below therefore set weights explicitly, from an uncertainty of two metres for the range observations
 * and of half a unit for both components of the angular position observations.
 */
const double testRangeUncertainty = 2.0;
const double testAngularPositionUncertainty = 0.5;
const double testRangeWeight = 1.0 / ( testRangeUncertainty * testRangeUncertainty );
const double testAngularPositionWeight = 1.0 / ( testAngularPositionUncertainty * testAngularPositionUncertainty );

//! Create a dataset with three one-way range observations and two angular position observations.
/*!
 * The dataset contains both a single-component and a two-component observable, so that the bookkeeping between
 * observations and the rows they occupy is exercised. Observation ids 0 to 2 are the range observations, ids 3 and 4
 * are the angular position observations, which occupy two rows each. The residuals grow with the observation id, so
 * that observation 4 is by far the largest outlier of the dataset.
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
                                std::vector< Eigen::VectorXd >( ),
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
                                std::vector< Eigen::VectorXd >( ),
                                { ( Eigen::Vector2d( ) << 2.0, 10.0 ).finished( ), ( Eigen::Vector2d( ) << 20.0, 40.0 ).finished( ) } );

    dataset->setConstantSingleObservationScalarWeight( ObservationSelectionCondition< double, double >::observableType( one_way_range ),
                                                       testRangeWeight );
    dataset->setConstantSingleObservationDiagonalWeight(
            ObservationSelectionCondition< double, double >::observableType( angular_position ),
            ( Eigen::Vector2d( ) << testAngularPositionWeight, testAngularPositionWeight ).finished( ) );

    return dataset;
}

//! Create a dataset with a requested number of one-way range observations.
/*!
 * Test cases that check a step of an algorithm in isolation provide the residuals or chi-squared values that the step
 * is to work on themselves, and only need the dataset to define how many observations exist. This function creates
 * such a dataset of any size, with one observation per observation id.
 */
std::shared_ptr< ObservationDataset< double, double > > createRangeDataset( const unsigned int numberOfObservations )
{
    std::vector< Eigen::VectorXd > observations;
    std::vector< double > times;
    for( unsigned int observationIndex = 0; observationIndex < numberOfObservations; observationIndex++ )
    {
        observations.push_back( ( Eigen::VectorXd( 1 ) << 10.0 + static_cast< double >( observationIndex ) ).finished( ) );
        times.push_back( 1.0 + static_cast< double >( observationIndex ) );
    }

    std::shared_ptr< ObservationDataset< double, double > > dataset = std::make_shared< ObservationDataset< double, double > >( );
    dataset->addObservationSet( one_way_range, createTestLinkDefinition( "StationA" ), observations, times, receiver );
    dataset->setConstantSingleObservationScalarWeight( ObservationSelectionCondition< double, double >::all( ), testRangeWeight );
    return dataset;
}

//! Create a dataset of one-way range observations that have the requested residuals.
/*!
 * Test cases of the simple algorithm compare the residual of an observation against a threshold, and are most readable
 * when the residuals themselves are written out. Observation ids follow the order of the residuals, and the
 * observations are one second apart, starting at a time of one second.
 */
std::shared_ptr< ObservationDataset< double, double > > createResidualDataset( const std::vector< double >& residualValues )
{
    std::vector< Eigen::VectorXd > observations;
    std::vector< double > times;
    std::vector< Eigen::VectorXd > residuals;
    for( unsigned int observationIndex = 0; observationIndex < residualValues.size( ); observationIndex++ )
    {
        observations.push_back( ( Eigen::VectorXd( 1 ) << 10.0 + static_cast< double >( observationIndex ) ).finished( ) );
        times.push_back( 1.0 + static_cast< double >( observationIndex ) );
        residuals.push_back( ( Eigen::VectorXd( 1 ) << residualValues.at( observationIndex ) ).finished( ) );
    }

    std::shared_ptr< ObservationDataset< double, double > > dataset = std::make_shared< ObservationDataset< double, double > >( );
    dataset->addObservationSet( one_way_range,
                                createTestLinkDefinition( "StationA" ),
                                observations,
                                times,
                                receiver,
                                std::vector< Eigen::VectorXd >( ),
                                nullptr,
                                nullptr,
                                std::vector< Eigen::VectorXd >( ),
                                residuals );
    dataset->setConstantSingleObservationScalarWeight( ObservationSelectionCondition< double, double >::all( ), testRangeWeight );
    return dataset;
}

//! Compute the observation covariance the way the estimation does: by inverting the complete weight matrix.
Eigen::MatrixXd createObservationCovariance( const FlattenedObservationData< double, double >& flattenedObservationData )
{
    return Eigen::MatrixXd( flattenedObservationData.getSparseWeightMatrix( ) ).inverse( );
}

//! Data of one estimation iteration, of the size that a dataset requires.
/*!
 * An OutlierRejectionInput stores references only, so everything it refers to must stay alive for as long as the
 * input is used. This object owns that data, and hands out inputs that refer to it. The design matrix, the parameter
 * covariance and the parameter correction are all zero, which is enough for algorithms that only inspect residuals.
 */
struct TestIterationData
{
    TestIterationData( const std::shared_ptr< ObservationDataset< double, double > >& dataset, const int numberOfParameters = 2 ):
        flattenedData_( dataset->createOrderedFlattenedObservationData( true ) ),
        observationCovariance_( createObservationCovariance( flattenedData_ ) ),
        residuals_( flattenedData_.getResidualVector( ) ),
        designMatrix_( Eigen::MatrixXd::Zero( residuals_.rows( ), numberOfParameters ) ),
        parameterCovariance_( Eigen::MatrixXd::Zero( numberOfParameters, numberOfParameters ) ),
        parameterCorrection_( Eigen::VectorXd::Zero( numberOfParameters ) )
    { }

    OutlierRejectionInput< double, double > getInput( const int iterationNumber ) const
    {
        return OutlierRejectionInput< double, double >( iterationNumber,
                                                        flattenedData_,
                                                        observationCovariance_,
                                                        residuals_,
                                                        designMatrix_,
                                                        parameterCovariance_,
                                                        parameterCorrection_ );
    }

    FlattenedObservationData< double, double > flattenedData_;
    Eigen::MatrixXd observationCovariance_;
    Eigen::VectorXd residuals_;
    Eigen::MatrixXd designMatrix_;
    Eigen::MatrixXd parameterCovariance_;
    Eigen::VectorXd parameterCorrection_;
};

//! Return the ids of the observations that are rejected, in increasing order, to be compared against an expected list.
std::vector< unsigned int > getRejectedObservationIds( const std::vector< bool >& isRejected )
{
    std::vector< unsigned int > rejectedObservationIds;
    for( unsigned int observationId = 0; observationId < isRejected.size( ); observationId++ )
    {
        if( isRejected.at( observationId ) )
        {
            rejectedObservationIds.push_back( observationId );
        }
    }
    return rejectedObservationIds;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////  Tests of the outlier rejection framework, independent of any algorithm       ////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE( outlier_rejection_framework )

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
}

//! Check that the rejection status computed by an algorithm is applied to the observation dataset.
BOOST_AUTO_TEST_CASE( test_OutlierRejectionStatusIsAppliedToDataset )
{
    const std::shared_ptr< ObservationDataset< double, double > > dataset = createTestDataset( );

    const FlattenedObservationData< double, double > flattenedData = dataset->createOrderedFlattenedObservationData( true );
    const Eigen::MatrixXd observationCovariance = createObservationCovariance( flattenedData );
    const Eigen::VectorXd residuals = flattenedData.getResidualVector( );
    const Eigen::MatrixXd designMatrix = Eigen::MatrixXd::Ones( residuals.rows( ), 2 );
    const Eigen::MatrixXd parameterCovariance = Eigen::MatrixXd::Identity( 2, 2 );
    const Eigen::VectorXd parameterCorrection = Eigen::VectorXd::Zero( 2 );

    const OutlierRejectionInput< double, double > outlierRejectionInput(
            0, flattenedData, observationCovariance, residuals, designMatrix, parameterCovariance, parameterCorrection );

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

//! Check that the per-observation covariance is the inverse of the weight matrix, also with off-diagonal weights.
BOOST_AUTO_TEST_CASE( test_ObservationCovariance )
{
    const Eigen::MatrixXd emptyMatrix = Eigen::MatrixXd::Zero( 0, 0 );
    const Eigen::VectorXd emptyVector = Eigen::VectorXd::Zero( 0 );

    const std::shared_ptr< ObservationDataset< double, double > > dataset = createTestDataset( );
    dataset->setConstantSingleObservationScalarWeight( ObservationSelectionCondition< double, double >::observableType( one_way_range ),
                                                       4.0 );
    dataset->setConstantSingleObservationDiagonalWeight(
            ObservationSelectionCondition< double, double >::observableType( angular_position ),
            ( Eigen::Vector2d( ) << 2.0, 8.0 ).finished( ) );

    {
        const FlattenedObservationData< double, double > flattenedData = dataset->createOrderedFlattenedObservationData( true );

        // The input stores references only, so the data it refers to must outlive it and cannot be a temporary
        const Eigen::MatrixXd observationCovariance = createObservationCovariance( flattenedData );
        const OutlierRejectionInput< double, double > input(
                0, flattenedData, observationCovariance, flattenedData.getResidualVector( ), emptyMatrix, emptyMatrix, emptyVector );

        // For diagonal weights, the covariance of an observation is the inverse of its own weights
        BOOST_CHECK_CLOSE_FRACTION( input.getObservationCovariance( 0 )( 0, 0 ), 1.0 / 4.0, 1.0E-12 );

        const Eigen::MatrixXd angularCovariance = input.getObservationCovariance( 3 );
        BOOST_CHECK_EQUAL( angularCovariance.rows( ), 2 );
        BOOST_CHECK_EQUAL( angularCovariance.cols( ), 2 );
        BOOST_CHECK_CLOSE_FRACTION( angularCovariance( 0, 0 ), 1.0 / 2.0, 1.0E-12 );
        BOOST_CHECK_CLOSE_FRACTION( angularCovariance( 1, 1 ), 1.0 / 8.0, 1.0E-12 );
        BOOST_CHECK_SMALL( angularCovariance( 0, 1 ), 1.0E-12 );
    }

    {
        // Correlate the first two range observations. The covariance of a single observation is then a block of the
        // inverse of the complete weight matrix, and not the inverse of that observation's own weight block.
        Eigen::MatrixXd weightBlock = Eigen::MatrixXd::Zero( 2, 2 );
        weightBlock << 4.0, 1.0, 1.0, 4.0;
        dataset->setWeightBlock( { 0, 1 }, { 0, 1 }, weightBlock );

        const FlattenedObservationData< double, double > flattenedData = dataset->createOrderedFlattenedObservationData( true );
        const Eigen::MatrixXd observationCovariance = createObservationCovariance( flattenedData );
        const OutlierRejectionInput< double, double > input(
                0, flattenedData, observationCovariance, flattenedData.getResidualVector( ), emptyMatrix, emptyMatrix, emptyVector );

        // The correlation must be accounted for: the covariance of the observation is a block of the inverse of the
        // complete weight matrix, which for a correlated pair differs from the inverse of its own weight entry
        const Eigen::MatrixXd expectedCovariance = weightBlock.inverse( );
        BOOST_CHECK_CLOSE_FRACTION( input.getObservationCovariance( 0 )( 0, 0 ), expectedCovariance( 0, 0 ), 1.0E-12 );
        BOOST_CHECK( std::fabs( expectedCovariance( 0, 0 ) - 1.0 / 4.0 ) > 1.0E-6 );

        // The uncorrelated observations are unaffected
        BOOST_CHECK_CLOSE_FRACTION( input.getObservationCovariance( 2 )( 0, 0 ), 1.0 / 4.0, 1.0E-12 );
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

BOOST_AUTO_TEST_SUITE_END( )  // outlier_rejection_framework

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////  Tests of the individual steps of the algorithm of Carpino et al. (2003)      ////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE( carpino_outlier_rejection )

//! Carpino algorithm with the individual steps of the algorithm made publicly accessible.
/*!
 * The steps of the algorithm are protected members of CarpinoOutlierRejection: they are implementation details, and
 * are not part of the interface that the estimation uses. Deriving a class from it is the standard way of reaching
 * such members in a test, since a derived class may access the protected members of its base class. The 'using'
 * declarations below re-declare those members as public members of this test-only class, without changing anything
 * about the functions themselves.
 */
class TestableCarpinoOutlierRejection : public CarpinoOutlierRejection< double, double >
{
public:
    TestableCarpinoOutlierRejection( const std::shared_ptr< CarpinoOutlierRejectionSettings >& outlierRejectionSettings,
                                     const std::shared_ptr< ObservationDataset< double, double > >& observationDataset ):
        CarpinoOutlierRejection< double, double >( outlierRejectionSettings, observationDataset )
    { }

    //! Set the rejection status that the steps of the algorithm are to treat as the status of the current iteration.
    void setRejectionStatus( const std::vector< bool >& isRejected )
    {
        isRejected_ = isRejected;
    }

    using CarpinoOutlierRejection< double, double >::applyMaximumRejectedFraction;
    using CarpinoOutlierRejection< double, double >::decideRejectionStatus;
    using CarpinoOutlierRejection< double, double >::getRejectionThreshold;
    using CarpinoOutlierRejection< double, double >::computeChiSquared;
};

//! Term that the algorithm adds to the rejection threshold, which raises the threshold (and thereby makes rejection
//! less likely) when only few observations are left in the fit.
double carpinoFudgeTerm( const int numberOfAcceptedObservations )
{
    return 400.0 * std::pow( 1.2, -numberOfAcceptedObservations );
}

//! Check that observations are only rejected from the iteration at which the settings allow rejection onwards.
/*!

 */
BOOST_AUTO_TEST_CASE( test_FirstIterationWithRejection )
{
    const int firstIterationWithRejection = 2;

    const std::shared_ptr< ObservationDataset< double, double > > dataset = createTestDataset( );

    const FlattenedObservationData< double, double > flattenedData = dataset->createOrderedFlattenedObservationData( true );
    const Eigen::MatrixXd observationCovariance = createObservationCovariance( flattenedData );
    const Eigen::VectorXd residuals = flattenedData.getResidualVector( );

    const int numberOfParameters = 2;
    const Eigen::MatrixXd designMatrix = Eigen::MatrixXd::Ones( residuals.rows( ), numberOfParameters );
    const Eigen::MatrixXd parameterCovariance = 1.0E-6 * Eigen::MatrixXd::Identity( numberOfParameters, numberOfParameters );
    const Eigen::VectorXd parameterCorrection = Eigen::VectorXd::Zero( numberOfParameters );

    CarpinoOutlierRejection< double, double > outlierRejection(
            std::make_shared< CarpinoOutlierRejectionSettings >( 9.0, 8.0, 0.25, firstIterationWithRejection ), dataset );

    // In every iteration before the first iteration with rejection, the status of all observations is left unchanged,
    // even though the last angular position observation has by far the largest residual of the dataset
    for( int iterationNumber = 0; iterationNumber < firstIterationWithRejection; iterationNumber++ )
    {
        const OutlierRejectionInput< double, double > outlierRejectionInput( iterationNumber,
                                                                             flattenedData,
                                                                             observationCovariance,
                                                                             residuals,
                                                                             designMatrix,
                                                                             parameterCovariance,
                                                                             parameterCorrection );
        outlierRejection.updateRejectionStatus( outlierRejectionInput );
        BOOST_CHECK_EQUAL( outlierRejection.getNumberOfRejectedObservations( ), 0 );
    }

    // From that iteration onwards, that observation, and only that observation, is rejected
    const OutlierRejectionInput< double, double > outlierRejectionInput( firstIterationWithRejection,
                                                                         flattenedData,
                                                                         observationCovariance,
                                                                         residuals,
                                                                         designMatrix,
                                                                         parameterCovariance,
                                                                         parameterCorrection );
    outlierRejection.updateRejectionStatus( outlierRejectionInput );

    BOOST_CHECK_EQUAL( outlierRejection.getNumberOfRejectedObservations( ), 1 );
    BOOST_CHECK_EQUAL( outlierRejection.getRejectionStatus( ).at( 4 ), true );
}

//! Check the decision that is taken for a single observation, from its chi-squared and its current status.
BOOST_AUTO_TEST_CASE( test_DecideRejectionStatus )
{
    const double settingsRejectionThreshold = 9.0;
    const double recoveryThreshold = 8.0;

    TestableCarpinoOutlierRejection outlierRejection(
            std::make_shared< CarpinoOutlierRejectionSettings >( settingsRejectionThreshold, recoveryThreshold ),
            createRangeDataset( 4 ) );

    // The threshold against which observations are rejected is recomputed in every iteration, and is passed to the
    // function as an argument. It is deliberately chosen larger than the threshold from the settings here, so that the
    // two cannot be confused.
    const double iterationRejectionThreshold = 25.0;

    // Named for readability: the first argument of the function is the status of the observation so far
    const bool isAccepted = false;
    const bool isRejected = true;

    // An observation that is in the fit is rejected when its chi-squared exceeds the threshold of this iteration
    BOOST_CHECK_EQUAL( outlierRejection.decideRejectionStatus( isAccepted, 30.0, iterationRejectionThreshold ), true );

    // A chi-squared above the threshold from the settings, but below the threshold of this iteration, is not enough:
    // only the threshold of this iteration is used for rejection
    BOOST_CHECK_EQUAL( outlierRejection.decideRejectionStatus( isAccepted, 20.0, iterationRejectionThreshold ), false );

    // The comparison is strict, so an observation exactly at the threshold is kept
    BOOST_CHECK_EQUAL( outlierRejection.decideRejectionStatus( isAccepted, iterationRejectionThreshold, iterationRejectionThreshold ),
                       false );

    // An observation that is out of the fit is recovered when its chi-squared drops below the recovery threshold, for
    // which the comparison is strict as well
    BOOST_CHECK_EQUAL( outlierRejection.decideRejectionStatus( isRejected, 7.0, iterationRejectionThreshold ), false );
    BOOST_CHECK_EQUAL( outlierRejection.decideRejectionStatus( isRejected, recoveryThreshold, iterationRejectionThreshold ), true );

    // The recovery threshold is lower than the rejection threshold, so an observation with a chi-squared in between
    // the two keeps whichever status it had. This hysteresis is what prevents observations from oscillating between
    // the two states in successive iterations.
    const double chiSquaredBetweenThresholds = 0.5 * ( recoveryThreshold + iterationRejectionThreshold );
    BOOST_CHECK_EQUAL( outlierRejection.decideRejectionStatus( isRejected, chiSquaredBetweenThresholds, iterationRejectionThreshold ),
                       isRejected );
    BOOST_CHECK_EQUAL( outlierRejection.decideRejectionStatus( isAccepted, chiSquaredBetweenThresholds, iterationRejectionThreshold ),
                       isAccepted );

    // A chi-squared that is not a usable number cannot be compared against a threshold, and leaves the status of the
    // observation unchanged. A negative value indicates a residual covariance that is not positive definite.
    const std::vector< double > invalidChiSquaredValues = { -1.0,
                                                            0.0,
                                                            std::numeric_limits< double >::quiet_NaN( ),
                                                            std::numeric_limits< double >::infinity( ) };
    for( const double invalidChiSquared : invalidChiSquaredValues )
    {
        BOOST_CHECK_EQUAL( outlierRejection.decideRejectionStatus( isAccepted, invalidChiSquared, iterationRejectionThreshold ),
                           isAccepted );
        BOOST_CHECK_EQUAL( outlierRejection.decideRejectionStatus( isRejected, invalidChiSquared, iterationRejectionThreshold ),
                           isRejected );
    }
}

//! Check that the rejection threshold of an iteration follows the observations that are in the fit, with the threshold
//! from the settings as a lower limit.
BOOST_AUTO_TEST_CASE( test_RejectionThreshold )
{
    const double settingsRejectionThreshold = 9.0;
    const unsigned int numberOfObservations = 6;

    TestableCarpinoOutlierRejection outlierRejection( std::make_shared< CarpinoOutlierRejectionSettings >( settingsRejectionThreshold, 8.0 ),
                                                      createRangeDataset( numberOfObservations ) );

    // The last two observations are out of the fit, so their (very large) chi-squared may not influence the threshold
    outlierRejection.setRejectionStatus( { false, false, false, false, true, true } );
    const double fudgeTerm = carpinoFudgeTerm( numberOfObservations - 2 );

    {
        // The largest chi-squared of the observations that are in the fit is 100. A quarter of that is larger than the
        // threshold from the settings, so the threshold follows the data.
        const std::vector< double > chiSquaredPerObservation = { 4.0, 100.0, 40.0, 12.0, 1.0E4, 2.0E4 };
        BOOST_CHECK_CLOSE_FRACTION(
                outlierRejection.getRejectionThreshold( chiSquaredPerObservation ), 0.25 * 100.0 + fudgeTerm, 1.0E-12 );
    }

    {
        // The largest chi-squared of the observations that are in the fit is now 20. A quarter of that is below the
        // threshold from the settings, so the threshold from the settings is used instead.
        const std::vector< double > chiSquaredPerObservation = { 4.0, 20.0, 10.0, 12.0, 1.0E4, 2.0E4 };
        BOOST_CHECK_CLOSE_FRACTION(
                outlierRejection.getRejectionThreshold( chiSquaredPerObservation ), settingsRejectionThreshold + fudgeTerm, 1.0E-12 );
    }

    {
        // The term that is added to the threshold shrinks as the number of observations in the fit grows: for the six
        // observations above it dominates the threshold, while for a large dataset it is negligible
        BOOST_CHECK( fudgeTerm > 100.0 );

        const unsigned int largeNumberOfObservations = 100;
        TestableCarpinoOutlierRejection largeDatasetOutlierRejection(
                std::make_shared< CarpinoOutlierRejectionSettings >( settingsRejectionThreshold, 8.0 ),
                createRangeDataset( largeNumberOfObservations ) );

        std::vector< double > chiSquaredPerObservation( largeNumberOfObservations, 4.0 );
        chiSquaredPerObservation.at( 1 ) = 100.0;
        BOOST_CHECK_CLOSE_FRACTION( largeDatasetOutlierRejection.getRejectionThreshold( chiSquaredPerObservation ), 0.25 * 100.0, 1.0E-6 );
    }
}

//! Check that no more observations are rejected than the settings allow, and that the worst outliers are the ones
//! that stay rejected.
BOOST_AUTO_TEST_CASE( test_MaximumRejectedFraction )
{
    const unsigned int numberOfObservations = 10;
    const double maximumRejectedFraction = 0.4;  // 4 of the 10 observations may be rejected

    TestableCarpinoOutlierRejection outlierRejection(
            std::make_shared< CarpinoOutlierRejectionSettings >( 9.0, 8.0, maximumRejectedFraction ),
            createRangeDataset( numberOfObservations ) );

    const std::vector< bool > noObservationRejected( numberOfObservations, false );

    // Chi-squared per observation id. The values of observations 8 and 9 are large, so that they would be the worst
    // outliers if they were to be considered here.
    const std::vector< double > chiSquaredPerObservation = { 10.0, 60.0, 30.0, 50.0, 20.0, 40.0, 15.0, 25.0, 1.0E3, 2.0E3 };

    {
        // Six observations meet the rejection criterion while four may be rejected, so the two with the smallest
        // chi-squared (observations 0 and 4) are put back into the fit
        outlierRejection.setRejectionStatus( noObservationRejected );
        std::vector< bool > newRejectionStatus = { true, true, true, true, true, true, false, false, false, false };
        outlierRejection.applyMaximumRejectedFraction( newRejectionStatus, chiSquaredPerObservation );

        const std::vector< unsigned int > expectedRejectedObservationIds = { 1, 2, 3, 5 };
        const std::vector< unsigned int > rejectedObservationIds = getRejectedObservationIds( newRejectionStatus );
        BOOST_CHECK_EQUAL_COLLECTIONS( rejectedObservationIds.begin( ),
                                       rejectedObservationIds.end( ),
                                       expectedRejectedObservationIds.begin( ),
                                       expectedRejectedObservationIds.end( ) );
    }

    {
        // Observations 8 and 9 were already out of the fit before this iteration. They take up two of the four allowed
        // rejections and are not reconsidered, so only the two worst of the four new outliers remain rejected.
        outlierRejection.setRejectionStatus( { false, false, false, false, false, false, false, false, true, true } );
        std::vector< bool > newRejectionStatus = { true, true, true, true, false, false, false, false, true, true };
        outlierRejection.applyMaximumRejectedFraction( newRejectionStatus, chiSquaredPerObservation );

        const std::vector< unsigned int > expectedRejectedObservationIds = { 1, 3, 8, 9 };
        const std::vector< unsigned int > rejectedObservationIds = getRejectedObservationIds( newRejectionStatus );
        BOOST_CHECK_EQUAL_COLLECTIONS( rejectedObservationIds.begin( ),
                                       rejectedObservationIds.end( ),
                                       expectedRejectedObservationIds.begin( ),
                                       expectedRejectedObservationIds.end( ) );
    }

    {
        // Five observations were already out of the fit, which is more than the maximum allows. No new observation may
        // then be rejected, but the observations that were already rejected are left alone.
        const std::vector< bool > previouslyRejected = { false, false, false, false, true, true, true, true, true, false };
        outlierRejection.setRejectionStatus( previouslyRejected );
        std::vector< bool > newRejectionStatus = { true, true, false, false, true, true, true, true, true, false };
        outlierRejection.applyMaximumRejectedFraction( newRejectionStatus, chiSquaredPerObservation );

        const std::vector< unsigned int > expectedRejectedObservationIds = { 4, 5, 6, 7, 8 };
        const std::vector< unsigned int > rejectedObservationIds = getRejectedObservationIds( newRejectionStatus );
        BOOST_CHECK_EQUAL_COLLECTIONS( rejectedObservationIds.begin( ),
                                       rejectedObservationIds.end( ),
                                       expectedRejectedObservationIds.begin( ),
                                       expectedRejectedObservationIds.end( ) );
    }

    {
        // Fewer rejections than the maximum are left untouched
        outlierRejection.setRejectionStatus( noObservationRejected );
        std::vector< bool > newRejectionStatus = { false, true, false, true, false, true, false, false, false, false };
        outlierRejection.applyMaximumRejectedFraction( newRejectionStatus, chiSquaredPerObservation );

        const std::vector< unsigned int > expectedRejectedObservationIds = { 1, 3, 5 };
        const std::vector< unsigned int > rejectedObservationIds = getRejectedObservationIds( newRejectionStatus );
        BOOST_CHECK_EQUAL_COLLECTIONS( rejectedObservationIds.begin( ),
                                       rejectedObservationIds.end( ),
                                       expectedRejectedObservationIds.begin( ),
                                       expectedRejectedObservationIds.end( ) );
    }

    {
        // Five observations with an equal chi-squared meet the criterion, of which one has to be put back into the
        // fit. The choice is made by observation id, so that it does not depend on the order in which the
        // observations happen to be stored.
        outlierRejection.setRejectionStatus( noObservationRejected );
        const std::vector< double > equalChiSquaredPerObservation( numberOfObservations, 50.0 );
        std::vector< bool > newRejectionStatus = { true, false, true, false, true, false, true, false, true, false };
        outlierRejection.applyMaximumRejectedFraction( newRejectionStatus, equalChiSquaredPerObservation );

        const std::vector< unsigned int > expectedRejectedObservationIds = { 0, 2, 4, 6 };
        const std::vector< unsigned int > rejectedObservationIds = getRejectedObservationIds( newRejectionStatus );
        BOOST_CHECK_EQUAL_COLLECTIONS( rejectedObservationIds.begin( ),
                                       rejectedObservationIds.end( ),
                                       expectedRejectedObservationIds.begin( ),
                                       expectedRejectedObservationIds.end( ) );
    }
}


//! Test computation of Chi-squared for an arbitrary case. Expected values computed seperately
//! using a Python implementation of Carpino outlier rejection
BOOST_AUTO_TEST_CASE( test_ChiSquaredCalculation )
{
    std::shared_ptr<CarpinoOutlierRejectionSettings> settings =
        std::make_shared<CarpinoOutlierRejectionSettings>(9.0,8.0);
    TestableCarpinoOutlierRejection outlierRejection(
        settings, createRangeDataset( 10 ));

    Eigen::MatrixXd covariance(2,2);
    covariance << 1E-9, 1E-10, 1E-10, 1E-9;

    Eigen::VectorXd residuals(2);
    residuals << 2E-9, 1E-9;

    Eigen::MatrixXd partialsMatrix = Eigen::MatrixXd::Ones(2,6) * 1E-5;
    partialsMatrix(0,0) = 2E-5;

    Eigen::MatrixXd parameterCovariance = Eigen::MatrixXd::Ones(6,6) * 1E-5;
    parameterCovariance(0,0) = 2E-5;

    Eigen::VectorXd parameterCorrection = Eigen::VectorXd::Ones(6)*1E-5;
    parameterCorrection(0) = 2E-5;

    // Check case where observation is rejected
    const double chiSquaredRejectedObservation = outlierRejection.computeChiSquared(
        partialsMatrix,
        residuals,
        parameterCorrection,
        parameterCovariance,
        covariance,
        true);
    const double expectedChiSquaredRejectedObservation = 1.246383124908181e-09;

    BOOST_CHECK_CLOSE_FRACTION(chiSquaredRejectedObservation, expectedChiSquaredRejectedObservation, 1E-15);

    // check case where observation is accepted
    const double chiSquaredAcceptedObservation = outlierRejection.computeChiSquared(
        partialsMatrix,
        residuals,
        parameterCorrection,
        parameterCovariance,
        covariance,
        false);

    const double expectedChiSquaredAcceptedObservation = 1.246546181332082e-09;
    BOOST_CHECK_CLOSE(chiSquaredAcceptedObservation, expectedChiSquaredAcceptedObservation, 1E-15);

}

BOOST_AUTO_TEST_SUITE_END( )  // carpino_outlier_rejection

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////  Tests of the simple algorithm, which compares residuals against a threshold  ////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE( simple_outlier_rejection )

//! Check that the settings are validated on creation, and that they return what was provided.
BOOST_AUTO_TEST_CASE( test_SimpleOutlierRejectionSettings )
{
    const std::shared_ptr< OutlierRejectionSettings > settings = simpleOutlierRejectionSettings( 10.0, 2, false );

    const std::shared_ptr< SimpleOutlierRejectionSettings > simpleSettings =
            std::dynamic_pointer_cast< SimpleOutlierRejectionSettings >( settings );
    BOOST_REQUIRE( simpleSettings != nullptr );
    BOOST_CHECK( ( simpleSettings->getOutlierRejectionType( ) == OutlierRejectionType::simple_outlier_rejection ) );
    BOOST_CHECK_CLOSE_FRACTION( simpleSettings->getMaximumAllowedResidualValue( ), 10.0, 1.0E-15 );
    BOOST_CHECK_EQUAL( simpleSettings->getFirstIterationWithRejection( ), 2 );
    BOOST_CHECK_EQUAL( simpleSettings->getAllowRestore( ), false );

    // Observations are restored and rejection starts in the second iteration, unless the user asks otherwise
    const std::shared_ptr< SimpleOutlierRejectionSettings > defaultSettings =
            std::dynamic_pointer_cast< SimpleOutlierRejectionSettings >( simpleOutlierRejectionSettings( 10.0 ) );
    BOOST_CHECK_EQUAL( defaultSettings->getFirstIterationWithRejection( ), 1 );
    BOOST_CHECK_EQUAL( defaultSettings->getAllowRestore( ), true );

    // A threshold that is not a positive value rejects every observation, and an iteration cannot be negative
    BOOST_CHECK_THROW( simpleOutlierRejectionSettings( 0.0 ), std::runtime_error );
    BOOST_CHECK_THROW( simpleOutlierRejectionSettings( -1.0 ), std::runtime_error );
    BOOST_CHECK_THROW( simpleOutlierRejectionSettings( 10.0, -1 ), std::runtime_error );

    // A threshold per observable type is validated in the same way, and may not be left empty
    const std::map< ObservableType, double > thresholdMap = { { one_way_range, 6.0 }, { angular_position, 15.0 } };
    const std::shared_ptr< SimpleOutlierRejectionSettings > mapSettings =
            std::dynamic_pointer_cast< SimpleOutlierRejectionSettings >( simpleOutlierRejectionSettings( thresholdMap ) );
    BOOST_REQUIRE( mapSettings != nullptr );
    BOOST_CHECK_EQUAL( mapSettings->getMaximumAllowedResidualValueMap( ).size( ), 2 );
    BOOST_CHECK_CLOSE_FRACTION( mapSettings->getMaximumAllowedResidualValueMap( ).at( one_way_range ), 6.0, 1.0E-15 );

    BOOST_CHECK_THROW( simpleOutlierRejectionSettings( std::map< ObservableType, double >( ) ), std::runtime_error );
    BOOST_CHECK_THROW( simpleOutlierRejectionSettings( std::map< ObservableType, double >( { { one_way_range, -1.0 } } ) ),
                       std::runtime_error );

    // The single threshold is not used when a threshold per observable type was provided
    BOOST_CHECK( mapSettings->getMaximumAllowedResidualValueMap( ).empty( ) == false );
    BOOST_CHECK( simpleSettings->getMaximumAllowedResidualValueMap( ).empty( ) == true );
}

//! Check that the outlier rejection object created from the settings is of the corresponding type.
BOOST_AUTO_TEST_CASE( test_SimpleOutlierRejectionCreation )
{
    const std::shared_ptr< ObservationDataset< double, double > > dataset = createResidualDataset( { 4.0, 5.0, 8.0 } );

    const std::shared_ptr< OutlierRejection< double, double > > outlierRejection =
            createOutlierRejection< double, double >( simpleOutlierRejectionSettings( 10.0 ), dataset );
    BOOST_REQUIRE( outlierRejection != nullptr );
    BOOST_CHECK( ( std::dynamic_pointer_cast< SimpleOutlierRejection< double, double > >( outlierRejection ) != nullptr ) );
    BOOST_CHECK( ( outlierRejection->getOutlierRejectionType( ) == OutlierRejectionType::simple_outlier_rejection ) );

    BOOST_CHECK_EQUAL( outlierRejection->getRejectionStatus( ).size( ), dataset->getNumberOfObservations( ) );
    BOOST_CHECK_EQUAL( outlierRejection->getNumberOfRejectedObservations( ), 0 );

    // Unlike the algorithm of Carpino et al. (2003), this algorithm compares residuals against a threshold that the
    // user provides, and therefore does not need weights that represent the uncertainty of the observations
    std::shared_ptr< ObservationDataset< double, double > > datasetWithoutWeights = std::make_shared< ObservationDataset< double, double > >( );
    datasetWithoutWeights->addObservationSet( one_way_range,
                                              createTestLinkDefinition( "StationA" ),
                                              { ( Eigen::VectorXd( 1 ) << 10.0 ).finished( ) },
                                              { 1.0 },
                                              receiver );
    BOOST_CHECK_NO_THROW( ( createOutlierRejection< double, double >( simpleOutlierRejectionSettings( 10.0 ), datasetWithoutWeights ) ) );
}

//! Check that observations are rejected on the size of their residual, whatever the sign of that residual is.
BOOST_AUTO_TEST_CASE( test_SimpleRejectionByResidualSize )
{
    const std::shared_ptr< ObservationDataset< double, double > > dataset = createResidualDataset( { 4.0, -20.0, 8.0, 30.0 } );
    const TestIterationData iterationData( dataset );

    // Rejection is allowed from the first iteration onwards, so that this test case is not concerned with iterations
    const std::shared_ptr< OutlierRejection< double, double > > outlierRejection =
            createOutlierRejection< double, double >( simpleOutlierRejectionSettings( 10.0, 0 ), dataset );
    outlierRejection->updateRejectionStatus( iterationData.getInput( 0 ) );

    // Observation 1 has a residual of -20, which is further from zero than the threshold of 10 and is therefore
    // rejected, even though the value itself is smaller than the threshold
    const std::vector< unsigned int > expectedRejectedObservationIds = { 1, 3 };
    const std::vector< unsigned int > rejectedObservationIds = getRejectedObservationIds( outlierRejection->getRejectionStatus( ) );
    BOOST_CHECK_EQUAL_COLLECTIONS( rejectedObservationIds.begin( ),
                                   rejectedObservationIds.end( ),
                                   expectedRejectedObservationIds.begin( ),
                                   expectedRejectedObservationIds.end( ) );

    // The decision must be applied to the dataset, which is what the estimation reads
    BOOST_CHECK_EQUAL( dataset->getObservationRow( 0 ).isActive_, true );
    BOOST_CHECK_EQUAL( dataset->getObservationRow( 1 ).isActive_, false );
    BOOST_CHECK_EQUAL( dataset->getObservationRow( 2 ).isActive_, true );
    BOOST_CHECK_EQUAL( dataset->getObservationRow( 3 ).isActive_, false );

    // An observation exactly at the threshold is kept, since the comparison is strict
    const std::shared_ptr< ObservationDataset< double, double > > thresholdDataset = createResidualDataset( { 10.0, -10.0 } );
    const TestIterationData thresholdIterationData( thresholdDataset );
    const std::shared_ptr< OutlierRejection< double, double > > thresholdOutlierRejection =
            createOutlierRejection< double, double >( simpleOutlierRejectionSettings( 10.0, 0 ), thresholdDataset );
    thresholdOutlierRejection->updateRejectionStatus( thresholdIterationData.getInput( 0 ) );
    BOOST_CHECK_EQUAL( thresholdOutlierRejection->getNumberOfRejectedObservations( ), 0 );
}

//! Check that observations are only rejected from the iteration at which the settings allow rejection onwards.
BOOST_AUTO_TEST_CASE( test_SimpleFirstIterationWithRejection )
{
    const int firstIterationWithRejection = 2;

    const std::shared_ptr< ObservationDataset< double, double > > dataset = createResidualDataset( { 4.0, -20.0, 8.0, 30.0 } );
    const TestIterationData iterationData( dataset );

    const std::shared_ptr< OutlierRejection< double, double > > outlierRejection = createOutlierRejection< double, double >(
            simpleOutlierRejectionSettings( 10.0, firstIterationWithRejection ), dataset );

    // The first iterations use the a priori parameter values and can have large residuals for every observation, so
    // the status of all observations is left unchanged until the requested iteration is reached
    for( int iterationNumber = 0; iterationNumber < firstIterationWithRejection; iterationNumber++ )
    {
        outlierRejection->updateRejectionStatus( iterationData.getInput( iterationNumber ) );
        BOOST_CHECK_EQUAL( outlierRejection->getNumberOfRejectedObservations( ), 0 );
    }

    outlierRejection->updateRejectionStatus( iterationData.getInput( firstIterationWithRejection ) );
    BOOST_CHECK_EQUAL( outlierRejection->getNumberOfRejectedObservations( ), 2 );
}

//! Check that an observation whose residual drops below the threshold is recovered only if the settings allow it.
BOOST_AUTO_TEST_CASE( test_SimpleRestore )
{
    // Residuals of the first iteration, of which the second observation is an outlier, and residuals of a later
    // iteration, in which that observation fits the estimated parameters again
    const std::vector< double > outlyingResiduals = { 4.0, 30.0, 8.0, 6.0 };
    const std::vector< Eigen::VectorXd > improvedResiduals = { ( Eigen::VectorXd( 1 ) << 4.0 ).finished( ),
                                                               ( Eigen::VectorXd( 1 ) << 2.0 ).finished( ),
                                                               ( Eigen::VectorXd( 1 ) << 8.0 ).finished( ),
                                                               ( Eigen::VectorXd( 1 ) << 6.0 ).finished( ) };

    {
        const bool allowRestore = true;
        const std::shared_ptr< ObservationDataset< double, double > > dataset = createResidualDataset( outlyingResiduals );
        const std::shared_ptr< OutlierRejection< double, double > > outlierRejection =
                createOutlierRejection< double, double >( simpleOutlierRejectionSettings( 10.0, 0, allowRestore ), dataset );

        const TestIterationData firstIterationData( dataset );
        outlierRejection->updateRejectionStatus( firstIterationData.getInput( 0 ) );
        BOOST_REQUIRE_EQUAL( outlierRejection->getRejectionStatus( ).at( 1 ), true );

        // The residual of the observation is now below the threshold, so it is put back into the fit
        dataset->setResidualsForSet( 0, improvedResiduals );
        const TestIterationData secondIterationData( dataset );
        outlierRejection->updateRejectionStatus( secondIterationData.getInput( 1 ) );

        BOOST_CHECK_EQUAL( outlierRejection->getNumberOfRejectedObservations( ), 0 );
        BOOST_CHECK_EQUAL( dataset->getObservationRow( 1 ).isActive_, true );
    }

    {
        const bool allowRestore = false;
        const std::shared_ptr< ObservationDataset< double, double > > dataset = createResidualDataset( outlyingResiduals );
        const std::shared_ptr< OutlierRejection< double, double > > outlierRejection =
                createOutlierRejection< double, double >( simpleOutlierRejectionSettings( 10.0, 0, allowRestore ), dataset );

        const TestIterationData firstIterationData( dataset );
        outlierRejection->updateRejectionStatus( firstIterationData.getInput( 0 ) );
        BOOST_REQUIRE_EQUAL( outlierRejection->getRejectionStatus( ).at( 1 ), true );

        // Once rejected, the observation stays out of the fit, however small its residual becomes
        dataset->setResidualsForSet( 0, improvedResiduals );
        const TestIterationData secondIterationData( dataset );
        outlierRejection->updateRejectionStatus( secondIterationData.getInput( 1 ) );

        BOOST_CHECK_EQUAL( outlierRejection->getNumberOfRejectedObservations( ), 1 );
        BOOST_CHECK_EQUAL( outlierRejection->getRejectionStatus( ).at( 1 ), true );
        BOOST_CHECK_EQUAL( dataset->getObservationRow( 1 ).isActive_, false );
    }
}

//! Check that a threshold can be provided per observable type, which is what a dataset of several observables needs.
BOOST_AUTO_TEST_CASE( test_SimpleThresholdPerObservableType )
{
    const std::shared_ptr< ObservationDataset< double, double > > dataset = createTestDataset( );
    const TestIterationData iterationData( dataset );

    // The residuals of the range observations are 4, 5 and 8, and those of the angular position observations are
    // (2, 10) and (20, 40). Each observable type is compared against its own threshold.
    const std::map< ObservableType, double > thresholdMap = { { one_way_range, 6.0 }, { angular_position, 15.0 } };
    const std::shared_ptr< OutlierRejection< double, double > > outlierRejection =
            createOutlierRejection< double, double >( simpleOutlierRejectionSettings( thresholdMap, 0 ), dataset );
    outlierRejection->updateRejectionStatus( iterationData.getInput( 0 ) );

    // Observation 2 exceeds the threshold of the range observations, and observation 4 that of the angular position
    // observations. Observation 3 has a component of 10, which would be an outlier under the range threshold but is
    // not under its own.
    const std::vector< unsigned int > expectedRejectedObservationIds = { 2, 4 };
    const std::vector< unsigned int > rejectedObservationIds = getRejectedObservationIds( outlierRejection->getRejectionStatus( ) );
    BOOST_CHECK_EQUAL_COLLECTIONS( rejectedObservationIds.begin( ),
                                   rejectedObservationIds.end( ),
                                   expectedRejectedObservationIds.begin( ),
                                   expectedRejectedObservationIds.end( ) );

    // An observable type of the dataset that has no threshold cannot be checked, which must be reported
    const std::map< ObservableType, double > incompleteThresholdMap = { { one_way_range, 6.0 } };
    const std::shared_ptr< OutlierRejection< double, double > > incompleteOutlierRejection =
            createOutlierRejection< double, double >( simpleOutlierRejectionSettings( incompleteThresholdMap, 0 ), createTestDataset( ) );
    BOOST_CHECK_THROW( incompleteOutlierRejection->updateRejectionStatus( iterationData.getInput( 0 ) ), std::runtime_error );
}

//! Check that observations that were excluded before the estimation are left out of the fit.
BOOST_AUTO_TEST_CASE( test_SimpleDisabledObservationsAreNotRestored )
{
    const std::shared_ptr< ObservationDataset< double, double > > dataset = createResidualDataset( { 4.0, 2.0, 8.0, 30.0 } );

    // The user excludes the first observation before the estimation starts, for a reason of their own that the
    // algorithm knows nothing about
    dataset->rejectObservations( ObservationSelectionCondition< double, double >::timeLessThan( 1.5 ) );
    BOOST_REQUIRE_EQUAL( dataset->getObservationRow( 0 ).isActive_, false );

    const TestIterationData iterationData( dataset );
    const bool allowRestore = true;
    const std::shared_ptr< OutlierRejection< double, double > > outlierRejection =
            createOutlierRejection< double, double >( simpleOutlierRejectionSettings( 10.0, 0, allowRestore ), dataset );
    BOOST_CHECK_EQUAL( outlierRejection->getNumberOfRejectedObservations( ), 1 );

    outlierRejection->updateRejectionStatus( iterationData.getInput( 0 ) );

    // The residual of that observation is below the threshold, but restoring it is not the algorithm's decision to
    // make. Only observation 3, which the algorithm rejected itself, is added to the excluded observations.
    const std::vector< unsigned int > expectedRejectedObservationIds = { 0, 3 };
    const std::vector< unsigned int > rejectedObservationIds = getRejectedObservationIds( outlierRejection->getRejectionStatus( ) );
    BOOST_CHECK_EQUAL_COLLECTIONS( rejectedObservationIds.begin( ),
                                   rejectedObservationIds.end( ),
                                   expectedRejectedObservationIds.begin( ),
                                   expectedRejectedObservationIds.end( ) );
    BOOST_CHECK_EQUAL( dataset->getObservationRow( 0 ).isActive_, false );
}

BOOST_AUTO_TEST_SUITE_END( )  // simple_outlier_rejection

BOOST_AUTO_TEST_SUITE_END( )  // test_outlier_rejection

}  // namespace unit_tests
}  // namespace tudat
