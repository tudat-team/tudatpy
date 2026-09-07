/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <memory>
#include <sstream>
#include <vector>

#include <boost/test/included/unit_test.hpp>
#include <Eigen/Core>

#include <cereal/archives/binary.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"

namespace tudat
{
namespace unit_tests
{

namespace
{

std::shared_ptr< observation_models::SingleObservationSet< double, double > > createSingleObservationSet(
        const observation_models::ObservableType observableType,
        const std::string& transmitterName,
        const std::string& receiverName,
        const std::vector< double >& observationTimes,
        const std::vector< double >& observationValues )
{
    std::map< observation_models::LinkEndType, std::pair< std::string, std::string > > linkEndsMap;
    linkEndsMap[ observation_models::transmitter ] = std::make_pair( transmitterName, std::string( ) );
    linkEndsMap[ observation_models::receiver ] = std::make_pair( receiverName, std::string( ) );
    observation_models::LinkDefinition linkEnds( linkEndsMap );

    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observations;
    observations.reserve( observationValues.size( ) );
    for( const double observationValue : observationValues )
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > observation( 1, 1 );
        observation( 0 ) = observationValue;
        observations.push_back( observation );
    }

    return std::make_shared< observation_models::SingleObservationSet< double, double > >(
            observableType, linkEnds, observations, observationTimes, observation_models::receiver );
}

void compareObservationCollections( observation_models::ObservationCollection< double, double >& left,
                                    observation_models::ObservationCollection< double, double >& right )
{
    BOOST_CHECK_EQUAL( left.getTotalObservableSize( ), right.getTotalObservableSize( ) );
    BOOST_CHECK_EQUAL( left.getConcatenatedObservationTimes( ).size( ), right.getConcatenatedObservationTimes( ).size( ) );
    BOOST_CHECK_EQUAL( left.getObservationsSets( ).size( ), right.getObservationsSets( ).size( ) );

    const Eigen::VectorXd leftObservationVector = left.getConcatenatedObservations( );
    const Eigen::VectorXd rightObservationVector = right.getConcatenatedObservations( );
    BOOST_REQUIRE_EQUAL( leftObservationVector.size( ), rightObservationVector.size( ) );
    for( int i = 0; i < leftObservationVector.size( ); ++i )
    {
        BOOST_CHECK_CLOSE_FRACTION(
                leftObservationVector( i ), rightObservationVector( i ), std::numeric_limits< double >::epsilon( ) * 10.0 );
    }

    const std::vector< double > leftTimes = left.getConcatenatedDoubleObservationTimes( );
    const std::vector< double > rightTimes = right.getConcatenatedDoubleObservationTimes( );
    BOOST_REQUIRE_EQUAL( leftTimes.size( ), rightTimes.size( ) );
    for( unsigned int i = 0; i < leftTimes.size( ); ++i )
    {
        BOOST_CHECK_CLOSE_FRACTION( leftTimes.at( i ), rightTimes.at( i ), std::numeric_limits< double >::epsilon( ) * 10.0 );
    }

    const std::vector< std::vector< double > > leftObservationTimes = left.getObservationTimesDouble( );
    const std::vector< std::vector< double > > rightObservationTimes = right.getObservationTimesDouble( );
    BOOST_REQUIRE_EQUAL( leftObservationTimes.size( ), rightObservationTimes.size( ) );
    for( unsigned int i = 0; i < leftObservationTimes.size( ); ++i )
    {
        BOOST_REQUIRE_EQUAL( leftObservationTimes.at( i ).size( ), rightObservationTimes.at( i ).size( ) );
        for( unsigned int j = 0; j < leftObservationTimes.at( i ).size( ); ++j )
        {
            BOOST_CHECK_CLOSE_FRACTION( leftObservationTimes.at( i ).at( j ),
                                        rightObservationTimes.at( i ).at( j ),
                                        std::numeric_limits< double >::epsilon( ) * 10.0 );
        }
    }

    const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > leftObservationSets = left.getObservations( );
    const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > rightObservationSets = right.getObservations( );
    BOOST_REQUIRE_EQUAL( leftObservationSets.size( ), rightObservationSets.size( ) );
    for( unsigned int i = 0; i < leftObservationSets.size( ); ++i )
    {
        BOOST_REQUIRE_EQUAL( leftObservationSets.at( i ).size( ), rightObservationSets.at( i ).size( ) );
        for( int j = 0; j < leftObservationSets.at( i ).size( ); ++j )
        {
            BOOST_CHECK_CLOSE_FRACTION(
                    leftObservationSets.at( i )( j ), rightObservationSets.at( i )( j ), std::numeric_limits< double >::epsilon( ) * 10.0 );
        }
    }

    BOOST_CHECK_EQUAL( left.getTimeBounds( ).first, right.getTimeBounds( ).first );
    BOOST_CHECK_EQUAL( left.getTimeBounds( ).second, right.getTimeBounds( ).second );
}

// Frozen writer layout from feature/data-refactor at 213fd175. This test-only
// record deliberately does not call the current facade serializer: it models
// files that users saved before the dataset backend existed.
struct BaseBranchObservationArchive {
    observation_models::ObservableType observable = observation_models::one_way_range;
    observation_models::LinkDefinition link;
    std::pair< double, double > bounds = { 1.0, 2.0 };
    std::vector< Eigen::VectorXd > observations = { Eigen::Vector1d::Constant( 10.0 ), Eigen::Vector1d::Constant( 20.0 ) };
    std::vector< double > times = { 1.0, 2.0 };
    observation_models::LinkEndType reference = observation_models::receiver;
    std::vector< Eigen::VectorXd > dependent = { Eigen::Vector1d::Constant( 100.0 ), Eigen::Vector1d::Constant( 200.0 ) };
    std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping > bookkeeping;
    std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillary;
    unsigned int count = 2;
    unsigned int dimension = 1;
    std::vector< Eigen::VectorXd > weights = { Eigen::Vector1d::Constant( 2.0 ), Eigen::Vector1d::Constant( 3.0 ) };
    std::vector< Eigen::VectorXd > residuals = { Eigen::Vector1d::Constant( 0.5 ), Eigen::Vector1d::Constant( -0.25 ) };
    std::shared_ptr< BaseBranchObservationArchive > filtered;

    BaseBranchObservationArchive( )
    {
        link[ observation_models::transmitter ] = observation_models::LinkEndId( "Earth", "" );
        link[ observation_models::receiver ] = observation_models::LinkEndId( "Vehicle", "" );
        ancillary = std::make_shared< observation_models::ObservationAncillarySimulationSettings >( );
        ancillary->setAncillaryDoubleData( observation_models::doppler_integration_time, 10.0 );
    }

    template< class Archive >
    void serialize( Archive& ar )
    {
        ar( observable,
            link,
            bounds,
            observations,
            times,
            reference,
            dependent,
            bookkeeping,
            ancillary,
            count,
            dimension,
            weights,
            residuals,
            filtered );
    }
};

}  // namespace

//! Test serialization and deserialization of SingleObservationSet
BOOST_AUTO_TEST_CASE( test_SingleObservationSetSerialization )
{
    using namespace observation_models;
    using namespace simulation_setup;

    // Create simple test data
    const ObservableType observableType = one_way_range;

    // Create link definition: transmitter -> receiver
    LinkDefinition linkEnds;
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Transmitter", "" );
    linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Receiver", "" );

    // Create observation data: 5 observations of range values
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observations;
    std::vector< double > observationTimes;

    const double startTime = 0.0;
    const double timeStep = 10.0;
    for( unsigned int i = 0; i < 5; i++ )
    {
        // Create a range observation (scalar)
        Eigen::Matrix< double, Eigen::Dynamic, 1 > obs( 1, 1 );
        obs( 0, 0 ) = 1.0e6 + i * 1000.0;  // Range values in meters
        observations.push_back( obs );
        observationTimes.push_back( startTime + i * timeStep );
    }

    // Create the original observation set
    auto dependentVariableBookkeeping = std::make_shared< ObservationDependentVariableBookkeeping >( observableType, linkEnds );
    auto ancillarySettings = std::make_shared< ObservationAncillarySimulationSettings >( );
    ancillarySettings->setAncillaryDoubleData( doppler_integration_time, 60.0 );
    std::shared_ptr< SingleObservationSet< double, double > > originalObservationSet =
            std::make_shared< SingleObservationSet< double, double > >( observableType,
                                                                        linkEnds,
                                                                        observations,
                                                                        observationTimes,
                                                                        receiver,
                                                                        std::vector< Eigen::VectorXd >( ),
                                                                        dependentVariableBookkeeping,
                                                                        ancillarySettings );

    // Serialize to binary stream
    std::stringstream serializationStream;
    {
        cereal::BinaryOutputArchive outputArchive( serializationStream );
        outputArchive( originalObservationSet );
    }

    // Deserialize from binary stream
    std::shared_ptr< SingleObservationSet< double, double > > deserializedObservationSet;
    {
        cereal::BinaryInputArchive inputArchive( serializationStream );
        inputArchive( deserializedObservationSet );
    }

    // Verify that the deserialized object is not null
    BOOST_REQUIRE( deserializedObservationSet != nullptr );
    BOOST_CHECK( originalObservationSet != deserializedObservationSet );
    BOOST_CHECK( *originalObservationSet == *deserializedObservationSet );

    // Verify observable type
    BOOST_CHECK_EQUAL( deserializedObservationSet->getObservableType( ), originalObservationSet->getObservableType( ) );

    // Verify link ends
    LinkDefinition deserializedLinkEnds = deserializedObservationSet->getLinkEnds( );
    BOOST_CHECK_EQUAL( deserializedLinkEnds[ transmitter ].getBodyName( ), linkEnds[ transmitter ].getBodyName( ) );
    BOOST_CHECK_EQUAL( deserializedLinkEnds[ receiver ].getBodyName( ), linkEnds[ receiver ].getBodyName( ) );

    // Verify reference link end
    BOOST_CHECK_EQUAL( deserializedObservationSet->getReferenceLinkEnd( ), originalObservationSet->getReferenceLinkEnd( ) );

    // Verify number of observations
    BOOST_CHECK_EQUAL( deserializedObservationSet->getNumberOfObservables( ), originalObservationSet->getNumberOfObservables( ) );

    // Verify observation times
    std::vector< double > deserializedTimes = deserializedObservationSet->getObservationTimes( );
    BOOST_CHECK_EQUAL( deserializedTimes.size( ), observationTimes.size( ) );
    for( unsigned int i = 0; i < observationTimes.size( ); i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( deserializedTimes[ i ], observationTimes[ i ], std::numeric_limits< double >::epsilon( ) * 10.0 );
    }

    // Verify observation values
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > deserializedObs = deserializedObservationSet->getObservations( );
    BOOST_CHECK_EQUAL( deserializedObs.size( ), observations.size( ) );
    for( unsigned int i = 0; i < observations.size( ); i++ )
    {
        BOOST_CHECK_EQUAL( deserializedObs[ i ].rows( ), observations[ i ].rows( ) );
        BOOST_CHECK_EQUAL( deserializedObs[ i ].cols( ), observations[ i ].cols( ) );
        for( int j = 0; j < observations[ i ].rows( ); j++ )
        {
            BOOST_CHECK_CLOSE_FRACTION(
                    deserializedObs[ i ]( j ), observations[ i ]( j ), std::numeric_limits< double >::epsilon( ) * 10.0 );
        }
    }

    // Verify single observation access (functional test)
    for( unsigned int i = 0; i < 5; i++ )
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > originalObs = originalObservationSet->getObservation( i );
        Eigen::Matrix< double, Eigen::Dynamic, 1 > deserializedObs = deserializedObservationSet->getObservation( i );

        BOOST_CHECK_CLOSE_FRACTION( deserializedObs( 0 ), originalObs( 0 ), std::numeric_limits< double >::epsilon( ) * 10.0 );
    }

    // Verify observation time access (functional test)
    for( unsigned int i = 0; i < 5; i++ )
    {
        double originalTime = originalObservationSet->getObservationTime( i );
        double deserializedTime = deserializedObservationSet->getObservationTime( i );

        BOOST_CHECK_CLOSE_FRACTION( deserializedTime, originalTime, std::numeric_limits< double >::epsilon( ) * 10.0 );
    }

    // Verify observations vector conversion
    Eigen::Matrix< double, Eigen::Dynamic, 1 > originalObsVector = originalObservationSet->getObservationsVector( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > deserializedObsVector = deserializedObservationSet->getObservationsVector( );

    BOOST_CHECK_EQUAL( originalObsVector.size( ), deserializedObsVector.size( ) );
    for( int i = 0; i < originalObsVector.size( ); i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( deserializedObsVector( i ), originalObsVector( i ), std::numeric_limits< double >::epsilon( ) * 10.0 );
    }

    // Verify time bounds
    std::pair< double, double > originalTimeBounds = originalObservationSet->getTimeBounds( );
    std::pair< double, double > deserializedTimeBounds = deserializedObservationSet->getTimeBounds( );

    BOOST_CHECK_CLOSE_FRACTION( deserializedTimeBounds.first, originalTimeBounds.first, std::numeric_limits< double >::epsilon( ) * 10.0 );
    BOOST_CHECK_CLOSE_FRACTION(
            deserializedTimeBounds.second, originalTimeBounds.second, std::numeric_limits< double >::epsilon( ) * 10.0 );
}

//! Test serialization with vector observations (3D)
BOOST_AUTO_TEST_CASE( test_SingleObservationSetSerializationVector )
{
    using namespace observation_models;
    using namespace simulation_setup;

    // Create test data for 3D position observations
    const ObservableType observableType = position_observable;

    LinkDefinition linkEnds;
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Sat", "" );
    linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Station", "" );

    // Create 3D position observations
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > observations;
    std::vector< double > observationTimes;

    for( unsigned int i = 0; i < 3; i++ )
    {
        // 3D position vector
        Eigen::Matrix< double, Eigen::Dynamic, 1 > obs( 3, 1 );
        obs( 0 ) = 1.0e6 + i * 100.0;
        obs( 1 ) = 2.0e6 + i * 200.0;
        obs( 2 ) = 3.0e6 + i * 300.0;
        observations.push_back( obs );
        observationTimes.push_back( static_cast< double >( i ) * 5.0 );
    }

    std::shared_ptr< SingleObservationSet< double, double > > originalObservationSet =
            std::make_shared< SingleObservationSet< double, double > >(
                    observableType, linkEnds, observations, observationTimes, receiver );

    // Serialize
    std::stringstream stream;
    {
        cereal::BinaryOutputArchive oa( stream );
        oa( originalObservationSet );
    }

    // Deserialize
    std::shared_ptr< SingleObservationSet< double, double > > deserializedObservationSet;
    {
        cereal::BinaryInputArchive ia( stream );
        ia( deserializedObservationSet );
    }

    // Verify
    BOOST_REQUIRE( deserializedObservationSet != nullptr );
    BOOST_CHECK_EQUAL( deserializedObservationSet->getSingleObservableSize( ), originalObservationSet->getSingleObservableSize( ) );
    BOOST_CHECK_EQUAL( deserializedObservationSet->getTotalObservationSetSize( ), originalObservationSet->getTotalObservationSetSize( ) );

    // Check individual values
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > deserializedObs = deserializedObservationSet->getObservations( );

    for( unsigned int i = 0; i < observations.size( ); i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_CLOSE_FRACTION(
                    deserializedObs[ i ]( j ), observations[ i ]( j ), std::numeric_limits< double >::epsilon( ) * 10.0 );
        }
    }
}

BOOST_AUTO_TEST_SUITE( test_ObservationCollection_serialization )

BOOST_AUTO_TEST_CASE( test_ObservationCollectionSerialization )
{
    using observation_models::ObservationCollection;
    using observation_models::one_way_range;

    const auto firstSet = createSingleObservationSet( one_way_range, "TxA", "RxA", { 0.0, 10.0, 20.0 }, { 1.0, 2.0, 3.0 } );
    const auto secondSet = createSingleObservationSet( one_way_range, "TxB", "RxB", { 5.0, 15.0 }, { 4.0, 5.0 } );

    const std::vector< std::shared_ptr< observation_models::SingleObservationSet< double, double > > > observationSets = { firstSet,
                                                                                                                           secondSet };

    ObservationCollection< double, double > originalCollection( observationSets );

    std::stringstream serializationStream;
    {
        cereal::BinaryOutputArchive outputArchive( serializationStream );
        outputArchive( originalCollection );
    }

    ObservationCollection< double, double > deserializedCollection;
    {
        cereal::BinaryInputArchive inputArchive( serializationStream );
        inputArchive( deserializedCollection );
    }

    compareObservationCollections( originalCollection, deserializedCollection );
    BOOST_CHECK( originalCollection == deserializedCollection );
    const auto& originalSets = originalCollection.getObservationsReference( );
    const auto& deserializedSets = deserializedCollection.getObservationsReference( );
    BOOST_REQUIRE( !originalSets.empty( ) );
    BOOST_REQUIRE( !deserializedSets.empty( ) );
    BOOST_REQUIRE( !originalSets.cbegin( )->second.empty( ) );
    BOOST_REQUIRE( !deserializedSets.cbegin( )->second.empty( ) );
    BOOST_REQUIRE( !originalSets.cbegin( )->second.cbegin( )->second.empty( ) );
    BOOST_REQUIRE( !deserializedSets.cbegin( )->second.cbegin( )->second.empty( ) );
    BOOST_CHECK( originalSets.cbegin( )->second.cbegin( )->second.front( ) !=
                 deserializedSets.cbegin( )->second.cbegin( )->second.front( ) );

    Eigen::VectorXd originalConcatenatedObservations = originalCollection.getConcatenatedObservations( );
    Eigen::VectorXd originalConcatenatedResiduals = originalCollection.getConcatenatedResiduals( );
    Eigen::VectorXd newResiduals =
            Eigen::VectorXd::LinSpaced( originalConcatenatedResiduals.size( ), 10.0, 10.0 + originalConcatenatedResiduals.size( ) - 1 );

    originalCollection.setResiduals( newResiduals );
    deserializedCollection.setResiduals( newResiduals );

    BOOST_CHECK_EQUAL( originalCollection.getConcatenatedResiduals( ).size( ), deserializedCollection.getConcatenatedResiduals( ).size( ) );
    for( int i = 0; i < newResiduals.size( ); ++i )
    {
        BOOST_CHECK_CLOSE_FRACTION( originalCollection.getConcatenatedResiduals( )( i ),
                                    deserializedCollection.getConcatenatedResiduals( )( i ),
                                    std::numeric_limits< double >::epsilon( ) * 10.0 );
    }

    BOOST_CHECK_EQUAL( originalCollection.getObservationVector( ).size( ), deserializedCollection.getObservationVector( ).size( ) );
    for( int i = 0; i < originalCollection.getObservationVector( ).size( ); ++i )
    {
        BOOST_CHECK_CLOSE_FRACTION( originalCollection.getObservationVector( )( i ),
                                    deserializedCollection.getObservationVector( )( i ),
                                    std::numeric_limits< double >::epsilon( ) * 10.0 );
    }
}

BOOST_AUTO_TEST_CASE( test_dataset_serialization_preserves_surviving_identity_and_sparse_weights )
{
    using namespace observation_models;
    auto dataset = std::make_shared< ObservationDataset<> >( );
    LinkDefinition link;
    link[ transmitter ] = LinkEndId( "Earth", "" );
    link[ receiver ] = LinkEndId( "Vehicle", "" );
    dataset->addObservationSet( one_way_range, link, { Eigen::Vector1d::Constant( 10.0 ) }, { 1.0 }, receiver );
    const auto angularSet = dataset->addObservationSet(
            angular_position, link, { Eigen::Vector2d( 20.0, 21.0 ), Eigen::Vector2d( 30.0, 31.0 ) }, { 2.0, 3.0 }, receiver );
    Eigen::Matrix4d weights;
    weights << 2.0, 0.1, 0.2, 0.0, 0.1, 3.0, 0.0, 0.3, 0.2, 0.0, 4.0, 0.4, 0.0, 0.3, 0.4, 5.0;
    dataset->setWeightMatrixForSet( angularSet, weights );
    dataset->removeObservations( ObservationSelectionCondition<>::timeLessThan( 2.0 ) );
    dataset->rejectObservations( ObservationSelectionCondition<>::timeBounds( 2.0, 2.0 ), "saved reason" );
    ObservationCollection<> original( dataset );
    std::stringstream stream;
    {
        cereal::BinaryOutputArchive archive( stream );
        archive( original );
    }
    ObservationCollection<> restored;
    {
        cereal::BinaryInputArchive archive( stream );
        archive( restored );
    }
    auto restoredDataset = restored.getObservationDataset( );
    BOOST_CHECK( *restoredDataset == *dataset );
    BOOST_CHECK_EQUAL( restoredDataset->getObservationRow( 1 ).rejectionReason_, "saved reason" );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            restoredDataset->createEstimationProjection( true ).getSparseWeightMatrix( ).toDense( ), weights, 1.0E-15 );
    restoredDataset->deleteRejectedObservations( );
    BOOST_CHECK_EQUAL( restoredDataset->getObservationRow( 2 ).setId_, angularSet );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( restoredDataset->createEstimationProjection( ).getSparseWeightMatrix( ).toDense( ),
                                       weights.bottomRightCorner( 2, 2 ),
                                       1.0E-15 );
    restoredDataset->addObservationsToSet( angularSet, { Eigen::Vector2d( 40.0, 41.0 ) }, { 4.0 } );
    BOOST_CHECK_EQUAL( restoredDataset->getObservationIdsForSet( angularSet ).back( ), 3 );
    BOOST_CHECK_EQUAL( dataset->getNumberOfObservations( ), 2 );
    BOOST_CHECK_EQUAL( restored.getTotalObservableSize( ), 4 );
}

BOOST_AUTO_TEST_CASE( test_serialized_overlapping_collections_preserve_shared_set_ownership )
{
    using namespace observation_models;
    const auto shared = createSingleObservationSet( one_way_range, "Earth", "Vehicle", { 1.0 }, { 10.0 } );
    const auto separate = createSingleObservationSet( one_way_range, "Mars", "Vehicle", { 2.0 }, { 20.0 } );
    using Sets = std::vector< std::shared_ptr< SingleObservationSet<> > >;
    ObservationCollection<> first( Sets{ shared } );
    ObservationCollection<> second( Sets{ separate, shared } );
    std::stringstream stream;
    {
        cereal::BinaryOutputArchive archive( stream );
        archive( first, second );
    }
    ObservationCollection<> restoredFirst, restoredSecond;
    {
        cereal::BinaryInputArchive archive( stream );
        archive( restoredFirst, restoredSecond );
    }
    restoredFirst.setConstantWeight( 7.0 );
    Eigen::Vector2d expected( 7.0, 1.0 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( restoredSecond.getConcatenatedWeights( ), expected, 1.0E-15 );
    BOOST_CHECK_EQUAL( first.getConcatenatedWeights( )( 0 ), 1.0 );
}

BOOST_AUTO_TEST_CASE( test_base_branch_binary_observation_archives_remain_readable )
{
    using namespace observation_models;
    auto legacy = std::make_shared< BaseBranchObservationArchive >( );
    legacy->filtered = std::make_shared< BaseBranchObservationArchive >( *legacy );
    legacy->filtered->count = 1;
    legacy->filtered->observations = { Eigen::Vector1d::Constant( 30.0 ) };
    legacy->filtered->times = { 3.0 };
    legacy->filtered->bounds = { 3.0, 3.0 };
    legacy->filtered->dependent = { Eigen::Vector1d::Constant( 300.0 ) };
    legacy->filtered->weights = { Eigen::Vector1d::Constant( 4.0 ) };
    legacy->filtered->residuals = { Eigen::Vector1d::Constant( -1.0 ) };

    std::stringstream singleStream;
    {
        cereal::BinaryOutputArchive archive( singleStream );
        archive( legacy );
    }
    std::shared_ptr< SingleObservationSet<> > restoredSingle;
    {
        cereal::BinaryInputArchive archive( singleStream );
        archive( restoredSingle );
    }
    Eigen::Vector2d values( 10.0, 20.0 ), weights( 2.0, 3.0 ), residuals( 0.5, -0.25 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( restoredSingle->getObservationDataset( )->getObservationVectorForSet( 0 ), values, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( restoredSingle->getWeightsVector( ), weights, 1.0E-15 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( restoredSingle->getResidualsVector( ), residuals, 1.0E-15 );
    BOOST_CHECK_EQUAL( restoredSingle->getNumberOfFilteredObservations( ), 1 );
    BOOST_CHECK_EQUAL( restoredSingle->getFilteredObservationSet( )->getWeightsVector( )( 0 ), 4.0 );
    BOOST_CHECK_EQUAL( restoredSingle->getAncillarySettings( )->getAncillaryDoubleData( doppler_integration_time ), 10.0 );
    BOOST_CHECK_EQUAL( restoredSingle->getObservationDataset( )->getDependentVariables( 1 )( 0 ), 200.0 );

    using LegacyMap = std::map< ObservableType, std::map< LinkEnds, std::vector< std::shared_ptr< BaseBranchObservationArchive > > > >;
    LegacyMap legacyCollection;
    legacyCollection[ one_way_range ][ legacy->link.linkEnds_ ] = { legacy, legacy };
    std::stringstream collectionStream;
    {
        cereal::BinaryOutputArchive archive( collectionStream );
        archive( legacyCollection );
    }
    ObservationCollection<> restoredCollection;
    {
        cereal::BinaryInputArchive archive( collectionStream );
        archive( restoredCollection );
    }
    const auto sets = restoredCollection.getSingleObservationSets( );
    BOOST_REQUIRE_EQUAL( sets.size( ), 2 );
    BOOST_CHECK( sets.front( ) == sets.back( ) );
    BOOST_CHECK_EQUAL( restoredCollection.getTotalObservableSize( ), 4 );
    Eigen::Vector4d repeatedWeights( 2.0, 3.0, 2.0, 3.0 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            restoredCollection.getObservationDataset( )->createEstimationProjection( ).getWeightVector( ), repeatedWeights, 1.0E-15 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
