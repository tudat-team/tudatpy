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
#include <memory>
#include <sstream>
#include <vector>

#include <boost/test/unit_test.hpp>
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

void compareObservationCollections(
    observation_models::ObservationCollection< double, double >& left,
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
        BOOST_CHECK_CLOSE_FRACTION( leftObservationVector( i ), rightObservationVector( i ), std::numeric_limits< double >::epsilon( ) * 10.0 );
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
            BOOST_CHECK_CLOSE_FRACTION( leftObservationTimes.at( i ).at( j ), rightObservationTimes.at( i ).at( j ),
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
            BOOST_CHECK_CLOSE_FRACTION( leftObservationSets.at( i )( j ), rightObservationSets.at( i )( j ),
                                        std::numeric_limits< double >::epsilon( ) * 10.0 );
        }
    }

    BOOST_CHECK_EQUAL( left.getTimeBounds( ).first, right.getTimeBounds( ).first );
    BOOST_CHECK_EQUAL( left.getTimeBounds( ).second, right.getTimeBounds( ).second );
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_ObservationCollection_serialization )

BOOST_AUTO_TEST_CASE( test_ObservationCollectionSerialization )
{
    using observation_models::ObservationCollection;
    using observation_models::one_way_range;

    const auto firstSet = createSingleObservationSet( one_way_range, "TxA", "RxA", { 0.0, 10.0, 20.0 }, { 1.0, 2.0, 3.0 } );
    const auto secondSet = createSingleObservationSet( one_way_range, "TxB", "RxB", { 5.0, 15.0 }, { 4.0, 5.0 } );

    const std::vector< std::shared_ptr< observation_models::SingleObservationSet< double, double > > > observationSets =
            { firstSet, secondSet };

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

    Eigen::VectorXd originalConcatenatedObservations = originalCollection.getConcatenatedObservations( );
    Eigen::VectorXd originalConcatenatedResiduals = originalCollection.getConcatenatedResiduals( );
    Eigen::VectorXd newResiduals = Eigen::VectorXd::LinSpaced( originalConcatenatedResiduals.size( ), 10.0, 10.0 + originalConcatenatedResiduals.size( ) - 1 );

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

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat