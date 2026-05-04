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

#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <Eigen/Core>

#include <cereal/archives/binary.hpp>

#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation_setup/singleObservationSet.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_SingleObservationSet_serialization )

//! Test serialization and deserialization of SingleObservationSet
BOOST_AUTO_TEST_CASE( test_SingleObservationSetSerialization )
{
    using namespace observation_models;
    using namespace simulation_setup;

    // Create simple test data
    const ObservableType observableType = range;
    
    // Create link definition: transmitter -> receiver
    LinkDefinition linkEnds;
    linkEnds[ transmission_end ] = std::make_pair< std::string, std::string >( "Transmitter", "" );
    linkEnds[ reception_end ] = std::make_pair< std::string, std::string >( "Receiver", "" );

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
    std::shared_ptr< SingleObservationSet< double, double > > originalObservationSet =
            std::make_shared< SingleObservationSet< double, double > >(
                observableType,
                linkEnds,
                observations,
                observationTimes,
                reception_end  // reference link end
            );

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

    // Verify observable type
    BOOST_CHECK_EQUAL( 
        deserializedObservationSet->getObservableType( ),
        originalObservationSet->getObservableType( ) );

    // Verify link ends
    LinkDefinition deserializedLinkEnds = deserializedObservationSet->getLinkEnds( );
    BOOST_CHECK_EQUAL( 
        deserializedLinkEnds[ transmission_end ].first,
        linkEnds[ transmission_end ].first );
    BOOST_CHECK_EQUAL( 
        deserializedLinkEnds[ reception_end ].first,
        linkEnds[ reception_end ].first );

    // Verify reference link end
    BOOST_CHECK_EQUAL( 
        deserializedObservationSet->getReferenceLinkEnd( ),
        originalObservationSet->getReferenceLinkEnd( ) );

    // Verify number of observations
    BOOST_CHECK_EQUAL( 
        deserializedObservationSet->getNumberOfObservables( ),
        originalObservationSet->getNumberOfObservables( ) );

    // Verify observation times
    std::vector< double > deserializedTimes = deserializedObservationSet->getObservationTimes( );
    BOOST_CHECK_EQUAL( deserializedTimes.size( ), observationTimes.size( ) );
    for( unsigned int i = 0; i < observationTimes.size( ); i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( 
            deserializedTimes[ i ],
            observationTimes[ i ],
            std::numeric_limits< double >::epsilon( ) * 10.0 );
    }

    // Verify observation values
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > deserializedObs = 
        deserializedObservationSet->getObservations( );
    BOOST_CHECK_EQUAL( deserializedObs.size( ), observations.size( ) );
    for( unsigned int i = 0; i < observations.size( ); i++ )
    {
        BOOST_CHECK_EQUAL( deserializedObs[ i ].rows( ), observations[ i ].rows( ) );
        BOOST_CHECK_EQUAL( deserializedObs[ i ].cols( ), observations[ i ].cols( ) );
        for( int j = 0; j < observations[ i ].rows( ); j++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( 
                deserializedObs[ i ]( j ),
                observations[ i ]( j ),
                std::numeric_limits< double >::epsilon( ) * 10.0 );
        }
    }

    // Verify single observation access (functional test)
    for( unsigned int i = 0; i < 5; i++ )
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > originalObs = originalObservationSet->getObservation( i );
        Eigen::Matrix< double, Eigen::Dynamic, 1 > deserializedObs = deserializedObservationSet->getObservation( i );

        BOOST_CHECK_CLOSE_FRACTION( 
            deserializedObs( 0 ),
            originalObs( 0 ),
            std::numeric_limits< double >::epsilon( ) * 10.0 );
    }

    // Verify observation time access (functional test)
    for( unsigned int i = 0; i < 5; i++ )
    {
        double originalTime = originalObservationSet->getObservationTime( i );
        double deserializedTime = deserializedObservationSet->getObservationTime( i );

        BOOST_CHECK_CLOSE_FRACTION( 
            deserializedTime,
            originalTime,
            std::numeric_limits< double >::epsilon( ) * 10.0 );
    }

    // Verify observations vector conversion
    Eigen::Matrix< double, Eigen::Dynamic, 1 > originalObsVector = originalObservationSet->getObservationsVector( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > deserializedObsVector = deserializedObservationSet->getObservationsVector( );

    BOOST_CHECK_EQUAL( originalObsVector.size( ), deserializedObsVector.size( ) );
    for( int i = 0; i < originalObsVector.size( ); i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( 
            deserializedObsVector( i ),
            originalObsVector( i ),
            std::numeric_limits< double >::epsilon( ) * 10.0 );
    }

    // Verify time bounds
    std::pair< double, double > originalTimeBounds = originalObservationSet->getTimeBounds( );
    std::pair< double, double > deserializedTimeBounds = deserializedObservationSet->getTimeBounds( );

    BOOST_CHECK_CLOSE_FRACTION( 
        deserializedTimeBounds.first,
        originalTimeBounds.first,
        std::numeric_limits< double >::epsilon( ) * 10.0 );
    BOOST_CHECK_CLOSE_FRACTION( 
        deserializedTimeBounds.second,
        originalTimeBounds.second,
        std::numeric_limits< double >::epsilon( ) * 10.0 );
}

//! Test serialization with vector observations (3D)
BOOST_AUTO_TEST_CASE( test_SingleObservationSetSerializationVector )
{
    using namespace observation_models;
    using namespace simulation_setup;

    // Create test data for 3D position observations
    const ObservableType observableType = position_observable;
    
    LinkDefinition linkEnds;
    linkEnds[ transmission_end ] = std::make_pair< std::string, std::string >( "Sat", "" );
    linkEnds[ reception_end ] = std::make_pair< std::string, std::string >( "Station", "" );

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
                observableType,
                linkEnds,
                observations,
                observationTimes,
                reception_end
            );

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
    BOOST_CHECK_EQUAL( 
        deserializedObservationSet->getSingleObservableSize( ),
        originalObservationSet->getSingleObservableSize( ) );
    BOOST_CHECK_EQUAL( 
        deserializedObservationSet->getTotalObservationSetSize( ),
        originalObservationSet->getTotalObservationSetSize( ) );

    // Check individual values
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > deserializedObs = 
        deserializedObservationSet->getObservations( );
    
    for( unsigned int i = 0; i < observations.size( ); i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( 
                deserializedObs[ i ]( j ),
                observations[ i ]( j ),
                std::numeric_limits< double >::epsilon( ) * 10.0 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
