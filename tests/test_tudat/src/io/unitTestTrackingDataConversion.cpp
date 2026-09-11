#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/simulation/estimation_setup/createObservationCollection.h"

namespace tudat
{
namespace unit_tests
{
using namespace observation_models;

BOOST_AUTO_TEST_SUITE( test_tracking_data_conversion )

BOOST_AUTO_TEST_CASE( testVectorWeightSetterShapeAndExceptionSafety )
{
    const LinkEnds linkEnds = { { transmitter, LinkEndId( "433", "" ) }, { receiver, LinkEndId( "Earth", "500" ) } };
    for( const auto type : { one_way_range, angular_position, position_observable } )
    {
        const int size = getObservableSize( type );
        for( const unsigned int count : { 1, 2 } )
        {
            std::vector< Eigen::VectorXd > observations( count, Eigen::VectorXd::Ones( size ) );
            std::vector< double > epochs( count );
            for( unsigned int i = 0; i < count; ++i )
            {
                epochs[ i ] = i;
            }
            SingleObservationSet< double, double > observationSet( type, LinkDefinition( linkEnds ), observations, epochs, receiver );
            const std::vector< Eigen::VectorXd > weights( count, Eigen::VectorXd::Constant( size, 3.0 ) );
            BOOST_CHECK_NO_THROW( observationSet.setWeights( weights ) );
            BOOST_CHECK( observationSet.getWeightsVector( ).isApprox( Eigen::VectorXd::Constant( count * size, 3.0 ) ) );
            auto invalid = weights;
            invalid.pop_back( );
            BOOST_CHECK_THROW( observationSet.setWeights( invalid ), std::runtime_error );
            invalid = weights;
            invalid.push_back( weights.front( ) );
            BOOST_CHECK_THROW( observationSet.setWeights( invalid ), std::runtime_error );
            invalid = weights;
            invalid.front( ).setConstant( 7.0 );
            invalid.back( ) = Eigen::VectorXd::Ones( size + 1 );
            BOOST_CHECK_THROW( observationSet.setWeights( invalid ), std::runtime_error );
            BOOST_CHECK( observationSet.getWeightsVector( ).isApprox( Eigen::VectorXd::Constant( count * size, 3.0 ) ) );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )
}  // namespace unit_tests
}  // namespace tudat
