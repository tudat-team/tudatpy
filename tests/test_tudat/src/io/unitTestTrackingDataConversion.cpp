#define BOOST_TEST_MAIN
#include <boost/test/included/unit_test.hpp>
#include "tudat/simulation/estimation_setup/createObservationCollection.h"
#include "tudat/astro/ground_stations/transmittingFrequencies.h"
#include "tudat/astro/system_models/vehicleSystems.h"

namespace tudat
{
namespace unit_tests
{
using namespace observation_models;

std::shared_ptr< data::TrackingData<> > angularTracking( const double epoch = 10.0 )
{
    return std::make_shared< data::TrackingData<> >(
            "AngularPosition",
            data::PlainLinkDefinition{ { { "433", "" }, "transmitter" }, { { "Earth", "500" }, "receiver" } },
            std::vector< Eigen::VectorXd >{ Eigen::Vector2d( 1.0, 2.0 ) },
            std::vector< double >{ epoch },
            "receiver",
            "TDB" );
}

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

BOOST_AUTO_TEST_CASE( testOpticalMetadataIsNotSimulationAncillaryData )
{
    auto tracking = angularTracking( );
    for( const auto key : { "band", "catalog", "note2", "custom_name", "discovery", "mag", "phottype", "number" } )
    {
        tracking->addAncillarySettings( key, std::vector< std::string >{ "metadata" } );
    }
    tracking->addObservationMetadata( "observer", { "Alice" } );
    simulation_setup::SystemOfBodies bodies;
    BOOST_CHECK_NO_THROW( createSingleObservationSetFromTrackingData( tracking, bodies ) );
    tracking->addAncillarySettings( "frequency bands", std::vector< std::string >{ "X-band", "S-band" } );
    auto ancillary = getAncillarySettingsFromTrackingData( tracking );
    BOOST_CHECK( ancillary->getAncillaryDoubleVectorData( frequency_bands ) ==
                 std::vector< double >( { convertFrequencyBandToDouble( x_band ), convertFrequencyBandToDouble( s_band ) } ) );
    tracking->addAncillarySettings( "frequency bands", std::vector< std::string >{ "invalid-band" } );
    BOOST_CHECK_THROW( getAncillarySettingsFromTrackingData( tracking ), std::runtime_error );
    tracking = angularTracking( );
    tracking->addAncillarySettings( "unknown simulation setting", std::vector< std::string >{ "X-band" } );
    BOOST_CHECK_THROW( getAncillarySettingsFromTrackingData( tracking ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( testMixedCorrectedAndUncorrectedInputs )
{
    auto corrected = angularTracking( );
    corrected->setObservationCorrections( { Eigen::Vector2d( -0.1, -0.2 ) } );
    auto uncorrected = std::make_shared< data::TrackingData<> >( "OneWayRange",
                                                                 corrected->getLinkEnds( ),
                                                                 std::vector< Eigen::VectorXd >{ Eigen::VectorXd::Constant( 1, 1000.0 ) },
                                                                 std::vector< double >{ 20.0 },
                                                                 "receiver",
                                                                 "TDB" );
    simulation_setup::SystemOfBodies bodies;
    const std::vector< std::shared_ptr< data::TrackingData<> > > input = { corrected, uncorrected };
    BOOST_CHECK_NO_THROW( createObservationCollection( input, bodies, true ) );
    auto correctedSet = createSingleObservationSetFromTrackingData( corrected, bodies, true );
    auto uncorrectedSet = createSingleObservationSetFromTrackingData( uncorrected, bodies, true );
    BOOST_CHECK( correctedSet->getObservations( ).at( 0 ).isApprox( Eigen::Vector2d( 0.9, 1.8 ) ) );
    BOOST_CHECK_EQUAL( uncorrectedSet->getObservations( ).at( 0 )( 0 ), 1000.0 );
    BOOST_CHECK( corrected->getObservations( ).at( 0 ).isApprox( Eigen::Vector2d( 1.0, 2.0 ) ) );
}

BOOST_AUTO_TEST_CASE( testUnsupportedFrequencyHistoryFailsWithoutChangingEnvironment )
{
    simulation_setup::SystemOfBodies bodies;
    bodies.createEmptyBody( "Probe" );
    auto vehicle = std::make_shared< system_models::VehicleSystems >( );
    bodies.at( "Probe" )->setVehicleSystems( vehicle );
    auto history = std::make_shared< data::PiecewiseConstantFrequencySupplementaryData >( );
    history->setFrequency( 0.0, 8.4E9 );
    history->setFrequency( 60.0, 8.5E9 );
    const std::map< std::pair< std::string, std::string >, std::vector< std::shared_ptr< data::FrequencySupplementaryData > > > input = {
        { { "Probe", "" }, { history } }
    };
    BOOST_CHECK_THROW( setFrequencySupplementaryDataInBodies( bodies, input ), std::runtime_error );
    BOOST_CHECK( vehicle->getTransmittedFrequencyCalculator( ) == nullptr );
    auto existing = std::make_shared< ground_stations::ConstantFrequencyInterpolator >( 1.0E9 );
    vehicle->setTransmittedFrequencyCalculator( existing );
    BOOST_CHECK_THROW( setFrequencySupplementaryDataInBodies( bodies, input ), std::runtime_error );
    BOOST_CHECK( vehicle->getTransmittedFrequencyCalculator( ) == existing );
}

BOOST_AUTO_TEST_SUITE_END( )
}  // namespace unit_tests
}  // namespace tudat
