#define BOOST_TEST_MAIN

#include <string>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Basics/testMacros.h"

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"

#include "Astrodynamics/Bodies/vehicle.h"
#include "Astrodynamics/ObservationModels/nWayRangeRateObservationModel.h"
#include "Astrodynamics/ObservationModels/oneWayRangeObservationModel.h"
#include "Astrodynamics/ObservationModels/oneWayRangeRateObservationModel.h"
#include "Astrodynamics/ObservationModels/twoWayRangeObservationModel.h"
#include "Astrodynamics/ObservationModels/twoWayRangeRateObservationModel.h"
#include "Astrodynamics/ObservationModels/createObservationModel.h"

#include "SimulationSetup/defaultBodies.h"
#include "SimulationSetup/createBodies.h"
#include "SimulationSetup/createGroundStations.h"
#include "InputOutput/writeDataToFile.h"


namespace tudat
{
namespace unit_tests
{

//Using declarations.
using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::bodies;
using namespace tudat::simulation_setup;
using namespace tudat::ephemerides;
using namespace tudat::hardware_models;

BOOST_AUTO_TEST_SUITE( testNWayObservationModels )

NamedBodyMap creatEnvironment( )
{
    //Load spice kernels.
    std::string kernelsPath = input_output::getDataFilesRootPath( ) + "SpiceKernels/";
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    //Define setting for total number of bodies and those which need to be integrated numerically.
    //The first numberOfNumericalBodies from the bodyNames vector will be integrated numerically.

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );

    // Specify initial time
    double initialEphemerisTime = double( 1.0E7 );
    double finalEphemerisTime = double( 1.2E7 );
    double maximumTimeStep = 3600.0;

    double buffer = 10.0 * maximumTimeStep;

    // Create bodies needed in simulation
    std::map< std::string, boost::shared_ptr< Body > > bodyMap = createCelestialBodies(
                getDefaultBodySettings< double, double >(
                    bodyNames,initialEphemerisTime - buffer, finalEphemerisTime + buffer, simple ) );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Create ground stations
    std::pair< std::string, std::string > grazStation = std::pair< std::string, std::string >( "Earth", "Graz" );
    std::pair< std::string, std::string > mcDonaldStation = std::pair< std::string, std::string >( "Earth", "McDonald Observatory" );
    std::pair< std::string, std::string > yarragadeeStation = std::pair< std::string, std::string >( "Earth", "Yarragadee" );
    std::pair< std::string, std::string > hartebeesthoekStation = std::pair< std::string, std::string >( "Earth", "Hartebeesthoek" );
    std::pair< std::string, std::string > mslStation = std::pair< std::string, std::string >( "Mars", "MSL" );
    std::pair< std::string, std::string > vikingStation = std::pair< std::string, std::string >( "Mars", "Viking" );

    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.push_back( grazStation );
    groundStations.push_back( hartebeesthoekStation );
    groundStations.push_back( yarragadeeStation );
    groundStations.push_back( mcDonaldStation );
    groundStations.push_back( mslStation );
    groundStations.push_back( vikingStation );
    createGroundStations( bodyMap, groundStations );

    return bodyMap;
}

BOOST_AUTO_TEST_CASE( test_n_way_observation_models )
{
    NamedBodyMap bodyMap = creatEnvironment( );

    LinkEnds oneWayLinkEndMap;
    oneWayLinkEndMap[ transmitter ] = std::pair< std::string, std::string >( "Earth", "Graz" );
    oneWayLinkEndMap[ receiver ] = std::pair< std::string, std::string >( "Mars", "Viking" );

    boost::shared_ptr< NWayRangeObservationModel< > > emulatedOneWayRangeModel =
            boost::make_shared< NWayRangeObservationModel< > >(
                oneWayLinkEndMap, bodyMap );
    boost::shared_ptr< OneWayRangeObservationModel< > > oneWayRangeObservationModel =
            boost::make_shared< OneWayRangeObservationModel< > >(
                oneWayLinkEndMap.at( receiver ), oneWayLinkEndMap.at( transmitter ), bodyMap );

    double testTime = 1.1E7;

    double emulatedOneWayRange = emulatedOneWayRangeModel->computeObservation( testTime, transmitter );
    double directOneWayRange = oneWayRangeObservationModel->computeObservation( testTime, transmitter );

    BOOST_CHECK_CLOSE_FRACTION( emulatedOneWayRange, directOneWayRange, std::numeric_limits< double >::epsilon( ) );

    emulatedOneWayRange = emulatedOneWayRangeModel->computeObservation( testTime, receiver );

    directOneWayRange = oneWayRangeObservationModel->computeObservation( testTime, receiver );

    BOOST_CHECK_CLOSE_FRACTION( emulatedOneWayRange, directOneWayRange, std::numeric_limits< double >::epsilon( ) );

    boost::function< double( ) > countIntervalFunction = boost::lambda::constant( 1200.0 );

    boost::shared_ptr< NWayRangeRateObservationModel< > > emulatedOneWayRangeRateModel =
            boost::make_shared< NWayRangeRateObservationModel< > >(
                oneWayLinkEndMap, bodyMap, countIntervalFunction );
    boost::shared_ptr< OneWayRangeRateObservationModel< > > oneWayRangeRateObservationModel =
            boost::make_shared< OneWayRangeRateObservationModel< > >(
                oneWayLinkEndMap.at( receiver ), oneWayLinkEndMap.at( transmitter ), bodyMap, countIntervalFunction );

    double emulatedOneWayRangeRate = emulatedOneWayRangeRateModel->computeObservation( testTime, transmitter );
    double directOneWayRangeRate = oneWayRangeRateObservationModel->computeObservation( testTime, transmitter );

    BOOST_CHECK_CLOSE_FRACTION( emulatedOneWayRangeRate, directOneWayRangeRate, 10E-11 );

    emulatedOneWayRangeRate = emulatedOneWayRangeRateModel->computeObservation( testTime, receiver );
    directOneWayRangeRate = oneWayRangeRateObservationModel->computeObservation( testTime, receiver );

    BOOST_CHECK_CLOSE_FRACTION( emulatedOneWayRangeRate, directOneWayRangeRate, 10E-11 );


    LinkEnds twoWayDirectLinkEndMap, twoWayEmulatedLinkEndMap;

    twoWayDirectLinkEndMap[ transmitter ] = std::pair< std::string, std::string >( "Earth", "Graz" );
    twoWayDirectLinkEndMap[ reflector ] = std::pair< std::string, std::string >( "Mars", "Viking" );
    twoWayDirectLinkEndMap[ receiver ] = std::pair< std::string, std::string >( "Earth", "Yarragadee" );

    twoWayEmulatedLinkEndMap[ transmitter ] = std::pair< std::string, std::string >( "Earth", "Graz" );
    twoWayEmulatedLinkEndMap[ reflector1 ] = std::pair< std::string, std::string >( "Mars", "Viking" );
    twoWayEmulatedLinkEndMap[ receiver ] = std::pair< std::string, std::string >( "Earth", "Yarragadee" );

    boost::shared_ptr< NWayRangeObservationModel< > > emulatedTwoWayRangeModel =
            boost::make_shared< NWayRangeObservationModel< > >(
                twoWayEmulatedLinkEndMap, bodyMap );
    boost::shared_ptr< TwoWayRangeObservationModel< > > twoWayRangeObservationModel =
            boost::make_shared< TwoWayRangeObservationModel< > >(
                twoWayDirectLinkEndMap.at( transmitter ), twoWayDirectLinkEndMap.at( reflector ),
                twoWayDirectLinkEndMap.at( receiver ), bodyMap );

    double emulatedTwoWayRange = emulatedTwoWayRangeModel->computeObservation( testTime, transmitter );
    double directTwoWayRange = twoWayRangeObservationModel->computeObservation( testTime, transmitter );

    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayRange, directTwoWayRange, std::numeric_limits< double >::epsilon( ) );

    emulatedTwoWayRange = emulatedTwoWayRangeModel->computeObservation( testTime, receiver );
    directTwoWayRange = twoWayRangeObservationModel->computeObservation( testTime, receiver );

    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayRange, directTwoWayRange, std::numeric_limits< double >::epsilon( ) );

    emulatedTwoWayRange = emulatedTwoWayRangeModel->computeObservation( testTime, reflector1 );
    directTwoWayRange = twoWayRangeObservationModel->computeObservation( testTime, reflector );

    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayRange, directTwoWayRange, std::numeric_limits< double >::epsilon( ) );

    boost::shared_ptr< NWayRangeRateObservationModel< > > emulatedTwoWayRangeRateModel =
            boost::make_shared< NWayRangeRateObservationModel< > >(
                twoWayEmulatedLinkEndMap, bodyMap, countIntervalFunction );
    boost::shared_ptr< TwoWayRangeRateObservationModel< > > twoWayRangeRateObservationModel =
            boost::make_shared< TwoWayRangeRateObservationModel< > >(
                twoWayDirectLinkEndMap.at( transmitter ), twoWayDirectLinkEndMap.at( reflector ),
                twoWayDirectLinkEndMap.at( receiver ), bodyMap, countIntervalFunction );

    std::vector< double > emulatedTwoWayrangeRateLinkEndTimes;
    std::vector< basic_mathematics::Vector6d > emulatedTwoWayrangeRateLinkEndStates;

    std::vector< double > directTwoWayrangeRateLinkEndTimes;
    std::vector< basic_mathematics::Vector6d > directTwoWayrangeRateLinkEndStates;

    double emulatedTwoWayRangeRate = emulatedTwoWayRangeRateModel->computeObservationAndLinkEndData(
                testTime, transmitter, emulatedTwoWayrangeRateLinkEndTimes, emulatedTwoWayrangeRateLinkEndStates );

    double directTwoWayRangeRate = twoWayRangeRateObservationModel->computeObservationAndLinkEndData(
                testTime, transmitter, directTwoWayrangeRateLinkEndTimes, directTwoWayrangeRateLinkEndStates );

    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayRangeRate, directTwoWayRangeRate, 2.0E-11 );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 0 ], directTwoWayrangeRateLinkEndTimes[ 0 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 1 ], directTwoWayrangeRateLinkEndTimes[ 1 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 2 ], directTwoWayrangeRateLinkEndTimes[ 1 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 3 ], directTwoWayrangeRateLinkEndTimes[ 2 ],
                                std::numeric_limits< double >::epsilon( ) );


    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 4 ], directTwoWayrangeRateLinkEndTimes[ 3 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 5 ], directTwoWayrangeRateLinkEndTimes[ 4 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 6 ], directTwoWayrangeRateLinkEndTimes[ 4 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 7 ], directTwoWayrangeRateLinkEndTimes[ 5 ],
                                std::numeric_limits< double >::epsilon( ) );

    emulatedTwoWayRangeRate = emulatedTwoWayRangeRateModel->computeObservationAndLinkEndData(
                testTime, receiver, emulatedTwoWayrangeRateLinkEndTimes, emulatedTwoWayrangeRateLinkEndStates );
    directTwoWayRangeRate = twoWayRangeRateObservationModel->computeObservationAndLinkEndData(
                testTime, receiver, directTwoWayrangeRateLinkEndTimes, directTwoWayrangeRateLinkEndStates );

    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayRangeRate, directTwoWayRangeRate, 1.0E-11 );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 0 ], directTwoWayrangeRateLinkEndTimes[ 0 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 1 ], directTwoWayrangeRateLinkEndTimes[ 1 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 2 ], directTwoWayrangeRateLinkEndTimes[ 1 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 3 ], directTwoWayrangeRateLinkEndTimes[ 2 ],
                                std::numeric_limits< double >::epsilon( ) );

    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 4 ], directTwoWayrangeRateLinkEndTimes[ 3 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 5 ], directTwoWayrangeRateLinkEndTimes[ 4 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 6 ], directTwoWayrangeRateLinkEndTimes[ 4 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 7 ], directTwoWayrangeRateLinkEndTimes[ 5 ],
                                std::numeric_limits< double >::epsilon( ) );


    emulatedTwoWayRangeRate = emulatedTwoWayRangeRateModel->computeObservationAndLinkEndData(
                testTime, reflector1, emulatedTwoWayrangeRateLinkEndTimes, emulatedTwoWayrangeRateLinkEndStates  );
    directTwoWayRangeRate = twoWayRangeRateObservationModel->computeObservationAndLinkEndData(
                testTime, reflector, directTwoWayrangeRateLinkEndTimes, directTwoWayrangeRateLinkEndStates );

    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayRangeRate, directTwoWayRangeRate, 1.0E-11 );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 0 ], directTwoWayrangeRateLinkEndTimes[ 0 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 1 ], directTwoWayrangeRateLinkEndTimes[ 1 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 2 ], directTwoWayrangeRateLinkEndTimes[ 1 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 3 ], directTwoWayrangeRateLinkEndTimes[ 2 ],
                                std::numeric_limits< double >::epsilon( ) );

    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 4 ], directTwoWayrangeRateLinkEndTimes[ 3 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 5 ], directTwoWayrangeRateLinkEndTimes[ 4 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 6 ], directTwoWayrangeRateLinkEndTimes[ 4 ],
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( emulatedTwoWayrangeRateLinkEndTimes[ 7 ], directTwoWayrangeRateLinkEndTimes[ 5 ],
                                std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_CASE( test_delayed_n_way_observation_models )
{
    NamedBodyMap bodyMap = creatEnvironment( );

    LinkEnds nWayLinkEndMap;
    std::vector< LinkEnds > oneWayLinkEnds;
    oneWayLinkEnds.resize( 3 );

    nWayLinkEndMap[ transmitter ] = std::pair< std::string, std::string >( "Earth", "Graz" );
    nWayLinkEndMap[ reflector1 ] = std::pair< std::string, std::string >( "Mars", "Viking" );
    nWayLinkEndMap[ reflector2 ] = std::pair< std::string, std::string >( "Earth", "Yarragadee" );
    nWayLinkEndMap[ receiver ] = std::pair< std::string, std::string >( "Mars", "MSL" );

    oneWayLinkEnds[ 0 ][ transmitter ] = std::pair< std::string, std::string >( "Earth", "Graz" );
    oneWayLinkEnds[ 0 ][ receiver ] = std::pair< std::string, std::string >( "Mars", "Viking" );

    oneWayLinkEnds[ 1 ][ transmitter ] = std::pair< std::string, std::string >( "Mars", "Viking"  );
    oneWayLinkEnds[ 1 ][ receiver ] = std::pair< std::string, std::string >( "Earth", "Yarragadee" );

    oneWayLinkEnds[ 2 ][ transmitter ] = std::pair< std::string, std::string >( "Earth", "Yarragadee" );
    oneWayLinkEnds[ 2 ][ receiver ] = std::pair< std::string, std::string >( "Mars", "MSL" );

    boost::shared_ptr< NWayRangeObservationModel< > > nWayRangeModel =
            boost::make_shared< NWayRangeObservationModel< > >(
                nWayLinkEndMap, bodyMap );

    std::vector< boost::shared_ptr< OneWayRangeObservationModel< > > > oneWayRangeObservationModels;

    for( unsigned int i = 0; i < oneWayLinkEnds.size( ); i++ )
    {
        oneWayRangeObservationModels.push_back( boost::make_shared< OneWayRangeObservationModel< > >(
                                                    oneWayLinkEnds.at( i ).at( receiver ), oneWayLinkEnds.at( i ).at( transmitter ), bodyMap ) );
    }

    double observationTime = 1.1E7;

    boost::function< std::vector< double >( ) > retransmissionDelayFunctions;
    retransmissionDelayFunctions = boost::lambda::constant( boost::assign::list_of( 31.2 )( -12.4 ) );
    nWayRangeModel->resetRetransmissionDelaysFunction( retransmissionDelayFunctions );

    std::vector< double > oneWayRanges;
    std::vector< std::vector< basic_mathematics::Vector6d > > oneWayVectorOfStates;
    std::vector< std::vector< double > > oneWayVectorOfTimes;
    oneWayRanges.resize( 3 );
    oneWayVectorOfStates.resize( 3 );
    oneWayVectorOfTimes.resize( 3 );

    std::vector< basic_mathematics::Vector6d > nWayVectorOfStates;
    std::vector< double > nWayVectorOfTimes;


    for( LinkEnds::iterator linkEndIterator = nWayLinkEndMap.begin( ); linkEndIterator != nWayLinkEndMap.end( ); linkEndIterator++ )
    {
        double nWayRangeObservation = nWayRangeModel->computeObservationsAndLinkEndData(
                    observationTime, linkEndIterator->first, nWayVectorOfTimes, nWayVectorOfStates ).x( );


        for( unsigned int i = 0; i < 3; i++ )
        {
            oneWayRanges[ i ] = oneWayRangeObservationModels.at( i )->computeObservationsAndLinkEndData(
                        nWayVectorOfTimes.at( 2 * i ), transmitter, oneWayVectorOfTimes[ i ], oneWayVectorOfStates[ i ] ).x( );
            BOOST_CHECK_CLOSE_FRACTION( oneWayVectorOfTimes[ i ][ 0 ], nWayVectorOfTimes[ 2 * i ], std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( oneWayVectorOfTimes[ i ][ 1 ], nWayVectorOfTimes[ 2 * i + 1 ], std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( oneWayVectorOfStates[ i ][ 0 ], nWayVectorOfStates[ 2 * i ], std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( oneWayVectorOfStates[ i ][ 1 ], nWayVectorOfStates[ 2 * i + 1 ], std::numeric_limits< double >::epsilon( ) )
        }
        double reconstructredNWayRange =
                oneWayRanges[ 0 ] + oneWayRanges[ 1 ] + oneWayRanges[ 2 ] +
                ( retransmissionDelayFunctions( )[ 0 ] + retransmissionDelayFunctions( )[ 1 ] ) * physical_constants::SPEED_OF_LIGHT;


        BOOST_CHECK_CLOSE_FRACTION( nWayRangeObservation, reconstructredNWayRange, std::numeric_limits< double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_CASE( test_biased_n_way_observation_models )
{
    NamedBodyMap bodyMap = creatEnvironment( );

    LinkEnds nWayLinkEndMap;

    nWayLinkEndMap[ transmitter ] = std::pair< std::string, std::string >( "Earth", "Graz" );
    nWayLinkEndMap[ reflector1 ] = std::pair< std::string, std::string >( "Mars", "Viking" );
    nWayLinkEndMap[ reflector2 ] = std::pair< std::string, std::string >( "Earth", "Yarragadee" );
    nWayLinkEndMap[ reflector3 ] = std::pair< std::string, std::string >( "Mars", "MSL" );
    nWayLinkEndMap[ receiver ] = std::pair< std::string, std::string >( "Earth", "Graz" );


    std::vector< Time > clockErrorArcTimes;
    double singleClockArcLength = 1.0E6;
    Time currentTime = Time( 1.05E7 );
    clockErrorArcTimes.push_back( 1.05E7 );
    clockErrorArcTimes.push_back( 1.05E7 + singleClockArcLength );
    clockErrorArcTimes.push_back( 1.05E7 + 2.0 * singleClockArcLength );
    clockErrorArcTimes.push_back( 1.05E7 + 3.0 * singleClockArcLength );
    std::vector< double > clockArcPolynomialErrors;
    clockArcPolynomialErrors.push_back( 0.0 );
    clockArcPolynomialErrors.push_back( 0.0 );
    clockArcPolynomialErrors.push_back( 0.0 );


    boost::shared_ptr< GroundStation > grazStation = boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Earth"] )->
            getGroundStation( "Graz" );
    boost::shared_ptr< TimingSystem > grazTimingSystem = boost::make_shared< TimingSystem >(
                clockErrorArcTimes, clockArcPolynomialErrors );
    boost::shared_ptr< SystemHardware > grazTransmitterSystem = boost::make_shared< SystemHardware >( );
    grazTransmitterSystem->setTimingSystem( grazTimingSystem );
    grazStation->setSystemHardware( grazTransmitterSystem );

    boost::shared_ptr< GroundStation > yarragadeeStation = boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Earth"] )->
            getGroundStation( "Yarragadee" );
    boost::shared_ptr< TimingSystem > yarragadeeTimingSystem = boost::make_shared< TimingSystem >(
                clockErrorArcTimes, clockArcPolynomialErrors );
    boost::shared_ptr< SystemHardware > yarragadeeTransmitterSystem = boost::make_shared< SystemHardware >( );
    yarragadeeTransmitterSystem->setTimingSystem( yarragadeeTimingSystem );
    yarragadeeStation->setSystemHardware( yarragadeeTransmitterSystem );

    boost::shared_ptr< GroundStation > vikingStation = boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Mars"] )->
            getGroundStation( "Viking" );
    boost::shared_ptr< TimingSystem > vikingTimingSystem = boost::make_shared< TimingSystem >(
                clockErrorArcTimes, clockArcPolynomialErrors );
    boost::shared_ptr< SystemHardware > vikingTransmitterSystem = boost::make_shared< SystemHardware >( );
    vikingTransmitterSystem->setTimingSystem( vikingTimingSystem );
    vikingStation->setSystemHardware( vikingTransmitterSystem );

    boost::shared_ptr< GroundStation > mslStation = boost::dynamic_pointer_cast< CelestialBody >( bodyMap[ "Mars"] )->
            getGroundStation( "MSL" );
    boost::shared_ptr< TimingSystem > mslTimingSystem = boost::make_shared< TimingSystem >(
                clockErrorArcTimes, clockArcPolynomialErrors );
    boost::shared_ptr< SystemHardware > mslTransmitterSystem = boost::make_shared< SystemHardware >( );
    mslTransmitterSystem->setTimingSystem( mslTimingSystem );
    mslStation->setSystemHardware( mslTransmitterSystem );

    std::map< ObservableType, std::vector< LinkEnds > > biasList;
    biasList[ nWayRange ].push_back( nWayLinkEndMap );
    setObservationBiasesWithTimingErrors( biasList, bodyMap, 0 );

    boost::shared_ptr< NWayRangeObservationModel< > > nWayRangeModel =
            boost::dynamic_pointer_cast< NWayRangeObservationModel< > >(
                ObservationModelCreator< 1, double, double, double >::createObservationModel(
                    nWayRange, nWayLinkEndMap, bodyMap, std::vector< boost::shared_ptr< LightTimeCorrectionSettings > >( ),
                    AuxiliaryTypedefs< double >::AuxiliaryObservationSettingsInterfaceMap( ) ,
                    observationBiasInterfaceList.at( nWayRange ).at( nWayLinkEndMap ) ) );

    double observationTime = 1.1E7;

    boost::function< std::vector< double >( ) > retransmissionDelayFunctions;
    retransmissionDelayFunctions = boost::lambda::constant( boost::assign::list_of( 31.2 )( -12.4 )( 68.7 ) );
    nWayRangeModel->resetRetransmissionDelaysFunction( retransmissionDelayFunctions );

    std::vector< double > oneWayRanges;
    std::vector< std::vector< basic_mathematics::Vector6d > > oneWayVectorOfStates;
    std::vector< std::vector< double > > oneWayVectorOfTimes;
    oneWayRanges.resize( 3 );
    oneWayVectorOfStates.resize( 3 );
    oneWayVectorOfTimes.resize( 3 );

    std::vector< basic_mathematics::Vector6d > nWayVectorOfStates;
    std::vector< double > nWayVectorOfTimes;

    double nominalNWayRangeObservation = nWayRangeModel->computeObservationsAndLinkEndData(
                observationTime, transmitter, nWayVectorOfTimes, nWayVectorOfStates ).x( );

    std::vector< boost::shared_ptr< TimingSystem > > reflectorTimingSystems;
    reflectorTimingSystems.push_back( vikingTimingSystem );
    reflectorTimingSystems.push_back( yarragadeeTimingSystem );
    reflectorTimingSystems.push_back( mslTimingSystem );

    double currentObservation;

    std::map< int, double > clockArcPolynomialErrorsMap;
    clockArcPolynomialErrorsMap[ 0 ] = 0.0;
    clockArcPolynomialErrorsMap[ 1 ] = 0.0;
    clockArcPolynomialErrorsMap[ 2 ] = 0.0;


    for( unsigned int i = 0; i < reflectorTimingSystems.size( ); i++ )
    {

        clockArcPolynomialErrorsMap[ 0 ] = 1.0E-3;
        reflectorTimingSystems.at( i )->setGlobalPolynomialClockCorrections( clockArcPolynomialErrorsMap );
        currentObservation = nWayRangeModel->computeObservationsAndLinkEndData(
                    observationTime, transmitter, nWayVectorOfTimes, nWayVectorOfStates ).x( );
        BOOST_CHECK_SMALL( currentObservation - nominalNWayRangeObservation, std::numeric_limits< double >::epsilon( ) );
        clockArcPolynomialErrorsMap[ 0 ] = 0;
        reflectorTimingSystems.at( i )->setGlobalPolynomialClockCorrections( clockArcPolynomialErrorsMap );

        clockArcPolynomialErrorsMap[ 1 ] = 1.0E-5;
        reflectorTimingSystems.at( i )->setGlobalPolynomialClockCorrections( clockArcPolynomialErrorsMap );
        currentObservation = nWayRangeModel->computeObservationsAndLinkEndData(
                    observationTime, transmitter, nWayVectorOfTimes, nWayVectorOfStates ).x( );
        BOOST_CHECK_SMALL( currentObservation - nominalNWayRangeObservation - (
                   physical_constants::SPEED_OF_LIGHT * clockArcPolynomialErrorsMap.at( 1 ) *
                                    retransmissionDelayFunctions( ).at( i ) ), 1.0E-3 );
        clockArcPolynomialErrorsMap[ 1 ] = 0;
        reflectorTimingSystems.at( i )->setGlobalPolynomialClockCorrections( clockArcPolynomialErrorsMap );

        clockArcPolynomialErrorsMap[ 2 ] = 1.0E-7;
        reflectorTimingSystems.at( i )->setGlobalPolynomialClockCorrections( clockArcPolynomialErrorsMap );
        currentObservation = nWayRangeModel->computeObservationsAndLinkEndData(
                    observationTime, transmitter, nWayVectorOfTimes, nWayVectorOfStates ).x( );
        BOOST_CHECK_SMALL( currentObservation - nominalNWayRangeObservation - (
                   physical_constants::SPEED_OF_LIGHT * clockArcPolynomialErrorsMap.at( 2 ) * (
                       std::pow( nWayVectorOfTimes.at( 2 * i + 2 )  - 1.05E7, 2.0 ) -
                       std::pow( nWayVectorOfTimes.at( 2 * i + 1 )  - 1.05E7, 2.0 ) ) ), 1.0E-3 );

        clockArcPolynomialErrorsMap[ 2 ] = 0;
        reflectorTimingSystems.at( i )->setGlobalPolynomialClockCorrections( clockArcPolynomialErrorsMap );
    }

    clockArcPolynomialErrorsMap[ 0 ] = 1.0E-3;
    grazTimingSystem->setGlobalPolynomialClockCorrections( clockArcPolynomialErrorsMap );
    currentObservation = nWayRangeModel->computeObservationsAndLinkEndData(
                observationTime, transmitter, nWayVectorOfTimes, nWayVectorOfStates ).x( );
    BOOST_CHECK_SMALL( currentObservation - nominalNWayRangeObservation, std::numeric_limits< double >::epsilon( ) );
    clockArcPolynomialErrorsMap[ 0 ] = 0;
    grazTimingSystem->setGlobalPolynomialClockCorrections( clockArcPolynomialErrorsMap );

    clockArcPolynomialErrorsMap[ 1 ] = 1.0E-5;
    grazTimingSystem->setGlobalPolynomialClockCorrections( clockArcPolynomialErrorsMap );
    currentObservation = nWayRangeModel->computeObservationsAndLinkEndData(
                observationTime, transmitter, nWayVectorOfTimes, nWayVectorOfStates ).x( );
    BOOST_CHECK_SMALL( currentObservation - nominalNWayRangeObservation - (
               physical_constants::SPEED_OF_LIGHT * clockArcPolynomialErrorsMap.at( 1 ) * (
                   nWayVectorOfTimes.at( nWayVectorOfTimes.size( ) - 1 ) - nWayVectorOfTimes.at( 0 ) ) ), 1.0E-3 );
    clockArcPolynomialErrorsMap[ 1 ] = 0;
    grazTimingSystem->setGlobalPolynomialClockCorrections( clockArcPolynomialErrorsMap );

    clockArcPolynomialErrorsMap[ 2 ] = 1.0E-7;
    grazTimingSystem->setGlobalPolynomialClockCorrections( clockArcPolynomialErrorsMap );
    currentObservation = nWayRangeModel->computeObservationsAndLinkEndData(
                observationTime, transmitter, nWayVectorOfTimes, nWayVectorOfStates ).x( );
    BOOST_CHECK_SMALL( currentObservation - nominalNWayRangeObservation - (
               physical_constants::SPEED_OF_LIGHT * clockArcPolynomialErrorsMap.at( 2 ) * (
                   std::pow( nWayVectorOfTimes.at( nWayVectorOfTimes.size( ) - 1 )  - 1.05E7, 2.0 ) -
                   std::pow( nWayVectorOfTimes.at( 0 )  - 1.05E7, 2.0 ) ) ), 1.0E-3 );
    clockArcPolynomialErrorsMap[ 2 ] = 0;
    grazTimingSystem->setGlobalPolynomialClockCorrections( clockArcPolynomialErrorsMap );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
