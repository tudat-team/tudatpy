
#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/InputOutput/basicInputOutput.h"

#include "Astrodynamics/Bodies/body.h"
#include "Astrodynamics/Bodies/celestialBody.h"
#include "Astrodynamics/ObservationModels/angularPositionObservationModel.h"
#include "SimulationSetup/defaultBodies.h"
#include "SimulationSetup/createBodies.h"
#include "SimulationSetup/createGroundStations.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;
using namespace tudat::bodies;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;


BOOST_AUTO_TEST_SUITE( test_angular_position_model )


BOOST_AUTO_TEST_CASE( testAngularPositionModel )
{
    std::string kernelsPath = input_output::getDataFilesRootPath( ) + "SpiceKernels/";
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    std::map< std::string, BodySettingLevel > bodySettings;
    bodySettings[ "Earth" ] = simple;
    bodySettings[ "Sun" ] = simple;
    bodySettings[ "Moon" ] = simple;
    bodySettings[ "Mars" ] = simple;

    // Specify initial time
    Time initialEphemerisTime = basic_astrodynamics::calculateJulianDaySinceEpoch(
                boost::gregorian::date( 2013, 3, 1 ), 0.0 ) * physical_constants::JULIAN_DAY;

    Time finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;

    double buffer = 10.0 * maximumTimeStep;

    // Create bodies needed in simulation
    std::map< std::string, boost::shared_ptr< BodySettings > > defaultBodySettings =
            getDefaultBodySettings< double, Time >(
                bodySettings, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    defaultBodySettings[ "Earth" ]->gravityFieldVariationSettings.clear( );

    std::map< std::string, boost::shared_ptr< Body > > bodyMap = createCelestialBodies( defaultBodySettings );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Create ground stations
    std::pair< std::string, std::string > mslStation = std::pair< std::string, std::string >( "Mars", "MSL" );

    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.push_back( mslStation );
    createGroundStations( bodyMap, groundStations );

    AngularPositionObservationModel< double, Time, double > observationModel =
            AngularPositionObservationModel< double,Time, double >(
                std::make_pair( "Earth", "" ), std::make_pair( "Mars", "MSL" ), bodyMap );

    Time receiverObservationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;

    std::vector< double > linkEndTimes;
    std::vector< basic_mathematics::Vector6d > linkEndStates;


    Eigen::Vector2d observationFromReceptionTime = observationModel.computeObservations( receiverObservationTime, receiver );
    Eigen::Vector2d observationFromReceptionTime2 = observationModel.computeObservationsAndLinkEndData(
                receiverObservationTime, receiver, linkEndTimes, linkEndStates );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( observationFromReceptionTime, observationFromReceptionTime2, std::numeric_limits< double >::epsilon( ) );

    BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( receiverObservationTime ), linkEndTimes[ 1 ], std::numeric_limits< double >::epsilon( ) );

    Eigen::Vector3d positionDifference = ( linkEndStates[ 0 ] - linkEndStates[ 1 ] ).segment( 0, 3 );

    Eigen::Vector3d sphericalRelativeCoordinates = coordinate_conversions::convertCartesianToSpherical(
                positionDifference );

    BOOST_CHECK_CLOSE_FRACTION( sphericalRelativeCoordinates.z( ), observationFromReceptionTime( 0 ), std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( mathematical_constants::PI / 2.0 - sphericalRelativeCoordinates.y( ), observationFromReceptionTime( 1 ),
                                std::numeric_limits< double >::epsilon( ) );

    Time transmitterObservationTime = receiverObservationTime - ( linkEndTimes[ 1 ] - linkEndTimes[ 0 ] );

    BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( transmitterObservationTime ), linkEndTimes[ 0 ], std::numeric_limits< double >::epsilon( ) );

    std::vector< double > linkEndTimes2;
    std::vector< basic_mathematics::Vector6d > linkEndStates2;

    Eigen::Vector2d observationFromTransmissionTime = observationModel.computeObservations( transmitterObservationTime, transmitter );
    Eigen::Vector2d observationFromTransmissionTime2 = observationModel.computeObservationsAndLinkEndData(
                transmitterObservationTime, transmitter, linkEndTimes2, linkEndStates2 );


    for( unsigned int i = 0; i < 2; i++ )
    {
        BOOST_CHECK_SMALL( observationFromTransmissionTime( i ) - observationFromTransmissionTime2( i ), 5.0E-15 );
        BOOST_CHECK_SMALL( observationFromTransmissionTime( i ) - observationFromReceptionTime( i ), 5.0E-15 );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}

