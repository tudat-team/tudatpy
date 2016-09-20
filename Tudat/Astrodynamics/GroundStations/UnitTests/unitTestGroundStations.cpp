#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/Astrodynamics/GroundStations/nominalGroundStationState.h"
#include "Tudat/Astrodynamics/GroundStations/groundStation.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/SimulationSetup/createGroundStations.h"
#include "Tudat/SimulationSetup/createBodies.h"
#include "Tudat/SimulationSetup/defaultBodies.h"


namespace tudat
{
namespace unit_tests
{


using namespace tudat::simulation_setup;
using namespace tudat::input_output;
using namespace tudat::basic_astrodynamics;
using namespace tudat::physical_constants;
using namespace tudat::ground_stations;

BOOST_AUTO_TEST_SUITE( test_ground_stations )

BOOST_AUTO_TEST_CASE( test_GroundStations )
{
    //Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( ) + "SpiceKernels/";

    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    //Define setting for total number of bodies and those which need to be integrated numerically.
    //The first numberOfNumericalBodies from the bodyNames vector will be integrated numerically.

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = 1.002E7;
    double buffer = 10000.0;


    // Create bodies needed in simulation
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames,initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings[ "Moon" ]->gravityFieldVariationSettings.clear( );
    NamedBodyMap bodyMap =
            createBodies( bodySettings );

    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.push_back( std::make_pair( "Earth", "Concepcion" ) );
    createGroundStations( bodyMap, groundStations );

    boost::shared_ptr< GroundStation > concepcionStation =
            bodyMap[ "Earth" ]->getGroundStation( "Concepcion" );


    Eigen::Vector3d eccentricity1;
    eccentricity1<< 0.3075, -1.0109, -1.0329;
    Eigen::Vector3d eccentricity2;
    eccentricity2<< 0.3077,-1.0107,-1.0360 ;

    boost::gregorian::date eccentricitySwitchTime( 2004, 10, 04 );

    Eigen::Vector3d nominalStationPosition = concepcionStation->getNominalStationState( )->getNominalCartesianPosition( );

    Eigen::Vector3d expectedStationState;
    expectedStationState<<0.149203303135301E+07, -.488794607556858E+07, -.380356589716459E+07;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedStationState, nominalStationPosition, 1.0E-15 );

    Eigen::Vector3d linearStationVelocity = concepcionStation->getNominalStationState( )->getLinearVelocityInMetersPerYear( );
    Eigen::Vector3d expectedStationVelocity;
    expectedStationVelocity<<0.348885761109463E-01, -.151253270336384E-02, 0.169053027347585E-01 ;

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( linearStationVelocity, expectedStationVelocity, 1.0E-15 );

    boost::gregorian::date testDate( 2005, 1, 1 );

    double testTime = calculateJulianDaySinceEpoch( testDate, 0.0, JULIAN_DAY_ON_J2000 ) *  JULIAN_DAY;

    Eigen::Vector3d stationPositionInTime = concepcionStation->getNominalStationState( )->getCartesianPositionInTime( testTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( stationPositionInTime, ( nominalStationPosition + eccentricity2 ), 1.0E-14 );

    testTime += 5.0 * physical_constants::JULIAN_YEAR;
    stationPositionInTime = concepcionStation->getNominalStationState( )->getCartesianPositionInTime( testTime );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( stationPositionInTime,
                                       ( nominalStationPosition + 5.0 * expectedStationVelocity + eccentricity2 ), 1.0E-15 );

    testTime = calculateJulianDaySinceEpoch( eccentricitySwitchTime, 0.0, JULIAN_DAY_ON_J2000 ) *  JULIAN_DAY + 0.001;
    Eigen::Vector3d stationPositionInTimePost = concepcionStation->getNominalStationState( )->getCartesianPositionInTime( testTime );
    testTime = calculateJulianDaySinceEpoch( eccentricitySwitchTime, 0.0, JULIAN_DAY_ON_J2000 ) *  JULIAN_DAY - 0.001;
    Eigen::Vector3d stationPositionInTimePre = concepcionStation->getNominalStationState( )->getCartesianPositionInTime( testTime );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( stationPositionInTimePost - stationPositionInTimePre ),
                                       ( eccentricity2 - eccentricity1 ), 2.0E-6 );

    //setGroundStationPositionVariationFunction( concepcionStation);

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
