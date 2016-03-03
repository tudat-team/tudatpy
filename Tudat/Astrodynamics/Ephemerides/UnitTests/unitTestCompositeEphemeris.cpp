#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Basics/testMacros.h"

#include "Tudat/InputOutput/basicInputOutput.h"

#include "Astrodynamics/Bodies/body.h"
#include "Astrodynamics/Ephemerides/createLinkEndEphemeris.h"
#include "SimulationSetup/defaultBodies.h"
#include "SimulationSetup/createBodies.h"
#include "SimulationSetup/createGroundStations.h"



namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_composite_ephemeris )

BOOST_AUTO_TEST_CASE( testCompositeEphemeris )
{

    using namespace tudat::observation_models;
    using namespace tudat::orbit_determination;
    using namespace tudat::estimatable_parameters;
    using namespace tudat::interpolators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::spice_interface;
    using namespace tudat::bodies;
    using namespace tudat::simulation_setup;
    using namespace tudat::ephemerides;

    //Load spice kernels.
    std::string kernelsPath = input_output::getDataFilesRootPath( ) + "SpiceKernels/";
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "jup291.bsp");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    //Define setting for total number of bodies and those which need to be integrated numerically.
    //The first numberOfNumericalBodies from the bodyNames vector will be integrated numerically.

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = 1.1E7;
    double maximumTimeStep = 3600.0;

    double buffer = 5.0 * maximumTimeStep;

    // Create bodies needed in simulation
    std::map< std::string, boost::shared_ptr< Body > > bodyMap = createCelestialBodies(
                getDefaultBodySettings< double, double >( bodyNames,initialEphemerisTime - buffer, finalEphemerisTime + buffer, simple ) );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Create ground stations
    std::pair< std::string, std::string > grazStation = std::pair< std::string, std::string >( "Earth", "Graz" );
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.push_back( grazStation );
    createGroundStations( bodyMap, groundStations );

    setSingleBodyGroundStationPositionVariationFunctions( boost::dynamic_pointer_cast< CelestialBody >( bodyMap["Earth"] ), basic_solid_body,
                                                          initialEphemerisTime, finalEphemerisTime );

    boost::shared_ptr< Body > earth = bodyMap.at( "Earth" );
    boost::shared_ptr< RotationalEphemeris > rotationModel = earth->getRotationalEphemeris( );
    boost::shared_ptr< GroundStation > groundStation = boost::dynamic_pointer_cast< CelestialBody >( earth )->getGroundStation( "Graz" );

    boost::shared_ptr< Ephemeris > ephemeris1 = createReferencePointEphemeris< double, double >(
                earth, rotationModel, boost::bind( &GroundStation::getStateInPlanetFixedFrame< double >, groundStation, _1 ) );
    boost::shared_ptr< Ephemeris > ephemeris2 = createReferencePointEphemeris< double, long double >(
                earth, rotationModel, boost::bind( &GroundStation::getStateInPlanetFixedFrame< double >, groundStation, _1 ) );
    boost::shared_ptr< Ephemeris > ephemeris3 = createReferencePointEphemeris< Time, double >(
                earth, rotationModel, boost::bind( &GroundStation::getStateInPlanetFixedFrame< Time >, groundStation, _1 ) );
    boost::shared_ptr< Ephemeris > ephemeris4 = createReferencePointEphemeris< Time, long double >(
                earth, rotationModel, boost::bind( &GroundStation::getStateInPlanetFixedFrame< Time >, groundStation, _1 ) );

    double testTime = 1.05E7;
    Time testTime2 = Time( testTime );

    basic_mathematics::Vector6d doubleStateFromDoubleTime;
    doubleStateFromDoubleTime.segment( 0, 3 ) = earth->getStateInBaseFrameFromEphemeris( testTime ).segment( 0, 3 ) +
            rotationModel->getRotationToBaseFrame( testTime ) * groundStation->getStateInPlanetFixedFrame( testTime ).segment( 0, 3 );
    doubleStateFromDoubleTime.segment( 3, 3 ) = earth->getStateInBaseFrameFromEphemeris( testTime ).segment( 3, 3 ) +
            rotationModel->getRotationToBaseFrame( testTime ) * groundStation->getStateInPlanetFixedFrame( testTime ).segment( 3, 3 );

    basic_mathematics::Vector6d doubleStateFromTimeTime;
    doubleStateFromTimeTime.segment( 0, 3 ) = earth->getStateInBaseFrameFromEphemerisTime( testTime2 ).segment( 0, 3 ) +
            rotationModel->getRotationToBaseFrame( testTime ) * groundStation->getStateInPlanetFixedFrame( testTime ).segment( 0, 3 );
    doubleStateFromTimeTime.segment( 3, 3 ) = earth->getStateInBaseFrameFromEphemerisTime( testTime2 ).segment( 3, 3 ) +
            rotationModel->getRotationToBaseFrame( testTime ) * groundStation->getStateInPlanetFixedFrame( testTime ).segment( 3, 3 );

    Eigen::Matrix< long double, 6, 1 > longDoubleStateFromDoubleTime;
    longDoubleStateFromDoubleTime.segment( 0, 3 ) = earth->getLongStateInBaseFrameFromEphemeris( testTime ).segment( 0, 3 ) +
            ( rotationModel->getRotationToBaseFrame( testTime ) * groundStation->getStateInPlanetFixedFrame( testTime ).segment( 0, 3 ) ).cast< long double >( );
    longDoubleStateFromDoubleTime.segment( 3, 3 ) = earth->getLongStateInBaseFrameFromEphemeris( testTime ).segment( 3, 3 ) +
            ( rotationModel->getRotationToBaseFrame( testTime ) * groundStation->getStateInPlanetFixedFrame( testTime ).segment( 3, 3 ) ).cast< long double >( );

    Eigen::Matrix< long double, 6, 1 > longDoubleStateFromTimeTime;
    longDoubleStateFromTimeTime.segment( 0, 3 ) = earth->getLongStateInBaseFrameFromEphemerisTime( testTime2 ).segment( 0, 3 ) +
            ( rotationModel->getRotationToBaseFrame( testTime ) * groundStation->getStateInPlanetFixedFrame( testTime ).segment( 0, 3 ) ).cast< long double >( );
    longDoubleStateFromTimeTime.segment( 3, 3 ) = earth->getLongStateInBaseFrameFromEphemerisTime( testTime2 ).segment( 3, 3 ) +
            ( rotationModel->getRotationToBaseFrame( testTime ) * groundStation->getStateInPlanetFixedFrame( testTime ).segment( 3, 3 ) ).cast< long double >( );

    double doubleTolerance = std::numeric_limits< double >::epsilon( );
    double longDoubleTolerance = std::numeric_limits< long double >::epsilon( );

    //NOTE: Make test robust for smaller entries 2 and 5

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( doubleStateFromDoubleTime, doubleStateFromTimeTime, doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( doubleStateFromDoubleTime, longDoubleStateFromDoubleTime.cast< double >( ), doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( doubleStateFromDoubleTime, longDoubleStateFromTimeTime.cast< double >( ), doubleTolerance );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( doubleStateFromTimeTime, doubleStateFromDoubleTime, doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( doubleStateFromTimeTime, longDoubleStateFromDoubleTime.cast< double >( ), doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( doubleStateFromTimeTime, longDoubleStateFromTimeTime.cast< double >( ), doubleTolerance );


    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( longDoubleStateFromDoubleTime.cast< double >( ), doubleStateFromTimeTime, doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( longDoubleStateFromDoubleTime.cast< double >( ), doubleStateFromDoubleTime, doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( longDoubleStateFromDoubleTime.cast< double >( ), longDoubleStateFromTimeTime.cast< double >( ), doubleTolerance );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris1->getCartesianStateFromEphemeris( testTime ), doubleStateFromDoubleTime, doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris1->getCartesianStateFromEphemeris( testTime2 ), doubleStateFromDoubleTime, doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris1->getCartesianLongStateFromEphemeris( testTime ), doubleStateFromDoubleTime.cast< long double >( ), doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris1->getCartesianLongStateFromEphemeris( testTime2 ), doubleStateFromDoubleTime.cast< long double >( ), doubleTolerance );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris2->getCartesianStateFromEphemeris( testTime ), doubleStateFromDoubleTime, doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris2->getCartesianStateFromEphemeris( testTime2 ), doubleStateFromTimeTime, doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris2->getCartesianLongStateFromEphemeris( testTime ).segment( 0, 2 ),
                                       longDoubleStateFromDoubleTime.segment( 0, 2 ), longDoubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris2->getCartesianLongStateFromEphemeris( testTime2 ).segment( 0, 2 ),
                                       longDoubleStateFromDoubleTime.segment( 0, 2 ), longDoubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris2->getCartesianLongStateFromEphemeris( testTime ).segment( 3, 3 ),
                                       longDoubleStateFromDoubleTime.segment( 3, 3 ), longDoubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris2->getCartesianLongStateFromEphemeris( testTime2 ).segment( 3, 3 ),
                                       longDoubleStateFromDoubleTime.segment( 3, 3 ), longDoubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris2->getCartesianLongStateFromEphemeris( testTime ).segment( 2, 1 ),
                                       longDoubleStateFromDoubleTime.segment( 2, 1 ), ( 1000.0 * longDoubleTolerance ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris2->getCartesianLongStateFromEphemeris( testTime2 ).segment( 2, 1 ),
                                       longDoubleStateFromDoubleTime.segment( 2, 1 ), ( 1000.0 * longDoubleTolerance ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris3->getCartesianStateFromEphemeris( testTime ), doubleStateFromDoubleTime, doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris3->getCartesianStateFromEphemeris( testTime2 ), doubleStateFromDoubleTime, doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris3->getCartesianLongStateFromEphemeris( testTime ), doubleStateFromDoubleTime.cast< long double >( ), doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris3->getCartesianLongStateFromEphemeris( testTime2 ), doubleStateFromDoubleTime.cast< long double >( ), doubleTolerance );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris4->getCartesianStateFromEphemeris( testTime ), doubleStateFromDoubleTime, doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris4->getCartesianStateFromEphemeris( testTime2 ), doubleStateFromDoubleTime, doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris4->getCartesianLongStateFromEphemeris( testTime ).segment( 0, 2 ),
                                       longDoubleStateFromDoubleTime.segment( 0, 2 ), longDoubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris4->getCartesianLongStateFromEphemeris( testTime2 ).segment( 0, 2 ),
                                       longDoubleStateFromTimeTime.segment( 0, 2 ), longDoubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris4->getCartesianLongStateFromEphemeris( testTime ).segment( 3, 3 ),
                                       longDoubleStateFromDoubleTime.segment( 3, 3 ), longDoubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris4->getCartesianLongStateFromEphemeris( testTime2 ).segment( 3, 3 ),
                                       longDoubleStateFromTimeTime.segment( 3, 3 ), longDoubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris4->getCartesianLongStateFromEphemeris( testTime ).segment( 2, 1 ),
                                       longDoubleStateFromDoubleTime.segment( 2, 1 ), ( 1000.0 * longDoubleTolerance ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris4->getCartesianLongStateFromEphemeris( testTime2 ).segment( 2, 1 ),
                                       longDoubleStateFromTimeTime.segment( 2, 1 ), ( 1000.0 * longDoubleTolerance ) );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
