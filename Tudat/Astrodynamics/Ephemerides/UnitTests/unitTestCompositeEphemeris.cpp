#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/Ephemerides/compositeEphemeris.h"
#include "Tudat/Basics/testMacros.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/SimulationSetup/body.h"
#include "Tudat/SimulationSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/createBodies.h"


namespace tudat
{

namespace unit_tests
{

basic_mathematics::Vector6d getGroundStationPosition(
        const double time )
{
    Eigen::Vector3d nominalPosition;
    nominalPosition << 4194511.7, 1162789.7, 4647362.5;
    Eigen::Vector3d nominalVelocity;
    nominalVelocity << 3.0E-3, -2.3E-3, 2.4E3;
    nominalVelocity /= physical_constants::JULIAN_YEAR;
    return ( basic_mathematics::Vector6d( )<< nominalPosition + time * nominalVelocity, nominalVelocity ).finished( );
}

BOOST_AUTO_TEST_SUITE( test_composite_ephemeris )

//! Test functionality of composite ephemeris. NOTE: Part of the functionality is also tested by test_FrameManager
BOOST_AUTO_TEST_CASE( testCompositeEphemeris )
{
    using namespace tudat::interpolators;
    using namespace tudat::simulation_setup;
    using namespace tudat::ephemerides;

    //Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
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
    NamedBodyMap bodyMap = createBodies(
                getDefaultBodySettings( bodyNames,initialEphemerisTime - buffer, finalEphemerisTime + buffer ) );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    boost::shared_ptr< Ephemeris > earthEphemeris = bodyMap.at( "Earth" )->getEphemeris( );
    boost::shared_ptr< RotationalEphemeris > rotationModel = bodyMap.at( "Earth" )->getRotationalEphemeris( );

    boost::shared_ptr< Ephemeris > ephemeris1 = createReferencePointEphemeris< double, double >(
                earthEphemeris, rotationModel, &getGroundStationPosition );
    boost::shared_ptr< Ephemeris > ephemeris2 = createReferencePointEphemeris< double, long double >(
                earthEphemeris, rotationModel, &getGroundStationPosition );
    double testTime = 1.05E7;

    basic_mathematics::Vector6d doubleStateFromDoubleTime;
    doubleStateFromDoubleTime.segment( 0, 3 ) = earthEphemeris->getCartesianStateFromEphemeris( testTime ).segment( 0, 3 ) +
            rotationModel->getRotationToBaseFrame( testTime ) * getGroundStationPosition( testTime ).segment( 0, 3 );
    doubleStateFromDoubleTime.segment( 3, 3 ) = earthEphemeris->getCartesianStateFromEphemeris( testTime ).segment( 3, 3 ) +
            rotationModel->getRotationToBaseFrame( testTime ) * getGroundStationPosition( testTime ).segment( 3, 3 ) +
            rotationModel->getDerivativeOfRotationToBaseFrame( testTime ) * getGroundStationPosition( testTime ).segment( 0, 3 );

    Eigen::Matrix< long double, 6, 1 > longDoubleStateFromDoubleTime;
    longDoubleStateFromDoubleTime.segment( 0, 3 ) = earthEphemeris->getCartesianLongStateFromEphemeris( testTime ).segment( 0, 3 ) +
            ( rotationModel->getRotationToBaseFrame( testTime ) * getGroundStationPosition( testTime ).segment( 0, 3 ) ).cast< long double >( );
    longDoubleStateFromDoubleTime.segment( 3, 3 ) = earthEphemeris->getCartesianLongStateFromEphemeris( testTime ).segment( 3, 3 ) +
            ( rotationModel->getRotationToBaseFrame( testTime ) * getGroundStationPosition( testTime ).segment( 3, 3 ) ).cast< long double >( ) +
            ( rotationModel->getDerivativeOfRotationToBaseFrame( testTime ) * getGroundStationPosition( testTime ).segment( 0, 3 ) ).cast< long double >( );
    double doubleTolerance = std::numeric_limits< double >::epsilon( );
    double longDoubleTolerance = std::numeric_limits< long double >::epsilon( );

    //NOTE: Make test robust for smaller entries 2 and 5

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( doubleStateFromDoubleTime, longDoubleStateFromDoubleTime.cast< double >( ), doubleTolerance );


    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( longDoubleStateFromDoubleTime.cast< double >( ), doubleStateFromDoubleTime, doubleTolerance );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris1->getCartesianStateFromEphemeris( testTime ), doubleStateFromDoubleTime, doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris1->getCartesianLongStateFromEphemeris( testTime ), doubleStateFromDoubleTime.cast< long double >( ), doubleTolerance );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris2->getCartesianStateFromEphemeris( testTime ), doubleStateFromDoubleTime, doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris2->getCartesianLongStateFromEphemeris( testTime ).segment( 0, 2 ),
                                       longDoubleStateFromDoubleTime.segment( 0, 2 ), longDoubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris2->getCartesianLongStateFromEphemeris( testTime ).segment( 2, 1 ),
                                       longDoubleStateFromDoubleTime.segment( 2, 1 ), ( 1000.0 * longDoubleTolerance ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris2->getCartesianLongStateFromEphemeris( testTime ).segment( 3, 2 ),
                                       longDoubleStateFromDoubleTime.segment( 3, 2 ), ( 10.0 * longDoubleTolerance ) );//Tolerance is not fully met because rotation is in double
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ephemeris2->getCartesianLongStateFromEphemeris( testTime ).segment( 5, 1 ),
                                       longDoubleStateFromDoubleTime.segment( 5, 1 ), ( 10000.0 * longDoubleTolerance ) );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
