/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

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

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"


namespace tudat
{

namespace unit_tests
{

Eigen::Vector6d getGroundStationPosition(
        const double time )
{
    Eigen::Vector3d nominalPosition;
    nominalPosition << 4194511.7, 1162789.7, 4647362.5;
    Eigen::Vector3d nominalVelocity;
    nominalVelocity << 3.0E-3, -2.3E-3, 2.4E3;
    nominalVelocity /= physical_constants::JULIAN_YEAR;
    return ( Eigen::Vector6d( ) << nominalPosition + time * nominalVelocity, nominalVelocity ).finished( );
}

BOOST_AUTO_TEST_SUITE( test_composite_ephemeris )

//! Test functionality of composite ephemeris. NOTE: Part of the functionality is also tested by test_FrameManager
BOOST_AUTO_TEST_CASE( testCompositeEphemeris )
{
    using namespace tudat::interpolators;
    using namespace tudat::simulation_setup;
    using namespace tudat::ephemerides;

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );


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

    // Retrieve Earth state/rotation objects
    boost::shared_ptr< Ephemeris > earthEphemeris = bodyMap.at( "Earth" )->getEphemeris( );
    boost::shared_ptr< RotationalEphemeris > rotationModel = bodyMap.at( "Earth" )->getRotationalEphemeris( );

    // Create reference point CompositeEphemeris objects (double and long double state scalars).
    boost::shared_ptr< Ephemeris > ephemeris1 = createReferencePointEphemeris< double, double >(
                earthEphemeris, rotationModel, &getGroundStationPosition );
    boost::shared_ptr< Ephemeris > ephemeris2 = createReferencePointEphemeris< double, long double >(
                earthEphemeris, rotationModel, &getGroundStationPosition );
    double testTime = 1.05E7;

    // Manually compute double state
    Eigen::Vector6d doubleStateFromDoubleTime;
    doubleStateFromDoubleTime.segment( 0, 3 ) = earthEphemeris->getCartesianState( testTime ).segment( 0, 3 ) +
            rotationModel->getRotationToBaseFrame( testTime ) * getGroundStationPosition( testTime ).segment( 0, 3 );
    doubleStateFromDoubleTime.segment( 3, 3 ) = earthEphemeris->getCartesianState( testTime ).segment( 3, 3 ) +
            rotationModel->getRotationToBaseFrame( testTime ) * getGroundStationPosition( testTime ).segment( 3, 3 ) +
            rotationModel->getDerivativeOfRotationToBaseFrame( testTime ) *
            getGroundStationPosition( testTime ).segment( 0, 3 );

    // Manually compute long double state
    Eigen::Matrix< long double, 6, 1 > longDoubleStateFromDoubleTime;
    longDoubleStateFromDoubleTime.segment( 0, 3 ) =
            earthEphemeris->getCartesianLongState( testTime ).segment( 0, 3 ) +
            ( rotationModel->getRotationToBaseFrame( testTime ) *
              getGroundStationPosition( testTime ).segment( 0, 3 ) ).cast< long double >( );
    longDoubleStateFromDoubleTime.segment( 3, 3 ) =
            earthEphemeris->getCartesianLongState( testTime ).segment( 3, 3 ) +
            ( rotationModel->getRotationToBaseFrame( testTime ) *
              getGroundStationPosition( testTime ).segment( 3, 3 ) ).cast< long double >( ) +
            ( rotationModel->getDerivativeOfRotationToBaseFrame( testTime ) *
              getGroundStationPosition( testTime ).segment( 0, 3 ) ).cast< long double >( );
    double doubleTolerance = std::numeric_limits< double >::epsilon( );
    double longDoubleTolerance = std::numeric_limits< long double >::epsilon( );


    // Test consistency of manual functions
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                doubleStateFromDoubleTime, longDoubleStateFromDoubleTime.cast< double >( ), doubleTolerance );

    // Test double composte ephemeris
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ephemeris1->getCartesianState( testTime ),
                doubleStateFromDoubleTime, doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ephemeris1->getCartesianLongState( testTime ),
                doubleStateFromDoubleTime.cast< long double >( ), doubleTolerance );

    // Test long double composte ephemeris, tolerances are not fully met as longDoubleTolerance because rotation is only
    // defined using double state scalars. This is especially influential for the z-components, where the nominal value
    // is much smaller.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ephemeris2->getCartesianState( testTime ),
                doubleStateFromDoubleTime, doubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ephemeris2->getCartesianLongState( testTime ).segment( 0, 2 ),
                longDoubleStateFromDoubleTime.segment( 0, 2 ), longDoubleTolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ephemeris2->getCartesianLongState( testTime ).segment( 2, 1 ),
                longDoubleStateFromDoubleTime.segment( 2, 1 ), ( 1000.0 * longDoubleTolerance ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ephemeris2->getCartesianLongState( testTime ).segment( 3, 2 ),
                longDoubleStateFromDoubleTime.segment( 3, 2 ), ( 10.0 * longDoubleTolerance ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ephemeris2->getCartesianLongState( testTime ).segment( 5, 1 ),
                longDoubleStateFromDoubleTime.segment( 5, 1 ), ( 10000.0 * longDoubleTolerance ) );
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
