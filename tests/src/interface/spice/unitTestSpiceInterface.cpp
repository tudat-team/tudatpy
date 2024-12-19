/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      To run this unit tests, a number of spice kernels need to be placed in the
 *      Spice kernel folder, by default external/SpiceInterface/kernels or the
 *      SPICE_KERNEL_CUSTOM_FOLDER folder set as an argument to cmake or in UserSetings.txt.
 *      The required kernels are:
 *           de421.bsp
 *           pck00009.tpc
 *           naif0012.tls
 *           de-403-masses.tpc
 *      They can be found in a single zip file on the wiki at
 *      http://tudat.tudelft.nl/projects/tudat/wiki/SpiceInterface/ on the Tudat website or,
 *      alternatively, on the NAIF server at ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>




#include <functional>
#include <memory>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/basics/testMacros.h"
#include "tudat/basics/utilityMacros.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"

#include <limits>
#include <stdexcept>

namespace tudat
{
namespace unit_tests
{

using Eigen::Vector6d;

BOOST_AUTO_TEST_SUITE( test_spice_wrappers )

//BOOST_AUTO_TEST_CASE( testMultiSpiceKernels )
//{
//    using namespace spice_interface;
//    using namespace physical_constants;

//    spice_interface::loadStandardSpiceKernels( );
//    spice_interface::loadSpiceKernelInTudat( paths::getSpiceKernelPath() + "/MEX_ROB_040101_041231_001.BSP" );
//    double testTime = 4.0 * JULIAN_YEAR + 40.0 * JULIAN_DAY;
//    Eigen::Vector6d marsCentricState = getBodyCartesianStateAtEpoch( "Mars Express", "Mars", "IAU_Mars", "NONE", testTime );
//}

// Test 1: Test Julian day <-> Ephemeris time conversions at J2000.
BOOST_AUTO_TEST_CASE( testNoKernelCrash )
{
    using namespace spice_interface;
    using namespace input_output;
    using namespace physical_constants;

    // Create settings at which states are to be evaluated.
    const std::string aberrationCorrections = "NONE";
    const std::string observer = "Solar System Barycenter";
    const std::string target = "Mars";
    const std::string referenceFrame = "J2000";
    const double ephemerisTime = 1.0E6;

    // Get state from wrapper for state:
    toggleErrorReturn( );
    suppressErrorOutput( );
    const Eigen::Vector6d wrapperState = getBodyCartesianStateAtEpoch(
        target, observer, referenceFrame, aberrationCorrections, ephemerisTime );
    TUDAT_UNUSED_PARAMETER( wrapperState );
}

// Test 1: Test Julian day <-> Ephemeris time conversions at J2000.
BOOST_AUTO_TEST_CASE( testSpiceWrappers_1 )
{
    std::cout<<"Test*******************"<<std::endl;
    using namespace spice_interface;
    using namespace input_output;
    using namespace physical_constants;

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Exact ephemeris time at J2000.
    const double ephemerisTimeOneYearAfterJ2000 = JULIAN_YEAR;

    // Convert to Julian day.
    const double julianDayOneYearAfterJ2000Spice = convertEphemerisTimeToJulianDate(
                ephemerisTimeOneYearAfterJ2000 );

    // Exact Julian day at J2000.
    const double julianDayOneYearAfterJ2000 = 2451545.0 + JULIAN_YEAR_IN_DAYS;

    // Compare exact and converted values of Julian dat at J2000.
    BOOST_CHECK_CLOSE_FRACTION( julianDayOneYearAfterJ2000Spice, julianDayOneYearAfterJ2000,
                                std::numeric_limits< double >::epsilon( ) );

    // Convert exact Julian day at J2000 to ephemeris time.
    const double ephemerisTimeOneYearAfterJ2000Spice = convertJulianDateToEphemerisTime(
                julianDayOneYearAfterJ2000 );

    // Compare exact and converted values of ephemeris time at J2000.
    BOOST_CHECK_CLOSE_FRACTION( ephemerisTimeOneYearAfterJ2000Spice,
                                ephemerisTimeOneYearAfterJ2000,
                                std::numeric_limits< double >::epsilon( ) );
}

// Test 2: Test retrieval position and state of bodies.
BOOST_AUTO_TEST_CASE( testSpiceWrappers_2 )
{
    using namespace spice_interface;
    using namespace input_output;
    using namespace physical_constants;

    // Create settings at which states are to be evaluated.
    const std::string aberrationCorrections = "NONE";
    const std::string observer = "Solar System Barycenter";
    const std::string target = "Mars";
    const std::string referenceFrame = "J2000";
    const double ephemerisTime = 1.0E6;

    // Get state from wrapper for state:
    const Eigen::Vector6d wrapperState = getBodyCartesianStateAtEpoch(
                target, observer, referenceFrame, aberrationCorrections, ephemerisTime );

    // Get position from wrapper for position:
    const Eigen::Vector3d wrapperPosition = getBodyCartesianPositionAtEpoch(
                target, observer, referenceFrame, aberrationCorrections, ephemerisTime );

    // Get state directly from spice.
    double directSpiceState[ 6 ];
    double directLightTime;
    spkezr_c( target.c_str( ), ephemerisTime, referenceFrame.c_str( ),
              aberrationCorrections.c_str( ), observer.c_str( ),
              directSpiceState, &directLightTime );

    // Compare direct and wrapped results for state.
    for ( unsigned int i = 0; i < 6; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( directSpiceState[ i ] * 1000.0, wrapperState[ i ],
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Compare direct and wrapped results for position.
    for ( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( directSpiceState[ i ] * 1000.0 , wrapperPosition[ i ],
                                    std::numeric_limits< double >::epsilon( ) );
    }
}

// Test 3: Test retrieval of rotational state from Spice.
BOOST_AUTO_TEST_CASE( testSpiceWrappers_3 )
{
    using namespace spice_interface;
    using namespace input_output;
    using namespace physical_constants;

    // Create settings at which orientation is to be evaluated.
    const std::string observer = "J2000";
    std::string target = "J2000";
    const double ephemerisTime = 1.0e6;

    // Check if identity matrix (no rotation) is obtained if observer and target frame are equal.
    Eigen::Quaterniond rotationQuaternion = computeRotationQuaternionBetweenFrames(
                observer, target, ephemerisTime );
    Eigen::Vector3d testVector = Eigen::Vector3d::Zero( );
    Eigen::Vector3d testVector2 = Eigen::Vector3d::Zero( );

    for ( unsigned int i = 0; i < 3; i++ )
    {
        testVector = Eigen::Vector3d::Zero( );
        testVector[ i ] = 1.0;
        testVector2 = rotationQuaternion * testVector;
        TUDAT_CHECK_MATRIX_CLOSE( testVector, testVector2,
                                  std::numeric_limits< double >::epsilon( ) );
    }

    // Check if direct and wrapped results of non-trivial rotation are equal.
    target = "IAU_EARTH";

    // Retrieve rotation from wrapper.
    rotationQuaternion = computeRotationQuaternionBetweenFrames(
                observer, target, ephemerisTime );
    Eigen::Matrix3d rotationMatrixDerivative = computeRotationMatrixDerivativeBetweenFrames(
                observer, target, ephemerisTime );

    // Create rotational ephemeris with Spice
    std::shared_ptr< ephemerides::SpiceRotationalEphemeris > spiceRotationalEphemeris =
            std::make_shared< ephemerides::SpiceRotationalEphemeris >(
                observer, target );
    Eigen::Quaterniond rotationQuaternionFromObject = spiceRotationalEphemeris->
            getRotationToTargetFrame( ephemerisTime );
    Eigen::Quaterniond inverseRotationQuaternionFromObject = spiceRotationalEphemeris->
            getRotationToBaseFrame( ephemerisTime );

    Eigen::Matrix3d rotationMatrixDerivativeFromObject = spiceRotationalEphemeris->
            getDerivativeOfRotationToTargetFrame( ephemerisTime );
    Eigen::Matrix3d inverseRotationMatrixDerivativeFromObject = spiceRotationalEphemeris->
            getDerivativeOfRotationToBaseFrame( ephemerisTime );

    // Convert result to Matrix3d for comparison.
    const Eigen::Matrix3d rotationMatrix = Eigen::Matrix3d( rotationQuaternion );
    const Eigen::Matrix3d rotationMatrixFromObject =
            Eigen::Matrix3d( rotationQuaternionFromObject );
    const Eigen::Matrix3d inverseRotationMatrixFromObject =
            Eigen::Matrix3d( inverseRotationQuaternionFromObject );


    // Retrieve rotation directly from spice.
    double stateTransitionMatrix[ 6 ][ 6 ];
    sxform_c( observer.c_str( ), target.c_str( ), ephemerisTime, stateTransitionMatrix );

    double inverseStateTransitionMatrix[ 6 ][ 6 ];
    sxform_c( target.c_str( ), observer.c_str( ), ephemerisTime, inverseStateTransitionMatrix );

    // Check equality of results.
    for ( unsigned int i = 0; i < 3; i++ )
    {
        for ( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( stateTransitionMatrix[ i ][ j ] - rotationMatrix( i, j ),
                               2.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( stateTransitionMatrix[ i ][ j ] -
                               rotationMatrixFromObject( i, j ),
                               5.0 * std::numeric_limits< double >::epsilon( ) );

            BOOST_CHECK_SMALL( stateTransitionMatrix[ i + 3 ][ j ] -
                               rotationMatrixDerivative( i, j ),
                               2.0E-4 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( stateTransitionMatrix[ i + 3 ][ j ] -
                               rotationMatrixDerivativeFromObject( i, j ),
                               2.0E-4 * std::numeric_limits< double >::epsilon( ) );

            BOOST_CHECK_SMALL( inverseStateTransitionMatrix[ i ][ j ] -
                               inverseRotationMatrixFromObject( i, j ),
                               2.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_SMALL( inverseStateTransitionMatrix[ i + 3 ][ j ] -
                               inverseRotationMatrixDerivativeFromObject( i, j ),
                               2.0E-4 * std::numeric_limits< double >::epsilon( ) );

        }
    }

    // Check whether concatenation of forward and backward rotation from object yield no rotation.
    Eigen::Quaterniond forwardBackwardRotation =
            rotationQuaternionFromObject * inverseRotationQuaternionFromObject;
    BOOST_CHECK_SMALL( forwardBackwardRotation.w( ) - 1.0,
                       std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( forwardBackwardRotation.x( ),
                       std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( forwardBackwardRotation.y( ),
                       std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( forwardBackwardRotation.z( ),
                       std::numeric_limits< double >::epsilon( ) );

    // Check automatic state conversion to/from rotating frame.
    Eigen::Vector6d stateOfMoonInJ2000 = spice_interface::getBodyCartesianStateAtEpoch(
                "Moon", "SSB", "J2000", "NONE", 1.0E7 );
    Eigen::Vector6d expectedStateOfMoonWrtEarth = spice_interface::getBodyCartesianStateAtEpoch(
                "Moon", "SSB", "IAU_Earth", "NONE", 1.0E7 );

    Eigen::Vector6d computedStateOfMoonWrtEarth =
            ephemerides::transformStateToTargetFrame(
                stateOfMoonInJ2000, 1.0E7,  spiceRotationalEphemeris );
    Eigen::Vector6d computedStateOfMoonInJ2000 =
            ephemerides::transformStateToInertialOrientation(
                computedStateOfMoonWrtEarth, 1.0E7,  spiceRotationalEphemeris );

    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL(
                    std::fabs( expectedStateOfMoonWrtEarth( i ) - computedStateOfMoonWrtEarth( i ) ), 1.0E-4 );
        BOOST_CHECK_SMALL(
                    std::fabs( expectedStateOfMoonWrtEarth( i + 3 ) - computedStateOfMoonWrtEarth( i + 3 ) ), 1.0E-8 );

        BOOST_CHECK_SMALL(
                    std::fabs( stateOfMoonInJ2000( i ) - computedStateOfMoonInJ2000( i ) ), 1.0E-4 );
        BOOST_CHECK_SMALL(
                    std::fabs( stateOfMoonInJ2000( i + 3 ) - computedStateOfMoonInJ2000( i + 3 ) ), 1.0E-8 );

    }
}

// Test 4: Test retrieval of body properties.
BOOST_AUTO_TEST_CASE( testSpiceWrappers_4 )
{
    using namespace spice_interface;
    using namespace input_output;
    using namespace physical_constants;

    // Retrieve Sun's gravitational parameter from Spice.
    double sunGravitationalParameterSpice = getBodyGravitationalParameter( "Sun" );

    // Set Sun's gravitational parameter as read manually from kernel.
    const double sunGravitationalParameter = 1.327124400420322e+20;

    // Check if results are the same.
    BOOST_CHECK_CLOSE_FRACTION( sunGravitationalParameterSpice, sunGravitationalParameter,
                                std::numeric_limits< double >::epsilon( )  );

    // Check if same result is obtained through general property retrieval function.
    sunGravitationalParameterSpice = getBodyProperties( "Sun", "GM", 1 )[ 0 ] * 1.0e9;
    BOOST_CHECK_CLOSE_FRACTION( sunGravitationalParameterSpice, sunGravitationalParameter,
                                std::numeric_limits< double >::epsilon( )  );

    // Retrieve average Earth radius from Spice through wrapper.
    const double averageEarthRadius = getAverageRadius( "Earth" );

    // Retrieve Earth radii directly from Spice and compute average
    double spiceEarthRadii[ 3 ];
    SpiceInt numberOfReturnParameters = 0;
    bodvrd_c( "Earth","RADII", 3, &numberOfReturnParameters, spiceEarthRadii );
    const double directAverageEarthRadius = 1000.0 * ( spiceEarthRadii[ 0 ] + spiceEarthRadii[ 1 ]
                                                       + spiceEarthRadii[ 2 ] ) / 3.0;

    // Compare average Earth radii.
    BOOST_CHECK_CLOSE_FRACTION( directAverageEarthRadius, averageEarthRadius,
                                std::numeric_limits< double >::epsilon( )  );

    // Check correct conversion of name to NAIF ID.
    const double naifSunId = convertBodyNameToNaifId( "Sun" );
    BOOST_CHECK_EQUAL( 10, naifSunId );

    const double naifMoonId = convertBodyNameToNaifId( "Moon" );
    BOOST_CHECK_EQUAL( 301, naifMoonId );
}

// Test 5: Test Spice Ephemeris class using difference combinations of aberration corrections.
BOOST_AUTO_TEST_CASE( testSpiceWrappers_5 )
{
    using namespace spice_interface;
    using namespace input_output;
    using namespace physical_constants;
    using namespace ephemerides;

    // Create settings at which states are to be evaluated.
    std::string aberrationCorrections = "NONE";
    const std::string observer = "Moon";
    const std::string target = "Mars";
    const std::string referenceFrame = "IAU_EARTH";
    const double ephemerisTime = JULIAN_YEAR;

    SpiceEphemeris spiceEphemeris = SpiceEphemeris( target, observer, 0, 0, 0, referenceFrame );

    Eigen::Vector6d directState = Eigen::Vector6d( );
    Eigen::Vector6d ephemerisState = Eigen::Vector6d( );

    // Check calculated state with no aberration corrections.
    directState = getBodyCartesianStateAtEpoch( target, observer, referenceFrame,
                                                aberrationCorrections, ephemerisTime );
    ephemerisState = spiceEphemeris.getCartesianState( ephemerisTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directState, ephemerisState,
                                       std::numeric_limits< double >::epsilon( ) );

    // Check calculated state with light time correction.
    spiceEphemeris = SpiceEphemeris( target, observer, 0, 1, 0, referenceFrame );
    aberrationCorrections = "LT";
    directState = getBodyCartesianStateAtEpoch( target, observer, referenceFrame,
                                                aberrationCorrections, ephemerisTime );
    ephemerisState = spiceEphemeris.getCartesianState( ephemerisTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directState, ephemerisState,
                                       std::numeric_limits< double >::epsilon( ) );

    // Check calculated state with converged light time correction.
    spiceEphemeris = SpiceEphemeris( target, observer, 0, 1, 1, referenceFrame );
    aberrationCorrections = "CN";
    directState = getBodyCartesianStateAtEpoch( target, observer, referenceFrame,
                                                aberrationCorrections, ephemerisTime );
    ephemerisState = spiceEphemeris.getCartesianState( ephemerisTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directState, ephemerisState,
                                       std::numeric_limits< double >::epsilon( ) );

    // Check calculated state with light time correction and stellar aberration.
    spiceEphemeris = SpiceEphemeris( target, observer, 1, 1, 0, referenceFrame );
    aberrationCorrections = "LT+S";
    directState = getBodyCartesianStateAtEpoch( target, observer, referenceFrame,
                                                aberrationCorrections, ephemerisTime );
    ephemerisState = spiceEphemeris.getCartesianState( ephemerisTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directState, ephemerisState,
                                       std::numeric_limits< double >::epsilon( ) );

    // Check whether exceptions are handled correctly when giving inconsistent input to
    // SpiceEphemeris constrcutor.

    // Test exception when requesting stellar, but not light time aberration.
    bool areExceptionsHandledCorrectly = false;

    try
    {
        spiceEphemeris = SpiceEphemeris( target, observer, 1, 0, 0, referenceFrame );
    }
    catch( std::runtime_error const& )
    {
        areExceptionsHandledCorrectly = true;
    }

    BOOST_CHECK_EQUAL( areExceptionsHandledCorrectly, true );

    // Test exception when requesting stellar, but not light time aberration, with request for
    // converged light-time.
    areExceptionsHandledCorrectly = false;

    try
    {
        spiceEphemeris = SpiceEphemeris( target, observer, 1, 0, 1, referenceFrame );
    }
    catch( std::runtime_error const& )
    {
        areExceptionsHandledCorrectly = true;
    }

    BOOST_CHECK_EQUAL( areExceptionsHandledCorrectly, true );

    // Test exception when not requesting light time aberration, with request for converged
    // light-time.
    areExceptionsHandledCorrectly = false;

    try
    {
        spiceEphemeris = SpiceEphemeris( target, observer, 0, 0, 1, referenceFrame );
    }
    catch( std::runtime_error const& )
    {
        areExceptionsHandledCorrectly = true;
    }

    BOOST_CHECK_EQUAL( areExceptionsHandledCorrectly, true );
}

// Test 6: Compare Spice data with Horizons data.
BOOST_AUTO_TEST_CASE( testSpiceWrappers_6 )
{
    using namespace spice_interface;
    using namespace input_output;
    using namespace physical_constants;

    // Create settings at which states are to be evaluated.
    const std::string aberrationCorrections = "NONE";
    const std::string observer = "Solar System Barycenter";
    const std::string target = "Mars";
    const std::string referenceFrame = "ECLIPJ2000";
    const double julianDay = 2451556.500000000;

    // Get state from wrapper for state:
    const Eigen::Vector6d wrapperState = getBodyCartesianStateAtEpoch(
                target, observer, referenceFrame, aberrationCorrections,
                convertJulianDateToEphemerisTime( julianDay ) );

    // Set state as retrieved from Horizons (see Issue wiki-knowledgebase-spice interface)
    Eigen::Vector6d horizonsState;
    horizonsState << 2.066392047883538e8,
            2.364158324807732e7,
            -4.570656418319555e6,
            -1.850837582360033,
            2.612355357135549e1,
            5.930879066959573e-1;
    horizonsState *= 1000.0;

    // Compare direct and wrapped results for state.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( horizonsState, wrapperState, 1.0e-5 );
}

// Test 7: Loading and clearing kernels.
BOOST_AUTO_TEST_CASE( testSpiceWrappers_7 )
{
    using namespace spice_interface;
    using namespace input_output;

    // Initially clear all Spice kernels.
    clearSpiceKernels( );

    // Get ammount of loaded Spice kernels.
    int spiceKernelsLoaded = getTotalCountOfKernelsLoaded( );

    // Loaded kernels should be 0.
    BOOST_CHECK_EQUAL( spiceKernelsLoaded, 0 );

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Get ammount of loaded Spice kernels.
    spiceKernelsLoaded = getTotalCountOfKernelsLoaded( );

    // Loaded kernels should be 4.
    BOOST_CHECK_EQUAL( spiceKernelsLoaded, 13 );

    // Clear all Spice kernels.
    clearSpiceKernels( );

    // Get ammount of loaded Spice kernels.
    spiceKernelsLoaded = getTotalCountOfKernelsLoaded( );

    // Loaded kernels should be 0.
    BOOST_CHECK_EQUAL( spiceKernelsLoaded, 0 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
