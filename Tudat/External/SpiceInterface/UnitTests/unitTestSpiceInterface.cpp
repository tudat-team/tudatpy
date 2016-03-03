/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120717    D. Dirkx          Creation of file.
 *      121001    M. Ganeff         Added unit test for clearing and loading spice kernels.
 *      140127    D. Dirkx          Adapted for custom Spice kernel folder.
 *
 *    References
 *
 *    Notes
 *      To run this unit tests, a number of spice kernels need to be placed in the
 *      Spice kernel folder, by default External/SpiceInterface/Kernels or the
 *      SPICE_KERNEL_CUSTOM_FOLDER folder set as an argument to CMake or in UserSetings.txt.
 *      The required kernels are:
 *           de421.bsp
 *           pck00009.tpc
 *           naif0009.tls
 *           de-403-masses.tpc
 *      They can be found in a single zip file on the wiki at
 *      http://tudat.tudelft.nl/projects/tudat/wiki/SpiceInterface/ on the Tudat website or,
 *      alternatively, on the NAIF server at ftp://naif.jpl.nasa.gov/pub/naif/generic_kernels/.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/bind.hpp>
#include <boost/exception/all.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Basics/testMacros.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/External/SpiceInterface/spiceRotationalEphemeris.h"

#include <limits>
#include <iostream>
#include <stdexcept>

namespace tudat
{
namespace unit_tests
{

using basic_mathematics::Vector6d;

BOOST_AUTO_TEST_SUITE( test_spice_wrappers )

// Test 1: Test Julian day <-> Ephemeris time conversions at J2000.
BOOST_AUTO_TEST_CASE( testSpiceWrappers_1 )
{
    using namespace spice_interface;
    using namespace input_output;
    using namespace physical_constants;

    // Check if spice kernels exist.
    if ( !boost::filesystem::exists( getSpiceKernelPath( ) + "de421.bsp" ) )
    {
        std::cerr << "SPICE kernel de421.bsp not found in "<<getSpiceKernelPath( )
                  << std::endl;
    }

    if ( !boost::filesystem::exists( getSpiceKernelPath( ) + "pck00009.tpc" ) )
    {
        std::cerr << "SPICE kernel pck00009.tpc not found in "<<getSpiceKernelPath( )
                  << std::endl;
    }

    if ( !boost::filesystem::exists( getSpiceKernelPath( ) + "naif0009.tls" ) )
    {
        std::cerr << "SPICE kernel naif0009.tls not found in "<<getSpiceKernelPath( )
                  << std::endl;
    }

    if ( !boost::filesystem::exists( getSpiceKernelPath( ) + "de-403-masses.tpc" ) )
    {
        std::cerr << "SPICE kernel de-403-masses.tpc not found in "<<
                     getSpiceKernelPath( )<< std::endl;
    }

    // Load spice kernels.
    loadSpiceKernelInTudat( getSpiceKernelPath( ) + "pck00009.tpc" );
    loadSpiceKernelInTudat( getSpiceKernelPath( ) + "de-403-masses.tpc" );
    loadSpiceKernelInTudat( getSpiceKernelPath( ) + "de421.bsp" );
    loadSpiceKernelInTudat( getSpiceKernelPath( ) + "naif0009.tls" );

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
    const std::string abberationCorrections = "NONE";
    const std::string observer = "Solar System Barycenter";
    const std::string target = "Mars";
    const std::string referenceFrame = "J2000";
    const double ephemerisTime = 1.0E6;

    // Get state from wrapper for state:
    const Vector6d wrapperState = getBodyCartesianStateAtEpoch(
                target, observer, referenceFrame, abberationCorrections, ephemerisTime );

    // Get position from wrapper for position:
    const Eigen::Vector3d wrapperPosition = getBodyCartesianPositionAtEpoch(
                target, observer, referenceFrame, abberationCorrections, ephemerisTime );

    // Get state directly from spice.
    double directSpiceState[ 6 ];
    double directLightTime;
    spkezr_c( target.c_str( ), ephemerisTime, referenceFrame.c_str( ),
              abberationCorrections.c_str( ), observer.c_str( ),
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
    boost::shared_ptr< ephemerides::SpiceRotationalEphemeris > spiceRotationalEphemeris =
            boost::make_shared< ephemerides::SpiceRotationalEphemeris >(
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
                               4.0 * std::numeric_limits< double >::epsilon( ) );

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
    const double sunGravitationalParameter = 132712440023.310 * 1.0e9;

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
    std::string abberationCorrections = "NONE";
    const std::string observer = "Moon";
    const std::string target = "Mars";
    const std::string referenceFrame = "IAU_EARTH";
    const double ephemerisTime = JULIAN_YEAR;

    SpiceEphemeris spiceEphemeris = SpiceEphemeris( target, observer, 0, 0, 0, referenceFrame );

    basic_mathematics::Vector6d directState = basic_mathematics::Vector6d( );
    basic_mathematics::Vector6d ephemerisState = basic_mathematics::Vector6d( );

    // Check calculated state with no aberration corrections.
    directState = getBodyCartesianStateAtEpoch( target, observer, referenceFrame,
                                                abberationCorrections, ephemerisTime );
    ephemerisState = spiceEphemeris.getCartesianStateFromEphemeris( ephemerisTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directState, ephemerisState,
                                       std::numeric_limits< double >::epsilon( ) );

    // Check calculated state with light time correction.
    spiceEphemeris = SpiceEphemeris( target, observer, 0, 1, 0, referenceFrame );
    abberationCorrections = "LT";
    directState = getBodyCartesianStateAtEpoch( target, observer, referenceFrame,
                                                abberationCorrections, ephemerisTime );
    ephemerisState = spiceEphemeris.getCartesianStateFromEphemeris( ephemerisTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directState, ephemerisState,
                                       std::numeric_limits< double >::epsilon( ) );

    // Check calculated state with converged light time correction.
    spiceEphemeris = SpiceEphemeris( target, observer, 0, 1, 1, referenceFrame );
    abberationCorrections = "CN";
    directState = getBodyCartesianStateAtEpoch( target, observer, referenceFrame,
                                                abberationCorrections, ephemerisTime );
    ephemerisState = spiceEphemeris.getCartesianStateFromEphemeris( ephemerisTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( directState, ephemerisState,
                                       std::numeric_limits< double >::epsilon( ) );

    // Check calculated state with light time correction and stellar aberration.
    spiceEphemeris = SpiceEphemeris( target, observer, 1, 1, 0, referenceFrame );
    abberationCorrections = "LT+S";
    directState = getBodyCartesianStateAtEpoch( target, observer, referenceFrame,
                                                abberationCorrections, ephemerisTime );
    ephemerisState = spiceEphemeris.getCartesianStateFromEphemeris( ephemerisTime );
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

    catch( std::runtime_error )
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
    catch( std::runtime_error )
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

    catch( std::runtime_error )
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
    const std::string abberationCorrections = "NONE";
    const std::string observer = "Solar System Barycenter";
    const std::string target = "Mars";
    const std::string referenceFrame = "ECLIPJ2000";
    const double julianDay = 2451556.500000000;

    // Get state from wrapper for state:
    const Vector6d wrapperState = getBodyCartesianStateAtEpoch(
                target, observer, referenceFrame, abberationCorrections,
                convertJulianDateToEphemerisTime( julianDay ) );

    // Set state as retrieved from Horizons (see Issue wiki-knowledgebase-spice interface)
    Vector6d horizonsState;
    horizonsState << 2.066392047883538e8,
            2.364158324807732e7,
            -4.570656418319555e6,
            -1.850837582360033,
            2.612355357135549e1,
            5.930879066959573e-1;
    horizonsState *= 1000.0;

    // Compare direct and wrapped results for state.
    for ( int i = 0; i < 6; i++ )
    {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( horizonsState, wrapperState, 5.0e-7 );
    }
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
    // Load spice kernels.
    loadSpiceKernelInTudat( getSpiceKernelPath( ) + "de421.bsp" );
    loadSpiceKernelInTudat( getSpiceKernelPath( ) + "pck00009.tpc" );
    loadSpiceKernelInTudat( getSpiceKernelPath( ) + "naif0009.tls" );
    loadSpiceKernelInTudat( getSpiceKernelPath( ) + "de-403-masses.tpc" );

    // Get ammount of loaded Spice kernels.
    spiceKernelsLoaded = getTotalCountOfKernelsLoaded( );

    // Loaded kernels should be 4.
    BOOST_CHECK_EQUAL( spiceKernelsLoaded, 4 );

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
