/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Easy calculation. Gravitational Acceleration Tutorial,
 *          http://easycalculation.com/physics/classical-physics
 *          /learn-gravitational-acceleration.php, last accessed: 26th February, 2012.
 *      MathWorks. gravityzonal, MATLAB 2012b, 2012.
 *      Melman, J. Propagate software, J.C.P.Melman@tudelft.nl, 2012.
 *      Ronse, A. A parametric study of space debris impact footprints, MSc thesis, Delft
 *          University of Technlogy, Delft, The Netherlands, in preparation.
 *
 */

#define BOOST_TEST_MAIN

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/make_shared.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2GravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2J3GravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2J3J4GravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/UnitTests/planetTestData.h"

namespace tudat
{
namespace unit_tests
{

using namespace gravitation;

typedef std::map< int, double > KeyIntValueDoubleMap;

BOOST_AUTO_TEST_SUITE( test_gravitational_acceleration )

//! Test if gravitational acceleration is computed correctly.
BOOST_AUTO_TEST_CASE( testGravitationalAcceleration )
{
    // Test 1: compute gravitational acceleration exerted on surface of Earth
    //         (Easy calculation, 2012).
    {
        // Set gravitational parameter of Earth [m^3 s^-2].
        const double gravitationalParameterOfEarth = 6.6726e-11 * 5.9742e24;

        // Set position vector of Earth [m].
        const Eigen::Vector3d positionOfEarth = Eigen::Vector3d::Zero( );

        // Set position vector on Earth surface [m].
        const Eigen::Vector3d positionOnEarthSurface( 6.3781e6, 0.0, 0.0 );

        // Compute gravitational accelerating acting on Earth's surface [N].
        const Eigen::Vector3d gravitationalAccelerationExertedAtEarthSurface
                = computeGravitationalAcceleration(
                    positionOnEarthSurface, gravitationalParameterOfEarth, positionOfEarth );

        // Check if computed gravitational force matches expected value.
        BOOST_CHECK_CLOSE_FRACTION( 9.8, gravitationalAccelerationExertedAtEarthSurface.norm( ),
                                    1.0e-4 );
    }

    // Test 2: compute gravitational acceleration exerted on Lunar surface
    //         (Easy calculation, 2012).
    {
        // Set universal gravitational constant [m^3 kg^-1 s^-2].
        const double universalGravitationalConstant = 6.6726e-11;

        // Set mass of Moon [kg].
        const double massOfMoon = 7.36e22;

        // Set position vector of Moon [m].
        const Eigen::Vector3d positionOfMoon( 12.65, 0.23, -45.78 );

        // Set position vector on surface of Moon [m].
        const Eigen::Vector3d positionOfLunarSurface( 0.0, 1735771.89, 0.0 );

        // Compute gravitational accelerating acting on Lunar surface [N].
        const Eigen::Vector3d gravitationalAccelerationExertedAtLunarSurface
                = computeGravitationalAcceleration(
                    universalGravitationalConstant, positionOfLunarSurface,
                    massOfMoon, positionOfMoon );

        // Check if computed gravitational force matches expected value.
        BOOST_CHECK_CLOSE_FRACTION( 1.63, gravitationalAccelerationExertedAtLunarSurface.norm( ),
                                    1.0e-6 );
    }
}

//! Test if gravitational acceleration sum due to zonal terms is computed correctly using MATLAB.
BOOST_AUTO_TEST_CASE( testGravitationalAccelarationSumZonalMatlab )
{
    // These tests check if total acceleration due to zonal terms is computed correctly by
    // comparing to output generated using gravityzonal() function in MATLAB (Mathworks, 2012).
    // The planet data used is obtained from the documentation for the gravityzonal() function.

    // Get planet test data.
    std::vector< PlanetTestData > planetData = getPlanetMatlabTestData( );

    // Loop over all planet test data and recompute the results using Tudat code. Check that the
    // values computed match MATLAB's output (Mathworks, 2012).
    for ( unsigned int planet = 0; planet < planetData.size( ); planet++ )
    {
        for ( unsigned int body1 = 0; body1 < planetData.at( planet ).body1Positions.size( );
              body1++ )
        {
            for ( unsigned int body2 = 0; body2 < planetData.at( planet ).body2Positions.size( );
                  body2++ )
            {
                // Compute central gravitational acceleration term [m s^-2].
                const Eigen::Vector3d computedCentralAcceleration
                        = computeGravitationalAcceleration(
                            planetData.at( planet ).body2Positions.at( body2 ),
                            planetData.at( planet ).gravitationalParameter,
                            planetData.at( planet ).body1Positions.at( body1 ) );

                // Check that the computed central gravitational acceleration matches the expected
                // values.
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                            planetData.at( planet ).expectedAcceleration[ body1 ][ body2 ][
                            central ],
                            computedCentralAcceleration,
                            1.0e-15 );

                // Declare zonal coefficients used.
                KeyIntValueDoubleMap zonalCoefficients;

                // Loop over all available zonal gravity field coefficients.
                for ( KeyIntValueDoubleMap::iterator zonalCoefficientIterator
                      = planetData.at( planet ).zonalCoefficients.begin( );
                      zonalCoefficientIterator
                      != planetData.at( planet ).zonalCoefficients.end( );
                      zonalCoefficientIterator++ )
                {
                    // Add current zonal coefficient to local list.
                    zonalCoefficients[ zonalCoefficientIterator->first ]
                            = zonalCoefficientIterator->second;

                    // Compute gravitational acceleration sum [m s^-2].
                    const Eigen::Vector3d computedAccelerationSum
                            = computeGravitationalAccelerationZonalSum(
                                planetData.at( planet ).body2Positions.at( body2 ),
                                planetData.at( planet ).gravitationalParameter,
                                planetData.at( planet ).equatorialRadius,
                                zonalCoefficients,
                                planetData.at( planet ).body1Positions.at( body1 ) );

                    // Check that computed gravitational acceleration sum matches expected values.
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                planetData.at( planet ).expectedAcceleration[ body1 ][ body2 ][
                                zonalCoefficientIterator->first ],
                            computedAccelerationSum,
                            1.0e-15 );
                }
            }
        }
    }
}

//! Test if wrapper classes compute gravitational acceleration due to zonal terms correctly using
//! MATLAB.
BOOST_AUTO_TEST_CASE( testGravitationalAccelerationZonalSumWrapperClassesMatlab )
{
    // These tests check if total acceleration due to zonal terms is computed correctly by wrapper
    // class comparing to output generated using gravityzonal() function in MATLAB
    // (Mathworks, 2012). The planet data used is obtained from the documentation for the
    // gravityzonal() function. This check is for consistency purposes, since the wrapper class
    // wraps the free functions that are all tested too.

    // Short-cuts.
    using namespace gravitation;

    // Get planet test data.
    std::vector< PlanetTestData > planetData = getPlanetMatlabTestData( );

    // Loop over all planet test data and recompute the results using Tudat code. Check that the
    // values computed match MATLAB's output (Mathworks, 2012).
    for ( unsigned int planet = 0; planet < planetData.size( ); planet++ )
    {
        for ( unsigned int body1 = 0; body1 < planetData.at( planet ).body1Positions.size( );
              body1++ )
        {
            for ( unsigned int body2 = 0; body2 < planetData.at( planet ).body2Positions.size( );
                  body2++ )
            {
                // Test central gravitational acceleration wrapper class.
                {
                    // Declare central acceleration wrapper class object.
                    CentralGravitationalAccelerationModel3dPointer centralGravity
                            = std::make_shared< CentralGravitationalAccelerationModel3d >(
                                [ & ]( ){ return
                                    planetData.at( planet ).body2Positions.at( body2 ); },
                                planetData.at( planet ).gravitationalParameter,
                                [ & ]( ){ return
                                    planetData.at( planet ).body1Positions.at( body1 ); } );

                    // Compute central gravitational acceleration term [m s^-2].
                    const Eigen::Vector3d computedCentralAcceleration
                            = centralGravity->getAcceleration( );

                    // Check that the computed central gravitational acceleration matches the
                    // expected values.
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                planetData.at( planet ).expectedAcceleration[ body1 ][ body2 ][
                                central ],
                                computedCentralAcceleration,
                                1.0e-15 );
                }

                // Test central + J2 gravitational acceleration wrapper class.
                if ( planetData.at( planet ).zonalCoefficients.rbegin( )->first == 2 )
                {
                    // Declare central + J2 acceleration wrapper class object.
                    CentralJ2GravitationalAccelerationModelPointer centralJ2Gravity
                            = std::make_shared< CentralJ2GravitationalAccelerationModel >(
                                [ & ]( ){ return
                                    planetData.at( planet ).body2Positions.at( body2 ); },
                                planetData.at( planet ).gravitationalParameter,
                                planetData.at( planet ).equatorialRadius,
                                planetData.at( planet ).zonalCoefficients[ 2 ],
                            [ & ]( ){ return
                                planetData.at( planet ).body1Positions.at( body1 ); } );

                    // Compute gravitational acceleration sum [m s^-2].
                    const Eigen::Vector3d computedCentralJ2AccelerationSum
                            = centralJ2Gravity->getAcceleration( );

                    // Check that computed central + J2 gravitational acceleration matches
                    // expected values.
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                planetData.at( planet ).expectedAcceleration[
                                body1 ][ body2 ][ 2 ],
                            computedCentralJ2AccelerationSum,
                            1.0e-15 );
                }

                // Test central + J2 + J3 gravitational acceleration wrapper class.
                if ( planetData.at( planet ).zonalCoefficients.rbegin( )->first == 3 )
                {
                    // Declare central + J2 + J3 acceleration wrapper class object.
                    CentralJ2J3GravitationalAccelerationModelPointer centralJ2J3Gravity
                            = std::make_shared< CentralJ2J3GravitationalAccelerationModel >(
                                [ & ]( ){ return
                                    planetData.at( planet ).body2Positions.at( body2 ); },
                                planetData.at( planet ).gravitationalParameter,
                                planetData.at( planet ).equatorialRadius,
                                planetData.at( planet ).zonalCoefficients[ 2 ],
                            planetData.at( planet ).zonalCoefficients[ 3 ],
                            [ & ]( ){ return
                                planetData.at( planet ).body1Positions.at( body1 ); } );

                    // Compute gravitational acceleration sum [m s^-2].
                    const Eigen::Vector3d computedCentralJ2J3AccelerationSum
                            = centralJ2J3Gravity->getAcceleration( );

                    // Check that computed central + J2 + J3 gravitational acceleration matches
                    // expected values.
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                planetData.at( planet ).expectedAcceleration[
                                body1 ][ body2 ][ 3 ],
                            computedCentralJ2J3AccelerationSum,
                            1.0e-15 );
                }

                // Test central + J2 + J3 + J4 gravitational acceleration wrapper class.
                if ( planetData.at( planet ).zonalCoefficients.rbegin( )->first == 4 )
                {
                    // Declare pointer to wrapper class object.
                    CentralJ2J3J4GravitationalAccelerationModelPointer centralJ2J3J4Gravity;

                    // Check if only J2 and J4 are given; if so set J3 to 0.0.
                    if ( planetData.at( planet ).zonalCoefficients.size( ) == 2 )
                    {
                        // Declare central + J2 + J3 + J4 acceleration wrapper class object.
                        centralJ2J3J4Gravity = std::make_shared<
                                CentralJ2J3J4GravitationalAccelerationModel >(
                                    [ & ]( ){ return
                                        planetData.at( planet ).body2Positions.at( body2 ); },
                                    planetData.at( planet ).gravitationalParameter,
                                    planetData.at( planet ).equatorialRadius,
                                    planetData.at( planet ).zonalCoefficients[ 2 ],
                                0.0,
                                planetData.at( planet ).zonalCoefficients[ 4 ],
                                [ & ]( ){ return
                                    planetData.at( planet ).body1Positions.at( body1 ); } );
                    }

                    // Else, include given J3.
                    else
                    {
                        // Declare central + J2 + J3 + J4 acceleration wrapper class object.
                        centralJ2J3J4Gravity = std::make_shared<
                                CentralJ2J3J4GravitationalAccelerationModel >(
                                    [ & ]( ){ return
                                        planetData.at( planet ).body2Positions.at( body2 ); },
                                    planetData.at( planet ).gravitationalParameter,
                                    planetData.at( planet ).equatorialRadius,
                                    planetData.at( planet ).zonalCoefficients[ 2 ],
                                planetData.at( planet ).zonalCoefficients[ 3 ],
                                planetData.at( planet ).zonalCoefficients[ 4 ],
                                [ & ]( ){ return
                                    planetData.at( planet ).body1Positions.at( body1 ); } );
                    }

                    // Compute gravitational acceleration sum [m s^-2].
                    const Eigen::Vector3d computedCentralJ2J3J4AccelerationSum
                            = centralJ2J3J4Gravity->getAcceleration( );

                    // Check that computed central + J2 gravitational acceleration matches
                    // expected values.
                    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                                planetData.at( planet ).expectedAcceleration[
                                body1 ][ body2 ][ 4 ],
                            computedCentralJ2J3J4AccelerationSum,
                            1.0e-15 );
                }
            }
        }
    }
}

//! Test if gravitational acceleration due to zonal terms is computed correctly (Melman, 2012).
BOOST_AUTO_TEST_CASE( testGravitationalAccelarationZonalMelman )
{
    typedef Eigen::Vector3d ( *GravitationalAccelerationPointer )(
                const Eigen::Vector3d&, const double, const double, const double,
                const Eigen::Vector3d& );

    // These tests check if acceleration due to zonal terms is computed correctly by comparing to
    // output generated by (Melman, 2012).

    // Get planet test data.
    PlanetTestData earthData = getEarthMelmanTestData( );

    // Set map of function pointers for zonal coefficients.
    std::map< int, GravitationalAccelerationPointer > zonalGravitationalAccelerationPointers
            = { { 2, &computeGravitationalAccelerationDueToJ2 }, { 3, &computeGravitationalAccelerationDueToJ3 },
                { 4, &computeGravitationalAccelerationDueToJ4 } };

    // Loop over all planet test data and recompute the results using Tudat code. Check that the
    // values computed match results obtained by (Melman, 2012).
    for ( unsigned int body2 = 0; body2 < earthData.body2Positions.size( ); body2++ )
    {
        // Loop over all available zonal gravity field coefficients.
        for ( KeyIntValueDoubleMap::iterator zonalCoefficientIterator
              = earthData.zonalCoefficients.begin( );
              zonalCoefficientIterator != earthData.zonalCoefficients.end( );
              zonalCoefficientIterator++ )
        {
            // Compute gravitational acceleration due to given zonal term [m s^-2].
            const Eigen::Vector3d computedZonalGravitationalAcceleration
                    = zonalGravitationalAccelerationPointers[ zonalCoefficientIterator->first ](
                        earthData.body2Positions.at( body2 ),
                        earthData.gravitationalParameter,
                        earthData.equatorialRadius,
                        zonalCoefficientIterator->second,
                        earthData.body1Positions.at( 0 ) );

            // Check that computed gravitational acceleration sum matches expected values.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        earthData.expectedAcceleration[ 0 ][ body2 ][
                    zonalCoefficientIterator->first ],
                    computedZonalGravitationalAcceleration,
                    1.0e-14 );
        }
    }
}

//! Test if gravitational acceleration due to zonal terms is computed correctly (Ronse, 2012).
BOOST_AUTO_TEST_CASE( testGravitationalAccelarationZonalRonse )
{
    typedef Eigen::Vector3d ( *GravitationalAccelerationPointer )(
                const Eigen::Vector3d&, const double, const double, const double,
                const Eigen::Vector3d& );

    // These tests check if acceleration due to zonal terms is computed correctly by comparing to
    // output generated by (Ronse, 2012).

    // Get planet test data.
    PlanetTestData earthData = getEarthRonseTestData( );

    // Set map of function pointers for zonal coefficients.
    std::map< int, GravitationalAccelerationPointer > zonalGravitationalAccelerationPointers
            = { { 2, &computeGravitationalAccelerationDueToJ2 }, { 3, &computeGravitationalAccelerationDueToJ3 },
                { 4, &computeGravitationalAccelerationDueToJ4 } };

    // Loop over all planet test data and recompute the results using Tudat code. Check that the
    // values computed match results obtained by (Ronse, 2012).
    for ( unsigned int body2 = 0; body2 < earthData.body2Positions.size( ); body2++ )
    {
        // Compute central gravitational acceleration term [m s^-2].
        const Eigen::Vector3d computedCentralAcceleration
                = computeGravitationalAcceleration(
                    earthData.body2Positions.at( body2 ),
                    earthData.gravitationalParameter,
                    earthData.body1Positions.at( 0 ) );

        // Check that the computed central gravitational acceleration matches the expected
        // values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    earthData.expectedAcceleration[ 0 ][ body2 ][ central ],
                computedCentralAcceleration,
                1.0e-15 );

        // Loop over all available zonal gravity field coefficients.
        for ( KeyIntValueDoubleMap::iterator zonalCoefficientIterator
              = earthData.zonalCoefficients.begin( );
              zonalCoefficientIterator != earthData.zonalCoefficients.end( );
              zonalCoefficientIterator++ )
        {
            // Compute gravitational acceleration due to given zonal term [m s^-2].
            const Eigen::Vector3d computedZonalGravitationalAcceleration
                    = zonalGravitationalAccelerationPointers[ zonalCoefficientIterator->first ](
                        earthData.body2Positions.at( body2 ),
                        earthData.gravitationalParameter,
                        earthData.equatorialRadius,
                        zonalCoefficientIterator->second,
                        earthData.body1Positions.at( 0 ) );

            // Check that computed gravitational acceleration sum matches expected values.
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( earthData.expectedAcceleration[ 0 ][ body2 ][
                    zonalCoefficientIterator->first ],
                    computedZonalGravitationalAcceleration,
                    1.0e-13 );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
