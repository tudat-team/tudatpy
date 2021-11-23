/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      NASA, Goddard Spaceflight Center. Orbit Determination Toolbox (ODTBX), NASA - GSFC Open
 *          Source Software, http://opensource.gsfc.nasa.gov/projects/ODTBX/, last accessed:
 *          31st January, 2012.
 *      Fortescue, P. W., et al. Spacecraft systems engineering, Third Edition,
 *          Wiley, England, 2003.
 *      Vallado, D. A., McClain, W. D. Fundamentals of astro and applications, 2nd Edition,
 *          Kluwer Academic Publishers, The Netherlands, 2004.
 *      Harvard. Minor Planet Center.
 *          http://scully.cfa.harvard.edu/cgi-bin/returnprepeph.cgi?d=d&o=02060, last accessed:
 *          1st February, 2012.
 *      Wikipedia. Geostationary orbit, http://en.wikipedia.org/wiki/Geostationary_orbit, last
 *      accessed: 1st February, 2012, last modified: 29th January, 2012.
 *      Rocket and Space Technology. Example problems, http://www.braeunig.us/space/problem.htm,
 *          last accessed: 4th February, 2012.
 *      Jenab. http://jenab6.livejournal.com/15054.html, last accessed: 4th February, 2012, last
 *          modified: 6th August, 2008.
 *      Advanced Concepts Team, ESA. Keplerian Toolbox, http://sourceforge.net/projects/keptoolbox,
 *          last accessed: 21st April, 2012.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>
#include <iostream>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/basics/basicTypedefs.h"

namespace tudat
{
namespace unit_tests
{

using namespace mathematical_constants;

BOOST_AUTO_TEST_SUITE( test_orbital_element_conversions )

//! Test if conversion from Keplerian elements to Cartesian elements is working correctly.
BOOST_AUTO_TEST_CASE( testKeplerianToCartesianElementConversion )
{
    // Case 1: Elliptical orbit around the Earth.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Using declarations.
        using namespace orbital_element_conversions;

        // Set Earth gravitational parameter [m^3/s^2] .
        const double earthGravitationalParameter = 3.986004415e14;

        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        Eigen::Matrix< double, 6, 1 > keplerianElements ;
        keplerianElements( semiMajorAxisIndex ) = 8000.0 * 1000.0;
        keplerianElements( eccentricityIndex ) = 0.23;
        keplerianElements( inclinationIndex ) = 20.6 / 180.0 * PI;
        keplerianElements( argumentOfPeriapsisIndex ) = 274.78 / 180.0 * PI;
        keplerianElements( longitudeOfAscendingNodeIndex ) = 108.77 / 180.0 * PI;
        keplerianElements( trueAnomalyIndex ) = 46.11 / 180.0 * PI;

        // Set expected Cartesian elements [m,m,m,m/s,m/s,m/s].
        Eigen::Matrix< double, 6, 1 > expectedCartesianElements;
        expectedCartesianElements( xCartesianPositionIndex ) = 2.021874804243437e6;
        expectedCartesianElements( yCartesianPositionIndex ) = 6.042523817035284e6;
        expectedCartesianElements( zCartesianPositionIndex ) = -1.450371183512575e6;
        expectedCartesianElements( xCartesianVelocityIndex ) = -7.118283509842652e3;
        expectedCartesianElements( yCartesianVelocityIndex ) = 4.169050171542199e3;
        expectedCartesianElements( zCartesianVelocityIndex ) = 2.029066072016241e3;

        // Compute Cartesian elements.
        Eigen::Matrix< double, 6, 1 > computedCartesianElements;
        computedCartesianElements = orbital_element_conversions::
                convertKeplerianToCartesianElements( keplerianElements,
                                                     earthGravitationalParameter );

        // Check if computed Cartesian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCartesianElements,
                                           computedCartesianElements, 1.0e-15 );
    }

    // Case 2: Equatorial, circular orbit around Mars.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Using declarations.
        using namespace orbital_element_conversions;

        // Set Mars gravitational parameter [m^3/s^2].
        const double marsGravitationalParameter = 4.2828018915e13;

        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        Eigen::Matrix< double, 6, 1 > keplerianElements;
        keplerianElements( semiMajorAxisIndex ) = 9201.61 * 1000.0;
        keplerianElements( eccentricityIndex ) = 0.0;
        keplerianElements( inclinationIndex ) = 0.0;
        keplerianElements( argumentOfPeriapsisIndex ) = 12.54 / 180.0 * PI;
        keplerianElements( longitudeOfAscendingNodeIndex ) = 201.55 / 180.0 * PI;
        keplerianElements( trueAnomalyIndex ) = -244.09 / 180.0 * PI;

        // Set expected Cartesian elements [m,m,m,m/s,m/s,m/s].
        Eigen::Matrix< double, 6, 1 > expectedCartesianElements;
        expectedCartesianElements( xCartesianPositionIndex ) = 7.968828015716932e6;
        expectedCartesianElements( yCartesianPositionIndex ) = -4.600804999999997e6;
        expectedCartesianElements( zCartesianPositionIndex ) = 0.0;
        expectedCartesianElements( xCartesianVelocityIndex ) = 1.078703495685965e3;
        expectedCartesianElements( yCartesianVelocityIndex ) = 1.868369260830248e3;
        expectedCartesianElements( zCartesianVelocityIndex ) = 0.0;

        // Compute Cartesian elements.
        Eigen::Matrix< double, 6, 1 > computedCartesianElements;
        computedCartesianElements = orbital_element_conversions::
                convertKeplerianToCartesianElements( keplerianElements,
                                                     marsGravitationalParameter );

        // Check if computed Cartesian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCartesianElements,
                                           computedCartesianElements, 1.0e-15 );
    }

    // Case 3: Hyperbolic orbit around the Sun.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Using declarations.
        using namespace orbital_element_conversions;

        // Set Sun gravitational parameter [m^3/s^2].
        const double sunGravitationalParameter = 1.32712440018e20;

        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        Eigen::Matrix< double, 6, 1 > keplerianElements;
        keplerianElements( semiMajorAxisIndex ) = -4.5e11;
        keplerianElements( eccentricityIndex ) = 2.3;
        keplerianElements( inclinationIndex ) = 25.5 / 180.0 * PI;
        keplerianElements( argumentOfPeriapsisIndex ) = 156.11 / 180.0 * PI;
        keplerianElements( longitudeOfAscendingNodeIndex ) = -215.03 / 180.0 * PI;
        keplerianElements( trueAnomalyIndex ) = 123.29 / 180.0 * PI;

        // Set expected Cartesian elements [m,m,m,m/s,m/s,m/s].
        Eigen::Matrix< double, 6, 1 > expectedCartesianElements;
        expectedCartesianElements( xCartesianPositionIndex ) = -2.776328224174438e12;
        expectedCartesianElements( yCartesianPositionIndex ) = -6.053823869632723e12;
        expectedCartesianElements( zCartesianPositionIndex ) = 3.124576293512172e12;
        expectedCartesianElements( xCartesianVelocityIndex ) = 7.957674684798018e3;
        expectedCartesianElements( yCartesianVelocityIndex ) = 1.214817382001788e4;
        expectedCartesianElements( zCartesianVelocityIndex ) = -6.923442392618828e3;

        // Compute Cartesian elements.
        Eigen::Matrix< double, 6, 1 > computedCartesianElements;
        computedCartesianElements = orbital_element_conversions::
                convertKeplerianToCartesianElements( keplerianElements,
                                                     sunGravitationalParameter );

        // Check if computed Cartesian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCartesianElements,
                                           computedCartesianElements, 1.0e-15 );
    }

    // Case 4: Parabolic orbit around the Earth.
    // Earth-orbiting satellite example (Rocket and Space Technology, 2012).
    {
        // Using declarations.
        using namespace orbital_element_conversions;

        // Set Earth gravitational parameter [m^3/s^2].
        const double earthGravitationalParameter = 3.986005e14;

        // Set Keplerian elements [m,-,rad,rad,rad,rad].
        Eigen::Matrix< double, 6, 1 > keplerianElements;
        keplerianElements( semiLatusRectumIndex ) = 2.0 * 6678140.0;
        keplerianElements( eccentricityIndex ) = 1.0;
        keplerianElements( inclinationIndex ) = 45.0 / 180.0 * PI;
        keplerianElements( argumentOfPeriapsisIndex ) = 0.0;
        keplerianElements( longitudeOfAscendingNodeIndex ) = 0.0;
        keplerianElements( trueAnomalyIndex ) = 0.0;

        // Set expected escape velocity [m/s].
        const double expectedEscapeVelocity = 10926.0;

        // Compute Cartesian elements.
        Eigen::Matrix< double, 6, 1 > computedCartesianElements;
        computedCartesianElements = orbital_element_conversions::
                convertKeplerianToCartesianElements( keplerianElements,
                                                     earthGravitationalParameter );

        // Check if computed escape veloicty matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedEscapeVelocity,
                                    computedCartesianElements.segment( 3, 3 ).norm( ), 1.0e-4 );
    }
}

//! Test if conversion from Cartesian elements to Keplerian elements is working correctly, using
//! benchmark data.
BOOST_AUTO_TEST_CASE( testCartesianToKeplerianElementConversionBenchmark )
{
    // Case 1: Elliptical orbit around the Earth.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Using declarations.
        using namespace orbital_element_conversions;

        // Earth gravitational parameter.
        const double earthGravitationalParameter = 3.986004415e14;

        // Set Cartesian elements.
        Eigen::Matrix< double, 6, 1 > cartesianElements;
        cartesianElements( xCartesianPositionIndex ) = 3.75e6;
        cartesianElements( yCartesianPositionIndex ) = 4.24e6;
        cartesianElements( zCartesianPositionIndex ) = -1.39e6;
        cartesianElements( xCartesianVelocityIndex ) = -4.65e3;
        cartesianElements( yCartesianVelocityIndex ) = -2.21e3;
        cartesianElements( zCartesianVelocityIndex ) = 1.66e3;

        // Set expected Keplerian elements.
        Eigen::Matrix< double, 6, 1 > expectedKeplerianElements;
        expectedKeplerianElements( semiMajorAxisIndex ) = 3.707478199246163e6;
        expectedKeplerianElements( eccentricityIndex ) = 0.949175203660321;
        expectedKeplerianElements( inclinationIndex ) = 0.334622356632438;
        expectedKeplerianElements( argumentOfPeriapsisIndex ) = 2.168430616511167;
        expectedKeplerianElements( longitudeOfAscendingNodeIndex ) = 1.630852596545341;
        expectedKeplerianElements( trueAnomalyIndex ) = 3.302032232567084;

        // Compute Keplerian elements.
        Eigen::Matrix< double, 6, 1 > computedKeplerianElements;
        computedKeplerianElements = orbital_element_conversions::
                convertCartesianToKeplerianElements( cartesianElements,
                                                     earthGravitationalParameter );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, 1.0e-14 );
    }

    // Case 2: Equatorial, circular orbit around Venus.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Using declarations.
        using namespace orbital_element_conversions;

        // Venus gravitational parameter.
        const double venusGravitationalParameter = 3.2485504415e14;

        // Set Cartesian elements.
        Eigen::Matrix< double, 6, 1 > cartesianElements;
        cartesianElements( xCartesianPositionIndex ) = 5.580537430785387e6;
        cartesianElements( yCartesianPositionIndex ) = 2.816487703435473e6;
        cartesianElements( zCartesianPositionIndex ) = 0.0;
        cartesianElements( xCartesianVelocityIndex ) = -3.248092722413634e3;
        cartesianElements( yCartesianVelocityIndex ) = 6.435711753323540e3;
        cartesianElements( zCartesianVelocityIndex ) = 0.0;

        // Set expected Keplerian elements.
        Eigen::Matrix< double, 6, 1 > expectedKeplerianElements;
        expectedKeplerianElements( semiMajorAxisIndex ) = 6.251e6;
        expectedKeplerianElements( eccentricityIndex ) = 0.0;
        expectedKeplerianElements( inclinationIndex ) = 0.0;
        expectedKeplerianElements( argumentOfPeriapsisIndex ) = 0.0;
        expectedKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.0;
        expectedKeplerianElements( trueAnomalyIndex ) = 26.78 / 180.0 * PI;

        // Declare and compute converted Keplerian elements.
        Eigen::Matrix< double, 6, 1 > computedKeplerianElements;
        computedKeplerianElements = orbital_element_conversions::
                convertCartesianToKeplerianElements( cartesianElements,
                                                     venusGravitationalParameter );

        // Check if computed semi-major axis matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedKeplerianElements( semiMajorAxisIndex ),
                                    computedKeplerianElements( semiMajorAxisIndex ), 1.0e-15 );

        // Check if computed eccentricity matches the expected value.
        BOOST_CHECK_SMALL( computedKeplerianElements( eccentricityIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        // Check if computed inclination matches the expected value.
        BOOST_CHECK_SMALL( computedKeplerianElements( inclinationIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        // Check if computed argument of periapsis matches the expected value.
        BOOST_CHECK_SMALL( std::fmod( computedKeplerianElements( argumentOfPeriapsisIndex ),
                                      2.0 * PI ),
                           std::numeric_limits< double >::epsilon( ) );

        // Check if computed right ascension of ascending node matches the expected value.
        BOOST_CHECK_SMALL( computedKeplerianElements( longitudeOfAscendingNodeIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        // Check if computed true anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedKeplerianElements( trueAnomalyIndex ),
                                    computedKeplerianElements( trueAnomalyIndex ), 1.0e-15 );
    }

    // Case 3: Hyperbolic orbit around the Sun.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Using declarations.
        using namespace orbital_element_conversions;

        // Sun gravitational parameter.
        const double sunGravitationalParameter = 1.32712440018e20;

        // Declare and set Cartesian elements.
        Eigen::Matrix< double, 6, 1 > cartesianElements;
        cartesianElements( xCartesianPositionIndex ) = 7.035635643405699e11;
        cartesianElements( yCartesianPositionIndex ) = -2.351218213055550e11;
        cartesianElements( zCartesianPositionIndex ) = 0.037960971564309e11;
        cartesianElements( xCartesianVelocityIndex ) = -1.731375459746510e4;
        cartesianElements( yCartesianVelocityIndex ) = -1.535713656317794e4;
        cartesianElements( zCartesianVelocityIndex ) = 0.423498718768347e4;

        // Set expected Keplerian elements.
        Eigen::Matrix< double, 6, 1 > expectedKeplerianElements;
        expectedKeplerianElements( semiMajorAxisIndex ) = -6.78e11;
        expectedKeplerianElements( eccentricityIndex ) = 1.89;
        expectedKeplerianElements( inclinationIndex ) = 167.91 / 180 * PI;
        expectedKeplerianElements( argumentOfPeriapsisIndex ) = 45.78 / 180.0 * PI;
        expectedKeplerianElements( longitudeOfAscendingNodeIndex ) = 342.89 / 180.0 * PI;
        expectedKeplerianElements( trueAnomalyIndex ) = 315.62 / 180.0 * PI;

        // Compute Keplerian elements.
        Eigen::Matrix< double, 6, 1 > computedKeplerianElements;
        computedKeplerianElements = orbital_element_conversions::
                convertCartesianToKeplerianElements( cartesianElements,
                                                     sunGravitationalParameter );

        // Check if computed Keplerian elements match the expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedKeplerianElements,
                                           computedKeplerianElements, 1.0e-15 );
    }

    // Case 4: Low-eccentricity, low-inclination orbit around Uranus.
    // The benchmark data is obtained by running Keplerian Toolbox (ESA, 2012). It is important to
    // note that Keplerian Toolbox uses a different order of elements than Tudat (argument of
    // periapsis and longitude of ascending node positions in vector are switched). Not all values
    // are tested to the same precision, as there are slight differences in the computations using
    // Keplerian toolbox and Tudat. This is likely down to the use of functions in Eigen, such as
    // .cross, which are written out in full in Keplerian toolbox. If this turns out to be an
    // issue, the use of Eigen functions will need to be evaluated.
    {
        // Using declarations.
        using namespace orbital_element_conversions;

        // Uranus gravitational parameter [m^3 s^-2].
        const double uranusGravitationalParameter = 5.793943348799999e15;

        // Declare and set Cartesian elements.
        Eigen::Matrix< double, 6, 1 > cartesianElements;
        cartesianElements( xCartesianPositionIndex ) = -33760437.1526459;
        cartesianElements( yCartesianPositionIndex ) = -91719029.3283878;
        cartesianElements( zCartesianPositionIndex ) = -757.744826269064;
        cartesianElements( xCartesianVelocityIndex ) = 7225.59158151846;
        cartesianElements( yCartesianVelocityIndex ) = -2659.60535464239;
        cartesianElements( zCartesianVelocityIndex ) = -0.0486289999748532;

        // Set expected Keplerian elements.
        Eigen::Matrix< double, 6, 1 > expectedKeplerianElements;
        expectedKeplerianElements( semiMajorAxisIndex ) = 97736000.0;
        expectedKeplerianElements( eccentricityIndex ) = 1.0e-5;
        expectedKeplerianElements( inclinationIndex ) = 1.0e-5;
        expectedKeplerianElements( argumentOfPeriapsisIndex ) = 4.39704289113435;
        expectedKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.330903954202613;
        expectedKeplerianElements( trueAnomalyIndex ) = 5.914936209560839;

        // Compute Keplerian elements.
        Eigen::Matrix< double, 6, 1 > computedKeplerianElements;
        computedKeplerianElements = orbital_element_conversions::
                convertCartesianToKeplerianElements( cartesianElements,
                                                     uranusGravitationalParameter );

        // Check if computed Keplerian elements match the expected values.
        BOOST_CHECK_CLOSE_FRACTION( expectedKeplerianElements( semiMajorAxisIndex ),
                                    computedKeplerianElements( semiMajorAxisIndex ),
                                    1.0e-14 );

        BOOST_CHECK_SMALL( std::fabs( expectedKeplerianElements( eccentricityIndex )
                                      - computedKeplerianElements( eccentricityIndex ) ),
                           1.0e-15 );

        BOOST_CHECK_SMALL( std::fabs( expectedKeplerianElements( inclinationIndex )
                                      - computedKeplerianElements( inclinationIndex ) ),
                           1.0e-10 );

        BOOST_CHECK_CLOSE_FRACTION( expectedKeplerianElements( argumentOfPeriapsisIndex ),
                                    computedKeplerianElements( argumentOfPeriapsisIndex ),
                                    1.0e-10 );

        BOOST_CHECK_CLOSE_FRACTION( expectedKeplerianElements( longitudeOfAscendingNodeIndex ),
                                    computedKeplerianElements( longitudeOfAscendingNodeIndex ),
                                    1.0e-14 );

        BOOST_CHECK_CLOSE_FRACTION( expectedKeplerianElements( trueAnomalyIndex ),
                                    computedKeplerianElements( trueAnomalyIndex ),
                                    1.0e-6 );
    }

}

//! Test back and forth Kepler <-> Cartesian conversion for parabolic orbit
template< typename ScalarType >
void convertParabolicOrbitBackAndForth(
        const ScalarType tolerance )
{
    // Using declarations.
    using namespace orbital_element_conversions;

    // Set Earth gravitational parameter [m^3/s^2].
    const ScalarType earthGravitationalParameter = 3.986005e14;

    // Set Keplerian elements [m,-,rad,rad,rad,rad].
    Eigen::Matrix< ScalarType, 6, 1 > keplerianElements;
    keplerianElements( semiLatusRectumIndex ) = static_cast< ScalarType >(
                2.0 * 6678140.0 );
    keplerianElements( eccentricityIndex ) = getFloatingInteger< ScalarType >( 1 );
    keplerianElements( inclinationIndex ) = static_cast< ScalarType >(
                45.0 / 180.0 * PI );
    keplerianElements( argumentOfPeriapsisIndex ) = getFloatingInteger< ScalarType >( 0 );
    keplerianElements( longitudeOfAscendingNodeIndex ) = getFloatingInteger< ScalarType >( 0 );
    keplerianElements( trueAnomalyIndex ) = getFloatingInteger< ScalarType >( 0 );

    // Compute Cartesian elements.
    Eigen::Matrix< ScalarType, 6, 1 > computedCartesianElements =
            convertKeplerianToCartesianElements< ScalarType >(
                keplerianElements, earthGravitationalParameter );

    // Recompute Keplerian elements.
    Eigen::Matrix< ScalarType, 6, 1 > recomputedKeplerianElements =
            convertCartesianToKeplerianElements< ScalarType >(
                computedCartesianElements, earthGravitationalParameter );

    // Check if recomputed Keplerian elements match the expected values.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( keplerianElements,
                                       recomputedKeplerianElements, tolerance );

    std::cout<<"Keplerian "<<std::setprecision( 16 )<<keplerianElements.transpose( )<<std::endl;
    std::cout<<"Keplerian "<<recomputedKeplerianElements.transpose( )<<std::endl;
    std::cout<<"Keplerian "<<( keplerianElements - recomputedKeplerianElements ).transpose( )<<std::endl;
    recomputedKeplerianElements( eccentricityIndex ) = getFloatingInteger< ScalarType >( 1 );

    // Convert recomputed Keplerian elements to Cartesian elements.
    Eigen::Matrix< ScalarType, 6, 1 > recomputedCartesianElements =
            convertKeplerianToCartesianElements(
                recomputedKeplerianElements, earthGravitationalParameter );

    // Check that computed Cartesian elements match.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                computedCartesianElements, recomputedCartesianElements, ( 10.0 * tolerance ) );

    std::cout<<"Cartesian "<<std::setprecision( 16 )<<computedCartesianElements.transpose( )<<std::endl;
    std::cout<<"Cartesian "<<recomputedCartesianElements.transpose( )<<std::endl;
}

//! Test back and forth Kepler <-> Cartesian conversion for circular equatorial orbit
template< typename ScalarType >
void convertCircularEquatorialOrbitBackAndForth(
        const ScalarType tolerance )
{
    // Using declarations.
    using namespace orbital_element_conversions;

    // Earth gravitational parameter [m^3 s^-2].
    const ScalarType earthGravitationalParameter = 3.9859383624e14;

    // Set Keplerian elements.
    Eigen::Matrix< ScalarType, 6, 1 > keplerianElements;
    keplerianElements( semiMajorAxisIndex ) =
            static_cast< ScalarType >( 8000.0 );
    keplerianElements( eccentricityIndex ) = getFloatingInteger< ScalarType >( 0 );
    keplerianElements( inclinationIndex ) = getFloatingInteger< ScalarType >( 0 );
    keplerianElements( argumentOfPeriapsisIndex ) =
            static_cast< ScalarType >( 243.0 / 180.0 * PI );
    keplerianElements( longitudeOfAscendingNodeIndex ) =
            static_cast< ScalarType >( -79.6 / 180.0 * PI );
    keplerianElements( trueAnomalyIndex ) =
            static_cast< ScalarType >( 126.45 / 180.0 * PI );

    // Convert Keplerian elements to Cartesian elements.
    Eigen::Matrix< ScalarType, 6, 1 > computedCartesianElements;
    computedCartesianElements = convertKeplerianToCartesianElements< ScalarType >(
                keplerianElements, earthGravitationalParameter );

    // Convert Cartesian elements back to Keplerian elements.
    Eigen::Matrix< ScalarType, 6, 1 > recomputedKeplerianElements;
    recomputedKeplerianElements = convertCartesianToKeplerianElements< ScalarType >(
                computedCartesianElements, earthGravitationalParameter );

    // Check that recomputed Keplerian elements match the input values. (In this limit case,
    // the argument of periapsis and longitude of ascending node are set to zero, and the
    // true anomaly "absorbs" everything).
    BOOST_CHECK_CLOSE_FRACTION( keplerianElements( semiMajorAxisIndex ),
                                recomputedKeplerianElements( semiMajorAxisIndex ),
                                tolerance );

    BOOST_CHECK_SMALL( recomputedKeplerianElements( eccentricityIndex ), tolerance );

    BOOST_CHECK_SMALL( recomputedKeplerianElements( inclinationIndex ),
                       std::numeric_limits< ScalarType >::epsilon( ) );

    BOOST_CHECK_SMALL( recomputedKeplerianElements( argumentOfPeriapsisIndex ),
                       std::numeric_limits< ScalarType >::epsilon( ) );

    BOOST_CHECK_SMALL( recomputedKeplerianElements( longitudeOfAscendingNodeIndex ),
                       std::numeric_limits< ScalarType >::epsilon( ) );

    BOOST_CHECK_CLOSE_FRACTION( std::fmod(
                                    keplerianElements( argumentOfPeriapsisIndex )
                                    + keplerianElements( longitudeOfAscendingNodeIndex )
                                    + keplerianElements( trueAnomalyIndex ), 2.0 * PI ),
                                recomputedKeplerianElements( trueAnomalyIndex ),
                                std::numeric_limits< ScalarType >::epsilon( ) );

    // Convert recomputed Keplerian elements to Cartesian elements.
    Eigen::Matrix< ScalarType, 6, 1 > recomputedCartesianElements =
            convertKeplerianToCartesianElements(
                recomputedKeplerianElements, earthGravitationalParameter );

    // Check that computed Cartesian elements match.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                computedCartesianElements, recomputedCartesianElements, ( 10.0 * tolerance ) );
}

//! Test back and forth Kepler <-> Cartesian conversion for circular non-equatorial orbit
template< typename ScalarType >
void convertCircularNonEquatorialOrbitBackAndForth(
        const ScalarType tolerance )
{
    // Using declarations.
    using namespace orbital_element_conversions;

    // Earth gravitational parameter [m^3 s^-2].
    const ScalarType earthGravitationalParameter = 3.9859383624e14;

    // Set Keplerian elements.
    Eigen::Matrix< ScalarType, 6, 1 > keplerianElements;
    keplerianElements( semiMajorAxisIndex ) = static_cast< ScalarType >( 8000.0 );
    keplerianElements( eccentricityIndex ) =  getFloatingInteger< ScalarType >( 0 );
    keplerianElements( inclinationIndex ) = static_cast< ScalarType >( 176.11 / 180.0 * PI );
    keplerianElements( argumentOfPeriapsisIndex ) =
            static_cast< ScalarType >(  243.0 / 180.0 * PI );
    keplerianElements( longitudeOfAscendingNodeIndex ) =
            static_cast< ScalarType >( -79.6 / 180.0 * PI );
    keplerianElements( trueAnomalyIndex ) =
            static_cast< ScalarType >( 126.45 / 180.0 * PI );

    // Convert Keplerian elements to Cartesian elements.
    Eigen::Matrix< ScalarType, 6, 1 > computedCartesianElements =
            convertKeplerianToCartesianElements< ScalarType >(
                keplerianElements, earthGravitationalParameter );

    // Convert Cartesian elements back to Keplerian elements.
    Eigen::Matrix< ScalarType, 6, 1 > recomputedKeplerianElements =
            convertCartesianToKeplerianElements< ScalarType >(
                computedCartesianElements, earthGravitationalParameter );

    // Check that recomputed Keplerian elements match the input values. (In this limit case,
    // the argument of periapsis is set to zero, and the true anomaly "absorbs" this).
    BOOST_CHECK_CLOSE_FRACTION(
                keplerianElements( semiMajorAxisIndex ),
                recomputedKeplerianElements( semiMajorAxisIndex ), tolerance );

    BOOST_CHECK_SMALL( recomputedKeplerianElements( eccentricityIndex ), tolerance );

    BOOST_CHECK_CLOSE_FRACTION(
                keplerianElements( inclinationIndex ),
                recomputedKeplerianElements( inclinationIndex ), tolerance );

    BOOST_CHECK_SMALL( recomputedKeplerianElements( argumentOfPeriapsisIndex ),
                       std::numeric_limits< ScalarType >::epsilon( ) );

    BOOST_CHECK_CLOSE_FRACTION(
                keplerianElements( longitudeOfAscendingNodeIndex ) +
                getFloatingInteger< ScalarType >( 2 ) * getPi< ScalarType >( ),
                recomputedKeplerianElements( longitudeOfAscendingNodeIndex ), tolerance );

    BOOST_CHECK_CLOSE_FRACTION(
                std::fmod( keplerianElements( argumentOfPeriapsisIndex )
                           + keplerianElements( trueAnomalyIndex ),
                           getFloatingInteger< ScalarType >( 2 ) * getPi< ScalarType >( ) ),
                recomputedKeplerianElements( trueAnomalyIndex ), 10.0 * tolerance );

    // Convert recomputed Keplerian elements to Cartesian elements.
    Eigen::Matrix< ScalarType, 6, 1 > recomputedCartesianElements =
            convertKeplerianToCartesianElements< ScalarType >(
                recomputedKeplerianElements, earthGravitationalParameter );

    // Check that computed Cartesian elements match.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                computedCartesianElements, recomputedCartesianElements, ( 100.0 * tolerance ) );
}

//! Test back and forth Kepler <-> Cartesian conversion for non-circular equatorial orbit
template< typename ScalarType >
void convertNonCircularEquatorialOrbitBackAndForth(
        const ScalarType tolerance )
{
    // Using declarations.
    using namespace orbital_element_conversions;

    // Earth gravitational parameter [m^3 s^-2].
    const ScalarType earthGravitationalParameter = 3.9859383624e14;

    // Set Keplerian elements.
    Eigen::Matrix< ScalarType, 6, 1 > keplerianElements;
    keplerianElements( semiMajorAxisIndex ) = static_cast< ScalarType >( 8000.0 );
    keplerianElements( eccentricityIndex ) = static_cast< ScalarType >( 0.2 );
    keplerianElements( inclinationIndex ) = getFloatingInteger< ScalarType >( 0 );
    keplerianElements( argumentOfPeriapsisIndex ) = static_cast< ScalarType >( 243.0 / 180.0 * PI );
    keplerianElements( longitudeOfAscendingNodeIndex ) = static_cast< ScalarType >(
                -79.6 / 180.0 * PI );
    keplerianElements( trueAnomalyIndex ) = static_cast< ScalarType >( 126.45 / 180.0 * PI );

    // Convert Keplerian elements to Cartesian elements.
    Eigen::Matrix< ScalarType, 6, 1 > computedCartesianElements =
            convertKeplerianToCartesianElements( keplerianElements, earthGravitationalParameter );

    // Convert Cartesian elements back to Keplerian elements.
    Eigen::Matrix< ScalarType, 6, 1 > recomputedKeplerianElements =
            convertCartesianToKeplerianElements(
                computedCartesianElements, earthGravitationalParameter );

    // Check that recomputed Keplerian elements match the input values. (In this limit case,
    // the longitude of the ascending node is set to zero, and the argument of periapsis and
    // true anomaly "absorb" this).
    BOOST_CHECK_CLOSE_FRACTION( keplerianElements( semiMajorAxisIndex ),
                                recomputedKeplerianElements( semiMajorAxisIndex ),
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( keplerianElements( eccentricityIndex ),
                                recomputedKeplerianElements( eccentricityIndex ),
                                tolerance );

    BOOST_CHECK_CLOSE_FRACTION( keplerianElements( inclinationIndex ),
                                recomputedKeplerianElements( inclinationIndex ),
                                std::numeric_limits< ScalarType >::epsilon( ) );

    BOOST_CHECK_SMALL( recomputedKeplerianElements( longitudeOfAscendingNodeIndex ),
                       std::numeric_limits< ScalarType >::epsilon( ) );

    BOOST_CHECK_CLOSE_FRACTION( keplerianElements( longitudeOfAscendingNodeIndex )
                                + keplerianElements( argumentOfPeriapsisIndex )
                                + keplerianElements( trueAnomalyIndex ),
                                recomputedKeplerianElements( argumentOfPeriapsisIndex )
                                + recomputedKeplerianElements( trueAnomalyIndex ),
                                std::numeric_limits< ScalarType >::epsilon( ) );

    // Convert recomputed Keplerian elements to Cartesian elements.
    Eigen::Matrix< ScalarType, 6, 1 > recomputedCartesianElements;
    recomputedCartesianElements = orbital_element_conversions::
            convertKeplerianToCartesianElements(
                recomputedKeplerianElements, earthGravitationalParameter );

    // Check that computed Cartesian elements match.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                computedCartesianElements, recomputedCartesianElements, ( 10.0 * tolerance ) );
}

//! Test back and forth Kepler <-> Cartesian conversion for non-circular non-equatorial orbit
template< typename ScalarType >
void convertNonCircularNonEquatorialOrbitBackAndForth(
        const ScalarType tolerance )
{
    // Using declarations.
    using namespace orbital_element_conversions;

    // Earth gravitational parameter [m^3 s^-2].
    const ScalarType earthGravitationalParameter = 3.9859383624e14;

    // Set Keplerian elements.
    Eigen::Matrix< ScalarType, 6, 1 > keplerianElements;
    keplerianElements( semiMajorAxisIndex ) = static_cast< ScalarType >( 8000.0 );
    keplerianElements( eccentricityIndex ) = static_cast< ScalarType >( 0.2 );
    keplerianElements( inclinationIndex ) = getFloatingInteger< ScalarType >(
                176.11 / 180.0 * PI );
    keplerianElements( argumentOfPeriapsisIndex ) = static_cast< ScalarType >( 243.0 / 180.0 * PI );
    keplerianElements( longitudeOfAscendingNodeIndex ) = static_cast< ScalarType >(
                -79.6 / 180.0 * PI );
    keplerianElements( trueAnomalyIndex ) = static_cast< ScalarType >( 126.45 / 180.0 * PI );

    // Convert Keplerian elements to Cartesian elements.
    Eigen::Matrix< ScalarType, 6, 1 > computedCartesianElements =
            convertKeplerianToCartesianElements( keplerianElements, earthGravitationalParameter );

    // Convert Cartesian elements back to Keplerian elements.
    Eigen::Matrix< ScalarType, 6, 1 > recomputedKeplerianElements =
            convertCartesianToKeplerianElements(
                computedCartesianElements, earthGravitationalParameter );

    recomputedKeplerianElements( longitudeOfAscendingNodeIndex ) =
            recomputedKeplerianElements( longitudeOfAscendingNodeIndex ) -
            getFloatingInteger< ScalarType >( 2 ) * getPi< ScalarType >( );

    // Check that recomputed Keplerian elements match the input values. (In this limit case,
    // the longitude of the ascending node is set to zero, and the argument of periapsis and
    // true anomaly "absorb" this).
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                keplerianElements, recomputedKeplerianElements, ( tolerance ) );

    // Convert recomputed Keplerian elements to Cartesian elements.
    Eigen::Matrix< ScalarType, 6, 1 > recomputedCartesianElements;
    recomputedCartesianElements = orbital_element_conversions::
            convertKeplerianToCartesianElements(
                recomputedKeplerianElements, earthGravitationalParameter );

    // Check that computed Cartesian elements match.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                computedCartesianElements, recomputedCartesianElements, ( 20.0 * tolerance ) );
}

//! Test if conversion from Cartesian elements to Keplerian elements is working correctly, using
//! back and forth conversion, for both double and long double precision.
BOOST_AUTO_TEST_CASE( testCartesianToKeplerianElementConversionBackAndForth )
{
    double ratioOfPrecision = std::numeric_limits< long double >::epsilon( ) /
            std::numeric_limits< double >::epsilon( );

    // Case 1: Parabolic orbit around the Sun.
    // This test is based on converting Keplerian elements to Cartesian element and then
    // recomputing the input Keplerian element values. Ideally, an independent check will replace
    // this test in future.
    {
        convertParabolicOrbitBackAndForth< double >( 1.0E-15 );
        convertParabolicOrbitBackAndForth< long double >( 1.0E-15L * ratioOfPrecision );
    }
    // Case 2: Converting to and fro between Keplerian and Cartesian elements for a
    // zero-eccentricity, zero-inclination orbit (circular, equatorial) around the Earth. This
    // test ensures internal consistency within Tudat between the Cartesian <-> Keplerian element
    // converters.
    {
        convertCircularEquatorialOrbitBackAndForth< double >( 1.0E-15 );
        convertCircularEquatorialOrbitBackAndForth< long double >( 1.0E-15L * ratioOfPrecision );
    }

    // Case 3: Converting to and fro between Keplerian and Cartesian elements for a
    // zero-eccentricity, non-zero-inclination orbit (circular, inclined) around the Earth. This
    // test ensures internal consistency within Tudat between the Cartesian <-> Keplerian element
    // converters.
    {
        convertCircularNonEquatorialOrbitBackAndForth< double >( 1.0E-15 );
        convertCircularNonEquatorialOrbitBackAndForth< long double >( 1.0E-15L * ratioOfPrecision );
    }

    // Case 4: Converting to and fro between Keplerian and Cartesian elements for a
    // non-zero-eccentricity, zero-inclination orbit (non-circular, equatorial) around the Earth.
    // This test ensures internal consistency within Tudat between the Cartesian <-> Keplerian
    // element converters.
    {
        convertNonCircularEquatorialOrbitBackAndForth< double >( 1.0E-15 );
        convertNonCircularEquatorialOrbitBackAndForth< long double >( 1.0E-15L * ratioOfPrecision );

    }

    // Case 5: Converting to and fro between Keplerian and Cartesian elements for a
    // non-zero-eccentricity, non-zero-inclination orbit (non-circular, non-equatorial).
    // This test ensures internal consistency within Tudat between the Cartesian <-> Keplerian
    // element converters.
    {
        convertNonCircularNonEquatorialOrbitBackAndForth< double >( 5.0E-15 );
        convertNonCircularNonEquatorialOrbitBackAndForth< long double >( 5.0E-15L * ratioOfPrecision );

    }
}

//! Test if conversion from true anomaly to eccentric anomaly is working correctly.
BOOST_AUTO_TEST_CASE( testTrueAnomalyToEccentricAnomalyConversion )
{
    // Case 1: General elliptical orbit.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Set eccentricity.
        const double eccentricity = 0.146;

        // Set true anomaly.
        const double trueAnomaly = 82.16 / 180.0 * PI;

        // Set expected elliptic eccentric anomaly.
        const double expectedEllipticEccentricAnomaly = 1.290237398010989;

        // Compute elliptic eccentric anomaly.
        const double computedEllipticEccentricAnomaly
                = orbital_element_conversions::
                convertTrueAnomalyToEllipticalEccentricAnomaly( trueAnomaly, eccentricity );

        // Check if computed elliptic eccentric anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedEllipticEccentricAnomaly,
                                    computedEllipticEccentricAnomaly,
                                    2.0 * std::numeric_limits< double >::epsilon( ) );
    }

    // Case 2: Circular orbit.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Set eccentricity.
        const double eccentricity = 0.0;

        // Set true anomaly.
        const double trueAnomaly = 160.43 / 180.0 * PI;

        // Set expected elliptic eccentric anomaly.
        const double expectedEllipticEccentricAnomaly = 2.800031718974503;

        // Compute elliptic eccentric anomaly.
        const double computedEllipticEccentricAnomaly
                = orbital_element_conversions::
                convertTrueAnomalyToEllipticalEccentricAnomaly( trueAnomaly, eccentricity );

        // Check if computed elliptic eccentric anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedEllipticEccentricAnomaly,
                                    computedEllipticEccentricAnomaly,
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Case 3: Circular orbit at periapsis.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Set eccentricity.
        const double eccentricity = 0.0;

        // Set true anomaly.
        const double trueAnomaly = 0.0;

        // Set expected elliptic eccentric anomaly.
        const double expectedEllipticEccentricAnomaly = 0.0;

        // Compute elliptic eccentric anomaly.
        const double computedEllipticEccentricAnomaly
                = orbital_element_conversions::
                convertTrueAnomalyToEllipticalEccentricAnomaly( trueAnomaly, eccentricity );

        // Check if computed elliptic eccentric anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedEllipticEccentricAnomaly,
                                    computedEllipticEccentricAnomaly,
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Case 4: General hyperbolic orbit.
    // The benchmark data is obtained from (Fortescue, 2003).
    {
        // Set eccentricity.
        const double eccentricity = 3.0;

        // Set true anomaly.
        const double trueAnomaly = 0.5291;

        // Set expected hyperbolic eccentric anomaly.
        const double expectedHyperbolicEccentricAnomaly = 0.3879;

        // Compute hyperbolic eccentric anomaly.
        const double convertedHyperbolicEccentricAnomaly
                = orbital_element_conversions
                ::convertTrueAnomalyToHyperbolicEccentricAnomaly( trueAnomaly, eccentricity );

        // Check if computed hyperbolic eccentric anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedHyperbolicEccentricAnomaly,
                                    convertedHyperbolicEccentricAnomaly, 1.0e-5 );
    }

    // Case 5: General elliptical orbit (test for wrapper function).
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Set eccentricity.
        const double eccentricity = 0.146;

        // Set true anomaly.
        const double trueAnomaly = 82.16 / 180.0 * PI;

        // Set expected elliptic eccentric anomaly.
        const double expectedEllipticEccentricAnomaly = 1.290237398010989;

        // Compute elliptic eccentric anomaly.
        const double computedEllipticEccentricAnomaly
                = orbital_element_conversions::
                convertTrueAnomalyToEccentricAnomaly( trueAnomaly, eccentricity );

        // Check if computed elliptic eccentric anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedEllipticEccentricAnomaly,
                                    computedEllipticEccentricAnomaly,
                                    2.0 * std::numeric_limits< double >::epsilon( ) );
    }

    // Case 6: General hyperbolic orbit (test for wrapper function).
    // The benchmark data is obtained from (Fortescue, 2003).
    {
        // Set eccentricity.
        const double eccentricity = 3.0;

        // Set true anomaly.
        const double trueAnomaly = 0.5291;

        // Set expected hyperbolic eccentric anomaly.
        const double expectedHyperbolicEccentricAnomaly = 0.3879;

        // Compute hyperbolic eccentric anomaly.
        const double convertedHyperbolicEccentricAnomaly
                = orbital_element_conversions
                ::convertTrueAnomalyToEccentricAnomaly( trueAnomaly, eccentricity );

        // Check if computed hyperbolic eccentric anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedHyperbolicEccentricAnomaly,
                                    convertedHyperbolicEccentricAnomaly, 1.0e-5 );
    }
}

//! Test if conversion from eccentric anomaly to true anomaly is working correctly.
BOOST_AUTO_TEST_CASE( testEccentricAnomalyToTrueAnomalyConversion )
{
    // Case 1: General elliptical orbit.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Set eccentricity.
        const double eccentricity = 0.639;

        // Set elliptic eccentric anomaly.
        const double ellipticEccentricAnomaly = 239.45 / 180.0 * PI;

        // Set expected true anomaly.
        const double expectedTrueAnomaly = 3.665218735816221;

        // Compute true anomaly, modulo 2*pi.
        const double convertedTrueAnomaly
                = orbital_element_conversions::
                convertEllipticalEccentricAnomalyToTrueAnomaly( ellipticEccentricAnomaly,
                                                                eccentricity ) + 2.0 * PI;

        // Check if computed true anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedTrueAnomaly, convertedTrueAnomaly,
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Case 2: Circular orbit.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Set eccentricity.
        const double eccentricity = 0.0;

        // Set elliptic eccentric anomaly.
        const double ellipticEccentricAnomaly = -99.54 / 180.0 * PI;

        // Set expected true anomaly.
        const double expectedTrueAnomaly = 4.545884569744431;

        // Compute true anomaly.
        const double convertedTrueAnomaly
                = orbital_element_conversions::
                convertEllipticalEccentricAnomalyToTrueAnomaly( ellipticEccentricAnomaly,
                                                                eccentricity ) + 2.0 * PI;

        // Check if computed true anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedTrueAnomaly, convertedTrueAnomaly,
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Case 3: Circular orbit at periapsis.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Set eccentricity.
        const double eccentricity = 0.0;

        // Set elliptic eccentric anomaly.
        const double ellipticEccentricAnomaly = 0.0;

        // Set expected true anomaly.
        const double expectedTrueAnomaly = 0.0;

        // Compute true anomaly.
        const double convertedTrueAnomaly
                = orbital_element_conversions::
                convertEllipticalEccentricAnomalyToTrueAnomaly( ellipticEccentricAnomaly,
                                                                eccentricity );

        // Check if computed true anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedTrueAnomaly, convertedTrueAnomaly,
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Case 4: General hyperbolic orbit.
    // The benchmark data is obtained from (Fortescue, 2003).
    {
        // Set eccentricity.
        const double eccentricity = 3.0;

        // Set hyperbolic eccentric anomaly.
        const double hyperbolicEccentricAnomaly = 0.3879;

        // Set expected true anomaly.
        const double expectedTrueAnomaly = 0.5291;

        // Compute true anomaly.
        const double convertedTrueAnomaly
                = orbital_element_conversions
                ::convertHyperbolicEccentricAnomalyToTrueAnomaly(
                    hyperbolicEccentricAnomaly, eccentricity );

        // Check if computed true anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedTrueAnomaly, convertedTrueAnomaly, 1.0e-5 );
    }

    // Case 5: General elliptical orbit (test for wrapper function).
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Set eccentricity.
        const double eccentricity = 0.639;

        // Set elliptic eccentric anomaly.
        const double ellipticEccentricAnomaly = 239.45 / 180.0 * PI;

        // Set expected true anomaly.
        const double expectedTrueAnomaly = 3.665218735816221;

        // Compute true anomaly, modulo 2*pi.
        const double convertedTrueAnomaly
                = orbital_element_conversions::
                convertEccentricAnomalyToTrueAnomaly( ellipticEccentricAnomaly,
                                                      eccentricity ) + 2.0 * PI;

        // Check if computed true anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedTrueAnomaly, convertedTrueAnomaly,
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Case 6: General hyperbolic orbit (test for wrapper function).
    // The benchmark data is obtained from (Fortescue, 2003).
    {
        // Set eccentricity.
        const double eccentricity = 3.0;

        // Set hyperbolic eccentric anomaly.
        const double hyperbolicEccentricAnomaly = 0.3879;

        // Set expected true anomaly.
        const double expectedTrueAnomaly = 0.5291;

        // Compute true anomaly.
        const double convertedTrueAnomaly
                = orbital_element_conversions
                ::convertEccentricAnomalyToTrueAnomaly(
                    hyperbolicEccentricAnomaly, eccentricity );

        // Check if computed true anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedTrueAnomaly, convertedTrueAnomaly, 1.0e-5 );
    }
}

//! Test if conversion from eccentric anomaly to mean anomaly is working correctly.
BOOST_AUTO_TEST_CASE( testEccentricAnomalyToMeanAnomalyConversion )
{
    // Case 1: General elliptical orbit.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Set eccentricity.
        const double eccentricity = 0.541;

        // Set elliptical eccentric anomaly.
        const double ellipticalEccentricAnomaly = 176.09 / 180.0 * PI;

        // Set expected mean anomaly.
        const double expectedMeanAnomaly = 3.036459804491048;

        // Compute mean anomaly.
        const double computedMeanAnomaly = orbital_element_conversions
                ::convertEllipticalEccentricAnomalyToMeanAnomaly(
                    ellipticalEccentricAnomaly, eccentricity );

        // Check if computed mean anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedMeanAnomaly, computedMeanAnomaly,
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Case 2: Circular orbit.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Set eccentricity.
        const double eccentricity = 0.0;

        // Set elliptical eccentric anomaly.
        const double ellipticalEccentricAnomaly = 320.12 / 180.0 * PI;

        // Set expected mean anomaly.
        const double expectedMeanAnomaly = 5.587148001484247;

        // Compute mean anomaly.
        const double computedMeanAnomaly = orbital_element_conversions
                ::convertEllipticalEccentricAnomalyToMeanAnomaly(
                    ellipticalEccentricAnomaly, eccentricity );

        // Check if computed mean anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedMeanAnomaly, computedMeanAnomaly,
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Case 3: Circular orbit at periapsis.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Set eccentricity.
        const double eccentricity = 0.0;

        // Set elliptical eccentric anomaly.
        const double ellipticalEccentricAnomaly = 0.0;

        // Set expected mean anomaly.
        const double expectedMeanAnomaly = 0.0;

        // Compute mean anomaly.
        const double computedMeanAnomaly = orbital_element_conversions
                ::convertEllipticalEccentricAnomalyToMeanAnomaly(
                    ellipticalEccentricAnomaly, eccentricity );

        // Check if computed mean anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedMeanAnomaly, computedMeanAnomaly,
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Case 4: General hyperbolic orbit.
    // The benchmark data is obtained from (Vallado, 2004).
    {
        // Set eccentricity.
        const double eccentricity = 2.4;

        // Set hyperbolic eccentric anomaly.
        const double hyperbolicEccentricAnomaly = 1.6013761449;

        // Set expected mean anomaly.
        const double expectedMeanAnomaly = 235.4 / 180.0 * PI;

        // Compute mean anomaly.
        const double computedMeanAnomaly = orbital_element_conversions
                ::convertHyperbolicEccentricAnomalyToMeanAnomaly( hyperbolicEccentricAnomaly,
                                                                  eccentricity );

        // Check if computed mean anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedMeanAnomaly, computedMeanAnomaly, 1.0e-8 );
    }

    // Case 5: General elliptical orbit (test for wrapper function).
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    {
        // Set eccentricity.
        const double eccentricity = 0.541;

        // Set elliptical eccentric anomaly.
        const double ellipticalEccentricAnomaly = 176.09 / 180.0 * PI;

        // Set expected mean anomaly.
        const double expectedMeanAnomaly = 3.036459804491048;

        // Compute mean anomaly.
        const double computedMeanAnomaly = orbital_element_conversions
                ::convertEccentricAnomalyToMeanAnomaly(
                    ellipticalEccentricAnomaly, eccentricity );

        // Check if computed mean anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedMeanAnomaly, computedMeanAnomaly,
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Case 6: General hyperbolic orbit.
    // The benchmark data is obtained from (Vallado, 2004).
    {
        // Set eccentricity.
        const double eccentricity = 2.4;

        // Set hyperbolic eccentric anomaly.
        const double hyperbolicEccentricAnomaly = 1.6013761449;

        // Set expected mean anomaly.
        const double expectedMeanAnomaly = 235.4 / 180.0 * PI;

        // Compute mean anomaly.
        const double computedMeanAnomaly = orbital_element_conversions
                ::convertEccentricAnomalyToMeanAnomaly( hyperbolicEccentricAnomaly, eccentricity );

        // Check if computed mean anomaly matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedMeanAnomaly, computedMeanAnomaly, 1.0e-8 );
    }
}

//! Test if conversion from elapsed time to mean anomaly change is working correctly.
BOOST_AUTO_TEST_CASE( testElapsedTimeToMeanAnomalyConversion )
{
    // Case 1: Earth-orbiting satellite.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    // The data is produced by running the kepprop2b_test() function through
    // kepprop2b('ValidationTest'). To see the mean anomaly values computed, in the MATLAB file,
    // remove the semi-colon at the end of line 68 in kepprop2b.m.
    {
        // Set elapsed time [s].
        const double elapsedTime = 8640.0;

        // Set Earth gravitational parameter.
        const double earthGravitationalParameter = 398600.4415;

        // Set semi-major axis [m].
        const double semiMajorAxis = 42165.3431351313;

        // Set expected elliptical mean anomaly change [rad].
        const double expectedEllipticalMeanAnomalyChange = 2.580579656848906 - 1.950567148859647;

        // Compute elliptical mean anomaly change [rad].
        const double computedEllipticalMeanAnomalyChange
                = orbital_element_conversions
                ::convertElapsedTimeToEllipticalMeanAnomalyChange(
                    elapsedTime, earthGravitationalParameter, semiMajorAxis );

        // Check if computed elliptical mean anomaly change matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedEllipticalMeanAnomalyChange,
                                    computedEllipticalMeanAnomalyChange, 1.0e-14 );
    }

    // Case 2: Earth-orbiting satellite with no elapsed time.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    // The data is produced by running the kepprop2b_test() function through
    // kepprop2b('ValidationTest'). To see the mean anomaly values computed, in the MATLAB file,
    // remove the semi-colon at the end of line 68 in kepprop2b.m.
    {
        // Set elapsed time [s].
        const double elapsedTime = 0.0;

        // Set Earth gravitational parameter.
        const double earthGravitationalParameter = 398600.4415;

        // Set semi-major axis [m].
        const double semiMajorAxis = 42165.3431351313;

        // Set expected elliptical mean anomaly change [rad].
        const double expectedEllipticalMeanAnomalyChange = 0.0;

        // Compute elliptical mean anomaly change [rad].
        const double computedEllipticalMeanAnomalyChange
                = orbital_element_conversions
                ::convertElapsedTimeToEllipticalMeanAnomalyChange(
                    elapsedTime, earthGravitationalParameter, semiMajorAxis );

        // Check if computed elliptical mean anomaly change matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedEllipticalMeanAnomalyChange,
                                    computedEllipticalMeanAnomalyChange,
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Case 3: Hyperbolic orbit around the Sun.
    // The benchmark data is obtained by simply repeating the equations in MATLAB. This test needs
    // to be replaced with real "external" benchmark data.
    {
        // Set elapsed time [s].
        const double elapsedTime = 1000.0;

        // Set Sun gravitational parameter [m^3/s^-2].
        const double sunGravitationalParameter = 3.9859383624e14;

        // Set semi-major axis [m].
        const double semiMajorAxis = -40000.0;

        // Set expected mean anomaly change [rad].
        const double expectedMeanAnomalyChange = 2.495601869539691e3;

        // Compute mean anomaly change [rad].
        const double computedMeanAnomalyChange
                = orbital_element_conversions
                ::convertElapsedTimeToHyperbolicMeanAnomalyChange(
                    elapsedTime, sunGravitationalParameter, semiMajorAxis );

        // Check if computed mean anomaly change matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedMeanAnomalyChange, computedMeanAnomalyChange,
                                    1.0e-15 );
    }

    // Case 4: Earth-orbiting satellite (test for wrapper function).
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    // The data is produced by running the kepprop2b_test() function through
    // kepprop2b('ValidationTest'). To see the mean anomaly values computed, in the MATLAB file,
    // remove the semi-colon at the end of line 68 in kepprop2b.m.
    {
        // Set elapsed time [s].
        const double elapsedTime = 8640.0;

        // Set Earth gravitational parameter.
        const double earthGravitationalParameter = 398600.4415;

        // Set semi-major axis [m].
        const double semiMajorAxis = 42165.3431351313;

        // Set expected elliptical mean anomaly change [rad].
        const double expectedEllipticalMeanAnomalyChange = 2.580579656848906 - 1.950567148859647;

        // Compute elliptical mean anomaly change [rad].
        const double computedEllipticalMeanAnomalyChange
                = orbital_element_conversions
                ::convertElapsedTimeToMeanAnomalyChange(
                    elapsedTime, earthGravitationalParameter, semiMajorAxis );

        // Check if computed elliptical mean anomaly change matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedEllipticalMeanAnomalyChange,
                                    computedEllipticalMeanAnomalyChange, 1.0e-14 );
    }

    // Case 5: Hyperbolic orbit around the Sun (test for wrapper function).
    // The benchmark data is obtained by simply repeating the equations in MATLAB. This test needs
    // to be replaced with real "external" benchmark data.
    {
        // Set elapsed time [s].
        const double elapsedTime = 1000.0;

        // Set Sun gravitational parameter [m^3/s^-2].
        const double sunGravitationalParameter = 3.9859383624e14;

        // Set semi-major axis [m].
        const double semiMajorAxis = -40000.0;

        // Set expected mean anomaly change [rad].
        const double expectedMeanAnomalyChange = 2.495601869539691e3;

        // Compute mean anomaly change [rad].
        const double computedMeanAnomalyChange
                = orbital_element_conversions
                ::convertElapsedTimeToMeanAnomalyChange(
                    elapsedTime, sunGravitationalParameter, semiMajorAxis );

        // Check if computed mean anomaly change matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedMeanAnomalyChange, computedMeanAnomalyChange,
                                    1.0e-15 );
    }
}

//! Test if conversion from mean anomaly change to elapsed time is working correctly.
BOOST_AUTO_TEST_CASE( testMeanAnomalyToElaspedTimeConversion )
{
    // Case 1: Earth-orbiting satellite.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    // The data is produced by running the kepprop2b_test() function through
    // kepprop2b('ValidationTest'). To see the mean anomaly values computed, in the MATLAB file,
    // remove the semi-colon at the end of line 68 in kepprop2b.m.
    {
        // Set elliptical mean anomaly change [rad].
        const double ellipticalMeanAnomalyChange = 3.210592164838165 - 1.950567148859647;

        // Set Earth gravitational parameter [m^3/s^2].
        const double earthGravitationalParameter = 398600.4415;

        // Set semi-major axis [m].
        const double semiMajorAxis = 42165.3431351313;

        // Set expected elapsed time [s].
        const double expectedElapsedTime = 17280.0;

        // Compute elapsed time [s].
        const double computedElapsedTime = orbital_element_conversions
                ::convertEllipticalMeanAnomalyChangeToElapsedTime(
                    ellipticalMeanAnomalyChange, earthGravitationalParameter, semiMajorAxis );

        // Check if computed elapsed time matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedElapsedTime, computedElapsedTime, 1.0e-15 );
    }

    // Case 2: Earth-orbiting satellite with no mean anomaly change.
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    // The data is produced by running the kepprop2b_test() function through
    // kepprop2b('ValidationTest'). To see the mean anomaly values computed, in the MATLAB file,
    // remove the semi-colon at the end of line 68 in kepprop2b.m.
    {
        // Set elliptical mean anomaly change [rad].
        const double ellipticalMeanAnomalyChange = 0.0;

        // Set Earth gravitational parameter [m^3/s^2].
        const double earthGravitationalParameter = 398600.4415;

        // Set semi-major axis [m].
        const double semiMajorAxis = 42165.3431351313;

        // Set expected elapsed time [s].
        const double expectedElapsedTime = 0.0;

        // Compute elapsed time [s].
        const double computedElapsedTime = orbital_element_conversions
                ::convertEllipticalMeanAnomalyChangeToElapsedTime(
                    ellipticalMeanAnomalyChange, earthGravitationalParameter, semiMajorAxis );

        // Check if computed elapsed time matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedElapsedTime, computedElapsedTime,
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Case 3: Hyperbolic orbit around the Sun.
    // The benchmark data is obtained by simply repeating the equations in MATLAB. This test needs
    // to be replaced with real "external" benchmark data.
    {
        // Set hyperbolic mean anomaly [rad].
        const double hyperbolicMeanAnomaly = 2.495601869539691e3;

        // Set Sun gravitational parameter [m^3/s^-2].
        const double sunGravitationalParameter = 3.9859383624e14;

        // Set semi-major axis [m].
        const double semiMajorAxis = -40000.0;

        // Set expected elapsed time [s].
        const double expectedElapsedTime = 1000.0;

        // Compute elapsed time [s].
        const double computedElapsedTime = orbital_element_conversions
                ::convertHyperbolicMeanAnomalyChangeToElapsedTime(
                    hyperbolicMeanAnomaly, sunGravitationalParameter, semiMajorAxis );

        // Check if computed elapsed time matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedElapsedTime, computedElapsedTime, 1.0e-15 );
    }

    // Case 4: Earth-orbiting satellite (test wrapper function).
    // The benchmark data is obtained by running ODTBX (NASA, 2012).
    // The data is produced by running the kepprop2b_test() function through
    // kepprop2b('ValidationTest'). To see the mean anomaly values computed, in the MATLAB file,
    // remove the semi-colon at the end of line 68 in kepprop2b.m.
    {
        // Set elliptical mean anomaly change [rad].
        const double ellipticalMeanAnomalyChange = 3.210592164838165 - 1.950567148859647;

        // Set Earth gravitational parameter [m^3/s^2].
        const double earthGravitationalParameter = 398600.4415;

        // Set semi-major axis [m].
        const double semiMajorAxis = 42165.3431351313;

        // Set expected elapsed time [s].
        const double expectedElapsedTime = 17280.0;

        // Compute elapsed time [s].
        const double computedElapsedTime = orbital_element_conversions
                ::convertMeanAnomalyChangeToElapsedTime(
                    ellipticalMeanAnomalyChange, earthGravitationalParameter, semiMajorAxis );

        // Check if computed elapsed time matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedElapsedTime, computedElapsedTime, 1.0e-15 );
    }

    // Case 5: Hyperbolic orbit around the Sun.
    // The benchmark data is obtained by simply repeating the equations in MATLAB. This test needs
    // to be replaced with real "external" benchmark data.
    {
        // Set hyperbolic mean anomaly [rad].
        const double hyperbolicMeanAnomaly = 2.495601869539691e3;

        // Set Sun gravitational parameter [m^3/s^-2].
        const double sunGravitationalParameter = 3.9859383624e14;

        // Set semi-major axis [m].
        const double semiMajorAxis = -40000.0;

        // Set expected elapsed time [s].
        const double expectedElapsedTime = 1000.0;

        // Compute elapsed time [s].
        const double computedElapsedTime = orbital_element_conversions
                ::convertMeanAnomalyChangeToElapsedTime(
                    hyperbolicMeanAnomaly, sunGravitationalParameter, semiMajorAxis );

        // Check if computed elapsed time matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedElapsedTime, computedElapsedTime, 1.0e-15 );
    }
}

//! Test if conversion from mean motion to semi-major axis is working correctly.
BOOST_AUTO_TEST_CASE( testMeanMotionToSemiMajorAxisConversion )
{
    // Case 1: Geostationary satellite around the Earth.
    // The benchmark data is obtained from (Wikipedia, 2012).
    {
        // Set elliptical mean motion [rad/s].
        const double ellipticalMeanMotion = 7.2921e-5;

        // Set Earth gravitational parameter [m^3/s^2].
        const double earthGravitationalParameter = 5.9736e24 * 6.67428e-11;

        // Set expected semi-major axis [m].
        const double expectedSemiMajorAxis = 42164.0e3;

        // Compute semi-major axis [s].
        const double computedSemiMajorAxis
                = orbital_element_conversions
                ::convertEllipticalMeanMotionToSemiMajorAxis( ellipticalMeanMotion,
                                                              earthGravitationalParameter );

        // Check if computed semi-major axis matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedSemiMajorAxis, computedSemiMajorAxis, 1.0e-4 );
    }

    // Case 2: Geostationary satellite around Mars.
    // The benchmark data is obtained from (Wikipedia, 2012).
    {
        // Set elliptical mean motion [rad/s].
        const double ellipticalMeanMotion = 7.088218e-5;

        // Set Mars gravitational parameter [m^3/s^2].
        const double marsGravitationalParameter = 42828.0e9;

        // Set expected semi-major axis [m].
        const double expectedSemiMajorAxis = 20427.0e3;

        // Compute semi-major axis [m].
        const double computedSemiMajorAxis
                = orbital_element_conversions
                ::convertEllipticalMeanMotionToSemiMajorAxis( ellipticalMeanMotion,
                                                              marsGravitationalParameter );

        // Check if computed semi-major axis matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedSemiMajorAxis, computedSemiMajorAxis, 1.0e-4 );
    }
}

//! Test if conversion from semi-major axis to mean motion is working correctly.
BOOST_AUTO_TEST_CASE( testSemiMajorAxisToMeanMotionConversion )
{
    // Case 1: Geostationary satellite around the Earth.
    // The benchmark data is obtained from (Wikipedia, 2012).
    {
        // Set semi-major axis [m].
        const double semiMajorAxis = 42164.0e3;

        // Set Earth gravitational parameter [m^3/s^2].
        const double earthGravitationalParameter = 5.9736e24 * 6.67428e-11;

        // Set expected elliptical mean motion [rad/s].
        const double expectedEllipticalMeanMotion = 7.2921e-5;

        // Compute elliptical mean motion [rad/s].
        const double computedEllipticalMeanMotion
                = orbital_element_conversions
                ::convertSemiMajorAxisToEllipticalMeanMotion( semiMajorAxis,
                                                              earthGravitationalParameter );

        // Check if computed elliptical mean motion matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedEllipticalMeanMotion,
                                    computedEllipticalMeanMotion, 1.0e-3 );
    }

    // Case 2: Geostationary satellite around Mars.
    // The benchmark data is obtained from (Wikipedia, 2012).
    {
        // Set semi-major axis [m].
        const double semiMajorAxis = 20427.0e3;

        // Set Mars gravitational parameter [m^3/s^2].
        const double marsGravitationalParameter = 42828.0e9;

        // Set expected elliptical mean motion [rad/s].
        const double expectedEllipticalMeanMotion = 7.088218e-5;

        // Compute mean motion [rad/s].
        const double computedEllipticalMeanMotion
                = orbital_element_conversions
                ::convertSemiMajorAxisToEllipticalMeanMotion( semiMajorAxis,
                                                              marsGravitationalParameter );

        // Check if computed elliptical mean motion matches the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedEllipticalMeanMotion,
                                    computedEllipticalMeanMotion, 1.0e-4 );
    }

}

//! Test bug-fix to orbital element conversions for arguments of periapsis close to zero, pi or two pi.
//! Thanks to Andreas Kleinschneider for finding the bug.
BOOST_AUTO_TEST_CASE( test_ArgumentOfPeriapsisBugfix )
{

    // Define previously offending states and gravitational parameters.
    Eigen::Vector6d centralBodyCartesianState;
    centralBodyCartesianState << -521142852074.35858154296875, 595141535550.8497314453125, 8056093690.3882598876953125,
            3904.21508802941389149054884911, -18262.3110776023386279121041298, 103.914916159730324807242141105;

    Eigen::Vector6d orbitingBodyCartesianState;
    orbitingBodyCartesianState << -520894562964.4881591796875, 595481881288.00537109375, 8071968111.222164154052734375,
            -10065.8995464500585512723773718, -7984.05278915743656398262828588, 261.213167975313467650266829878;

    Eigen::Vector6d relativeCartesianState = orbitingBodyCartesianState - centralBodyCartesianState;

    double gravitationalParameterOfCentralBody = 126686534921800800;
    double gravitationalParameterOfOrbitingBody = 5959916033410.404296875;

    // Convert to Keplerian state
    Eigen::Vector6d keplerianState =
            orbital_element_conversions::convertCartesianToKeplerianElements(
                relativeCartesianState, gravitationalParameterOfCentralBody + gravitationalParameterOfOrbitingBody );

    // Check whether argument of periapsis is defined.
    BOOST_CHECK_EQUAL( ( keplerianState( orbital_element_conversions::argumentOfPeriapsisIndex ) ==
                         keplerianState( orbital_element_conversions::argumentOfPeriapsisIndex ) ), true );

}

BOOST_AUTO_TEST_CASE( test_LongitudeOfNodeBugfix )
{
    std::cout << std::endl << std::endl;
    using namespace tudat;

    // Define previously offending states and gravitational parameters.
    Eigen::Vector6d cartesianState;
    cartesianState << 146378739288.0336,  -32851886854.1209, -14241055658.24648,  6603.183449760695,
            26444.27948911581,  11463.63694918501;
    double gravitationalParameterOfCentralBody = 1.327128386237518e+20;

    // Convert to Keplerian state
    Eigen::Vector6d keplerianState =
            orbital_element_conversions::convertCartesianToKeplerianElements(
                cartesianState, gravitationalParameterOfCentralBody );

    Eigen::Vector6d recomputedCartesianState =
            orbital_element_conversions::convertKeplerianToCartesianElements(
                keplerianState, gravitationalParameterOfCentralBody );

    Eigen::Vector6d  recomputedKeplerianState =
            orbital_element_conversions::convertCartesianToKeplerianElements(
                recomputedCartesianState, gravitationalParameterOfCentralBody );

//    std::cout << std::setprecision( 16 ) << cartesianState.transpose( ) << std::endl;
//    std::cout << std::setprecision( 16 ) << recomputedCartesianState.transpose( ) << std::endl << std::endl;

//    std::cout << std::setprecision( 16 ) << keplerianState.transpose( ) << std::endl;
//    std::cout << std::setprecision( 16 ) << recomputedKeplerianState.transpose( ) << std::endl;
}


//! Test if conversion from eccentric anomaly to mean anomaly is working correctly.
BOOST_AUTO_TEST_CASE( testMeanToTrueAnomalyConversion )
{
    using namespace orbital_element_conversions;
    using namespace basic_mathematics;

    std::vector< double > meanAnomalies =
    { 1.0E-12, 1.0, mathematical_constants::PI, 5.0, 2.0 * mathematical_constants::PI - 1.0E-12 };
    std::vector< double > eccentricities =
    { 0.0, 0.1, 0.5, 0.9, 0.99, 0.99999999, 1.00000001, 1.001, 2.0, 5.0, 10.0 };
    for( unsigned int i = 0; i < meanAnomalies.size( ); i++ )
    {
        for( unsigned int j = 0; j < eccentricities.size( ); j++ )
        {
            double eccentricity = eccentricities.at( j );
            double meanAnomaly = meanAnomalies.at( i );

            double trueAnomaly = convertMeanAnomalyToTrueAnomaly( eccentricity, meanAnomaly );

            double recomputedMeanAnomaly = convertTrueAnomalyToMeanAnomaly( eccentricity, trueAnomaly );

            double anomalyDifference = meanAnomaly - recomputedMeanAnomaly;
            if( anomalyDifference > 1.0 )
            {
                anomalyDifference -= 2.0 * mathematical_constants::PI;
            }

            double tolerance = 1.0E-13;
            if( std::fabs( eccentricity - 1 ) < 0.1 )
            {
                tolerance *= 1.0E3;
            }

            if( std::fabs( eccentricity - 1 ) < 0.00001 )
            {
                tolerance *= 1.0E6;
            }
            BOOST_CHECK_SMALL( std::fabs( anomalyDifference ), tolerance );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
