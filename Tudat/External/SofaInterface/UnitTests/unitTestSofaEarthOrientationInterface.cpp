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

#include <limits>
#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/External/SofaInterface/earthOrientation.h"
#include "Tudat/External/SofaInterface/sofaTimeConversions.h"
#include "Tudat/Basics/timeType.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::sofa_interface;

BOOST_AUTO_TEST_SUITE( test_sofa_earth_orientation )

//! Test precession/nutation calculations from Sofa, using Sofa Earth orientation cookbook results.
BOOST_AUTO_TEST_CASE( testSofaPrecessionNutation )
{
    // Two-part Terrestrial time at which rotation is calculated.
    double testJulianDay1 = 2400000.5;
    double testJulianDay2 = 54195.500754444444444 * physical_constants::JULIAN_DAY;

    double arcSecondToRadian = 4.848136811095359935899141E-6;

    // Perform test for IAU_2000 and IAU_2006 cnventions
    for( unsigned int i = 0; i < 2; i++ )
    {
        // Get output data from SOFA cookbook.
        Eigen::Vector2d expectedPolePosition;
        double expectedCioLocator = 0.0;
        double dXTest = 0.0;
        double dYTest = 0.0;
        basic_astrodynamics::IAUConventions iauConventions;
        if( i == 0 )
        {
            iauConventions = basic_astrodynamics::iau_2000_a;

            expectedPolePosition.x( ) = 0.000712264729599;
            expectedPolePosition.y( ) = 0.000044385250426;
            expectedCioLocator = -0.002200475 * arcSecondToRadian;
            dXTest = 0.0001750 * arcSecondToRadian;
            dYTest = -0.000265 * arcSecondToRadian;
        }
        else
        {
            iauConventions = basic_astrodynamics::iau_2006;

            expectedPolePosition.x( ) = 0.000712264729525;
            expectedPolePosition.y( ) = 0.000044385248875;
            expectedCioLocator = -0.002200475 * arcSecondToRadian;
            dXTest = 0.0001750 * arcSecondToRadian;
            dYTest = -0.0002259 * arcSecondToRadian;

        }

        // Calculate X, Y pole position and CIO locator and compare against cookbook results.
        std::pair< Eigen::Vector2d, double > cipInGcrs =
                getPositionOfCipInGcrs( testJulianDay2, testJulianDay1, iauConventions );
        BOOST_CHECK_SMALL( expectedPolePosition.x( ) - ( cipInGcrs.first.x( ) + dXTest ), 5.0E-11 );
        BOOST_CHECK_SMALL( expectedPolePosition.y( ) - ( cipInGcrs.first.y( ) + dYTest ), 5.0E-11 );
        BOOST_CHECK_SMALL( expectedCioLocator - cipInGcrs.second, 2.0E-12 );

    }

}

//! Test GMST and ERA functions from Sofa.
BOOST_AUTO_TEST_CASE( testSofaEarthRotation )
{

    double testJulianDay1 = 2400000.5 + 54195;
    double testJulianDay2 = 0.500754444444444;

    double testUt1 = 2400000.5 + 54195;
    double testUt2 = 0.499999165813831;

    // Test correct copy/translation of Sofa fortran test from cookbook (GMST not directly given in cookbook).
    {
        double deltaPsi, deltaEpsilon;
        double meanObliquity;
        double rb[3][3], rp[3][3], rbp[3][3];
        double rn[3][3], rbpn[3][3];
        iauNut00a( testJulianDay1, testJulianDay2, &deltaPsi, &deltaEpsilon );
        iauPn00( testJulianDay1, testJulianDay2, deltaPsi, deltaEpsilon, &meanObliquity, rb, rp, rbp, rn, rbpn );
        // Transform dX,dY corrections from GCRS to mean of date.

        double arcSecondToRadian = 4.848136811095359935899141E-6;

        double dXTest = 0.17250 * arcSecondToRadian / 1000.0;
        double dYTest = -0.26500 * arcSecondToRadian / 1000.0;
        double angleVector1[3], angleVector2[3];
        angleVector1[0] = dXTest;
        angleVector1[1] = dYTest;
        angleVector1[2] = 0.0;
        iauRxp( rbpn, angleVector1, angleVector2 );
        double deltaPsiCorrection = angleVector2[0] / std::sin( meanObliquity );

        // Corrected nutation.;
        deltaPsi += deltaPsiCorrection;


        // Greenwich apparent sidereal time (IAU 1982/1994).
        double gst = ( iauGmst00( testUt1, testUt2, testJulianDay1, testJulianDay2 ) +
                       iauEe00( testJulianDay1, testJulianDay2, meanObliquity, deltaPsi ) )
                * 180.0 / mathematical_constants::PI;

        // Check if test data is correctly reprodcued
        double expectedAngle = 13.412417084674;
        BOOST_CHECK_SMALL( std::fabs( gst - expectedAngle ), 1.0E-12 );
    }

    // Compare direct against indirect GMST calculation.
    double expectedGmst = iauGmst00( testUt1, testUt2, testJulianDay1, testJulianDay2 );
    double calculatedGmst = calculateGreenwichMeanSiderealTime(
                testJulianDay2 * physical_constants::JULIAN_DAY, testUt2 * physical_constants::JULIAN_DAY,
                testJulianDay1, basic_astrodynamics::iau_2000_a );

    BOOST_CHECK_CLOSE_FRACTION( expectedGmst, calculatedGmst, std::numeric_limits< double >::epsilon( ) );

    // Calculate Earth rotation angle and comapre against result in cookbook.
    double earthRotationAngle = calculateEarthRotationAngle(
                testUt2 * physical_constants::JULIAN_DAY, testUt1 ) * 180.0 / mathematical_constants::PI;
    double expectedEarthRotationAngle = 13.318492966097;

    BOOST_CHECK_SMALL( std::fabs( expectedEarthRotationAngle - earthRotationAngle ), 1.0E-12 );

    earthRotationAngle = calculateEarthRotationAngleTemplated< double >(
                    ( testUt2 + testUt1 - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
                physical_constants::JULIAN_DAY ) * 180.0 / mathematical_constants::PI;
    BOOST_CHECK_SMALL( std::fabs( expectedEarthRotationAngle - earthRotationAngle ), 1.0E-7 );

    earthRotationAngle = calculateEarthRotationAngleTemplated< Time >(
                    tudat::Time( ( testUt1 - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * 24,
                                 static_cast< long double >( testUt2 ) * physical_constants::JULIAN_DAY_LONG ) ) *
            180.0 / mathematical_constants::PI;
    BOOST_CHECK_SMALL( std::fabs( expectedEarthRotationAngle - earthRotationAngle ), 1.0E-12 );


    // Test ERA for negative time since J2000
    testUt1 = basic_astrodynamics::JULIAN_DAY_ON_J2000 - 200.0;
    testUt2 = 0.499999165813831;
    double directEarthRotationAngle = calculateEarthRotationAngle(
                testUt2 * physical_constants::JULIAN_DAY, testUt1 ) * 180.0 / mathematical_constants::PI;

    earthRotationAngle = calculateEarthRotationAngleTemplated< double >(
                    ( testUt2 + testUt1 - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
                physical_constants::JULIAN_DAY ) * 180.0 / mathematical_constants::PI;
    BOOST_CHECK_SMALL( std::fabs( directEarthRotationAngle - earthRotationAngle ), 1.0E-7 );

    earthRotationAngle = calculateEarthRotationAngleTemplated< Time >(
                    tudat::Time( ( testUt1 - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * 24,
                                 static_cast< long double >( testUt2 ) * physical_constants::JULIAN_DAY_LONG ) ) *
            180.0 / mathematical_constants::PI;
    BOOST_CHECK_SMALL( std::fabs( directEarthRotationAngle - earthRotationAngle ), 1.0E-12 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat


