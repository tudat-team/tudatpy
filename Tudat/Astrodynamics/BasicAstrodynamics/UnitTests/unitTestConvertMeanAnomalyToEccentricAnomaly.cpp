/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      111215    T. Secretin       First creation of the code.
 *      111221    T. Secretin       Removed memory leaks. Added test for circular and
 *                                  near-parabolic orbits, as well as for negative eccentricities.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120421    K. Kumar          Updated test fixtures and cases to use updated conversion
 *                                  object.
 *
 *    References
 *      http://www.esa.int/gsp/ACT/doc/INF/Code/globopt/GTOPtoolbox.rar
 *
 */

#define BOOST_TEST_MAIN

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToEccentricAnomaly.h"

namespace tudat
{
namespace unit_tests
{

//! Conversion test fixture.
/*!
 * Conversion test fixture used by the Boost unit test framework. This code is executed before each
 * test.
 */
struct conversion_test_fixture
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    conversion_test_fixture( )
    {
        toleranceOrbitalElementConversion = 1e-8;
    }

    //! Conversion tolerance to test against.
    /*!
     * Conversion tolerance to test against.
     */
    double toleranceOrbitalElementConversion;

    //! Convert mean anomaly to eccentric anomaly.
    /*!
     * Converts mean anomaly to eccentric anomaly by creating a conversion object to test and
     * executing the conversion.
     * \param eccentricity Eccentricity [-].
     * \param meanAnomaly Mean anomaly [rad].
     * \return eccentricAnomaly Eccentric anomaly [rad].
     */
    double convertMeanAnomalyToEccentricAnomaly( const double eccentricity,
                                                 const double meanAnomaly )
    {
        // Conversion object to test; mean anomaly to eccentric anomaly conversion.
        tudat::orbital_element_conversions::ConvertMeanAnomalyToEccentricAnomaly
                meanToEccentricAnomaly( eccentricity, meanAnomaly,
                                        boost::make_shared< NewtonRaphson >( ) );

        // Convert to eccentric anomaly and return.
        return meanToEccentricAnomaly.convert( );
    }

protected:

private:
};

// Declare Boost fixture test suite.
BOOST_FIXTURE_TEST_SUITE( testsuite_convertMeanAnomalyToEccentricAnomaly, conversion_test_fixture )

// Test 1: Test conversion for circular orbits.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToEccentricAnomaly_circular )
{
    // Set test value for eccentricity.
    double testEccentricity = 0.0;

    // Set test value for mean anomaly.
    double testMeanAnomaly = 1.0472;

    // Set reference value for eccentric anomaly;
    double referenceEccentricAnomaly = 1.0472;

    // Compute eccentric anomaly.
    double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly(
                testEccentricity, testMeanAnomaly );

    // Check if computed eccentric anomaly is less than error tolerance.
    BOOST_CHECK_CLOSE( eccentricAnomaly, 
                       referenceEccentricAnomaly, 
                       toleranceOrbitalElementConversion );
}

// Test 2: Test conversion for valid range of eccentricities.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToEccentricAnomaly_range )
{
    // Set test values for eccentricity.
    double arrayOfTestEccentricities [4] = { 0.01671, 0.43582, 0.78514, 0.91525 };

    // Set test values for mean anomaly.
    double arrayOfTestMeanAnomalies [4] = { 
        tudat::unit_conversions::convertDegreesToRadians( 60.0 ),
        tudat::unit_conversions::convertDegreesToRadians( 90.0 ),
        tudat::unit_conversions::convertDegreesToRadians( 120.0 ),
        tudat::unit_conversions::convertDegreesToRadians( 220.0 ) };

    // Set reference values for eccentric anomaly;
    double arrayOfReferenceEccentricAnomalies [4] = { 1.06178920406832,
                                                      1.97200731113253,
                                                      2.5392410896466,
                                                      3.51006218528448 };

    // Loop over sets of data.
    for ( int i = 0; i < 4; i++ )
    {
        // Compute eccentric anomaly.
        double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly(
                    arrayOfTestEccentricities[ i ], arrayOfTestMeanAnomalies[ i ] );

        // Check if computed eccentric anomaly is less than error tolerance.
        BOOST_CHECK_CLOSE( eccentricAnomaly,
                           arrayOfReferenceEccentricAnomalies[ i ],
                           toleranceOrbitalElementConversion );
    }
}

// Test 3: Test conversion for negative eccentricities.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToEccentricAnomaly_negative )
{
    // Set test value for eccentricity.
    double testEccentricity = -0.5;

    // Set test value for mean anomaly.
    double testMeanAnomaly = 1.0472;

    // Compute eccentric anomaly.
    double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly(
                testEccentricity, testMeanAnomaly );

    // Check if computed eccentric anomaly is NaN for negative eccentricity.
    BOOST_CHECK( boost::math::isnan( eccentricAnomaly ) );
}

// Test 4: Test conversion for near-parabolic orbits.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToEccentricAnomaly_nearParabolic )
{
    // Set test value for eccentricity.
    double testEccentricity = 0.99;

    // Set test value for mean anomaly.
    double testMeanAnomaly = 1.0472;

    // Compute eccentric anomaly.
    double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly(
                testEccentricity, testMeanAnomaly );

    // Check if computed eccentric anomaly is NaN for parabolic orbits.
    BOOST_CHECK( boost::math::isnan( eccentricAnomaly ) );
}

// End Boost test suite.
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
