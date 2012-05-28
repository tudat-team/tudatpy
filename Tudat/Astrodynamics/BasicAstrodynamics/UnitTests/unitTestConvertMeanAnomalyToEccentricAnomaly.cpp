/*    Copyright (c) 2010-2012, Delft University of Technology
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
