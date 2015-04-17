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
 *      120822    P. Musegaas       First creation of the code. (7 test cases, including various
 *                                  random test cases due to chaotic nature of root-finder for
 *                                  some starter values.)
 *      120903    P. Musegaas       Improved random test (does not fail on mean anomaly of 0.0).
 *                                  Decreased number of random values for random tests.
 *      121205    P. Musegaas       Updated code to final version of rootfinders.
 *      130123    K. Kumar          Added separated test tolerance for near-parabolic cases in
 *                                  Test 4 to deal with conversion failure on some systems.
 *
 *    References
 *      GTOP, http://www.esa.int/gsp/ACT/doc/INF/Code/globopt/GTOPtoolbox.rar.
 *
 *    Notes
 *      Note that for some of the near-parabolic cases, the tolerance used for the to-and-fro
 *      conversions (Test 4) is several order of magnitudes higher than used for the regular cases.
 *      This should be investigated further in the future to fully characterize the nature of the
 *      conversions in the near-parabolic cases.
 *
 */

#define BOOST_TEST_MAIN

#include <ctime>
#include <fstream>
#include <string>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToHyperbolicEccentricAnomaly.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

using namespace orbital_element_conversions;

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
        toleranceOrbitalElementConversion = 1.0e-14;
        toleranceOrbitalElementConversionNearParabolic = 1.0e-9;
    }

    //! Conversion tolerance to test against.
    /*!
     * Conversion tolerance to test against.
     */
    double toleranceOrbitalElementConversion;

    //! Conversion tolerance to test against for near-parabolic cases.
    /*!
     * Conversion tolerance to test against for near-parabolic cases.
     */
    double toleranceOrbitalElementConversionNearParabolic;

    //! Convert mean anomaly to eccentric anomaly.
    /*!
     * Converts mean anomaly to eccentric anomaly by creating a conversion object to test and
     * executing the conversion.
     * \param eccentricity Eccentricity [-].
     * \param meanAnomaly Mean anomaly [rad].
     * \param useDefaultInitialGuess Boolean specifying whether to use default initial guess [-].
     * \param initialGuess Initial guess for rootfinder [rad].
     * \return eccentricAnomaly Eccentric anomaly [rad].
     */
    double convertMeanAnomalyToHyperbolicEccentricAnomaly( const double eccentricity,
                                                           const double meanAnomaly,
                                                           const bool useDefaultInitialGuess = true,
                                                           const double initialGuess = TUDAT_NAN )
    {
        // Conversion object to test; mean anomaly to eccentric anomaly conversion.
        orbital_element_conversions::
                ConvertMeanAnomalyToHyperbolicEccentricAnomaly  meanToHyperbolicEccentricAnomaly(
                    eccentricity, meanAnomaly, useDefaultInitialGuess, initialGuess );

        // Convert to eccentric anomaly and return.
        return meanToHyperbolicEccentricAnomaly.convert( );
    }

protected:

private:
};

//! Error writing function.
/*!
 * This function writes the input values that led to errors to a unique file, if any errors occured
 * during the random tests. To make the file unique, the date and time of execution is added. An
 * error message also shows the location of the file.
 * \param eccentricities A vector containing the eccentricities that caused a failure.
 * \param meanAnomalies A vector containing the mean anomalies that caused a failure.
 * \param testName A string specifying the name of the test that failed.
 */
void writeErrorsToFile( std::vector< double > eccentricities, std::vector< double > meanAnomalies,
                        std::string testName )
{
    // Obtain the current time.
    const boost::posix_time::ptime now = boost::posix_time::second_clock::local_time( );

    // Make a string containing the output file name. This output file is tagged with the date and
    // time at which the code was executed. The default date format is: YYYYMMDDTHHMMSS, in which T
    // separates date and time.
    const std::string outputFileName = input_output::getTudatRootPath( ) +
                                       "Astrodynamics/BasicAstrodynamics/UnitTests/" +
                                       "ErrorReportConversionMeanToHyperbolicEccentricAnomaly" +
                                       testName + "RunAt" + boost::posix_time::to_iso_string( now )
                                       + ".txt";

    // Make a stream to a file.
    std::ofstream errorFile( outputFileName.c_str( ) );

    // Write an introduction in the file explaining what happened. 70 lines long.
    errorFile << "This error report was generated because the unit test for the" << std::endl
              << "conversion of mean to hyperbolic eccentric anomaly has failed in" << std::endl
              << "one of the random tests. To ensure the data for which it failed is" << std::endl
              << "not lost, the corresponding input variables for these cases are" << std::endl
              << "listed below. Please report a bug on the Tudat website "<< std::endl
              << "(tudat.tudelft.nl), with these values, so that someone will look" << std::endl
              << "into it and the code can be improved." << std::endl << std::endl
              << "Eccentricities:           Mean anomalies:" << std::endl;

    // Set the precision for the output of the variables to 16 digits.
    errorFile.precision( 16 );

    // Add the corresponding eccentricities and mean anomalies at neatly arranged positions.
    for ( unsigned int counter = 0; counter < eccentricities.size( ); counter++ )
    {
        errorFile << std::setw( 25 ) << eccentricities[ counter ]
                  << std::setw( 25 ) << meanAnomalies[ counter ] << std::endl;
    }
    errorFile.close( );

    // Add an error message specifying the file that the values have been written to.
    std::cerr << "One or multiple errors occurred during random sampling. " << std::endl
              << "The values leading to these errors have been written to the following file: "
              << std::endl << outputFileName;
}

// Declare Boost fixture test suite.
BOOST_FIXTURE_TEST_SUITE( testsuite_convertMeanAnomalyToHyperbolicEccentricAnomaly,
                          conversion_test_fixture )

//! Test 1: Test a range of values for the conversion.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToHyperbolicEccentricAnomaly_range )
{
    // Set array of test values for eccentricity.
    const double arrayOfTestEccentricities [ 6 ] = { 1.03, 1.28, 1.97, 2.56, 10.87, 99.72 };

    // Set array of test values for mean anomaly.
    const double arrayOfTestMeanAnomalies [ 6 ] = { 1.5, 6.0, 0.5, -4.0, 5.5, 2.5 };

    // Set array of expected values for hyperbolic eccentric anomaly. These values were converted
    // back and forth to verify their correctness. Also the conversion was verified by comparing
    // with GTOP. GTOP uses a different definition for the hyperbolic eccentric anomaly. Hence
    // the mean anomalies were converted to cartesian positions and compared to the same conversion
    // using Tudat methods. (GTOP also does not use true anomaly as explicit step.)
    const double arrayOfExpectedHyperbolicAnomalies [ 6 ] = { 1.9132897042137,
                                                              2.60400218106873,
                                                              0.478057581067141,
                                                              -1.50971422579796,
                                                              0.529595060186511,
                                                              0.0253214157050963 };

    // Loop over sets of data.
    for ( int counter = 0; counter < 6; counter++ )
    {
        // Compute the hyperbolic eccentric anomaly.
        const double hyperbolicEccentricAnomaly = convertMeanAnomalyToHyperbolicEccentricAnomaly(
                    arrayOfTestEccentricities[ counter ], arrayOfTestMeanAnomalies[ counter ] );

        // Check if computed eccentric anomaly is less than error tolerance.
        BOOST_CHECK_CLOSE_FRACTION( arrayOfExpectedHyperbolicAnomalies[ counter ],
                                    hyperbolicEccentricAnomaly,
                                    toleranceOrbitalElementConversion );
    }
}

//! Test 2: Test a value that is out of range.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToHyperbolicEccentricAnomaly_TooLow )
{
    // Set test value for eccentricity.
    const double testEccentricity = 0.5;

    // Set test value for mean anomaly.
    const double testMeanAnomaly = 0.5;

    // Check if a runtime error is thrown if the anomaly is converted for this eccentricity.
    BOOST_CHECK_THROW( convertMeanAnomalyToHyperbolicEccentricAnomaly(
                           testEccentricity, testMeanAnomaly ), std::runtime_error );
}

//! Test 3: Test conversion for near-parabolic orbits.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToHyperbolicEccentricAnomaly_nearParabolic )
{
    // Set test value for eccentricity.
    const double testEccentricity = 1.0 + 1.0e-10;

     // Set array of test values for mean anomaly.
    const double arrayOfTestMeanAnomalies[ 4 ] = { -10.0, -1.4, 0.5, 7.6 };

    // Set array of expected values for hyperbolic eccentric anomaly. These values were converted
    // back and forth to verify their correctness. Also the conversion was verified by comparing
    // with GTOP. GTOP uses a different definition for the hyperbolic eccentric anomaly. Hence
    // the mean anomalies were converted to cartesian positions and compared to the same conversion
    // using Tudat methods. (GTOP also does not use true anomaly as explicit step)
    const double arrayOfExpectedHyperbolicAnomalies [ 4 ] = { -3.280887528670698,
                                                              -1.913052492643601,
                                                              1.396250871565077,
                                                              3.062027761338891 };

    // Loop over sets of data.
    for ( int counter = 0; counter < 4; counter++ )
    {
        // Compute the hyperbolic eccentric anomaly.
        const double hyperbolicEccentricAnomaly = convertMeanAnomalyToHyperbolicEccentricAnomaly(
                    testEccentricity, arrayOfTestMeanAnomalies[ counter ] );

        // Check if computed eccentric anomaly is less than error tolerance.
        BOOST_CHECK_CLOSE_FRACTION( arrayOfExpectedHyperbolicAnomalies[ counter ],
                                    hyperbolicEccentricAnomaly,
                                    toleranceOrbitalElementConversion );
    }

}

//! Test 4: Test a large number of anomalies for near-parabolic orbits.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToHyperbolicEccentricAnomaly_randomNearParabolic )
{
    // Create vectors that will store the input variables of a test that resulted in an error, such
    // that the error scenario can be reproduced.
    std::vector< double > failedMeanAnomalies, failedEccentricities;

    // Boolean that will be set true if a runtime error occurred.
    bool aRuntimeErrorOccurred = false;

    // Set test value for eccentricity, which is just above 1.0.
    const double testEccentricity = 1.0 + 1.0e-15;

    // Initialize both test and reverse calculated hyperbolic mean anomaly and the hyperbolic
    // eccentric anomaly.
    double testMeanAnomaly, reverseCalculatedMeanAnomaly, hyperbolicEccentricAnomaly = TUDAT_NAN;

    // Instantiate a random number generator for the mean anomaly generation, from -20 to 20.
    boost::mt19937 randomNumbergenerator( time( 0 ) );
    boost::random::uniform_real_distribution< > distributionMinus20To20( -20.0, 20.0 );
    boost::variate_generator< boost::mt19937&, boost::random::uniform_real_distribution < > >
            generateRandomNumberMinus20To20( randomNumbergenerator, distributionMinus20To20 );

    // Specify the number of random samples should be taken. A test of 100,000,000 was performed
    // by the author before the code was submitted. This test remains included to verify that any
    // future method will not fail. The behaviour of the conversion is namely very sensitive and
    // non-converging cases are highly sensitive to input values and the initial guess that is used.
    const int numberOfSamples = 10000;

    // Perform the conversion for the specified number of samples and test whether the values that
    // are subsequently converted back match the initial values.
    for ( int counter = 0; counter < numberOfSamples; counter++ )
    {
        // Set random value in test mean anomaly and the test eccentricity.
        testMeanAnomaly = generateRandomNumberMinus20To20( );

        // If the Rootfinder does not converge, it will produce a runtime error. In order to make
        // sure that these values that led to the error will not be lost, they will be stored in
        // the failed input data vectors. To do so, a try-catch sequence is used.
        try
        {
            // Compute eccentric anomaly.
            hyperbolicEccentricAnomaly = convertMeanAnomalyToHyperbolicEccentricAnomaly(
                    testEccentricity, testMeanAnomaly );
        }
        catch( std::runtime_error )
        {
            // Store the fact that a runtime error occurred, such that the values will be stored.
            aRuntimeErrorOccurred = true;
        }

        // Calculate the mean anomaly from this eccentric anomaly.
        reverseCalculatedMeanAnomaly =
                convertHyperbolicEccentricAnomalyToMeanAnomaly( hyperbolicEccentricAnomaly,
                                                                testEccentricity );

        // Test whether the computed mean anomaly is equal to the mean anomaly from the input and
        // that no runtime errors occurred. If an error was found, store the values leading to this
        // error in a vector for later use. '!' operator is there to ensure that a NaN value will
        // result in the values being written away. It is also checked that the mean anomaly is
        // not equal to 0.0, because that would result in falsely writing an error.
        if ( ( ( !( std::abs( testMeanAnomaly - reverseCalculatedMeanAnomaly ) / testMeanAnomaly <
                toleranceOrbitalElementConversionNearParabolic ) ) && ( testMeanAnomaly != 0.0 ) )
             || aRuntimeErrorOccurred )
        {
            failedMeanAnomalies.push_back( testMeanAnomaly );
            failedEccentricities.push_back( testEccentricity );
        }

        // Reset boolean.
        aRuntimeErrorOccurred = false;
    }

    // Check that no values have been written to the failedMeanAnomalies vector.  If so, this test
    // is passed. Otherwisely these values will be written away and this test will fail.
    BOOST_CHECK( failedMeanAnomalies.empty( ) );

    // If the vector is not empty, write the failed cases of this test case to a file.
    if ( !( failedMeanAnomalies.empty( ) ) )
    {
        writeErrorsToFile( failedEccentricities, failedMeanAnomalies, "Test4" );
    }
}

//! Test 5: Test large number of anomalies with a large number of eccentricities, in common regime.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToHyperbolicEccentricAnomaly_random )
{
    // Create vectors that will store the input variables of a test that resulted in an error, such
    // that the error scenario can be reproduced.
    std::vector< double > failedMeanAnomalies, failedEccentricities;

    // Boolean that will be set true if a runtime error occurred.
    bool aRuntimeErrorOccurred = false;

    // Initialize both test and reverse calculated hyperbolic mean anomaly, the hyperbolic
    // eccentric anomaly and the eccentricity.
    double testEccentricity, testMeanAnomaly, reverseCalculatedMeanAnomaly,
           hyperbolicEccentricAnomaly = TUDAT_NAN;

    // Instantiate random number generators. One for the mean anomaly generation, from -20 to 20,
    // another one for the eccentricity generation, from 1 to 10.
    boost::mt19937 randomNumbergenerator( time( 0 ) );
    boost::random::uniform_real_distribution< > distributionMinus20To20( -20.0, 20.0 );
    boost::variate_generator< boost::mt19937&, boost::random::uniform_real_distribution < > >
            generateRandomNumberMinus20To20( randomNumbergenerator, distributionMinus20To20 );
    boost::random::uniform_real_distribution< > distribution1To10( 1.0 + 1.0e-15, 10.0 );
    boost::variate_generator< boost::mt19937&, boost::random::uniform_real_distribution < > >
            generateRandomNumber1To10( randomNumbergenerator, distribution1To10 );

    // Specify the number of random samples should be taken. A test of 100,000,000 was performed
    // by the author before the code was submitted. This test remains included to verify that any
    // future method will not fail. The behaviour of the conversion is namely very sensitive and
    // non-converging cases are highly sensitive to input values and the initial guess that is used.
    const int numberOfSamples = 10000;

    // Perform the conversion for the specified number of samples and test whether the values that
    // are subsequently converted back match the initial values.
    for ( int counter = 0; counter < numberOfSamples; counter++ )
    {
        // Set random value in test mean anomaly and the test eccentricity.
        testMeanAnomaly = generateRandomNumberMinus20To20( );
        testEccentricity = generateRandomNumber1To10( );

        // If the Rootfinder does not converge, it will produce a runtime error. In order to make
        // sure that these values that led to the error will not be lost, they will be stored in
        // the failed input data vectors. To do so, a try-catch sequence is used.
        try
        {
            // Compute eccentric anomaly.
            hyperbolicEccentricAnomaly = convertMeanAnomalyToHyperbolicEccentricAnomaly(
                    testEccentricity, testMeanAnomaly );
        }
        catch( std::runtime_error )
        {
            // Store the fact that a runtime error occurred, such that the values will be stored.
            aRuntimeErrorOccurred = true;
        }

        // Calculate the mean anomaly from this eccentric anomaly.
        reverseCalculatedMeanAnomaly =
                convertHyperbolicEccentricAnomalyToMeanAnomaly( hyperbolicEccentricAnomaly,
                                                                testEccentricity );

        // Test whether the computed mean anomaly is equal to the mean anomaly from the input and
        // that no runtime errors occurred. If an error was found, store the values leading to this
        // error in a vector for later use. '!' operator is there to ensure that a NaN value will
        // result in the values being written away. It is also checked that the mean anomaly is
        // not equal to 0.0, because that would result in falsely writing an error.
        if ( ( ( !( std::abs( testMeanAnomaly - reverseCalculatedMeanAnomaly ) / testMeanAnomaly <
                toleranceOrbitalElementConversion ) ) && ( testMeanAnomaly != 0.0 ) )
             || aRuntimeErrorOccurred )
        {
            failedMeanAnomalies.push_back( testMeanAnomaly );
            failedEccentricities.push_back( testEccentricity );
        }

        // Reset boolean.
        aRuntimeErrorOccurred = false;
    }

    // Check that no values have been written to the failedMeanAnomalies vector.  If so, this test
    // is passed. Otherwisely these values will be written away and this test will fail.
    BOOST_CHECK( failedMeanAnomalies.empty( ) );

    // If the vector is not empty, write the failed cases of this test case to a file.
    if ( !( failedMeanAnomalies.empty( ) ) )
    {
        writeErrorsToFile( failedEccentricities, failedMeanAnomalies, "Test5" );
    }
}

//! Test 6: Test large number of high anomalies with a large number of high eccentricities.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToHyperbolicEccentricAnomaly_random_high )
{
    // Create vectors that will store the input variables of a test that resulted in an error, such
    // that the error scenario can be reproduced.
    std::vector< double > failedMeanAnomalies, failedEccentricities;

    // Boolean that will be set true if a runtime error occurred.
    bool aRuntimeErrorOccurred = false;

    // Initialize both test and reverse calculated hyperbolic mean anomaly, the hyperbolic
    // eccentric anomaly and the eccentricity.
    double testEccentricity, testMeanAnomaly, reverseCalculatedMeanAnomaly,
           hyperbolicEccentricAnomaly = TUDAT_NAN;

    // Instantiate random number generators. One for the mean anomaly generation, from -1.2e12 to
    // 1.2e12, another one for the eccentricity generation, from 1.0 to 1.0e15.
    boost::mt19937 randomNumbergenerator( time( 0 ) );
    boost::random::uniform_real_distribution< > distributionMinus12To12( -12, 12 );
    boost::variate_generator< boost::mt19937&, boost::random::uniform_real_distribution < > >
            generateRandomNumberMinus12To12( randomNumbergenerator, distributionMinus12To12 );
    boost::random::uniform_real_distribution< > distribution0To15( 0, 15 );
    boost::variate_generator< boost::mt19937&, boost::random::uniform_real_distribution < > >
            generateRandomNumber0To15( randomNumbergenerator, distribution0To15 );

    // Specify the number of random samples should be taken. A test of 100,000,000 was performed
    // by the author before the code was submitted. This test remains included to verify that any
    // future method will not fail. The behaviour of the conversion is namely very sensitive and
    // non-converging cases are highly sensitive to input values and the initial guess that is used.
    const int numberOfSamples = 10000;

    // Perform the conversion for the specified number of samples and test whether the values that
    // are subsequently converted back match the initial values.
    for ( int counter = 0; counter < numberOfSamples; counter++ )
    {
        // Set random value in test mean anomaly and the test eccentricity.
        testMeanAnomaly = generateRandomNumberMinus12To12( ) *
                          std::pow( 10, generateRandomNumberMinus12To12( ) );
        testEccentricity = 1.0 + 1.0 * std::pow( 10, generateRandomNumber0To15( ) );

        // If the Rootfinder does not converge, it will produce a runtime error. In order to make
        // sure that these values that led to the error will not be lost, they will be stored in
        // the failed input data vectors. To do so, a try-catch sequence is used.
        try
        {
            // Compute eccentric anomaly.
            hyperbolicEccentricAnomaly = convertMeanAnomalyToHyperbolicEccentricAnomaly(
                    testEccentricity, testMeanAnomaly );
        }
        catch( std::runtime_error )
        {
            // Store the fact that a runtime error occurred, such that the values will be stored.
            aRuntimeErrorOccurred = true;
        }

        // Calculate the mean anomaly from this eccentric anomaly.
        reverseCalculatedMeanAnomaly =
                convertHyperbolicEccentricAnomalyToMeanAnomaly( hyperbolicEccentricAnomaly,
                                                                testEccentricity );

        // Test whether the computed mean anomaly is equal to the mean anomaly from the input and
        // that no runtime errors occurred. If an error was found, store the values leading to this
        // error in a vector for later use. '!' operator is there to ensure that a NaN value will
        // result in the values being written away. It is also checked that the mean anomaly is
        // not equal to 0.0, because that would result in falsely writing an error.
        if ( ( ( !( std::abs( testMeanAnomaly - reverseCalculatedMeanAnomaly ) / testMeanAnomaly <
                toleranceOrbitalElementConversion ) ) && ( testMeanAnomaly != 0.0 ) )
             || aRuntimeErrorOccurred )
        {
            failedMeanAnomalies.push_back( testMeanAnomaly );
            failedEccentricities.push_back( testEccentricity );
        }

        // Reset boolean.
        aRuntimeErrorOccurred = false;
    }

    // Check that no values have been written to the failedMeanAnomalies vector.  If so, this test
    // is passed. Otherwisely these values will be written away and this test will fail.
    BOOST_CHECK( failedMeanAnomalies.empty( ) );

    // If the vector is not empty, write the failed cases of this test case to a file.
    if ( !( failedMeanAnomalies.empty( ) ) )
    {
        writeErrorsToFile( failedEccentricities, failedMeanAnomalies, "Test6" );
    }
}

//! Test 7: Test functionality of specifying the initial guess.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToHyperbolicEccentricAnomaly_specificInitialGuess )
{
    // Set test value for eccentricity.
    const double testEccentricity = 1.97;

    // Set test value for mean anomaly.
    const double testHyperbolicMeanAnomaly = 0.5;

    // Set expected hyperbolic eccentric anomaly. (Similar case as in Test 1.)
    const double expectedHyperbolicEccentricAnomaly = 0.478057581067141;

    // Set the initial guess.
    const double initialGuess = 2.0 * testHyperbolicMeanAnomaly / testEccentricity - 1.8;

    // Compute eccentric anomaly.
    const double hyperbolicEccentricAnomaly = convertMeanAnomalyToHyperbolicEccentricAnomaly(
                testEccentricity, testHyperbolicMeanAnomaly, false, initialGuess );

    // Check if computed eccentric anomaly is NaN for invalid eccentricity.
    BOOST_CHECK_CLOSE_FRACTION( expectedHyperbolicEccentricAnomaly, hyperbolicEccentricAnomaly,
                                toleranceOrbitalElementConversion );
}

// End Boost test suite.
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
