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

#include <boost/test/unit_test.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <fstream>

#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

using namespace orbital_element_conversions;
using namespace mathematical_constants;

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
              << "listed below. Please report a bug on the Tudat website " << std::endl
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

BOOST_AUTO_TEST_SUITE( test_mean_to_hyperbolic_eccentric_anomaly_conversion )

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
                                    1.0E-14 );
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
                                    1.0E-14 );
    }

}

//! Generalized function to test many to and fro mean to hyperbolic eccentric anomalies.
/*!
 *  Generalized function to test many to and fro mean to hyperbolic eccentric anomalies.
 *  The function allows random variations of both the mean anomlay and eccentricity. The
 *  eccentricity may be fixed, for instance to test near-parabolic orbits for many mean
 *  anomalies. Furthermore, the range of mean anomalies and eccentricities may be set
 *  as exponential, so that a random number of 4.5 generates an eccentricity of 10^(4.5), to allow
 *  testing of extreme orbits.
 *  \param caseId Name of test case, to be used for error output purposes.
 *  \param testTolerance Acceptance tolerance to be used for difference between original and
 *  reconverted mean anomaly. Tolerance is absolute if useExponentialValues is false and relative
 *  if useExponentialValues is true.
 *  \param meanAnomalyLimit Limit of mean anomaly values, range of mean anomalies is set as
 *  [-meanAnomalyLimit,meanAnomalyLimit].
 *  \param useConstantEccentricity Boolean determining if a constant eccentricity is used.
 *  \param useExponentialValues Boolean determining whether exponential values are used for
 *  random mean anomalies and eccentricities.
 *  \param minimumEccentricity Minimum value to be used for random eccentricities (only used if
 *  useConstantEccentricity is false).
 *  \param maximumEccentricity Maximum value to be used for random eccentricities (only used if
 *  useConstantEccentricity is false).
 *  \param constantEccentricity Constant value to be used for eccentricity (only used if
 *  useConstantEccentricity is true).
 *  \param numberOfSamples Number of random cases that are to be tested.
 */
template< typename ScalarType >
void testMeanToHyperbolicEccentricAnomalyConversions(
        const std::string& caseId,
        const ScalarType testTolerance,
        const ScalarType meanAnomalyLimit,
        const bool useConstantEccentricity = 0,
        const bool useExponentialValues = 0,
        const ScalarType minimumEccentricity = TUDAT_NAN,
        const ScalarType maximumEccentricity = TUDAT_NAN,
        const ScalarType constantEccentricity = TUDAT_NAN,
        const int numberOfSamples = 1E5)
{
    // Create vectors that will store the input variables of a test that resulted in an error, such
    // that the error scenario can be reproduced.
    std::vector< double > failedMeanAnomalies, failedEccentricities;

    // Boolean that will be set true if a runtime error occurred.
    bool aRuntimeErrorOccurred = false;

    // Set test value for eccentricity.
    ScalarType testEccentricity = constantEccentricity;

    // Initialize both test and reverse calculated mean anomaly and the eccentric anomaly.
    ScalarType testMeanAnomaly, reverseCalculatedMeanAnomaly, eccentricAnomaly = 0.0;

    // Instantiate random number generator.
    boost::mt19937 randomNumbergenerator( time( 0 ) );

    // Create generator for eccentricity (only used if useConstantEccentricity is false).
    boost::random::uniform_real_distribution< > eccentricityDistribution;
    if( !useConstantEccentricity )
    {
        eccentricityDistribution =
                boost::random::uniform_real_distribution< >(
                    minimumEccentricity, maximumEccentricity );
    }


    boost::variate_generator< boost::mt19937&, boost::random::uniform_real_distribution < > >
            eccentricityGenerator(
                randomNumbergenerator, eccentricityDistribution );

    // Create generator for mean anomaly.
    boost::random::uniform_real_distribution< ScalarType > meanAnomalyDistibution(
                -meanAnomalyLimit, meanAnomalyLimit );
    boost::variate_generator<
            boost::mt19937&, boost::random::uniform_real_distribution < ScalarType > >
            generateMeanAnomaly( randomNumbergenerator, meanAnomalyDistibution );

    // Perform the conversion for the specified number of samples and test whether the values that
    // are subsequently converted back match the initial values.
    for ( int counter = 0; counter < numberOfSamples; counter++ )
    {
        // Set random value in test mean anomaly.
        testMeanAnomaly = generateMeanAnomaly( );

        if( useExponentialValues )
        {
            testMeanAnomaly = testMeanAnomaly *
                    std::pow( getFloatingInteger< ScalarType >( 10 ), generateMeanAnomaly( ) );
        }
        // If eccentricity is to be varied, generate random value
        if( !useConstantEccentricity )
        {
            testEccentricity = eccentricityGenerator( );

            if( useExponentialValues )
            {
                testEccentricity = getFloatingInteger< ScalarType >( 1 ) +
                        std::pow( getFloatingInteger< ScalarType >( 10 ), testEccentricity );
            }
        }

        // If the Rootfinder does not converge, it will produce a runtime error. In order to make
        // sure that these values that led to the error will not be lost, they will be stored in
        // the failed input data vectors. To do so, a try-catch sequence is used.
        try
        {
            // Compute eccentric anomaly.
            eccentricAnomaly = convertMeanAnomalyToHyperbolicEccentricAnomaly< ScalarType>(
                        testEccentricity, testMeanAnomaly );
        }
        catch( std::runtime_error )
        {
            // Store the fact that a runtime error occurred, such that the values will be stored.
            aRuntimeErrorOccurred = true;
        }

        // Calculate the mean anomaly from this eccentric anomaly.
        reverseCalculatedMeanAnomaly = convertHyperbolicEccentricAnomalyToMeanAnomaly< ScalarType>(
                    eccentricAnomaly, testEccentricity );

        // Test whether the computed mean anomaly is equal to the mean anomaly from the input and
        // that no runtime errors occurred. If an error was found, store the values leading to this
        // error in a vector for later use. '!' operator is there to ensure that a NaN value will
        // result in the values being written away. It is also checked that the mean anomaly is
        // not equal to 0.0, because that would result in falsely writing an error.
        if( !useExponentialValues )
        {
            if ( ( ( !( std::abs( testMeanAnomaly - reverseCalculatedMeanAnomaly ) <
                        testTolerance ) )
                   && !( testMeanAnomaly == getFloatingInteger< ScalarType >( 0 ) ||
                         reverseCalculatedMeanAnomaly == getFloatingInteger< ScalarType >( 0 ) ) )
                 && !aRuntimeErrorOccurred )
            {
                failedMeanAnomalies.push_back( testMeanAnomaly );
                failedEccentricities.push_back( testEccentricity );
            }
        }
        else
        {
            if ( ( ( !( std::abs( testMeanAnomaly - reverseCalculatedMeanAnomaly ) /
                        testMeanAnomaly < testTolerance ) )
                   && !( testMeanAnomaly == getFloatingInteger< ScalarType >( 0 ) ||
                         reverseCalculatedMeanAnomaly == getFloatingInteger< ScalarType >( 0 ) ) )
                 && !aRuntimeErrorOccurred )
            {
                failedMeanAnomalies.push_back( testMeanAnomaly );
                failedEccentricities.push_back( testEccentricity );
            }
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
        writeErrorsToFile( failedEccentricities, failedMeanAnomalies, caseId );
    }
}

//! Test 4: Test large number of anomalies and eccentricities
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToEccentricAnomaly_nearParabolic_random_double )
{
    long double ratioOfPrecision = std::numeric_limits< long double >::epsilon( ) /
            std::numeric_limits< double >::epsilon( );

    // Test random conversions for near-parabolic orbits
    testMeanToHyperbolicEccentricAnomalyConversions< double >(
                "DoubleParabolic", 1.0E-13, 20.0, 1, 0, TUDAT_NAN, TUDAT_NAN, 1.0 + 1.0e-15 );
    testMeanToHyperbolicEccentricAnomalyConversions< long double >(
                "LongDoubleParabolic", 1.0E-13L * ratioOfPrecision, 20.0L, 1, 0,
                TUDAT_NAN, TUDAT_NAN, 1.0L + 1.0e-15L * ratioOfPrecision );

    // Test random conversions for random orbits
    testMeanToHyperbolicEccentricAnomalyConversions< double >(
                "DoubleTypical", 1.0E-13, 20.0, 0, 0, 1.0 + 1.0e-15, 10.0, TUDAT_NAN );
    testMeanToHyperbolicEccentricAnomalyConversions< long double >(
                "LongDoubleTypical", 1.0E-13L, 20.0L, 0, 0, 1.0L + 1.0e-15L * ratioOfPrecision,
                10.0L, TUDAT_NAN );

    // Test random conversions for extreme values of eccentricity and mean anomaly.
    testMeanToHyperbolicEccentricAnomalyConversions< double >(
                "DoubleHighlyEccentric", 1.0E-14, 12.0, 0, 1, 0.0, 15.0, TUDAT_NAN );
    testMeanToHyperbolicEccentricAnomalyConversions< long double >(
                "LongDoubleHighlyEccentric", 1.0E-14L, 12.0L, 0, 1, 0.0L, 15.0L, TUDAT_NAN );

}

//! Test 5: Test functionality of specifying the initial guess.
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
                                1.0E-14 );
}

// End Boost test suite.
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
