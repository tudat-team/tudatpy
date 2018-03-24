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
 *      GTOP, http://www.esa.int/gsp/ACT/doc/INF/Code/globopt/GTOPtoolbox.rar
 *      ODTBX tolbox: http://sourceforge.net/projects/odtbx/
 *
 *    Notes
 *      Note that for near-parabolic cases, the tolerance used for the to-and-fro conversions
 *      (Tests 5 & 6) is several order of magnitudes higher than used for the regular cases. This
 *      should be investigated further in the future to fully characterize the nature of the
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


using namespace mathematical_constants;
using namespace unit_conversions;
using namespace orbital_element_conversions;

BOOST_AUTO_TEST_SUITE( test_mean_to_eccentric_anomaly_conversion )

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
            "ErrorReportConversionMeanToEccentricAnomaly" +
            testName + "RunAt" + boost::posix_time::to_iso_string( now )
            + ".txt";

    // Make a stream to a file.
    std::ofstream errorFile( outputFileName.c_str( ) );

    // Write an introduction in the file explaining what happened. 70 lines long.
    errorFile << "This error report was generated because the unit test for the" << std::endl
              << "conversion of mean to eccentric anomaly has failed in one of the" << std::endl
              << "random tests. To ensure the data for which it failed is not lost," << std::endl
              << "the corresponding input variables for these cases are listed below." << std::endl
              << "Please report a bug on the Tudat website (tudat.tudelft.nl), with" << std::endl
              << "these values, so that someone will look into it and the code can be" << std::endl
              << "improved." << std::endl << std::endl
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

//! Test 1: Test conversion for circular orbits.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToEccentricAnomaly_circular )
{
    // Set test value for eccentricity.
    const double testEccentricity = 0.0;

    // Set test value for mean anomaly.
    const double testMeanAnomaly = 1.0472;

    // Set reference value for eccentric anomaly;
    const double referenceEccentricAnomaly = 1.0472;

    // Compute eccentric anomaly.
    const double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly( testEccentricity,
                                                                          testMeanAnomaly );

    // Check if computed eccentric anomaly is less than error tolerance.
    BOOST_CHECK_CLOSE_FRACTION( eccentricAnomaly, referenceEccentricAnomaly,
                                1.0E-13 );
}

//! Test 2: Test conversion for valid range of eccentricities.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToEccentricAnomaly_range )
{
    // Set test values for eccentricity.
    const double arrayOfTestEccentricities [4] = { 0.01671, 0.43582, 0.78514, 0.91525 };

    // Set test values for mean anomaly.
    const double arrayOfTestMeanAnomalies [4] = {
        convertDegreesToRadians( 60.0 ),
        convertDegreesToRadians( 90.0 ),
        convertDegreesToRadians( 120.0 ),
        convertDegreesToRadians( 220.0 ) };

    // Set reference values for eccentric anomaly. These were obtained using GTOP.
    const double arrayOfReferenceEccentricAnomalies [4] = { 1.06178920406832,
                                                            1.97200731113253,
                                                            2.5392410896466,
                                                            3.51006218528448 };

    // Loop over sets of data.
    for ( int i = 0; i < 4; i++ )
    {
        // Compute eccentric anomaly.
        const double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly(
                    arrayOfTestEccentricities[ i ], arrayOfTestMeanAnomalies[ i ] );

        // Check if computed eccentric anomaly is less than error tolerance.
        BOOST_CHECK_CLOSE_FRACTION( eccentricAnomaly, arrayOfReferenceEccentricAnomalies[ i ],
                                    1.0E-13 );
    }
}

//! Test 3: Test conversion for negative eccentricities.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToEccentricAnomaly_negative )
{
    // Set test value for eccentricity.
    const double testEccentricity = -0.5;

    // Set test value for mean anomaly.
    const double testMeanAnomaly = 1.0472;

    // Check if a runtime error is thrown if the anomaly is converted for this eccentricity.
    BOOST_CHECK_THROW( convertMeanAnomalyToEccentricAnomaly( testEccentricity, testMeanAnomaly ),
                       std::runtime_error );
}

//! Test 4: Test conversion for eccentricity above 1.0.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToEccentricAnomaly_TooHigh )
{
    // Set test value for eccentricity.
    const double testEccentricity = 2.0;

    // Set test value for mean anomaly.
    const double testMeanAnomaly = 1.0472;

    // Check if a runtime error is thrown if the anomaly is converted for this eccentricity.
    BOOST_CHECK_THROW( convertMeanAnomalyToEccentricAnomaly( testEccentricity, testMeanAnomaly ),
                       std::runtime_error );
}

//! Test 5: Test conversion for near-parabolic orbits.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToEccentricAnomaly_nearParabolic )
{
    // Set test value for eccentricity, which is just below  1.0.
    const double testEccentricity = 1.0 - 1.0e-15;

    // Set test values for mean anomaly, to verify that a wide range can be handled. The set was
    // selected such that it contains all possible limit cases, and a wide enough variety of random
    // values.
    const double arrayOfTestMeanAnomalies [ 17 ] = { 0.0,
                                                     1.0e-10,
                                                     0.5,
                                                     PI / 2 - 1.0e-10,
                                                     PI / 2,
                                                     PI / 2 + 1.0e-10,
                                                     2.5,
                                                     PI - 1.0e-10,
                                                     PI,
                                                     PI + 1.0e-10,
                                                     4.0,
                                                     3.0 / 2.0 * PI - 1.0e-10,
                                                     3.0 / 2.0 * PI,
                                                     3.0 / 2.0 * PI + 1.0e-10,
                                                     5.5,
                                                     2.0 * PI - 1.0e-10,
                                                     2.0 * PI };

    // Set the expected eccentric anomalies corresponding to the corresponding test mean anomaly
    // array. These values were obtained using the convertMeanAnomalyToEccentricAnomaly method,
    // and subsequently verified using the eccentric -> mean anomaly method, by running GTOP code
    // with the same values and by running ODTBX code with the same values. The iteration scheme
    // converges very rapidly except for values close to zero or 2*n*pi. Still Newton Raphson
    // converges in 20-40 iterations for these cases.
    const double arrayOfExpectedEccentricAnomalies [ 17 ] = { 0.0,
                                                              0.000843432672832182,
                                                              1.49730038909589,
                                                              2.30988145995031,
                                                              2.30988146001006,
                                                              2.30988146006981,
                                                              2.81798706288006,
                                                              3.14159265353979,
                                                              3.14159265358979,
                                                              3.14159265363979,
                                                              3.57764001198758,
                                                              3.97330384710978,
                                                              3.97330384722928,
                                                              3.97330384722972,
                                                              4.51869928040234,
                                                              6.28234187379524,
                                                              0.0 };

    // Test the values which are supposed to be 0.0.
    double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly(
                testEccentricity, arrayOfTestMeanAnomalies[ 0 ] );
    BOOST_CHECK_SMALL( eccentricAnomaly, 1.0E-9 );
    eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly(
                testEccentricity, arrayOfTestMeanAnomalies[ 16 ] );
    BOOST_CHECK_SMALL( eccentricAnomaly, 1.0E-9 );

    // Test the values that are supposed to be equal to certain other values.
    for ( int counter = 1; counter < 16; counter++ )
    {
        // Compute eccentric anomaly.
        eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly(
                    testEccentricity, arrayOfTestMeanAnomalies[ counter ] );

        // Check whether the expected value is set to the eccentric anomaly.
        BOOST_CHECK_CLOSE_FRACTION( eccentricAnomaly, arrayOfExpectedEccentricAnomalies[ counter ],
                                    1.0E-9 );
    }
}

//! Generalized function to test many to and fro mean to eccentric anomalies.
/*!
 *  Generalized function to test many to and fro mean to eccentric anomalies.
 *  The function allows random variations of both the mean anomlay and eccentricity. The
 *  eccentricity may be fixed, for instance to test near-parabolic orbits for many mean
 *  anomalies.
 *  \param caseId Name of test case, to be used for error output purposes.
 *  \param testTolerance Absolute acceptance tolerance to be used for difference between original
 *  and reconverted mean anomaly.
 *  \param useConstantEccentricity Boolean determining if a constant eccentricity is used.
 *  \param minimumEccentricity Minimum value to be used for random eccentricities (only used if
 *  useConstantEccentricity is false).
 *  \param maximumEccentricity Maximum value to be used for random eccentricities (only used if
 *  useConstantEccentricity is false).
 *  \param constantEccentricity Constant value to be used for eccentricity (only used if
 *  useConstantEccentricity is true).
 *  \param numberOfSamples Number of random cases that are to be tested.
 */
template< typename ScalarType >
void testMeanToEccentricAnomalyConversions(
        const std::string& caseId,
        const ScalarType testTolerance,
        const bool useConstantEccentricity = 0,
        const ScalarType minimumEccentricity = TUDAT_NAN,
        const ScalarType maximumEccentricity = TUDAT_NAN,
        const ScalarType constantEccentricity = TUDAT_NAN,
        const int numberOfSamples = 1E5 )
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
    boost::random::uniform_real_distribution< ScalarType > distribution0To2Pi(
                getFloatingInteger< ScalarType >( 0 ), getFloatingInteger< ScalarType >( 2 ) *
                getPi< ScalarType >( ) );
    boost::variate_generator<
            boost::mt19937&, boost::random::uniform_real_distribution < ScalarType > >
            generateRandomNumber0To2Pi( randomNumbergenerator, distribution0To2Pi );

    // Perform the conversion for the specified number of samples and test whether the values that
    // are subsequently converted back match the initial values.
    for ( int counter = 0; counter < numberOfSamples; counter++ )
    {
        // Set random value in test mean anomaly.
        testMeanAnomaly = generateRandomNumber0To2Pi( );

        // If eccentricity is to be varied, generate random value
        if( !useConstantEccentricity )
        {
            testEccentricity = eccentricityGenerator( );
        }

        // If the Rootfinder does not converge, it will produce a runtime error. In order to make
        // sure that these values that led to the error will not be lost, they will be stored in
        // the failed input data vectors. To do so, a try-catch sequence is used.
        try
        {
            // Compute eccentric anomaly.
            eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly< ScalarType>(
                        testEccentricity, testMeanAnomaly );
        }
        catch( std::runtime_error )
        {
            // Store the fact that a runtime error occurred, such that the values will be stored.
            aRuntimeErrorOccurred = true;
        }

        // Calculate the mean anomaly from this eccentric anomaly.
        reverseCalculatedMeanAnomaly = convertEccentricAnomalyToMeanAnomaly< ScalarType>(
                    eccentricAnomaly, testEccentricity );

        // Test whether the computed mean anomaly is equal to the mean anomaly from the input and
        // that no runtime errors occurred. If an error was found, store the values leading to this
        // error in a vector for later use. '!' operator is there to ensure that a NaN value will
        // result in the values being written away. It is also checked that the mean anomaly is
        // not equal to 0.0, because that would result in falsely writing an error.
        if ( ( ( !( std::abs( testMeanAnomaly - reverseCalculatedMeanAnomaly ) < testTolerance ) )
               && !( testMeanAnomaly == getFloatingInteger< ScalarType >( 0 ) ||
                  reverseCalculatedMeanAnomaly == getFloatingInteger< ScalarType >( 0 ) ) ) &&
             !aRuntimeErrorOccurred )
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
        writeErrorsToFile( failedEccentricities, failedMeanAnomalies, caseId );
    }
}

//! Test 6: Test large number of anomalies for near-parabolic orbits (double).
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToEccentricAnomaly_nearParabolic_random_double )
{
    long double ratioOfPrecision = std::numeric_limits< long double >::epsilon( ) /
            std::numeric_limits< double >::epsilon( );

    // Test random conversions for near-parabolic orbits
    testMeanToEccentricAnomalyConversions< double >(
                "DoubleParabolic", 1.0E-13, 1, TUDAT_NAN, TUDAT_NAN, 1.0 - 1.0e-15 );
    testMeanToEccentricAnomalyConversions< long double >(
                "LongDoubleParabolic", 1.0E-13L * ratioOfPrecision, 1, TUDAT_NAN, TUDAT_NAN,
                1.0L - 1.0e-15L * ratioOfPrecision );

    // Test random conversions for random orbits
    testMeanToEccentricAnomalyConversions< double >(
                "DoubleEccentric", 1.0E-15, 0, 0.0, 1.0 - 1.0e-11 );
    testMeanToEccentricAnomalyConversions< long double >(
                "LongDoubleEccentric", 1.0E-15L * ratioOfPrecision, 0, 0.0,
                1.0L - 1.0e-11L * ratioOfPrecision );
}


//! Test 7: Test functionality of specifying the initial guess. Case contained in Test 2 also.
BOOST_AUTO_TEST_CASE( test_convertMeanAnomalyToEccentricAnomaly_specificInitialGuess )
{
    // Set test value for eccentricity.
    const double testEccentricity = 0.78514;

    // Set test value for mean anomaly.
    const double testMeanAnomaly = convertDegreesToRadians( 120.0 );

    // Set reference value for eccentric anomaly;
    const double referenceEccentricAnomaly = 2.5392410896466027;

    // Specify initial guess.
    const bool useDefaultGuess = false;
    const double initialGuess = PI;

    // Compute eccentric anomaly.
    const double eccentricAnomaly = convertMeanAnomalyToEccentricAnomaly(
                testEccentricity, testMeanAnomaly, useDefaultGuess, initialGuess );

    // Check if computed eccentric anomaly is less than error tolerance.
    BOOST_CHECK_CLOSE( eccentricAnomaly,
                       referenceEccentricAnomaly,
                       1.0E-13 );
}

// End Boost test suite.
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
