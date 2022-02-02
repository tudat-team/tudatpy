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
 *      Introduction to Flight, Fifth edition, Appendix A, John D. Anderson Jr., McGraw Hill, 2005.
 *      US Standard Atmosphere 1976,
 *          http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770009539_1977009539.pdf.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <algorithm>
#include <vector>
#include <utility>

#include <boost/bind/bind.hpp>
using namespace boost::placeholders;

#include <boost/test/unit_test.hpp>

#include "tudat/astro/aerodynamics/nrlmsise00Atmosphere.h"
#include "tudat/astro/aerodynamics/tabulatedAtmosphere.h"
#include "tudat/io/basicInputOutput.h"

#include "tudat/astro/aerodynamics/nrlmsise00InputFunctions.h"

#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace unit_tests
{

//! Test NRLMSISE00 atmosphere model.
BOOST_AUTO_TEST_SUITE( test_nrlmsise00_atmosphere )

using tudat::aerodynamics::NRLMSISE00Input;
using tudat::aerodynamics::NRLMSISE00Atmosphere;

using tudat::mathematical_constants::PI;

// Global variable to be changed by tests and function.
NRLMSISE00Input data;

// Define input data generic (or almost completely) for all tests.
NRLMSISE00Input gen_data(0, 172, 29000.0, 16.0, 150.0, 150.0, 4.0);
std::vector< double > gen_input = { 400.0E3, -70.0*PI/180.0, 60.0*PI/180.0, 0.0 };

//! Test function to update the NRLMSISE00 iunput data as a function of time and geometry.
NRLMSISE00Input nrlmsiseTestFunction( double altitude, double longitude,
                                      double latitude, double time,
                                      bool computeLocalSolarTime,
                                      bool invariableLower )
{
    // Functionality encountered in the original wrapper class, these
    // have been moved out of the Tudat space and now into the
    // application space. This is given here as an example of how to
    // solve these problems and to keep the code for future generations :).
    // NOTE: both these functions are switched of for testing.
    if( computeLocalSolarTime )
    {
        data.localSolarTime = data.secondOfTheDay / 3600.0 + longitude / 15.0;
    }

    if( invariableLower && altitude < 80000.0 )
    {
        data.f107  = 150.0;
        data.f107a = 150.0;
        data.apDaily  = 4.0;
        data.switches[ 9 ] = 1;
    }

    return data;
}

//! Perform NRLMSISE-00 test of get functions.
//  Check the consistency between full output and get parameter output functions.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTestFunctions )
{
    // Define tolerance for equality
    double tolerance = 1.0E-18;

    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;

    // Define variations from the standard case
    // First case is the default case

    // Get full output and extract density and temperature
    std::pair< std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ], input[ 2 ], input[ 3 ] );
    double density1 = output.first[ 5 ] * 1000.0;
    double temperature1 = output.second[ 1 ];

    // Get density and temperature from functions
    double density2 = model.getDensity( input[ 0 ], input[ 1 ], input[ 2 ], input[ 3 ] );
    double temperature2 = model.getTemperature( input[ 0 ], input[ 1 ], input[ 2 ], input[ 3 ] );

    // Check densities and temperatures.
    BOOST_CHECK_CLOSE_FRACTION( density1, density2, tolerance );
    BOOST_CHECK_CLOSE_FRACTION( temperature1, temperature2, tolerance );

}

//! Perform test of hashing function.
//  Hashing is important to speed up the model and prevent double
//  calculations, however it can be tricky too, so make sure to reset the hash when in doubt.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTestHashing )
{
    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;
    // Define variations from the standard case
    // First case is the default case

    double density1 = model.getDensity( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Increase the daily F10.7 average (from 140.0)
    data.f107 = 180.0;
    double density2 = model.getDensity( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // See that although the F10.7 has changed, no recalculation was
    // triggered. The densities remain the same!
    BOOST_CHECK_EQUAL(density1, density2);

    // In order to trigger a recalculation we can
    // (A) reset the hash
    model.resetHashKey( );
    double density3 = model.getDensity( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // And indeed observe a new density is calculated and is not equal
    // anymore to the old one.
    BOOST_CHECK_PREDICATE( std::not_equal_to< double >( ), ( density1 )( density3 ) );

    // Recalculate density1, to "reset" the internal state of the
    // model , for the next demo
    data.f107 = 150.0;
    model.resetHashKey( );
    double density1new = model.getDensity( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );
    BOOST_CHECK_EQUAL( density1new, density1 );
    data.f107 = 180.0;

    // (B) change any of the four indepedent parameters, namely
    // altitude, latitude, longitude, and/or time. Since time is not
    // used in our current NRLMSISE00InputFunction of this test, we
    // can change time without changing anything.
    double density4 = model.getDensity( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] + 42.0);

    // The densities are not equal (ergo recalculated). Furthermore
    // the density using method (B) is exactly equal to the density
    // obtained using solution (A), this shows that both are equal in
    // triggering recalculations.
    BOOST_CHECK_PREDICATE(std::not_equal_to<double>( ), (density1)(density4));
    BOOST_CHECK_EQUAL(density3, density4);
}

//! Perform NRLMSISE-00 test 1
//  Check values for 11 output parameters on against (hardcoded) values
//  obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest1 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 6.665176904952E+05, 1.138805559752E+08, 1.998210925573E+07,
            4.022763585713E+05, 3.557464994516E+03, 4.074713532757E-15,
            3.475312399717E+04, 4.095913268293E+06, 2.667273209336E+04,
            1.250539943561E+03, 1.241416130019E+03 };

    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;
    // Define variations from the standard case
    // First case is the default case

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 2
// Check values for 11 output parameters on against (hardcoded) values
// obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest2 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 3.407293223161E+06, 1.586333369569E+08, 1.391117365461E+07,
            3.262559509596E+05, 1.559618150501E+03, 5.001845729072E-15,
            4.854208463340E+04, 4.380966712899E+06, 6.956681955942E+03,
            1.166754383757E+03, 1.161710451887E+03 };

    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;
    // Define variations from the standard case
    data.dayOfTheYear   =    81;

    // We change parameters other than alt, lat, long, and/or time,
    // so it's best to reset the hash.
    model.resetHashKey( );

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 3
// Check values for 11 output parameters on against (hardcoded) values
// obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest3 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 1.123767244038E+05, 6.934130086761E+04, 4.247105217477E+01,
            1.322750141475E-01, 2.618848418232E-05, 2.756772319269E-18,
            2.016749854321E+04, 5.741255934147E+03, 2.374394151990E+04,
            1.239892111717E+03, 1.239890640133E+03 };

    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;
    // Define variations from the standard case
    data.secondOfTheDay = 75000.0;
    input[0]            =  1000.0E3;

    // We change parameters other than alt, lat, long, and/or time,
    // so it's best to reset the hash.
    model.resetHashKey( );

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 4
// Check values for 11 output parameters on against (hardcoded) values
// obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest4 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 5.411554379937E+07, 1.918893443939E+11, 6.115825598225E+12,
            1.225201051740E+12, 6.023211973085E+10, 3.584426304113E-10,
            1.059879697741E+07, 2.615736693705E+05, 2.819879355928E-42,
            1.027318464900E+03, 2.068877764036E+02 };

    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;

    // Define variations from the standard case
    input[0] = 100.0E3;

    // We change parameters other than alt, lat, long, and/or time,
    // so it's best to reset the hash.
    model.resetHashKey( );

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 5
// Check values for 11 output parameters on against (hardcoded) values
// obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest5 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 1.851122486193E+06, 1.476554837927E+08, 1.579356228264E+07,
            2.633794977312E+05, 1.588781398384E+03, 4.809630239407E-15,
            5.816166780787E+04, 5.478984479069E+06, 1.264445941761E+03,
            1.212396152121E+03, 1.208135425212E+03 };

    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;

    // Define variations from the standard case
    input[ 2 ] = 0.0;

    // We change parameters other than alt, lat, long, and/or time,
    // so it's best to reset the hash.
    model.resetHashKey( );

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 6
// Check values for 11 output parameters on against (hardcoded) values
// obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest6 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 8.673095233906E+05, 1.278861768014E+08, 1.822576627172E+07,
            2.922214190618E+05, 2.402962436424E+03, 4.355865642645E-15,
            3.686389243751E+04, 3.897275503727E+06, 2.667273209336E+04,
            1.220146417915E+03, 1.212712083212E+03 };

    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;

    // Define variations from the standard case
    input[ 1 ] = 0.0;

    // We change parameters other than alt, lat, long, and/or time,
    // so it's best to reset the hash.
    model.resetHashKey( );

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 7
// Check values for 11 output parameters on against (hardcoded) values
// obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest7 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 5.776251216023E+05, 6.979138693660E+07, 1.236813559822E+07,
            2.492867715429E+05, 1.405738674178E+03, 2.470651391663E-15,
            5.291985567067E+04, 1.069814109367E+06, 2.667273209336E+04,
            1.116385376043E+03, 1.112998568217E+03 };

    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;

    // Define variations from the standard case
    data.localSolarTime = 4.0;

    // We change parameters other than alt, lat, long, and/or time,
    // so it's best to reset the hash.
    model.resetHashKey( );

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 8
// Check values for 11 output parameters on against (hardcoded) values
// obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest8 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 3.740304105508E+05, 4.782720123611E+07, 5.240380033324E+06,
            1.759874640391E+05, 5.501648779570E+02, 1.571888739255E-15,
            8.896775722935E+04, 1.979740836233E+06, 9.121814875991E+03,
            1.031247440715E+03, 1.024848492213E+03 };
    
    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;

    // Define variations from the standard case
    data.f107a = 70.0;

    // We change parameters other than alt, lat, long, and/or time,
    // so it's best to reset the hash.
    model.resetHashKey( );

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 9
// Check values for 11 output parameters on against (hardcoded) values
// obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest9 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 6.748338766624E+05, 1.245315260444E+08, 2.369009541053E+07,
            4.911583154750E+05, 4.578781099054E+03, 4.564420245361E-15,
            3.244594775161E+04, 5.370833087086E+06, 2.667273209336E+04,
            1.306052042027E+03, 1.293374040390E+03 };

    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;

    // Define variations from the standard case
    data.f107 = 180.0;

    // We change parameters other than alt, lat, long, and/or time,
    // so it's best to reset the hash.
    model.resetHashKey( );

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 10
// Check values for 11 output parameters on against (hardcoded) values
// obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest10 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 5.528600841645E+05, 1.198041324041E+08, 3.495797764558E+07,
            9.339618355028E+05, 1.096254765493E+04, 4.974543110322E-15,
            2.686427856260E+04, 4.889974232971E+06, 2.805444837126E+04,
            1.361868020785E+03, 1.347389183730E+03 };

    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;

    // Define variations from the standard case
    data.apDaily = 40.0;

    // We change parameters other than alt, lat, long, and/or time,
    // so it's best to reset the hash.
    model.resetHashKey( );

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 11
// Check values for 11 output parameters on against (hardcoded) values
// obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest11 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 1.375487584186E+14, 0.000000000000E+00, 2.049687044291E+19,
            5.498695433719E+18, 2.451733158028E+17, 1.261065661119E-03,
            0.000000000000E+00, 0.000000000000E+00, 0.000000000000E+00,
            1.027318464900E+03, 2.814647576632E+02 };
    
    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;

    // Define variations from the standard case
    input[ 0 ] = 0.0;

    // We change parameters other than alt, lat, long, and/or time,
    // so it's best to reset the hash.
    model.resetHashKey( );

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 12
// Check values for 11 output parameters on against (hardcoded) values
// obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest12 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 4.427442587677E+13, 0.000000000000E+00, 6.597567157737E+18,
            1.769929341406E+18, 7.891679955727E+16, 4.059139375799E-04,
            0.000000000000E+00, 0.000000000000E+00, 0.000000000000E+00,
            1.027318464900E+03, 2.274179808273E+02 };

    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;

    // Define variations from the standard case
    input[0] = 10.0E3;

    // We change parameters other than alt, lat, long, and/or time,
    // so it's best to reset the hash.
    model.resetHashKey( );

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 13
// Check values for 11 output parameters on against (hardcoded) values
// obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest13 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 2.127828756207E+12, 0.000000000000E+00, 3.170790550354E+17,
            8.506279809435E+16, 3.792741116806E+15, 1.950822245176E-05,
            0.000000000000E+00, 0.000000000000E+00, 0.000000000000E+00,
            1.027318464900E+03, 2.374389145877E+02 };

    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;

    // Define variations from the standard case
    input[0] = 30.0E3;

    // We change parameters other than alt, lat, long, and/or time,
    // so it's best to reset the hash.
    model.resetHashKey( );

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 14
// Check values for 11 output parameters on against (hardcoded) values
// obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest14 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 1.412183545593E+11, 0.000000000000E+00, 2.104369643783E+16,
            5.645392443377E+15, 2.517141749411E+14, 1.294709015929E-06,
            0.000000000000E+00, 0.000000000000E+00, 0.000000000000E+00,
            1.027318464900E+03, 2.795551129541E+02 };
    
    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;

    // Define variations from the standard case
    input[0] = 50.0E3;

    // We change parameters other than alt, lat, long, and/or time,
    // so it's best to reset the hash.
    model.resetHashKey( );

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 15
// Check values for 11 output parameters on against (hardcoded) values
// obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest15 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 1.254884400273E+10, 0.000000000000E+00, 1.874532829219E+15,
            4.923050980785E+14, 2.239685413856E+13, 1.147667671512E-07,
            0.000000000000E+00, 0.000000000000E+00, 0.000000000000E+00,
            1.027318464900E+03, 2.190732313642E+02 };

    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;

    // Define variations from the standard case
    input[0] = 70.0E3;

    // We change parameters other than alt, lat, long, and/or time,
    // so it's best to reset the hash.
    model.resetHashKey( );

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 16
// Check values for 11 output parameters on against (hardcoded) values
// obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest16 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 5.196477402973E+05, 1.274494072960E+08, 4.850449869853E+07,
            1.720837982575E+06, 2.354486590544E+04, 5.881940448652E-15,
            2.500078391081E+04, 6.279209825019E+06, 2.667273209336E+04,
            1.426411662282E+03, 1.408607795553E+03 };

    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;

    // Define variations from the standard case
    data.apVector = std::vector< double >( 7, 100.0 );
    data.switches[ 9 ] = -1;

    // We change parameters other than alt, lat, long, and/or time,
    // so it's best to reset the hash.
    model.resetHashKey( );

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 17
// Check values for 11 output parameters on against (hardcoded) values
// obtained from the similar nrlmsise-test.c program.
BOOST_AUTO_TEST_CASE( testNRLMSISE00AtmosphereTest17 )
{
    // Define verification data for specific test and tolerance
    double tolerance = 1.0E-12;
    std::vector< double > verificationData = { 4.260859748794E+07, 1.241342015549E+11, 4.929561542488E+12,
            1.048406749093E+12, 4.993465083056E+10, 2.914303550309E-10,
            8.831228592572E+06, 2.252515508626E+05, 2.415245929649E-42,
            1.027318464900E+03, 1.934071062577E+02 };

    // Create the model
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    // Create local copy of input and define variations
    data = gen_data;
    std::vector< double > input = gen_input;

    // Define variations from the standard case
    input[0]         = 100.0E3;
    data.apVector    = std::vector< double >(7, 100.0);
    data.switches[9] = -1;

    // We change parameters other than alt, lat, long, and/or time,
    // so it's best to reset the hash.
    model.resetHashKey( );

    // Get full output for current test case
    std::pair<std::vector< double >, std::vector< double > > output
            = model.getFullOutput( input[ 0 ], input[ 1 ],
            input[ 2 ], input[ 3 ] );

    // Loop through output and perform tests
    for( unsigned i = 0; i < 9; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i ],
                                    output.first[ i ], tolerance );
    }
    for( unsigned i = 0; i < 2; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( verificationData[ i + 9 ],
                output.second[ i ], tolerance );
    }
}

//! Perform NRLMSISE-00 test 18
// Check computation of speed of sound
BOOST_AUTO_TEST_CASE( testSpeedOfSound )
{
    // Construct model with default properties
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    double altitude = 0.0 ;
    double longitude = 0.0 ;
    double latitude = 0.0 ;
    double time = 0.0 ;

    // Check speed of sound : (https://en.wikipedia.org/wiki/Speed_of_sound  (at 20 deg celcius) )
    // 340.9 at sealevel (anderson, fundamentals of aerodynamics)
    BOOST_CHECK_CLOSE_FRACTION( model.getSpeedOfSound(altitude,longitude,latitude,time) , 343.21 , 2E-2 ); // < 2%
    BOOST_CHECK_CLOSE_FRACTION( model.getSpeedOfSound(altitude,longitude,latitude,time) , 340.9 , 1.8E-2 ); // < 1.8%

    // US Standard Atmosphere Model 1976
    std::string atmosphereTableFile = tudat::paths::getAtmosphereTablesPath( ) + "/USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";
    tudat::aerodynamics::TabulatedAtmosphere US76model( atmosphereTableFile );

    // Test wist US76 table (NRLMSISE should provide approximately the same results)
    BOOST_CHECK_CLOSE_FRACTION( model.getSpeedOfSound(altitude,longitude,latitude,time) , US76model.getSpeedOfSound(altitude) , 2E-2 ); // < 2%

    altitude = 10.0E3 ;
    BOOST_CHECK_CLOSE_FRACTION( model.getSpeedOfSound(altitude,longitude,latitude,time) , US76model.getSpeedOfSound(altitude) , 3.5E-2 ); // < 3.5%

    altitude = 40.0E3 ;
    BOOST_CHECK_CLOSE_FRACTION( model.getSpeedOfSound(altitude,longitude,latitude,time) , US76model.getSpeedOfSound(altitude) , 6.5E-2 ); // < 6.5%

    altitude = 80.0E3 ;
    BOOST_CHECK_CLOSE_FRACTION( model.getSpeedOfSound(altitude,longitude,latitude,time) , US76model.getSpeedOfSound(altitude) , 5E-2 ); // < 5%

    // Verify calculation
    double speedOfSound = std::sqrt( 1.4 * 8.3144598 * model.getTemperature(altitude,longitude,latitude,time) / model.getMeanMolarMass(altitude,longitude,latitude,time) ) ;
    BOOST_CHECK_CLOSE_FRACTION( speedOfSound , model.getSpeedOfSound(altitude,longitude,latitude,time), 1E-15 );
}

//! Perform NRLMSISE-00 test 19
// Check computation of molar mass
BOOST_AUTO_TEST_CASE( testMolarMass )
{
    // Construct model with default properties
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    double altitude = 0.0 ;
    double longitude = 0.0 ;
    double latitude = 0.0 ;
    double time = 0.0 ;

    // Check using textbook example: Dynamics of atmospheric re-entry, Regan , 1993
    BOOST_CHECK_CLOSE_FRACTION( model.getMeanMolarMass(altitude,longitude,latitude,time), 28.96643E-3 , 3E-4 ); // < 0.03%

    altitude = 50.0E3;
    BOOST_CHECK_CLOSE_FRACTION( model.getMeanMolarMass(altitude,longitude,latitude,time), 28.96643E-3 , 3E-4 ); // < 0.03%

    altitude = 200.0E3;
    BOOST_CHECK_CLOSE_FRACTION( model.getMeanMolarMass(altitude,longitude,latitude,time), 21.0E-3 , 3E-2 ); // < 3%

    altitude = 300.0E3;
    BOOST_CHECK_CLOSE_FRACTION( model.getMeanMolarMass(altitude,longitude,latitude,time), 18.0E-3 , 3E-2 ); // < 3%
}

//! Perform NRLMSISE-00 test 20
// Check computation mean free path
BOOST_AUTO_TEST_CASE( testMeanFreePath )
{
    // Construct model with default properties
    NRLMSISE00Atmosphere model( std::bind( &nrlmsiseTestFunction, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, false, false ) );

    double altitude = 20.0E3 ;
    double longitude = 0.5 ;
    double latitude = 0.2 ;
    double time = 1.0E5 ;

    // Verify calculation
    double meanFreePath = std::sqrt(2.0) * PI * std::pow( model.getWeightedAverageCollisionDiameter(altitude,longitude,latitude,time) , 2.0)
            * model.getAverageNumberDensity(altitude,longitude,latitude,time) ;
    meanFreePath = std::pow( meanFreePath , (-1.0) );

    BOOST_CHECK_CLOSE_FRACTION( model.getMeanFreePath(altitude,longitude,latitude,time) , meanFreePath , 1.0E-15 );

    // Verify using data - Logarithmic plot: physics of the Earth's space environment, Gerd W. Prolls (page 29)
    // Test using approximate values obtained from figure
    altitude = 0.0E3 ;
    BOOST_CHECK_CLOSE_FRACTION( model.getMeanFreePath(altitude,longitude,latitude,time) , 1.0E-7 , 10.0 );

    altitude = 60.0E3 ;
    BOOST_CHECK_CLOSE_FRACTION( model.getMeanFreePath(altitude,longitude,latitude,time) , 1.0E-3 , 1.2 );

    altitude = 300.0E3 ;
    BOOST_CHECK_CLOSE_FRACTION( model.getMeanFreePath(altitude,longitude,latitude,time) , 1.0E4 , 0.9 );
}


//! Perform NRLMSISE-00 test 21 - Test input function
// Corresponds to NRLMSISE test 1 but space weather file is used to set solar activity data.
// No adjustment values in space weather file are used.
BOOST_AUTO_TEST_CASE( test_nrlmise_InputFunction_no_adjustment )
{

//    NRLMSISE00Input gen_data(2030, 172, 29000.0, 16.0, 150.0, 150.0, 4.0);
    // DD-MM-YY = 21-06-2030 , Hrs-Min-Sec = 8-3-20
    double julianDate = tudat::basic_astrodynamics::convertCalendarDateToJulianDay< double >( 2030, 6, 21, 8, 3, 20.0 );
    double time = tudat::basic_astrodynamics::convertJulianDayToSecondsSinceEpoch(
                    julianDate , tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000) ;

    // std::vector< double > gen_input = { 400.0, -70.0, 60.0, 0.0 };
    double altitude     = 400.0E3       ; // km
    double longitude    = -70.0 * PI / 180.0  ;
    double latitude     = 60.0 * PI / 180.0  ;

    // find space weather file
    std::string spaceWeatherFilePath = tudat::paths::getTudatTestDataPath( ) + "/swAtmosTestWithAdjust.txt";

    tudat::input_output::solar_activity::SolarActivityDataMap solarActivityData =
            tudat::input_output::solar_activity::readSolarActivityData(spaceWeatherFilePath) ;

    // Create atmosphere model using NRLMISE00 input function
    // Local solar time is set to 16.0 to correspond with testing data!
    std::function< tudat::aerodynamics::NRLMSISE00Input (double,double,double,double) > inputFunction =
            std::bind( &tudat::aerodynamics::nrlmsiseInputFunction, std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3, std::placeholders::_4, solarActivityData , true , 16.0 );

    // Create Pointer to NRLMSISE model
    tudat::aerodynamics::NRLMSISE00Atmosphere atmosphereModel( inputFunction );

    // Verification data
    std::vector<double> verificationData = { 6.665176904952E+05, 1.138805559752E+08, 1.998210925573E+07,
            4.022763585713E+05, 3.557464994516E+03, 4.074713532757E-15,
            3.475312399717E+04, 4.095913268293E+06, 2.667273209336E+04,
            1.250539943561E+03, 1.241416130019E+03 };

    double computedDensity = atmosphereModel.getDensity( altitude , longitude, latitude , time ) ;

    BOOST_CHECK_CLOSE_FRACTION(verificationData[5]*1000 , computedDensity , 1E-11);
}

//! Perform NRLMSISE-00 test 22 - Test input function
// Corresponds to NRLMSISE test 1 but space weather file is used to set solar activity data.
// Adjustment values (flux qualifier is 1) in space weather file are used.
BOOST_AUTO_TEST_CASE( test_nrlmise_InputFunction_with_adjustment )
{

//    NRLMSISE00Input gen_data(2030, 172, 29000.0, 16.0, 150.0, 150.0, 4.0);
    // DD-MM-YY = 21-06-2030 , Hrs-Min-Sec = 8-3-20
    double julianDate = tudat::basic_astrodynamics::convertCalendarDateToJulianDay< double >( 2030, 6, 21, 8, 3, 20.0 );
    double time = tudat::basic_astrodynamics::convertJulianDayToSecondsSinceEpoch(
                    julianDate , tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000) ;

    // std::vector< double > gen_input = { 400.0, -70.0, 60.0, 0.0 };
    double altitude     = 400.0E3       ; // km
    double longitude    = -70.0 * PI / 180.0  ;
    double latitude     = 60.0 * PI / 180.0  ;

    // find space weather file
    std::string spaceWeatherFilePath = tudat::paths::getTudatTestDataPath( ) + "/swAtmosTestWithAdjust.txt";


    tudat::input_output::solar_activity::SolarActivityDataMap solarActivityData =
            tudat::input_output::solar_activity::readSolarActivityData(spaceWeatherFilePath) ;

    // Create atmosphere model using NRLMISE00 input function
    // Local solar time is set to 16.0 to correspond with testing data!
    std::function< tudat::aerodynamics::NRLMSISE00Input (double,double,double,double) > inputFunction =
            std::bind(&tudat::aerodynamics::nrlmsiseInputFunction, std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3, std::placeholders::_4, solarActivityData , true , 16.0 );

    // Create Pointer to NRLMSISE model
    tudat::aerodynamics::NRLMSISE00Atmosphere atmosphereModel( inputFunction );

    // Verification data
    std::vector<double> verificationData = { 6.665176904952E+05, 1.138805559752E+08, 1.998210925573E+07,
            4.022763585713E+05, 3.557464994516E+03, 4.074713532757E-15,
            3.475312399717E+04, 4.095913268293E+06, 2.667273209336E+04,
            1.250539943561E+03, 1.241416130019E+03 };

    double computedDensity = atmosphereModel.getDensity( altitude , longitude, latitude , time ) ;

    BOOST_CHECK_CLOSE_FRACTION(verificationData[5]*1000 , computedDensity , 1E-11);
}

BOOST_AUTO_TEST_CASE( test_nrlmise_FullFileLoad )
{

    // DD-MM-YY = 21-06-2030 , Hrs-Min-Sec = 8-3-20
    double julianDate = tudat::basic_astrodynamics::convertCalendarDateToJulianDay< double >( 2026, 7, 1, 0, 1, 0.0 );
    double time = tudat::basic_astrodynamics::convertJulianDayToSecondsSinceEpoch(
                    julianDate , tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000) ;

    // std::vector< double > gen_input = { 400.0, -70.0, 60.0, 0.0 };
    double altitude     = 400.0E3       ; // km
    double longitude    = -70.0 * PI / 180.0  ;
    double latitude     = 60.0 * PI / 180.0  ;

    // find space weather file
    std::string spaceWeatherFilePath = tudat::paths::getTudatTestDataPath( ) + "/sw19571001.txt";
    tudat::input_output::solar_activity::SolarActivityDataMap solarActivityData =
            tudat::input_output::solar_activity::readSolarActivityData(spaceWeatherFilePath) ;

    std::function< tudat::aerodynamics::NRLMSISE00Input (double,double,double,double) > inputFunction =
            std::bind(&tudat::aerodynamics::nrlmsiseInputFunction, std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3, std::placeholders::_4, solarActivityData , true , 16.0 );

    // Create Pointer to NRLMSISE model
    tudat::aerodynamics::NRLMSISE00Atmosphere atmosphereModel( inputFunction );
    double computedDensity = atmosphereModel.getDensity( altitude , longitude, latitude , time ) ;
    std::cout<<computedDensity<<std::endl;

    julianDate = tudat::basic_astrodynamics::convertCalendarDateToJulianDay< double >( 2026, 6, 30, 23, 59, 0.0 );
    time = tudat::basic_astrodynamics::convertJulianDayToSecondsSinceEpoch(
                        julianDate , tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000) ;
    computedDensity = atmosphereModel.getDensity( altitude , longitude, latitude , time ) ;
    std::cout<<computedDensity<<std::endl;

}

}

} // namespace unit_tests

} // namespace tudat
