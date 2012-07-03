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
 *      120606    T. Secretin       First creation of code.
 *      120608    E. Heeren         Placed tests in unit_tests namespace.
 */

#define BOOST_TEST_MAIN

#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/InputOutput/fieldValue.h"
#include "Tudat/InputOutput/fieldTransform.h"

namespace tudat
{
namespace input_output
{

// Define a derived linear field transform class
class LinearFieldTransform : public FieldTransform
{

public:
    // Constructor. Parameters correspond to the terms a and b in the equation y = a * x + b.
    LinearFieldTransform( double a, double b ) : factor( a ), term( b ) { }

    // Default destructor.
    ~LinearFieldTransform( ) { }

    // Transform function. Transformation is performed according to y = a * x + b.
    boost::shared_ptr< std::string > transform( boost::shared_ptr< std::string > input )
    {
        // Convert input string to a double.
        double number = boost::lexical_cast< double >( *input );

        // Perform transformation.
        double result = ( factor * number + term );

        // Convert result to string.
        std::string resultAsString = boost::lexical_cast< std::string > ( result );

        // Return pointer to string.
        return boost::shared_ptr< std::string >( new std::string( resultAsString ) );
    }

private:
    // a term.
    double factor;

    // b term.
    double term;

};

} // namespace input_output

namespace unit_tests
{

// Define Boost test suite.
BOOST_AUTO_TEST_SUITE( test_field_value )

//! Test that the raw data is returned.
BOOST_AUTO_TEST_CASE( fieldValue_getRaw )
{
    // Using declaration.
    using namespace tudat::input_output;

    // Create raw string data.
    std::string testAngleInDegrees = "60.0";

    // Create a field value object.
    FieldValue testFieldValue( field_types::state::inclination, testAngleInDegrees );

    // Retrieve raw data.
    boost::shared_ptr< std::string > returnedData = testFieldValue.getRaw( );

    // Verify that returned data is correct.
    BOOST_CHECK_EQUAL( *returnedData, testAngleInDegrees );

}

//! Test that the trasnformed data is returned.
BOOST_AUTO_TEST_CASE( fieldValue_get )
{
    // Using declaration.
    using namespace tudat::input_output;

    // Create raw string data.
    std::string testAngleInDegrees = "60.0";

    // Create the degrees to radians transform.
    boost::shared_ptr< FieldTransform > degreesToRadiansTransform(
                new LinearFieldTransform( tudat::mathematics::PI / 180.0, 0.0 ) );

    // Test 1: Test the get function with a transform.
    // Create a field value object with transform.
    FieldValue testFieldValue( field_types::state::inclination, testAngleInDegrees,
                               degreesToRadiansTransform );

    // Retrieve (transformed) data.
    boost::shared_ptr< std::string > returnedData = testFieldValue.get( );

    // Expected data.
    std::string expectedData = boost::lexical_cast< std::string > ( tudat::mathematics::PI / 3.0 );

    // Verify that returned data is correct.
    BOOST_CHECK_EQUAL( *returnedData, expectedData );

    // Test 2: Test the get function without transform.
    // Create a field value object.
    FieldValue testFieldValue2( field_types::state::inclination, testAngleInDegrees );

    // Retrieve (raw) data.
    boost::shared_ptr< std::string > returnedRawData = testFieldValue2.get( );

    // Verify that returned data is correct.
    BOOST_CHECK_EQUAL( *returnedRawData, testAngleInDegrees );

}

//! Test that the trasnformed data is returned.
BOOST_AUTO_TEST_CASE( fieldValue_operator )
{
    // Using declaration.
    using namespace tudat::input_output;

    // Create raw string data.
    std::string testAngleInDegrees = "60.0";

    // Create the degrees to radians transform.
    boost::shared_ptr< FieldTransform > degreesToRadiansTransform(
                new LinearFieldTransform( tudat::mathematics::PI / 180.0, 0.0 ) );

    // Test 1: Test the get function with a transform.
    // Create a field value object with transform.
    FieldValue testFieldValue( field_types::state::inclination, testAngleInDegrees,
                               degreesToRadiansTransform );

    // Retrieve (transformed) data.
    boost::shared_ptr< std::string > returnedData = testFieldValue( );

    // Expected data.
    std::string expectedData = boost::lexical_cast< std::string > ( tudat::mathematics::PI / 3.0 );

    // Verify that returned data is correct.
    BOOST_CHECK_EQUAL( *returnedData, expectedData );

    // Test 2: Test the get function without transform.
    // Create a field value object.
    FieldValue testFieldValue2( field_types::state::inclination, testAngleInDegrees );

    // Retrieve (raw) data.
    boost::shared_ptr< std::string > returnedRawData = testFieldValue2( );

    // Verify that returned data is correct.
    BOOST_CHECK_EQUAL( *returnedRawData, testAngleInDegrees );

}

// Close Boost test suite.
BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
