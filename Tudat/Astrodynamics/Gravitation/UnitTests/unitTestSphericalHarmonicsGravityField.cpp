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
 *      101215    K. Kumar          Creation of code.
 *      101216    K. Kumar          Updated to include test of getPotential( ).
 *      101230    K. Kumar          Updated to include test of getGradient( ) and getLaplacian( ).
 *      110104    J. Melman         Some minor comment and layout changes.
 *      110106    K. Kumar          Updated test using predefined gravity field and added machine
 *                                  precision variable.
 *      110107    K. Kumar          Updated call to predefined gravity field. Updated unit test to
 *                                  new protocol, added namespace, new filename, new file location.
 *      110113    K. Kumar          Added cerr statements.
 *      110115    J. Melman         Changed the error messages.
 *      110128    K. Kumar          Updated code to work with pointers.
 *      110202    K. Kumar          Updated code to work with State.
 *      110204    K. Kumar          Removed "vector" from naming.
 *      110310    K. Kumar          Changed naming from Laplacian to gradient tensor.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120521    F. Belien         Boostified unit test.
 *
 *    References
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Basics/testMacros.h>

#include "Tudat/Astrodynamics/Gravitation/centralGravityField.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"

namespace tudat
{
namespace unit_tests
{

//! Test the functionality of the spherical harmonics gravity field class.
BOOST_AUTO_TEST_SUITE( test_spherical_harmonics_gravity_field )

//! Test setting and getting gravitational parameter.
BOOST_AUTO_TEST_CASE( testSettingAndGettingGravitationalParameter )
{
    // Create gravity field for myPlanet.
    using astrodynamics::gravitation::SphericalHarmonicsGravityField;
    SphericalHarmonicsGravityField myPlanetGravityField;

    // Set gravitational parameter of myPlanet.
    const double gravitationalParameterOfMyPlanet = 22032.00;
    myPlanetGravityField.setGravitationalParameter( gravitationalParameterOfMyPlanet );

    // Set origin of gravity field of myPlanet with respect to geometric center.
    myPlanetGravityField.setOrigin( Eigen::Vector3d::Zero( ) );

    // Get the value for the gravitational parameter
    const double computedResultForTest1 = myPlanetGravityField.getGravitationalParameter( );

    // Determine the expected result.
    const double expectedResultForTest1 = gravitationalParameterOfMyPlanet;

    // Check that expected and computed values are equal.
    BOOST_CHECK_EQUAL( expectedResultForTest1, computedResultForTest1 );
}

//! Test getting gravitational parameter for predefined Earth central body gravity field.
BOOST_AUTO_TEST_CASE( testGetGravitationalParameterForPredefinedEarth )
{
    // Create predefined Earth central gravity field.
    using astrodynamics::gravitation::CentralGravityField;
    CentralGravityField predefinedEarthCentralGravityField( CentralGravityField::earth );

    const double computedResultForTest2 =
            predefinedEarthCentralGravityField.getGravitationalParameter( );

    // Determine the expected result
    const double expectedResultForTest2 = 3.9859383624e14;

    // Check that expected and computed values are equal.
    BOOST_CHECK_EQUAL ( expectedResultForTest2, computedResultForTest2 );
}

//! Test getting potential, the gradient of the potential and the gradient tensor of a given state.
BOOST_AUTO_TEST_CASE( testGetPotential )
{
    // Set position with respect to geometric center.
    const Eigen::Vector3d cartesianPosition( 5.0e6, 3.0e6, 1.0e6 );

    // Create gravity field for myPlanet.
    using astrodynamics::gravitation::SphericalHarmonicsGravityField;
    SphericalHarmonicsGravityField myPlanetGravityField;

    // Set gravitational parameter of myPlanet.
    const double gravitationalParameterOfMyPlanet = 22032.00;
    myPlanetGravityField.setGravitationalParameter( gravitationalParameterOfMyPlanet );

    // Get the potential of myPlanet
    const double computedResultForTest3 = myPlanetGravityField.getPotential( cartesianPosition );

    // Generate expected result for the potential of myPlanet
    double expectedResultForTest3 = gravitationalParameterOfMyPlanet / cartesianPosition.norm( );

    // Get the gradient of the potential of myPlanet
    const Eigen::Vector3d computedResultForTest4 = myPlanetGravityField.getGradientOfPotential(
                cartesianPosition );

    // Determine the expected result for the gradient of the potential of myPlanet
    const Eigen::VectorXd expectedResultForTest4 = -gravitationalParameterOfMyPlanet
            / std::pow( cartesianPosition.norm( ), 3.0 ) * cartesianPosition;

    // Get the gradient tensor of the potential of myPlanet
    const Eigen::Matrix3d computedResultForTest5
            = myPlanetGravityField.getGradientTensorOfPotential( cartesianPosition );

    // Set identity matrix.
    Eigen::Matrix3d identityMatrix = Eigen::Matrix3d::Identity( 3, 3 );

    // Determine the expected result.
    Eigen::Matrix3d expectedResultForTest5 = gravitationalParameterOfMyPlanet
            / std::pow( cartesianPosition.norm( ), 5.0 )
            * ( ( 3.0 * cartesianPosition * cartesianPosition.transpose( ) )
                - ( cartesianPosition.squaredNorm( ) * identityMatrix ) );

    // Check that compute values match expected values.
    BOOST_CHECK_EQUAL( expectedResultForTest3, computedResultForTest3 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedResultForTest4, computedResultForTest4,
                                       std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedResultForTest5, computedResultForTest5,
                                       std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace tudat
} // namespace unit_tests
