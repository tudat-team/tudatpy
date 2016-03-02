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
 *      101215    K. Kumar          File created.
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
 *      15XXXX    D. Dirkx          Changed some stuff.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"
#include "Tudat/Basics/testMacros.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"

#include <cmath>
#include <limits>

#include <Eigen/Core>

namespace tudat
{
namespace unit_tests
{

//! Test the functionality of the spherical harmonics gravity field class.
BOOST_AUTO_TEST_SUITE( test_spherical_harmonics_gravity_field )

//! Test setting and getting gravitational parameter.
BOOST_AUTO_TEST_CASE( testSettingAndGettingGravitationalParameter )
{
    // Set gravitational parameter of myPlanet.
    const double gravitationalParameterOfMyPlanet = 22032.00;

    // Create gravity field for myPlanet.
    using gravitation::SphericalHarmonicsGravityField;
    SphericalHarmonicsGravityField myPlanetGravityField( gravitationalParameterOfMyPlanet, 1.0 );

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
    boost::shared_ptr< gravitation::GravityFieldModel > predefinedEarthCentralGravityField =
            gravitation::getPredefinedCentralGravityField( gravitation::earth );

    const double computedResultForTest2 =
            predefinedEarthCentralGravityField->getGravitationalParameter( );

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

    // Set gravitational parameter of myPlanet.
    const double gravitationalParameterOfMyPlanet = 22032.00;

    // Create gravity field for myPlanet.
    gravitation::SphericalHarmonicsGravityField myPlanetGravityField(
                gravitationalParameterOfMyPlanet, 1.0 );

    // Get the potential of myPlanet
    double computedResultForTest3 = myPlanetGravityField.getGravitationalPotential(
                cartesianPosition );

    // Generate expected result for the potential of myPlanet
    double expectedResultForTest3 = gravitationalParameterOfMyPlanet / cartesianPosition.norm( );

    // Get the gradient of the potential of myPlanet
    Eigen::Vector3d computedResultForTest4 = myPlanetGravityField.getGradientOfPotential(
                cartesianPosition );

    // Determine the expected result for the gradient of the potential of myPlanet
    const Eigen::VectorXd expectedResultForTest4 = -gravitationalParameterOfMyPlanet
            / std::pow( cartesianPosition.norm( ), 3.0 ) * cartesianPosition;

    // Check that compute values match expected values.
    BOOST_CHECK_EQUAL( expectedResultForTest3, computedResultForTest3 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedResultForTest4, computedResultForTest4,
                                       std::numeric_limits< double >::epsilon( ) );

    // Create central gravity field for myPlanet.
    gravitation::GravityFieldModel myPlanetCentralGravityField(
                gravitationalParameterOfMyPlanet );

    // Get the potential of myPlanet
    computedResultForTest3 = myPlanetGravityField.getGravitationalPotential(
                cartesianPosition );

    // Get the gradient of the potential of myPlanet
    computedResultForTest4 = myPlanetGravityField.getGradientOfPotential(
                cartesianPosition );

    // Check that compute values match expected values.
    BOOST_CHECK_EQUAL( expectedResultForTest3, computedResultForTest3 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedResultForTest4, computedResultForTest4,
                                       std::numeric_limits< double >::epsilon( ) );
}

// Check the sum of all harmonics terms up to degree = 5 and order = 5.
BOOST_AUTO_TEST_CASE( test_SphericalHarmonicsGravitationalFieldFromFiniteDifferenceAcceleration1 )
{
    // Define gravitational parameter of Earth [m^3 s^-2]. The value is obtained from the Earth
    // Gravitational Model 2008 as described by Mathworks [2012].
    const double gravitationalParameter = 3.986004418e14;

    // Define radius of Earth [m]. The value is obtained from the Earth Gravitational Model 2008 as
    // described by Mathworks [2012].
    const double planetaryRadius = 6378137.0;

    // Define geodesy-normalized coefficients up to degree 5 and order 5. The values are obtained
    // from the Earth Gravitational Model 2008 as described by Mathworks [2012].
    const Eigen::MatrixXd cosineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              -4.841651437908150e-4, -2.066155090741760e-10, 2.439383573283130e-6, 0.0, 0.0, 0.0,
              9.571612070934730e-7, 2.030462010478640e-6, 9.047878948095281e-7,
              7.213217571215680e-7, 0.0, 0.0, 5.399658666389910e-7, -5.361573893888670e-7,
              3.505016239626490e-7, 9.908567666723210e-7, -1.885196330230330e-7, 0.0,
              6.867029137366810e-8, -6.292119230425290e-8, 6.520780431761640e-7,
              -4.518471523288430e-7, -2.953287611756290e-7, 1.748117954960020e-7
              ).finished( );

    const Eigen::MatrixXd sineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 1.384413891379790e-9, -1.400273703859340e-6, 0.0, 0.0, 0.0,
              0.0, 2.482004158568720e-7, -6.190054751776180e-7, 1.414349261929410e-6, 0.0, 0.0,
              0.0, -4.735673465180860e-7, 6.624800262758290e-7, -2.009567235674520e-7,
              3.088038821491940e-7, 0.0, 0.0, -9.436980733957690e-8, -3.233531925405220e-7,
              -2.149554083060460e-7, 4.980705501023510e-8, -6.693799351801650e-7
              ).finished( );

    // Define spherical harmonic gravity field
    gravitation::SphericalHarmonicsGravityField myPlanetGravityField(
                gravitationalParameter, planetaryRadius, cosineCoefficients, sineCoefficients );

    // Define arbitrary Cartesian position [m].
    const Eigen::Vector3d nominalPosition( 7.0e6, 8.0e6, 9.0e6 );

    // Declare variables for finite difference acceleration calculation
    double positionPerturbation = 1.0;
    Eigen::Vector3d finiteDifferenceAcceleration = Eigen::Vector3d::Zero( );
    Eigen::Vector3d currentPerturbedPosition = nominalPosition;

    // Calculate acceleration as finite spatial difference of potential.
    for( unsigned int i = 0; i < 3; i++ )
    {
        // Calculate potential at up-perturbed position
        currentPerturbedPosition( i ) += positionPerturbation;
        finiteDifferenceAcceleration( i ) = myPlanetGravityField.getGravitationalPotential(
                    currentPerturbedPosition );
        currentPerturbedPosition = nominalPosition;

        // Subtract potential at down-perturbed position
        currentPerturbedPosition( i ) -= positionPerturbation;
        finiteDifferenceAcceleration( i ) -= myPlanetGravityField.getGravitationalPotential(
                    currentPerturbedPosition );
        currentPerturbedPosition = nominalPosition;

        // Finalize finite difference calculation of entry i (x=0, y=1, z=2).
        finiteDifferenceAcceleration( i ) = finiteDifferenceAcceleration( i ) /
                ( 2.0 * positionPerturbation );
    }

    // Define expected acceleration according to the MATLAB function 'gravitysphericalharmonic'
    // described by Mathworks [2012] [m s^-2]
    // (see test_SphericalHarmonicsGravitationalAcceleration_Demo4).
    const Eigen::Vector3d expectedAcceleration(
                -1.032215878106932, -1.179683946769393, -1.328040277155269 );

    // Check if expected result matches computed result.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                expectedAcceleration, finiteDifferenceAcceleration, 5.0e-8 );

    // Calculate potential gradient analytically and compare to expected result.
    Eigen::Vector3d calculatedAcceleration = myPlanetGravityField.getGradientOfPotential(
                nominalPosition );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedAcceleration, calculatedAcceleration, 1.0E-15 );

}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace tudat
} // namespace unit_tests
