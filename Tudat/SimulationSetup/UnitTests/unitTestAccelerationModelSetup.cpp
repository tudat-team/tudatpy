/*    Copyright (c) 2010-2014, Delft University of Technology
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
 *      110207    B. Romgens        File created.
 *      110215    K. Kumar          Minor modifications to layout, comments
 *                                  and variable-naming.
 *      110411    K. Kumar          Added unit test for
 *                                  convertCartesianToSpherical( ) function.
 *      110701    K. Kumar          Updated failing tests with relative errors.
 *      110708    K. Kumar          Added unit tests for computeSampleMean( )
 *                                  and computeSampleVariance( ) functions.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      111111    K. Kumar          Strange error with convertCylindricalToCartesian function;
 *                                  achieved precision of results is less than machine precision,
 *                                  fixed by using slightly larger precision tolerance.
 *      120202    K. Kumar          Separated from unitTestBasicMathematics.cpp into new
 *                                  Interpolators sub-directory.
 *      120529    E.A.G. Heeren     Boostified unit test.
 *      120615    T. Secretin       Minor layout changes.
 *      120716    D. Dirkx          Updated with interpolator architecture.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/Basics/testMacros.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Mathematics/Interpolators/linearInterpolator.h"
#include "Tudat/SimulationSetup/createAccelerationModels.h"
#include "Tudat/SimulationSetup/createBodies.h"
#include "Tudat/SimulationSetup/defaultBodies.h"

namespace tudat
{
namespace unit_tests
{

using namespace simulation_setup;
using namespace basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_acceleration_model_setup )

//! Test set up of point mass gravitational accelerations, both direct and third-body.
BOOST_AUTO_TEST_CASE( test_centralGravityModelSetup )
{
    // Load Spice kernel with gravitational parameters.
    const std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc" );

    // Create bodies with gravitational parameters from Spice and JPL approximane positions
    // as ephemerides
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Mars" ] = boost::make_shared< BodySettings >( );
    bodySettings[ "Jupiter" ] = boost::make_shared< BodySettings >( );
    bodySettings[ "Sun" ] = boost::make_shared< BodySettings >( );
    bodySettings[ "Mars" ]->ephemerisSettings = boost::make_shared< ApproximatePlanetPositionSettings >(
                ephemerides::ApproximatePlanetPositionsBase::mars, 0 );
    bodySettings[ "Jupiter" ]->ephemerisSettings = boost::make_shared< ApproximatePlanetPositionSettings >(
                ephemerides::ApproximatePlanetPositionsBase::jupiter, 0 );
    bodySettings[ "Mars" ]->gravityFieldSettings =
            boost::make_shared< GravityFieldSettings >( central_spice );
    bodySettings[ "Jupiter" ]->gravityFieldSettings =
            boost::make_shared< GravityFieldSettings >( central_spice );
    bodySettings[ "Sun" ]->gravityFieldSettings =
            boost::make_shared< GravityFieldSettings >( central_spice );
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Defins state of Sun to be all zero.
    std::map< double, basic_mathematics::Vector6d > sunStateHistory;
    sunStateHistory[ -1.0E9 ] = basic_mathematics::Vector6d::Zero( );
    sunStateHistory[ 0.0 ] = basic_mathematics::Vector6d::Zero( );
    sunStateHistory[ -1.0E9 ] = basic_mathematics::Vector6d::Zero( );
    boost::shared_ptr< interpolators::LinearInterpolator< double, basic_mathematics::Vector6d > >
            sunStateInterpolaotor = boost::make_shared<
            interpolators::LinearInterpolator< double, basic_mathematics::Vector6d > >(
                sunStateHistory );
    bodyMap[ "Sun" ] ->setEphemeris( boost::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                         sunStateInterpolaotor ) );

    // Update bodies to current state (normally done by numerical integrator).
    for( NamedBodyMap::const_iterator bodyIterator = bodyMap.begin( ); bodyIterator !=
         bodyMap.end( ); bodyIterator++ )
    {
        bodyIterator->second->updateStateFromEphemeris( 1.0E7 );
    }


    // Define settings for accelerations: point  mass atraction by Jupiter and Sun on Mars
    SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Mars" ][ "Sun" ].push_back(
                boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationSettingsMap[ "Mars" ][ "Jupiter" ].push_back(
                boost::make_shared< AccelerationSettings >( central_gravity ) );

    // Define origin of integration to be barycenter.
    std::map< std::string, std::string > centralBodies;
    centralBodies[ "Mars" ] = "SSB";

    // Create accelerations
    AccelerationMap accelerationsMap = createAccelerationModelsMap(
                bodyMap, accelerationSettingsMap, centralBodies );

    // Retrieve created accelerations.
    boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
            sunAcceleration = accelerationsMap[ "Mars" ][ "Sun" ][ 0 ];
    boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
            jupiterAcceleration = accelerationsMap[ "Mars" ][ "Jupiter" ][ 0 ];

    // Create accelerations manually (point mass inertial).
    boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
            manualSunAcceleration =
            boost::make_shared< gravitation::CentralGravitationalAccelerationModel< > >(
                boost::bind( &Body::getPosition, bodyMap[ "Mars" ] ),
                spice_interface::getBodyGravitationalParameter( "Sun" ),
                boost::bind( &Body::getPosition, bodyMap[ "Sun" ] ) );
    boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
            manualJupiterAcceleration =
            boost::make_shared< gravitation::CentralGravitationalAccelerationModel< > >(
                boost::bind( &Body::getPosition, bodyMap[ "Mars" ] ),
                spice_interface::getBodyGravitationalParameter( "Jupiter" ),
                boost::bind( &Body::getPosition, bodyMap[ "Jupiter" ] ) );

    // Test equivalence of two acceleration models.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( basic_astrodynamics::updateAndGetAcceleration( sunAcceleration ) ),
                ( basic_astrodynamics::updateAndGetAcceleration( manualSunAcceleration ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( basic_astrodynamics::updateAndGetAcceleration( jupiterAcceleration ) ),
                ( basic_astrodynamics::updateAndGetAcceleration( manualJupiterAcceleration ) ),
                std::numeric_limits< double >::epsilon( ) );

    // Change central body to Sun, which will result in modified accelerations.
    centralBodies[ "Mars" ] = "Sun";

    // Recreate and retrieve accelerations.
    accelerationsMap = createAccelerationModelsMap(
                bodyMap, accelerationSettingsMap, centralBodies );
    sunAcceleration = accelerationsMap[ "Mars" ][ "Sun" ][ 0 ];
    jupiterAcceleration = accelerationsMap[ "Mars" ][ "Jupiter" ][ 0 ];

    // Manually create Sun's acceleration on Mars, which now include's Mars'gravitational parameter,
    // since the integration is done w.r.t. the Sun, not the barycenter.
    manualSunAcceleration =
            boost::make_shared< gravitation::CentralGravitationalAccelerationModel< > >(
                boost::bind( &Body::getPosition, bodyMap[ "Mars" ] ),
                spice_interface::getBodyGravitationalParameter( "Sun" ) +
                spice_interface::getBodyGravitationalParameter( "Mars" ),
                boost::bind( &Body::getPosition, bodyMap[ "Sun" ] ) );

    // Manually create Jupiter's acceleration on Mars, which now a third body acceleration,
    // with the Sun the central body.
    manualJupiterAcceleration =
            boost::make_shared< gravitation::ThirdBodyAcceleration<
            gravitation::CentralGravitationalAccelerationModel< > > >(
                boost::make_shared< gravitation::CentralGravitationalAccelerationModel< > >(
                    boost::bind( &Body::getPosition, bodyMap[ "Mars" ] ),
                    spice_interface::getBodyGravitationalParameter( "Jupiter" ),
                    boost::bind( &Body::getPosition, bodyMap[ "Jupiter" ] ) ),
                boost::make_shared< gravitation::CentralGravitationalAccelerationModel< > >(
                    boost::bind( &Body::getPosition, bodyMap[ "Sun" ] ),
                    spice_interface::getBodyGravitationalParameter( "Jupiter" ),
                    boost::bind( &Body::getPosition, bodyMap[ "Jupiter" ] ) ) );

    // Test equivalence of two acceleration models.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( basic_astrodynamics::updateAndGetAcceleration( sunAcceleration ) ),
                ( basic_astrodynamics::updateAndGetAcceleration( manualSunAcceleration ) ),
                std::numeric_limits< double >::epsilon( ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( basic_astrodynamics::updateAndGetAcceleration( jupiterAcceleration ) ),
                ( basic_astrodynamics::updateAndGetAcceleration( manualJupiterAcceleration ) ),
                std::numeric_limits< double >::epsilon( ) );
}

//! Test set up of spherical harmonic gravitational accelerations.
BOOST_AUTO_TEST_CASE( test_shGravityModelSetup )
{
    // Load Spice kernel with gravitational parameters.
    const std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc" );

    // Create body map
    NamedBodyMap bodyMap;
    bodyMap[ "Earth" ] = boost::make_shared< Body >( );
    bodyMap[ "Vehicle" ] = boost::make_shared< Body >( );

    // Set constant state for Earth and Vehicle
    basic_mathematics::Vector6d dummyEarthState =
            ( basic_mathematics::Vector6d ( ) << 1.1E11, 0.5E11, 0.01E11, 0.0
              ).finished( );
    bodyMap[ "Earth" ]->setCurrentTimeAndState( 0.0, dummyEarthState );
    bodyMap[ "Vehicle" ]->setCurrentTimeAndState(
                0.0, ( basic_mathematics::Vector6d ( ) << 7.0e6, 8.0e6, 9.0e6, 0.0, 0.0, 0.0
                       ).finished( ) + dummyEarthState );

    // Define Earth gravity field.
    double gravitationalParameter = 3.986004418e14;
    double planetaryRadius = 6378137.0;
    Eigen::MatrixXd cosineCoefficients =
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
    Eigen::MatrixXd sineCoefficients =
            ( Eigen::MatrixXd( 6, 6 ) <<
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
              0.0, 1.384413891379790e-9, -1.400273703859340e-6, 0.0, 0.0, 0.0,
              0.0, 2.482004158568720e-7, -6.190054751776180e-7, 1.414349261929410e-6, 0.0, 0.0,
              0.0, -4.735673465180860e-7, 6.624800262758290e-7, -2.009567235674520e-7,
              3.088038821491940e-7, 0.0, 0.0, -9.436980733957690e-8, -3.233531925405220e-7,
              -2.149554083060460e-7, 4.980705501023510e-8, -6.693799351801650e-7
              ).finished( );
    bodyMap[ "Earth" ]->setGravityFieldModel(
                boost::make_shared< gravitation::SphericalHarmonicsGravityField >(
                    gravitationalParameter, planetaryRadius, cosineCoefficients,
                    sineCoefficients ) );

    // Define settings for acceleration model (spherical harmonic due to Earth up to degree and
    // order 5.
    SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Vehicle" ][ "Earth" ].push_back(
                boost::make_shared< SphericalHarmonicAccelerationSettings >( 5, 5 ) );

    // Set accelerations to be calculated w.r.t. the Earth.
    std::map< std::string, std::string > centralBodies;
    centralBodies[ "Vehicle" ] = "Earth";

    // Create and retrieve acceleration.
    AccelerationMap accelerationsMap = createAccelerationModelsMap(
                bodyMap, accelerationSettingsMap, centralBodies );
    boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
            directAcceleration = accelerationsMap[ "Vehicle" ][ "Earth" ][ 0 ];

    // Manually create acceleration model.
    boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > >
            manualAcceleration =
            boost::make_shared< gravitation::SphericalHarmonicsGravitationalAccelerationModel< > >(
                boost::bind( &Body::getPosition, bodyMap[ "Vehicle" ] ),
                gravitationalParameter,
                planetaryRadius, cosineCoefficients, sineCoefficients,
                boost::bind( &Body::getPosition, bodyMap[ "Earth" ] ) );

    // Test equivalence of two acceleration models.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( basic_astrodynamics::updateAndGetAcceleration( manualAcceleration ) ),
                ( basic_astrodynamics::updateAndGetAcceleration( directAcceleration ) ),
                std::numeric_limits< double >::epsilon( ) );

    // Set (unrealistically) a gravity field model on the Vehicle, to test its
    // influence on acceleration.
    bodyMap[ "Vehicle" ]->setGravityFieldModel( boost::make_shared< gravitation::GravityFieldModel >(
                                                    0.1 * gravitationalParameter ) );

    // Recreate and retrieve acceleration.
    accelerationsMap = createAccelerationModelsMap(
                bodyMap, accelerationSettingsMap, centralBodies );
    directAcceleration = accelerationsMap[ "Vehicle" ][ "Earth" ][ 0 ];

    // Manually create acceleration.
    manualAcceleration =
            boost::make_shared< gravitation::SphericalHarmonicsGravitationalAccelerationModel< > >(
                boost::bind( &Body::getPosition, bodyMap[ "Vehicle" ] ),
                gravitationalParameter * 1.1,
                planetaryRadius, cosineCoefficients, sineCoefficients,
                boost::bind( &Body::getPosition, bodyMap[ "Earth" ] ) );

    // Test equivalence of two acceleration models.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( basic_astrodynamics::updateAndGetAcceleration( manualAcceleration ) ),
                ( basic_astrodynamics::updateAndGetAcceleration( directAcceleration ) ),
                std::numeric_limits< double >::epsilon( ) );

}

//! Test set up of spherical harmonic gravitational accelerations.
BOOST_AUTO_TEST_CASE( test_aerodynamicAccelerationModelSetup )
{
    using namespace tudat::simulation_setup;
    using namespace tudat;

    // Load Spice kernel with gravitational parameters.
    const std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp" );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc" );

    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Earth" ] = getDefaultSingleBodySettings( "Earth", 0.0, 86400.0 );
    bodySettings[ "Vehicle" ] = boost::make_shared< BodySettings >( );
    bodySettings[ "Vehicle" ]->aerodynamicCoefficientSettings =
            boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                1.0, 1.0, 1.0, Eigen::Vector3d::Zero( ), Eigen::Vector3d::UnitX( ), Eigen::Vector3d::Zero( ),
            1, 1 );
    NamedBodyMap bodyMap = createBodies( bodySettings );

    // Define settings for accelerations: point  mass atraction by Jupiter and Sun on Mars
    SelectedAccelerationMap accelerationSettingsMap;
    accelerationSettingsMap[ "Vehicle" ][ "Earth" ].push_back(
                boost::make_shared< AccelerationSettings >( aerodynamic ) );

    // Define origin of integration to be barycenter.
    std::map< std::string, std::string > centralBodies;
    centralBodies[ "Vehicle" ] = "Earth";

    // Create accelerations
    AccelerationMap accelerationsMap = createAccelerationModelsMap(
                bodyMap, accelerationSettingsMap, centralBodies );


}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat


