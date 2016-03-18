
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
#include <boost/assign/list_of.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"
#include "Tudat/Astrodynamics/Propagators/dynamicsSimulator.h"
#include "Tudat/Basics/testMacros.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Mathematics/Interpolators/linearInterpolator.h"
#include "Tudat/SimulationSetup/createAccelerationModels.h"
#include "Tudat/SimulationSetup/createBodies.h"
#include "Tudat/SimulationSetup/defaultBodies.h"
#include "Tudat/Astrodynamics/Aerodynamics/UnitTests/testApolloCapsuleCoefficients.h"

namespace tudat
{
namespace unit_tests
{

using namespace simulation_setup;
using namespace basic_astrodynamics;
using namespace spice_interface;
using namespace input_output;
using namespace reference_frames;
using namespace propagators;

BOOST_AUTO_TEST_SUITE( test_environment_updater )

//! Test set up of point mass gravitational accelerations, both direct and third-body.
BOOST_AUTO_TEST_CASE( test_centralGravityEnvironmentUpdate )
{

    double initialTime = 86400.0;

    // Load Spice kernels
    const std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp" );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc" );

    // Get settings for celestial bodies
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings;
    bodySettings[ "Earth" ] = getDefaultSingleBodySettings( "Earth", 0.0, 10.0 * 86400.0 );
    bodySettings[ "Sun" ] = getDefaultSingleBodySettings( "Sun", 0.0,10.0 * 86400.0 );
    bodySettings[ "Moon" ] = getDefaultSingleBodySettings( "Moon", 0.0,10.0 * 86400.0 );
    bodySettings[ "Mars" ] = getDefaultSingleBodySettings( "Mars", 0.0,10.0 * 86400.0 );
    bodySettings[ "Venus" ] = getDefaultSingleBodySettings( "Venus", 0.0,10.0 * 86400.0 );

    // Create settings for sh gravity field.
    double gravitationalParameter = 398600.4418E9;
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
    bodySettings[ "Earth" ]->gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >(
                gravitationalParameter, 6378.0E3, cosineCoefficients, sineCoefficients,
                "IAU_Earth" );

    // Create bodies
    NamedBodyMap bodyMap = createBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    SelectedAccelerationMap accelerationSettingsMap;
    std::map< std::string, std::string > centralBodies;
    std::vector< std::string > propagatedBodyList;
    std::vector< std::string > centralBodyList;

    std::map< IntegratedStateType, Eigen::VectorXd > integratedStateToSet;
    double testTime = 0.0;

    centralBodyList.push_back( centralBodies[ "Moon" ] );
    {

        Eigen::VectorXd testState = ( Eigen::VectorXd( 6 ) << 1.44E6, 2.234E8, -3343.246E7, 1.2E4, 1.344E3, -22.343E3 ).finished( );
        integratedStateToSet[ transational_state ] = testState;
        testTime = 2.0 * 86400.0;

        {
            accelerationSettingsMap[ "Moon" ][ "Sun" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );
            accelerationSettingsMap[ "Moon" ][ "Earth" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );

            // Define origin of integration to be barycenter.
            centralBodies[ "Moon" ] = "SSB";
            propagatedBodyList.push_back( "Moon" );
            centralBodyList.push_back( centralBodies[ "Moon" ] );

            // Create accelerations
            AccelerationMap accelerationsMap = createAccelerationModelsMap(
                        bodyMap, accelerationSettingsMap, centralBodies );

            boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodyList, accelerationsMap, propagatedBodyList, getInitialStateOfBody(
                            "Moon", centralBodies[ "Moon" ], bodyMap, initialTime ) );
            std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate =
                    createEnvironmentUpdaterSettings< double >( bodyMap, propagatorSettings );

            BOOST_CHECK_EQUAL( environmentModelsToUpdate.size( ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_transational_state_update ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_transational_state_update ).size( ), 3 );
            boost::shared_ptr< propagators::EnvironmentUpdater< double, double > > updater =
                    createEnvironmentUpdaterForDynamicalEquations< double, double >(
                        propagatorSettings, bodyMap );

            updater->updateEnvironment( testTime, integratedStateToSet );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Earth" )->getState( ),
                        bodyMap.at( "Earth" )->getEphemeris( )->getCartesianStateFromEphemeris( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Sun" )->getState( ),
                        bodyMap.at( "Sun" )->getEphemeris( )->getCartesianStateFromEphemeris( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Moon" )->getState( ), testState,
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Mars" )->getState( ), basic_mathematics::Vector6d::Zero( ),
                        std::numeric_limits< double >::epsilon( ) );

            updater->updateEnvironment(
                        0.5 * testTime, std::map< IntegratedStateType, Eigen::VectorXd >( ),
                        boost::assign::list_of( transational_state ) );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Earth" )->getState( ),
                        bodyMap.at( "Earth" )->getEphemeris( )->getCartesianStateFromEphemeris( 0.5 * testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Sun" )->getState( ),
                        bodyMap.at( "Sun" )->getEphemeris( )->getCartesianStateFromEphemeris( 0.5 * testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Moon" )->getState( ),
                        bodyMap.at( "Moon" )->getEphemeris( )->getCartesianStateFromEphemeris( 0.5 * testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Mars" )->getState( ), basic_mathematics::Vector6d::Zero( ),
                        std::numeric_limits< double >::epsilon( ) );
        }

        {
            accelerationSettingsMap.clear( );
            centralBodies.clear( );
            propagatedBodyList.clear( );
            centralBodyList.clear( );

            accelerationSettingsMap[ "Moon" ][ "Sun" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );

            accelerationSettingsMap[ "Moon" ][ "Mars" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );

            // Define origin of integration to be barycenter.
            centralBodies[ "Moon" ] = "Earth";
            propagatedBodyList.push_back( "Moon" );
            centralBodyList.push_back( centralBodies[ "Moon" ] );

            // Create accelerations
            AccelerationMap accelerationsMap = createAccelerationModelsMap(
                        bodyMap, accelerationSettingsMap, centralBodies );

            boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodyList, accelerationsMap, propagatedBodyList, getInitialStateOfBody(
                            "Moon", centralBodies[ "Moon" ], bodyMap, initialTime ) );
            std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate =
                    createEnvironmentUpdaterSettings< double >( bodyMap, propagatorSettings );

            BOOST_CHECK_EQUAL( environmentModelsToUpdate.size( ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_transational_state_update ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_transational_state_update ).size( ), 4 );
            boost::shared_ptr< propagators::EnvironmentUpdater< double, double > > updater =
                    createEnvironmentUpdaterForDynamicalEquations< double, double >(
                        propagatorSettings, bodyMap );

            updater->updateEnvironment( testTime, integratedStateToSet );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Earth" )->getState( ),
                        bodyMap.at( "Earth" )->getEphemeris( )->getCartesianStateFromEphemeris( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Sun" )->getState( ),
                        bodyMap.at( "Sun" )->getEphemeris( )->getCartesianStateFromEphemeris( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Moon" )->getState( ), testState,
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Mars" )->getState( ),
                        bodyMap.at( "Mars" )->getEphemeris( )->getCartesianStateFromEphemeris( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Venus" )->getState( ), basic_mathematics::Vector6d::Zero( ),
                        std::numeric_limits< double >::epsilon( ) );

            updater->updateEnvironment(
                        0.5 * testTime, std::map< IntegratedStateType, Eigen::VectorXd >( ),
                        boost::assign::list_of( transational_state ) );

        }

        {
            accelerationSettingsMap.clear( );
            centralBodies.clear( );
            propagatedBodyList.clear( );
            centralBodyList.clear( );

            accelerationSettingsMap[ "Moon" ][ "Sun" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );
            accelerationSettingsMap[ "Moon" ][ "Earth" ].push_back(
                        boost::make_shared< SphericalHarmonicAccelerationSettings >( 6, 6 ) );
            accelerationSettingsMap[ "Moon" ][ "Mars" ].push_back(
                        boost::make_shared< AccelerationSettings >( central_gravity ) );

            // Define origin of integration to be barycenter.
            centralBodies[ "Moon" ] = "Earth";
            propagatedBodyList.push_back( "Moon" );
            centralBodyList.push_back( centralBodies[ "Moon" ] );

            // Create accelerations
            AccelerationMap accelerationsMap = createAccelerationModelsMap(
                        bodyMap, accelerationSettingsMap, centralBodies );

            boost::shared_ptr< PropagatorSettings< double > > propagatorSettings =
                    boost::make_shared< TranslationalStatePropagatorSettings< double > >(
                        centralBodyList, accelerationsMap, propagatedBodyList, getInitialStateOfBody(
                            "Moon", centralBodies[ "Moon" ], bodyMap, initialTime ) );
            std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > > environmentModelsToUpdate =
                    createEnvironmentUpdaterSettings< double >( bodyMap, propagatorSettings );

            BOOST_CHECK_EQUAL( environmentModelsToUpdate.size( ), 3 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_transational_state_update ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_transational_state_update ).size( ), 4 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( body_rotational_state_update ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.at( body_rotational_state_update ).size( ), 1 );
            BOOST_CHECK_EQUAL( environmentModelsToUpdate.count( spherical_harmonic_gravity_field_update ), 1 );

            boost::shared_ptr< propagators::EnvironmentUpdater< double, double > > updater =
                    createEnvironmentUpdaterForDynamicalEquations< double, double >(
                        propagatorSettings, bodyMap );

            updater->updateEnvironment( testTime, integratedStateToSet );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Earth" )->getState( ),
                        bodyMap.at( "Earth" )->getEphemeris( )->getCartesianStateFromEphemeris( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Sun" )->getState( ),
                        bodyMap.at( "Sun" )->getEphemeris( )->getCartesianStateFromEphemeris( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Moon" )->getState( ), testState,
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Mars" )->getState( ),
                        bodyMap.at( "Mars" )->getEphemeris( )->getCartesianStateFromEphemeris( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Venus" )->getState( ), basic_mathematics::Vector6d::Zero( ),
                        std::numeric_limits< double >::epsilon( ) );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Earth" )->getCurrentRotationToGlobalFrame( ).toRotationMatrix( ),
                        bodyMap.at( "Earth" )->getRotationalEphemeris( )->getRotationToBaseFrame( testTime ).toRotationMatrix( ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Earth" )->getCurrentRotationToLocalFrame( ).toRotationMatrix( ),
                        bodyMap.at( "Earth" )->getRotationalEphemeris( )->getRotationToTargetFrame( testTime ).toRotationMatrix( ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Earth" )->getCurrentRotationMatrixDerivativeToGlobalFrame( ),
                        bodyMap.at( "Earth" )->getRotationalEphemeris( )->getDerivativeOfRotationToBaseFrame( testTime ),
                        std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Earth" )->getCurrentRotationMatrixDerivativeToLocalFrame( ),
                        bodyMap.at( "Earth" )->getRotationalEphemeris( )->getDerivativeOfRotationToTargetFrame( testTime ),
                        std::numeric_limits< double >::epsilon( ) );

            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Mars" )->getCurrentRotationToLocalFrame( ).toRotationMatrix( ),
                        Eigen::Matrix3d::Identity( ), std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Mars" )->getCurrentRotationToGlobalFrame( ).toRotationMatrix( ),
                        Eigen::Matrix3d::Identity( ), std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Mars" )->getCurrentRotationMatrixDerivativeToGlobalFrame( ),
                        Eigen::Matrix3d::Zero( ), std::numeric_limits< double >::epsilon( ) );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        bodyMap.at( "Mars" )->getCurrentRotationMatrixDerivativeToLocalFrame( ),
                        Eigen::Matrix3d::Zero( ), std::numeric_limits< double >::epsilon( ) );

            updater->updateEnvironment(
                        0.5 * testTime, std::map< IntegratedStateType, Eigen::VectorXd >( ),
                        boost::assign::list_of( transational_state ) );

        }
    }


    }

    BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat


