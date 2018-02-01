/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <boost/test/unit_test.hpp>
#include <boost/random/uniform_01.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/Gravitation/thirdBodyPerturbation.h"
#include "Tudat/Astrodynamics/Gravitation/mutualSphericalHarmonicGravityModel.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGravityField.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"


namespace tudat
{

namespace unit_tests
{

using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::gravitation;
using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;

//! Generate (dummy_ spherical harmonic coefficients.
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > generateCosineSineCoefficients(
        const int maximumDegree, const int maximumOrder, const int bodyIndex )
{
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumOrder + 1 );
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumOrder + 1 );

    cosineCoefficients( 0, 0 ) = 1.0;

    basic_mathematics::GlobalRandomNumberGeneratorType randomNumberGenerator(
                static_cast< unsigned int >( bodyIndex ) );
    boost::uniform_01< boost::mt19937> distribution( randomNumberGenerator );

    for( int i = 1; i < maximumDegree + 1; i++ )
    {
        for( int j = 0; ( j < maximumOrder + 1 ) && ( j <= i ); j++ )
        {
            cosineCoefficients( i, j ) = ( ( distribution( ) > 0.5 ) ? ( 1.0 ): ( -1.0 ) ) * distribution( ) * 1.0E-2;
            if( j > 0 )
            {
                sineCoefficients( i, j ) =  ( ( distribution( ) > 0.5 ) ? ( 1.0 ): ( -1.0 ) ) * distribution( ) * 1.0E-2;
            }

        }
    }

    return std::make_pair( cosineCoefficients, sineCoefficients );
}

//! Generate gravity field object (with severely exagerated magnitude).
boost::shared_ptr< tudat::simulation_setup::GravityFieldSettings > getDummyJovianSystemGravityField(
        const std::string& bodyName )
{
    boost::shared_ptr< GravityFieldSettings > gravityFieldSettings;

    std::vector< double > randomNumberSettings;
    randomNumberSettings.push_back( 0.0 );
    randomNumberSettings.push_back( 1.0E-4 );

    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > coefficients;

    if( bodyName == "Jupiter" )
    {
        coefficients = generateCosineSineCoefficients( 10, 10, 0 );

        gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >
                ( getBodyGravitationalParameter( "Jupiter" ), getAverageRadius( "Jupiter" ),
                  coefficients.first, coefficients.second, "IAU_Jupiter" );
    }
    else if( bodyName == "Io" )
    {
        coefficients = generateCosineSineCoefficients( 10, 10, 1 );

        gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >
                ( 5.959916033410404E012, 200.0 * getAverageRadius( "Io" ),
                  coefficients.first, coefficients.second, "IAU_Io" );
    }
    else if( bodyName == "Europa" )
    {
        coefficients = generateCosineSineCoefficients( 10, 10, 2 );

        gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >
                ( 3.202738774922892E12, 200.0 * getAverageRadius( "Europa" ),
                  coefficients.first, coefficients.second, "IAU_Europa" );
    }

    return gravityFieldSettings;

}


BOOST_AUTO_TEST_SUITE( test_mutual_spherical_harmonic_gravity )

//! Test mutual spherical harmonic acceleration against manually combined spherical harmonic accelerations.
//! Note that the size of the Galilean moons, as well as the magnitude of the spherical harmonic coefficients has
//! been exaggerated to perform a more robust test (i.e ensure that typical errors in implementation are well above
//! numerical errors).
BOOST_AUTO_TEST_CASE( testMutualSphericalHarmonicGravity )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create list of bodies to create.
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Sun" );

    // Specify initial time
    double initialTime = 1.0E7;
    double finalTime = 1.2E7;

    // Get body settings.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialTime, finalTime );
    bodySettings[ "Jupiter" ]->gravityFieldSettings = getDummyJovianSystemGravityField( "Jupiter" );
    bodySettings[ "Io" ]->gravityFieldSettings = getDummyJovianSystemGravityField( "Io" );
    bodySettings[ "Europa" ]->gravityFieldSettings = getDummyJovianSystemGravityField( "Europa" );

    bodySettings[ "Jupiter" ]->ephemerisSettings = boost::make_shared< KeplerEphemerisSettings >(
                ( Eigen::Vector6d( ) << 778.57E9, 0.0489, 1.3 / 60.0, 0.0, 0.0, 0.0 ).finished( ), 0.0,
                  getBodyGravitationalParameter( "Sun" ), "Sun", "ECLIPJ2000" );
    bodySettings[ "Io" ]->ephemerisSettings = boost::make_shared< KeplerEphemerisSettings >(
                ( Eigen::Vector6d( ) << 421.8E6, 0.004, 0.04 / 60.0, 0.0, 0.0, 0.0 ).finished( ), 0.0,
                  getBodyGravitationalParameter( "Jupiter" ), "Sun", "ECLIPJ2000" );
    bodySettings[ "Europa" ]->ephemerisSettings = boost::make_shared< KeplerEphemerisSettings >(
                ( Eigen::Vector6d( ) << 671.1E6, 0.009, 0.47 / 60.0, 0.0, 0.0, 0.0 ).finished( ), 0.0,
                  getBodyGravitationalParameter( "Jupiter" ), "Sun", "ECLIPJ2000" );


    // Create bodies needed in simulation
    NamedBodyMap bodyMap = createBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Set current state and rotation of bodies.
    double currentTime = 1.1E7;
    bodyMap[ "Jupiter" ]->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
    bodyMap[ "Jupiter" ]->setStateFromEphemeris( currentTime );
    bodyMap[ "Io" ]->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
    bodyMap[ "Io" ]->setStateFromEphemeris( currentTime );
    bodyMap[ "Europa" ]->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
    bodyMap[ "Europa" ]->setStateFromEphemeris( currentTime );

    // Retrieve gravity fields.
    boost::shared_ptr< SphericalHarmonicsGravityField > jupiterGravityField = boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                ( bodyMap.at( "Jupiter" ) )->getGravityFieldModel( ) );
    boost::shared_ptr< SphericalHarmonicsGravityField > ioGravityField = boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                ( bodyMap.at( "Io" ) )->getGravityFieldModel( ) );
    boost::shared_ptr< SphericalHarmonicsGravityField > europaGravityField = boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                ( bodyMap.at( "Europa" ) )->getGravityFieldModel( ) );

    // Create central gravity acceleration (mu = Io + Jupiter)
    boost::shared_ptr< AccelerationSettings > centralGravitySettings = boost::make_shared< AccelerationSettings >( central_gravity );
    boost::shared_ptr< CentralGravitationalAccelerationModel3d > centralGravity =
            boost::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                createAccelerationModel( bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), centralGravitySettings, "Io", "Jupiter",
                                         bodyMap.at( "Jupiter" ), "Jupiter" ) );

    // Calculate central gravity acceleration.
    centralGravity->updateMembers( );
    Eigen::Vector3d centralGravityAcceleration = centralGravity->getAcceleration( );

    // Create spherical harmonic gravity of Jupiter on Io, Jupiter-fixed (mu = Io + Jupiter)
    boost::shared_ptr< AccelerationSettings > sphericalHarmonicGravityOnIoFromJupiterSettings =
            boost::make_shared< SphericalHarmonicAccelerationSettings >( 7, 7 );
    boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicGravityOnIoFromJupiter =
            boost::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                createAccelerationModel(  bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), sphericalHarmonicGravityOnIoFromJupiterSettings,
                                          "Io", "Jupiter", bodyMap.at( "Jupiter" ), "Jupiter" ) );

    // Calculate spherical harmonic gravity of Jupiter on Io.
    sphericalHarmonicGravityOnIoFromJupiter->updateMembers( );
    Eigen::Vector3d sphericalHarmonicGravityOnIoFromJupiterAcceleration = sphericalHarmonicGravityOnIoFromJupiter->getAcceleration( );

    // Create spherical harmonic gravity of Io on Jupiter, Io-fixed (mu = Io + Jupiter)
    boost::shared_ptr< AccelerationSettings > sphericalHarmonicGravityOnJupiterFromIoSettings =
            boost::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 );
    boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicGravityOnJupiterFromIo =
            boost::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                createAccelerationModel( bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ), sphericalHarmonicGravityOnJupiterFromIoSettings,
                                         "Jupiter", "Io", bodyMap.at( "Io" ), "Io" ) );

    // Calculate spherical harmonic gravity of Io on Jupiter.
    sphericalHarmonicGravityOnJupiterFromIo->updateMembers( );
    Eigen::Vector3d sphericalHarmonicGravityOnJupiterFromIoAcceleration = sphericalHarmonicGravityOnJupiterFromIo->getAcceleration( );

    // Create mutual spherical harmonic gravity between Io and Jupiter on Io, Jupiter fixed (mu = Io + Jupiter)
    boost::shared_ptr< AccelerationSettings > mutualDirectJupiterIoShGravitySettings =
            boost::make_shared< MutualSphericalHarmonicAccelerationSettings >( 7, 7, 2, 2 );
    boost::shared_ptr< MutualSphericalHarmonicsGravitationalAccelerationModel > mutualDirectJupiterIoShGravity =
            boost::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                createAccelerationModel( bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), mutualDirectJupiterIoShGravitySettings,
                                         "Io", "Jupiter", bodyMap.at( "Jupiter" ), "Jupiter" ) );

    // Calculate mutual spherical harmonic gravity between Io and Jupiter on Io.
    mutualDirectJupiterIoShGravity->updateMembers( );
    Eigen::Vector3d mutualDirectJupiterIoShGravityAcceleration = mutualDirectJupiterIoShGravity->getAcceleration( );

    // Calculate expected mutual spherical harmonic gravity from sub-accelerations.
    Eigen::Vector3d expectedAcceleration = -centralGravityAcceleration + sphericalHarmonicGravityOnIoFromJupiterAcceleration -
            sphericalHarmonicGravityOnJupiterFromIoAcceleration;

    // Test against directly calculated mutual spherical harmonic gravity.
    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( expectedAcceleration( i ) - mutualDirectJupiterIoShGravityAcceleration( i ) ),
                           12.0 * std::numeric_limits< double >::epsilon( ) * expectedAcceleration.norm( ) );
    }

    // Create mutual spherical harmonic gravity between Io and Jupiter on Jupiter, Io fixed (mu = Io + Jupiter)
    boost::shared_ptr< AccelerationSettings > mutualDirectJupiterIoShGravitySettings2 =
            boost::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 7, 7 );
    boost::shared_ptr< MutualSphericalHarmonicsGravitationalAccelerationModel > mutualDirectJupiterIoShGravity2 =
            boost::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                createAccelerationModel( bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ), mutualDirectJupiterIoShGravitySettings2,
                                         "Jupiter", "Io", bodyMap.at( "Io" ), "Io" ) );

    // Calculate mutual spherical harmonic gravity between Io and Jupiter on Jupiter.
    mutualDirectJupiterIoShGravity2->updateMembers( );
    Eigen::Vector3d mutualDirectJupiterIoShGravityAcceleration2 = mutualDirectJupiterIoShGravity2->getAcceleration( );

    expectedAcceleration = centralGravityAcceleration - sphericalHarmonicGravityOnIoFromJupiterAcceleration + sphericalHarmonicGravityOnJupiterFromIoAcceleration;


    // Test against directly calculated mutual spherical harmonic gravity.
    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( expectedAcceleration( i ) - mutualDirectJupiterIoShGravityAcceleration2( i ) ),
                           12.0 * std::numeric_limits< double >::epsilon( ) * expectedAcceleration.norm( ) );
    }

    // Test against directly calculated mutual spherical harmonic gravity.
    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( mutualDirectJupiterIoShGravityAcceleration( i ) + mutualDirectJupiterIoShGravityAcceleration2( i ) ),
                           12.0 * std::numeric_limits< double >::epsilon( ) * mutualDirectJupiterIoShGravityAcceleration.norm( ) );
    }


    // Create 3rd body mutual spherical harmonics between Io and Europa on Europa, Jupiter fixed (mu = Io)
    boost::shared_ptr< AccelerationSettings > mutualThirdBodyIoOnEuropaShGravitySettings =
            boost::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 4, 4, 7, 7 );
    boost::shared_ptr< ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel > mutualThirdBodyIoOnEuropaShGravity =
            boost::dynamic_pointer_cast< ThirdBodyMutualSphericalHarmonicsGravitationalAccelerationModel >(
                createAccelerationModel( bodyMap.at( "Europa" ), bodyMap.at( "Io" ), mutualThirdBodyIoOnEuropaShGravitySettings,
                                         "Europa", "Io", bodyMap.at( "Jupiter" ), "Jupiter" ) );

    // Calculate 3rd body mutual spherical harmonics between Io and Europa on Europa.
    mutualThirdBodyIoOnEuropaShGravity->updateMembers( );
    Eigen::Vector3d mutualThirdBodyIoOnEuropaShGravityAcceleration = mutualThirdBodyIoOnEuropaShGravity->getAcceleration( );


    // Create mutual spherical harmonics between Io and Europa on Europa, Io fixed (mu = Io + Europa)
    boost::shared_ptr< MutualSphericalHarmonicsGravitationalAccelerationModel > mutualDirectIoOnEuropaShGravity =
            boost::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                createAccelerationModel( bodyMap.at( "Europa" ), bodyMap.at( "Io" ), mutualThirdBodyIoOnEuropaShGravitySettings,
                                         "Europa", "Io", bodyMap.at( "Io" ), "Io" ) );
    mutualDirectIoOnEuropaShGravity->updateMembers( );
    Eigen::Vector3d mutualDirectIoOnEuropaShGravityAcceleration = mutualDirectIoOnEuropaShGravity->getAcceleration( );

    // Get sub accelerations from 3rd body acceleration.
    Eigen::Vector3d directAccelerationFromThirdBodyModel = mutualThirdBodyIoOnEuropaShGravity->getAccelerationModelForBodyUndergoingAcceleration( )->getAcceleration( );
    Eigen::Vector3d centralBodyAccelerationFromThirdBodyModel = mutualThirdBodyIoOnEuropaShGravity->getAccelerationModelForCentralBody( )->
            getAcceleration( );

    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL(
                    ( directAccelerationFromThirdBodyModel( i ) -
                      ( ioGravityField->getGravitationalParameter( ) /
                        ( ioGravityField->getGravitationalParameter( ) + europaGravityField->getGravitationalParameter( ) ) *
                        mutualDirectIoOnEuropaShGravityAcceleration( i ) ) ),
                    ( 12.0 * std::numeric_limits< double >::epsilon( ) * directAccelerationFromThirdBodyModel.norm( ) ) );

        BOOST_CHECK_SMALL(
                    centralBodyAccelerationFromThirdBodyModel( i ) -
                    ( ioGravityField->getGravitationalParameter( ) /
                      ( ioGravityField->getGravitationalParameter( ) + jupiterGravityField->getGravitationalParameter( ) ) *
                      mutualDirectJupiterIoShGravityAcceleration2( i ) ),
                    ( 12.0 * std::numeric_limits< double >::epsilon( ) * centralBodyAccelerationFromThirdBodyModel.norm( ) ) );

        BOOST_CHECK_SMALL(
                    mutualThirdBodyIoOnEuropaShGravityAcceleration( i ) -
                    ( directAccelerationFromThirdBodyModel( i ) - centralBodyAccelerationFromThirdBodyModel( i ) ),
                    ( std::numeric_limits< double >::epsilon( ) * mutualThirdBodyIoOnEuropaShGravityAcceleration.norm( ) ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
