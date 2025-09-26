/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"

#include "tudat/interface/spice/spiceRotationalEphemeris.h"
#include "tudat/astro/ephemerides/fullPlanetaryRotationModel.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;

BOOST_AUTO_TEST_SUITE( test_planetary_rotation_model )

BOOST_AUTO_TEST_CASE( testPlanetaryRotationModel )
{
    double initialTime = 0.0;
    double finalTime = 1.0E8;

    spice_interface::loadStandardSpiceKernels( );

    SystemOfBodies bodies;
    bodies.createEmptyBody( "Mars" );

    std::shared_ptr< RotationModelSettings > defaultMarsRotationSettings = getHighAccuracyMarsRotationModel( );

    std::shared_ptr< RotationalEphemeris > marsRotationModel = createRotationModel( defaultMarsRotationSettings, "Mars" );

    double timeStep = 3600.0;
    double currentTime = initialTime + timeStep;
    double dt = 0.9;

    std::map< double, Eigen::MatrixXd > DerivativeDifferenceMap;

    while( currentTime < finalTime - timeStep )
    {
        Eigen::MatrixXd numericalRotationMatrixderivative =
                ( ( marsRotationModel->getRotationToBaseFrame( currentTime + dt ) ).toRotationMatrix( ) -
                  ( marsRotationModel->getRotationToBaseFrame( currentTime - dt ) ).toRotationMatrix( ) ) /
                ( 2 * dt );

        DerivativeDifferenceMap[ currentTime ] =
                marsRotationModel->getDerivativeOfRotationToBaseFrame( currentTime ) - numericalRotationMatrixderivative;

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( DerivativeDifferenceMap[ currentTime ]( i, j ), 1.0E-10 );
            }
        }

        currentTime += timeStep;
    }
}

BOOST_AUTO_TEST_CASE( testPlanetaryRotationModelBaseFrameConsistency )
{
    spice_interface::loadStandardSpiceKernels( );

    // Create rotation models for Mars with J2000 and ECLIPJ2000 base frames
    std::shared_ptr< RotationModelSettings > marsRotationSettingsJ2000 = getHighAccuracyMarsRotationModel( "J2000" );
    std::shared_ptr< RotationalEphemeris > marsRotationModelJ2000 = createRotationModel( marsRotationSettingsJ2000, "Mars" );

    std::shared_ptr< RotationModelSettings > marsRotationSettingsEclipJ2000 = getHighAccuracyMarsRotationModel( "ECLIPJ2000" );
    std::shared_ptr< RotationalEphemeris > marsRotationModelEclipJ2000 = createRotationModel( marsRotationSettingsEclipJ2000, "Mars" );

    // Define test time
    double testTime = 1.0E7;

    // Get rotations from body-fixed frame to the respective base frames
    Eigen::Quaterniond rotationToJ2000 = marsRotationModelJ2000->getRotationToBaseFrame( testTime );
    Eigen::Quaterniond rotationToEclipJ2000 = marsRotationModelEclipJ2000->getRotationToBaseFrame( testTime );

    // Get the expected rotation between the two base frames from SPICE
    Eigen::Quaterniond rotationJ2000ToEclipJ2000 =
            spice_interface::computeRotationQuaternionBetweenFrames( "J2000", "ECLIPJ2000", testTime );

    // Verify consistency: R_IAU_MARS_to_ECLIPJ2000 should equal R_J2000_to_ECLIPJ2000 * R_IAU_MARS_to_J2000
    Eigen::Quaterniond computedRotationToEclipJ2000 = rotationJ2000ToEclipJ2000 * rotationToJ2000;

    // Check that the difference between the computed and actual rotation is negligible
    Eigen::Matrix3d difference = rotationToEclipJ2000.toRotationMatrix( ) - computedRotationToEclipJ2000.toRotationMatrix( );
    for( int i = 0; i < 3; i++ )
    {
        for( int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::abs( difference( i, j ) ), 1.0E-15 );
        }
    }

    // Test for unsupported base frame, expecting an exception
    {
        bool exceptionThrown = false;
        try
        {
            std::shared_ptr< RotationModelSettings > marsRotationSettingsUnsupported =
                    getHighAccuracyMarsRotationModel( "UnsupportedFrame" );
            std::shared_ptr< RotationalEphemeris > failModel = createRotationModel( marsRotationSettingsUnsupported, "Mars" );
        }
        catch( const std::runtime_error& e )
        {
            // Verify that the error message is as expected
            std::string errorMessage = e.what( );
            bool correctMessage =
                    errorMessage.find( "Error in PlanetaryRotationModel, base frame UnsupportedFrame not recognized" ) != std::string::npos;
            BOOST_CHECK( correctMessage );
            exceptionThrown = true;
        }
        BOOST_CHECK( exceptionThrown );
    }
}

BOOST_AUTO_TEST_CASE( testPlanetaryRotationModelCustomParameters )
{
    using namespace tudat::unit_conversions;
    using namespace tudat::physical_constants;

    spice_interface::loadStandardSpiceKernels( );

    // Define default parameters (from Konopliv et al., 2016)
    const double milliArcSecondToRadian = mathematical_constants::PI / ( 180.0 * 1000.0 * 3600.0 );
    const double defaultAngleN = convertDegreesToRadians( 3.37919183 );
    const double defaultAngleJ = convertDegreesToRadians( 24.67682669 );
    const double defaultAnglePsiAtEpoch = convertDegreesToRadians( 81.9683988 );
    const double defaultAnglePsiRateAtEpoch = ( -7608.3 * milliArcSecondToRadian ) / JULIAN_YEAR;

    // 1. Create model with default parameters (no overrides)
    std::shared_ptr< RotationModelSettings > defaultSettings = getHighAccuracyMarsRotationModel( );
    std::shared_ptr< RotationalEphemeris > defaultModel = createRotationModel( defaultSettings, "Mars" );

    // 2. Create model by explicitly passing the default parameters
    std::shared_ptr< RotationModelSettings > explicitDefaultSettings = getHighAccuracyMarsRotationModel(
            "J2000", "Mars_Fixed", defaultAngleN, defaultAngleJ, defaultAnglePsiAtEpoch, defaultAnglePsiRateAtEpoch );
    std::shared_ptr< RotationalEphemeris > explicitDefaultModel = createRotationModel( explicitDefaultSettings, "Mars" );

    // 3. Create model with custom parameters
    const double customAngleN = defaultAngleN + 0.1;
    const double customAnglePsiAtEpoch = defaultAnglePsiAtEpoch - 0.05;
    std::shared_ptr< RotationModelSettings > customSettings = getHighAccuracyMarsRotationModel(
            "J2000", "Mars_Fixed_Custom", customAngleN, defaultAngleJ, customAnglePsiAtEpoch, defaultAnglePsiRateAtEpoch );
    std::shared_ptr< RotationalEphemeris > customModel = createRotationModel( customSettings, "Mars" );

    // Define test time
    double testTime = 1.0E7;

    // Get rotations
    Eigen::Quaterniond defaultRotation = defaultModel->getRotationToTargetFrame( testTime );
    Eigen::Quaterniond explicitDefaultRotation = explicitDefaultModel->getRotationToTargetFrame( testTime );
    Eigen::Quaterniond customRotation = customModel->getRotationToTargetFrame( testTime );

    // Verify that default and explicit-default models are identical
    Eigen::Matrix3d differenceDefault = defaultRotation.toRotationMatrix( ) - explicitDefaultRotation.toRotationMatrix( );
    for( int i = 0; i < 3; i++ )
    {
        for( int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::abs( differenceDefault( i, j ) ), 1.0E-15 );
        }
    }

    // Verify that default and custom models are different
    Eigen::Matrix3d differenceCustom = defaultRotation.toRotationMatrix( ) - customRotation.toRotationMatrix( );
    double differenceNorm = differenceCustom.norm( );
    BOOST_CHECK_GT( differenceNorm, 1.0E-3 );
}

BOOST_AUTO_TEST_CASE( testPlanetaryRotationModelFullyCustomParameters )
{
    using namespace tudat::unit_conversions;
    using namespace tudat::physical_constants;

    spice_interface::loadStandardSpiceKernels( );

    // Define default parameters (from Konopliv et al., 2016)
    const double milliArcSecondToRadian = mathematical_constants::PI / ( 180.0 * 1000.0 * 3600.0 );
    const double defaultAngleN = convertDegreesToRadians( 3.37919183 );
    const double defaultAngleJ = convertDegreesToRadians( 24.67682669 );
    const double defaultAnglePsiAtEpoch = convertDegreesToRadians( 81.9683988 );
    const double defaultAnglePsiRateAtEpoch = ( -7608.3 * milliArcSecondToRadian ) / JULIAN_YEAR;

    // Create default nutation correction settings
    std::map< double, std::pair< double, double > > defaultNutationCorrectionSettings;
    defaultNutationCorrectionSettings[ 0.0 ] = std::make_pair( -1.4 * milliArcSecondToRadian, 0.0 );
    defaultNutationCorrectionSettings[ 1.0 ] = std::make_pair( -0.4 * milliArcSecondToRadian, -632.6 * milliArcSecondToRadian );
    defaultNutationCorrectionSettings[ 2.0 ] = std::make_pair( 0.0, -44.2 * milliArcSecondToRadian );
    defaultNutationCorrectionSettings[ 3.0 ] = std::make_pair( 0.0, -4.0 * milliArcSecondToRadian );

    // Create default mean motion time dependent phase nutation corrections
    std::vector< std::map< double, std::pair< double, double > > > defaultMeanMotionTimeDependentPhaseNutationCorrections;
    std::map< double, std::pair< double, double > > meanMotionTimeDependentPhaseNutationCorrection;
    meanMotionTimeDependentPhaseNutationCorrection[ 1.0 ] =
            std::make_pair( -49.1 * milliArcSecondToRadian, -104.5 * milliArcSecondToRadian );
    meanMotionTimeDependentPhaseNutationCorrection[ 2.0 ] =
            std::make_pair( 515.7 * milliArcSecondToRadian, 1097.0 * milliArcSecondToRadian );
    meanMotionTimeDependentPhaseNutationCorrection[ 3.0 ] =
            std::make_pair( 112.8 * milliArcSecondToRadian, 240.1 * milliArcSecondToRadian );
    meanMotionTimeDependentPhaseNutationCorrection[ 4.0 ] = std::make_pair( 19.2 * milliArcSecondToRadian, 40.9 * milliArcSecondToRadian );
    meanMotionTimeDependentPhaseNutationCorrection[ 5.0 ] = std::make_pair( 3.0 * milliArcSecondToRadian, 6.5 * milliArcSecondToRadian );
    meanMotionTimeDependentPhaseNutationCorrection[ 6.0 ] = std::make_pair( 0.4 * milliArcSecondToRadian, 1.0 * milliArcSecondToRadian );
    defaultMeanMotionTimeDependentPhaseNutationCorrections.push_back( meanMotionTimeDependentPhaseNutationCorrection );

    // Create default rotation rate corrections
    std::map< double, std::pair< double, double > > defaultRotationRateCorrections;
    defaultRotationRateCorrections[ 1.0 ] =
            std::make_pair( 481.0 * milliArcSecondToRadian, -155.0 * milliArcSecondToRadian - 176 * milliArcSecondToRadian );
    defaultRotationRateCorrections[ 2.0 ] =
            std::make_pair( -103.0 * milliArcSecondToRadian, -93.0 * milliArcSecondToRadian - 8 * milliArcSecondToRadian );
    defaultRotationRateCorrections[ 3.0 ] =
            std::make_pair( -35.0 * milliArcSecondToRadian, -3.0 * milliArcSecondToRadian - 1 * milliArcSecondToRadian );
    defaultRotationRateCorrections[ 4.0 ] = std::make_pair( -10.0 * milliArcSecondToRadian, -8.0 * milliArcSecondToRadian );

    // Create default polar motion coefficients
    std::map< double, std::pair< double, double > > defaultXPolarMotionCoefficients;
    defaultXPolarMotionCoefficients[ 1.0 ] = std::make_pair( 2.8 * milliArcSecondToRadian * std::sin( convertDegreesToRadians( 46.5 ) ),
                                                             2.8 * milliArcSecondToRadian * std::cos( convertDegreesToRadians( 46.5 ) ) );
    defaultXPolarMotionCoefficients[ 2.0 ] = std::make_pair( 8.9 * milliArcSecondToRadian * std::sin( convertDegreesToRadians( -150.1 ) ),
                                                             8.9 * milliArcSecondToRadian * std::cos( convertDegreesToRadians( -150.1 ) ) );
    defaultXPolarMotionCoefficients[ 3.0 ] = std::make_pair( 0.0, 0.0 );
    defaultXPolarMotionCoefficients[ 4.0 ] = std::make_pair( 0.0, 0.0 );
    defaultXPolarMotionCoefficients[ 3.34 ] = std::make_pair( 0.0, 50.0 * milliArcSecondToRadian );  // Mars's Chandler wobble T=205 dd

    std::map< double, std::pair< double, double > > defaultYPolarMotionCoefficients;
    defaultYPolarMotionCoefficients[ 1.0 ] = std::make_pair( 11.7 * milliArcSecondToRadian * std::sin( convertDegreesToRadians( 118.7 ) ),
                                                             11.7 * milliArcSecondToRadian * std::cos( convertDegreesToRadians( 118.7 ) ) );
    defaultYPolarMotionCoefficients[ 2.0 ] = std::make_pair( 3.9 * milliArcSecondToRadian * std::sin( convertDegreesToRadians( 172.5 ) ),
                                                             3.9 * milliArcSecondToRadian * std::cos( convertDegreesToRadians( 118.7 ) ) );
    defaultYPolarMotionCoefficients[ 3.0 ] = std::make_pair( 0.0, 0.0 );
    defaultYPolarMotionCoefficients[ 4.0 ] = std::make_pair( 0.0, 0.0 );
    defaultYPolarMotionCoefficients[ 3.34 ] = std::make_pair( 0.0, 50.0 * milliArcSecondToRadian );  // Mars's Chandler wobble T=205 dd

    // 1. Create model with default parameters (no overrides)
    std::shared_ptr< RotationModelSettings > defaultSettings = getHighAccuracyMarsRotationModel( );
    std::shared_ptr< RotationalEphemeris > defaultModel = createRotationModel( defaultSettings, "Mars" );

    // 2. Create model by explicitly passing all the default parameters
    std::shared_ptr< RotationModelSettings > explicitDefaultSettings =
            getHighAccuracyMarsRotationModel( "J2000",
                                              "Mars_Fixed",
                                              defaultAngleN,
                                              defaultAngleJ,
                                              defaultAnglePsiAtEpoch,
                                              defaultAnglePsiRateAtEpoch,
                                              defaultNutationCorrectionSettings,
                                              defaultMeanMotionTimeDependentPhaseNutationCorrections,
                                              defaultRotationRateCorrections,
                                              defaultXPolarMotionCoefficients,
                                              defaultYPolarMotionCoefficients );
    std::shared_ptr< RotationalEphemeris > explicitDefaultModel = createRotationModel( explicitDefaultSettings, "Mars" );

    // 3. Create model with custom nutation correction parameters
    std::map< double, std::pair< double, double > > customNutationCorrectionSettings = defaultNutationCorrectionSettings;
    customNutationCorrectionSettings[ 0.0 ] = std::make_pair( -2.0 * milliArcSecondToRadian, 0.0 );  // Modified from -1.4 to -2.0

    std::shared_ptr< RotationModelSettings > customSettings =
            getHighAccuracyMarsRotationModel( "J2000",
                                              "Mars_Fixed_Custom",
                                              defaultAngleN,
                                              defaultAngleJ,
                                              defaultAnglePsiAtEpoch,
                                              defaultAnglePsiRateAtEpoch,
                                              customNutationCorrectionSettings,
                                              defaultMeanMotionTimeDependentPhaseNutationCorrections,
                                              defaultRotationRateCorrections,
                                              defaultXPolarMotionCoefficients,
                                              defaultYPolarMotionCoefficients );
    std::shared_ptr< RotationalEphemeris > customModel = createRotationModel( customSettings, "Mars" );

    // 4. Create model with custom polar motion coefficients
    std::map< double, std::pair< double, double > > customXPolarMotionCoefficients = defaultXPolarMotionCoefficients;
    customXPolarMotionCoefficients[ 1.0 ] =
            std::make_pair( 5.0 * milliArcSecondToRadian * std::sin( convertDegreesToRadians( 46.5 ) ),
                            5.0 * milliArcSecondToRadian * std::cos( convertDegreesToRadians( 46.5 ) ) );  // Modified from 2.8 to 5.0

    std::shared_ptr< RotationModelSettings > customPolarSettings =
            getHighAccuracyMarsRotationModel( "J2000",
                                              "Mars_Fixed_Custom2",
                                              defaultAngleN,
                                              defaultAngleJ,
                                              defaultAnglePsiAtEpoch,
                                              defaultAnglePsiRateAtEpoch,
                                              defaultNutationCorrectionSettings,
                                              defaultMeanMotionTimeDependentPhaseNutationCorrections,
                                              defaultRotationRateCorrections,
                                              customXPolarMotionCoefficients,
                                              defaultYPolarMotionCoefficients );
    std::shared_ptr< RotationalEphemeris > customPolarModel = createRotationModel( customPolarSettings, "Mars" );

    // Define test time
    double testTime = 1.0E7;

    // Get rotations
    Eigen::Quaterniond defaultRotation = defaultModel->getRotationToTargetFrame( testTime );
    Eigen::Quaterniond explicitDefaultRotation = explicitDefaultModel->getRotationToTargetFrame( testTime );
    Eigen::Quaterniond customRotation = customModel->getRotationToTargetFrame( testTime );
    Eigen::Quaterniond customPolarRotation = customPolarModel->getRotationToTargetFrame( testTime );

    // Verify that default and explicit-default models are identical
    Eigen::Matrix3d differenceDefault = defaultRotation.toRotationMatrix( ) - explicitDefaultRotation.toRotationMatrix( );
    for( int i = 0; i < 3; i++ )
    {
        for( int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::abs( differenceDefault( i, j ) ), 1.0E-15 );
        }
    }

    // Verify that default and custom nutation models are different
    Eigen::Matrix3d differenceCustom = defaultRotation.toRotationMatrix( ) - customRotation.toRotationMatrix( );
    double differenceNorm = differenceCustom.norm( );
    BOOST_CHECK_GT( differenceNorm, 1.0E-10 );

    // Verify that default and custom polar motion models are different
    Eigen::Matrix3d differenceCustomPolar = defaultRotation.toRotationMatrix( ) - customPolarRotation.toRotationMatrix( );
    double differenceNormPolar = differenceCustomPolar.norm( );
    BOOST_CHECK_GT( differenceNormPolar, 1.0E-10 );
}
BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
