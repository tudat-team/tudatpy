/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/astro/system_models/camera.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/readSumLmkFiles.h"

namespace tudat
{
namespace unit_tests
{

namespace
{

void checkClose( const double actual, const double expected, const double tolerance = 1.0E-12 )
{
    BOOST_CHECK_SMALL( std::fabs( actual - expected ), tolerance );
}

void writeTextFile( const boost::filesystem::path& filePath, const std::string& contents )
{
    std::ofstream stream( filePath.string( ) );
    stream << contents;
}

boost::filesystem::path makeTemporaryPath( const std::string& suffix )
{
    return boost::filesystem::temp_directory_path( ) / boost::filesystem::unique_path( "tudat_sum_lmk_%%%%%%" + suffix );
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_sum_lmk_file_reader )

BOOST_AUTO_TEST_CASE( testReadSumAndLmkFiles )
{
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const input_output::sum_lmk::SumImageData image = input_output::sum_lmk::readSumFile( testDataPath + "/sample.sum" );

    BOOST_CHECK_EQUAL( image.imageId_, "IMG001" );
    BOOST_CHECK_EQUAL( image.utcEpochString_, "2015 JUN 05 07:24:42.053" );
    BOOST_CHECK_EQUAL( image.imageSize_( 0 ), 1024 );
    BOOST_CHECK_EQUAL( image.imageSize_( 1 ), 1024 );
    BOOST_CHECK_EQUAL( image.threshold_, 500 );
    BOOST_CHECK_EQUAL( image.maxDn_, 65535 );
    checkClose( image.focalLengthMm_, 100.0 );
    checkClose( image.opticalCenter_( 0 ), 512.0 );
    checkClose( image.opticalCenter_( 1 ), 512.0 );
    checkClose( image.spacecraftObjectVector_( 2 ), 10000.0 );  // SCOBJ = spacecraft-to-object (target ahead)
    checkClose( image.cameraAxes_( 0, 0 ), 1.0 );
    checkClose( image.kMatrix_( 0, 0 ), 10.0 );
    checkClose( image.kMatrix_( 1, 1 ), 10.0 );
    checkClose( image.pointingSigma_( 0 ), 1.0E-4 );
    BOOST_REQUIRE_EQUAL( image.landmarkObservations_.size( ), 2 );
    BOOST_REQUIRE_EQUAL( image.limbFitObservations_.size( ), 1 );
    BOOST_CHECK_EQUAL( image.limbFitObservations_.at( 0 ).featureId_, "LIMB01" );

    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks =
            input_output::sum_lmk::readLmkFiles( { testDataPath + "/LMK0001.lmk", testDataPath + "/LMK0002.lmk" } );
    BOOST_REQUIRE_EQUAL( landmarks.size( ), 2 );
    BOOST_REQUIRE( landmarks.count( "LMK0001" ) == 1 );
    BOOST_REQUIRE( landmarks.count( "LMK0002" ) == 1 );
    checkClose( landmarks.at( "LMK0002" ).bodyFixedPosition_( 0 ), 1000.0 );
    checkClose( landmarks.at( "LMK0001" ).landmarkPositionSigma_( 0 ), 1.0 );
    BOOST_REQUIRE_EQUAL( landmarks.at( "LMK0001" ).pictures_.size( ), 1 );
    BOOST_CHECK( landmarks.at( "LMK0001" ).pictures_.at( 0 ).flagged_ );

    const std::shared_ptr< system_models::PsfCameraProjectionModel > projectionModel =
            std::make_shared< system_models::PsfCameraProjectionModel >(
                    image.focalLengthMm_, image.opticalCenter_, image.kMatrix_, Eigen::Matrix< double, 6, 1 >::Zero( ) );
    // SCOBJ is spacecraft-to-object, so spacecraft body-fixed position is -SCOBJ and the landmark
    // position relative to the spacecraft is VLM - (-SCOBJ) = VLM + SCOBJ.
    const Eigen::Vector3d relativeBodyFixedPosition = landmarks.at( "LMK0002" ).bodyFixedPosition_ + image.spacecraftObjectVector_;
    const Eigen::Vector2d projectedPixelLine =
            projectionModel->projectUnitVectorToPixelLine( image.cameraAxes_ * relativeBodyFixedPosition );
    checkClose( projectedPixelLine( 0 ), 612.0 );
    checkClose( projectedPixelLine( 1 ), 512.0 );
}

BOOST_AUTO_TEST_CASE( testSumLmkDuplicateFailures )
{
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";

    const boost::filesystem::path duplicateSumPath = makeTemporaryPath( ".sum" );
    writeTextFile( duplicateSumPath,
                   "IMGDUP\n"
                   "2015 JUN 05 07:24:42.053\n"
                   "1024 1024 500 65535 NPX, NLN, THRSH\n"
                   "100.0 512.0 512.0 MMFL, CTR\n"
                   "0.0 0.0 -10.0 SCOBJ\n"
                   "1.0 0.0 0.0 CX\n"
                   "0.0 1.0 0.0 CY\n"
                   "0.0 0.0 1.0 CZ\n"
                   "10.0 0.0 0.0 0.0 10.0 0.0 K-MATRIX\n"
                   "LANDMARKS\n"
                   "LMK0001 512.0 512.0\n"
                   "LMK0001 513.0 512.0\n"
                   "END FILE\n" );
    BOOST_CHECK_THROW( input_output::sum_lmk::readSumFile( duplicateSumPath.string( ) ), std::runtime_error );

    const boost::filesystem::path conflictingLmkPath = makeTemporaryPath( ".lmk" );
    writeTextFile( conflictingLmkPath,
                   "LMK0001 T\n"
                   "1.0 0.0 0.0 VLM\n"
                   "END FILE\n" );
    BOOST_CHECK_THROW( input_output::sum_lmk::readLmkFiles( { testDataPath + "/LMK0001.lmk", conflictingLmkPath.string( ) } ),
                       std::runtime_error );

    boost::filesystem::remove( duplicateSumPath );
    boost::filesystem::remove( conflictingLmkPath );
}

//! Parse a real SPC SUM file (W48230079013.SUM, Rosetta/OSIRIS) and assert the v1-used fields.
BOOST_AUTO_TEST_CASE( testReadRealSumFile )
{
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const input_output::sum_lmk::SumImageData image = input_output::sum_lmk::readSumFile( testDataPath + "/W48230079013.SUM" );

    BOOST_CHECK_EQUAL( image.imageId_, "W48230079013" );
    BOOST_CHECK_EQUAL( image.utcEpochString_, "2015 APR 14 16:26:27.974" );
    BOOST_CHECK_EQUAL( image.imageSize_( 0 ), 1024 );
    BOOST_CHECK_EQUAL( image.imageSize_( 1 ), 1024 );
    BOOST_CHECK_EQUAL( image.threshold_, 500 );
    BOOST_CHECK_EQUAL( image.maxDn_, 65535 );
    checkClose( image.focalLengthMm_, 135.68, 1.0E-9 );
    checkClose( image.opticalCenter_( 0 ), 511.7, 1.0E-9 );
    checkClose( image.opticalCenter_( 1 ), 510.6, 1.0E-9 );

    // SCOBJ km -> m.
    checkClose( image.spacecraftObjectVector_( 0 ), -143058.9240, 1.0E-4 );
    checkClose( image.spacecraftObjectVector_( 1 ), -38931.15914, 1.0E-4 );
    checkClose( image.spacecraftObjectVector_( 2 ), 84588.63431, 1.0E-4 );

    // Camera axes stored as rows CX/CY/CZ.
    checkClose( image.cameraAxes_( 0, 0 ), 0.5399383659, 1.0E-9 );
    checkClose( image.cameraAxes_( 0, 2 ), 0.8126088065, 1.0E-9 );
    checkClose( image.cameraAxes_( 1, 1 ), 0.9503508287, 1.0E-9 );
    checkClose( image.cameraAxes_( 2, 0 ), -0.8387783408, 1.0E-9 );
    checkClose( image.cameraAxes_( 2, 2 ), 0.4977460006, 1.0E-9 );

    // K-MATRIX (6 plain-decimal tokens, no Fortran exponent), row-major into 2x3.
    checkClose( image.kMatrix_( 0, 0 ), 37.0371, 1.0E-9 );
    checkClose( image.kMatrix_( 0, 1 ), 0.1240, 1.0E-9 );
    checkClose( image.kMatrix_( 0, 2 ), 0.0, 1.0E-9 );
    checkClose( image.kMatrix_( 1, 0 ), -0.1240, 1.0E-9 );
    checkClose( image.kMatrix_( 1, 1 ), 37.0371, 1.0E-9 );

    // 3-component SIGMA_PTG (radians) - distinct components, unlike the scalar synthetic fixture.
    checkClose( image.pointingSigma_( 0 ), 0.2231524849E-3, 1.0E-12 );
    checkClose( image.pointingSigma_( 2 ), 0.5515340338E-3, 1.0E-12 );
    BOOST_CHECK( std::fabs( image.pointingSigma_( 0 ) - image.pointingSigma_( 2 ) ) > 1.0E-9 );

    // SIGMA_VSO km -> m.
    checkClose( image.spacecraftObjectSigma_( 0 ), 39.07025850, 1.0E-6 );

    // DISTORTION parsed (zeros here), stored but unused in v1.
    checkClose( image.distortionCoefficients_( 0 ), 0.0, 1.0E-12 );

    // Three landmark observations with exact pixel values; empty LIMB FITS section.
    BOOST_REQUIRE_EQUAL( image.landmarkObservations_.size( ), 3 );
    BOOST_CHECK_EQUAL( image.landmarkObservations_.at( 0 ).landmarkId_, "FA0115" );
    checkClose( image.landmarkObservations_.at( 0 ).pixelCoordinates_( 0 ), 533.53, 1.0E-9 );
    checkClose( image.landmarkObservations_.at( 0 ).pixelCoordinates_( 1 ), 475.22, 1.0E-9 );
    BOOST_CHECK_EQUAL( image.landmarkObservations_.at( 2 ).landmarkId_, "GD0036" );
    BOOST_CHECK_EQUAL( image.limbFitObservations_.size( ), 0 );
}

//! Parse a real SPC LMK file (FA0115.LMK), including the combined multi-token labels
//! ("SIZE, SCALE(KM)" and "SIGKM, RMSLMK") used by production files.
BOOST_AUTO_TEST_CASE( testReadRealLmkFile )
{
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks = input_output::sum_lmk::readLmkFiles(
            { testDataPath + "/FA0115.LMK", testDataPath + "/FC0098.LMK", testDataPath + "/GD0036.LMK" } );
    BOOST_REQUIRE_EQUAL( landmarks.size( ), 3 );
    BOOST_REQUIRE( landmarks.count( "FA0115" ) == 1 );

    const input_output::sum_lmk::LmkLandmarkData& fa = landmarks.at( "FA0115" );
    BOOST_CHECK_EQUAL( fa.landmarkId_, "FA0115" );
    BOOST_CHECK_EQUAL( fa.typeFlag_, 'T' );

    // VLM km -> m.
    checkClose( fa.bodyFixedPosition_( 0 ), 2182.154009, 1.0E-6 );
    checkClose( fa.bodyFixedPosition_( 1 ), 479.8117875, 1.0E-6 );
    checkClose( fa.bodyFixedPosition_( 2 ), -424.4737920, 1.0E-6 );

    // Local axes (UX) and HORIZON flags.
    checkClose( fa.localXAxis_( 0 ), -0.5593237074, 1.0E-9 );
    BOOST_CHECK_EQUAL( fa.horizonFlags_( 0 ), -1 );
    BOOST_CHECK_EQUAL( fa.horizonFlags_( 3 ), -1 );

    // SIGMA_LMK km -> m (3 components).
    checkClose( fa.landmarkPositionSigma_( 0 ), 0.5522198847, 1.0E-7 );
    checkClose( fa.landmarkPositionSigma_( 2 ), 0.6491137634, 1.0E-7 );

    // Combined-label fields are now captured: "SIZE, SCALE(KM)" and "SIGKM, RMSLMK".
    BOOST_CHECK_EQUAL( fa.patchSize_, 49 );
    checkClose( fa.patchScaleKm_, 0.0100000, 1.0E-9 );
    checkClose( fa.sigKm_, 0.5000000000E-2, 1.0E-12 );
    checkClose( fa.rmsLmk_, 1.458559361, 1.0E-9 );

    // PICTURES section parsed.
    BOOST_REQUIRE_EQUAL( fa.pictures_.size( ), 3 );
    BOOST_CHECK_EQUAL( fa.pictures_.at( 0 ).imageId_, "N46299874022" );
    checkClose( fa.pictures_.at( 0 ).pixelCoordinates_( 0 ), 204.85, 1.0E-9 );
    BOOST_CHECK( !fa.pictures_.at( 0 ).flagged_ );
}

//! Missing required fields must throw with context.
BOOST_AUTO_TEST_CASE( testSumLmkRequiredFieldValidation )
{
    const std::string header =
            "IMGREQ\n"
            "2015 JUN 05 07:24:42.053\n"
            "1024 1024 500 65535 NPX, NLN, THRSH\n"
            "100.0 512.0 512.0 MMFL, CTR\n"
            "1.0 0.0 0.0 CX\n"
            "0.0 1.0 0.0 CY\n"
            "0.0 0.0 1.0 CZ\n";
    const std::string kMatrix = "10.0 0.0 0.0 0.0 10.0 0.0 K-MATRIX\n";
    const std::string landmarks = "LANDMARKS\nLMK0001 512.0 512.0\nEND FILE\n";

    // Missing K-MATRIX.
    {
        const boost::filesystem::path path = makeTemporaryPath( ".sum" );
        writeTextFile( path, header + landmarks );
        BOOST_CHECK_THROW( input_output::sum_lmk::readSumFile( path.string( ) ), std::runtime_error );
        boost::filesystem::remove( path );
    }
    // Missing CZ row.
    {
        std::string headerNoCz =
                "IMGREQ\n2015 JUN 05 07:24:42.053\n1024 1024 500 65535 NPX, NLN, THRSH\n"
                "100.0 512.0 512.0 MMFL, CTR\n1.0 0.0 0.0 CX\n0.0 1.0 0.0 CY\n";
        const boost::filesystem::path path = makeTemporaryPath( ".sum" );
        writeTextFile( path, headerNoCz + kMatrix + landmarks );
        BOOST_CHECK_THROW( input_output::sum_lmk::readSumFile( path.string( ) ), std::runtime_error );
        boost::filesystem::remove( path );
    }
    // No landmark rows.
    {
        const boost::filesystem::path path = makeTemporaryPath( ".sum" );
        writeTextFile( path, header + kMatrix + "LANDMARKS\nEND FILE\n" );
        BOOST_CHECK_THROW( input_output::sum_lmk::readSumFile( path.string( ) ), std::runtime_error );
        boost::filesystem::remove( path );
    }
    // LMK missing VLM.
    {
        const boost::filesystem::path path = makeTemporaryPath( ".lmk" );
        writeTextFile( path, "LMKNOVLM T\n0.1D-02 SIGMA_LMK\nEND FILE\n" );
        BOOST_CHECK_THROW( input_output::sum_lmk::readLmkFile( path.string( ) ), std::runtime_error );
        boost::filesystem::remove( path );
    }
}

//! Fortran D-exponent notation including negative mantissas must parse.
BOOST_AUTO_TEST_CASE( testFortranExponentParsing )
{
    const boost::filesystem::path path = makeTemporaryPath( ".lmk" );
    writeTextFile( path,
                   "LMKEXP T\n"
                   "-0.1545700894D+01 0.8640125542D+00 0.1009105347D+01 VLM\n"
                   "END FILE\n" );
    const input_output::sum_lmk::LmkLandmarkData landmark = input_output::sum_lmk::readLmkFile( path.string( ) );
    // km -> m, negative mantissa preserved.
    checkClose( landmark.bodyFixedPosition_( 0 ), -1545.700894, 1.0E-6 );
    checkClose( landmark.bodyFixedPosition_( 1 ), 864.0125542, 1.0E-6 );
    checkClose( landmark.bodyFixedPosition_( 2 ), 1009.105347, 1.0E-6 );
    boost::filesystem::remove( path );
}

//! Cross-file duplicate handling: identical duplicates accepted; conflicting (image, landmark) throws.
BOOST_AUTO_TEST_CASE( testSumCrossFileDuplicateHandling )
{
    const std::string baseSum =
            "IMGDUPX\n2015 JUN 05 07:24:42.053\n1024 1024 500 65535 NPX, NLN, THRSH\n"
            "100.0 512.0 512.0 MMFL, CTR\n1.0 0.0 0.0 CX\n0.0 1.0 0.0 CY\n0.0 0.0 1.0 CZ\n"
            "10.0 0.0 0.0 0.0 10.0 0.0 K-MATRIX\nLANDMARKS\n";

    // Identical (image, landmark) across two files: accepted (no double counting).
    {
        const boost::filesystem::path a = makeTemporaryPath( ".sum" );
        const boost::filesystem::path b = makeTemporaryPath( ".sum" );
        writeTextFile( a, baseSum + "LMK0001 512.0 512.0\nEND FILE\n" );
        writeTextFile( b, baseSum + "LMK0001 512.0 512.0\nEND FILE\n" );
        BOOST_CHECK_NO_THROW( input_output::sum_lmk::readSumFiles( { a.string( ), b.string( ) } ) );
        boost::filesystem::remove( a );
        boost::filesystem::remove( b );
    }
    // Conflicting pixel coordinates for the same (image, landmark) across files: throws.
    {
        const boost::filesystem::path a = makeTemporaryPath( ".sum" );
        const boost::filesystem::path b = makeTemporaryPath( ".sum" );
        writeTextFile( a, baseSum + "LMK0001 512.0 512.0\nEND FILE\n" );
        writeTextFile( b, baseSum + "LMK0001 600.0 512.0\nEND FILE\n" );
        BOOST_CHECK_THROW( input_output::sum_lmk::readSumFiles( { a.string( ), b.string( ) } ), std::runtime_error );
        boost::filesystem::remove( a );
        boost::filesystem::remove( b );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
