/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/io/basicInputOutput.h"
#include "tudat/io/readSumLmkFiles.h"
#include "tudat/support/testFileUtilities.h"

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

void writeTextFile( const std::filesystem::path& filePath, const std::string& contents )
{
    std::ofstream stream( filePath.string( ) );
    stream << contents;
}

std::filesystem::path makeTemporaryPath( const std::string& suffix )
{
    return std::filesystem::path( createTemporaryFilePath( "tudat_sum_lmk", suffix ) );
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
}

BOOST_AUTO_TEST_CASE( testSumLmkDuplicateFailures )
{
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";

    const std::filesystem::path duplicateSumPath = makeTemporaryPath( ".sum" );
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

    const std::filesystem::path conflictingLmkPath = makeTemporaryPath( ".lmk" );
    writeTextFile( conflictingLmkPath,
                   "LMK0001 T\n"
                   "1.0 0.0 0.0 VLM\n"
                   "END FILE\n" );
    BOOST_CHECK_THROW( input_output::sum_lmk::readLmkFiles( { testDataPath + "/LMK0001.lmk", conflictingLmkPath.string( ) } ),
                       std::runtime_error );

    std::filesystem::remove( duplicateSumPath );
    std::filesystem::remove( conflictingLmkPath );
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

    // DISTORTION is retained from the file. SUM/LMK conversion accepts the SPC-standard zero
    // row and rejects non-zero rows, whose separate Owen model is not encoded in a SUM file.
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
    const auto checkInvalidSum = []( const std::string& contents ) {
        const std::filesystem::path path = makeTemporaryPath( ".sum" );
        writeTextFile( path, contents );
        BOOST_CHECK_THROW( input_output::sum_lmk::readSumFile( path.string( ) ), std::runtime_error );
        std::filesystem::remove( path );
    };

    // Missing K-MATRIX.
    {
        const std::filesystem::path path = makeTemporaryPath( ".sum" );
        writeTextFile( path, header + landmarks );
        BOOST_CHECK_THROW( input_output::sum_lmk::readSumFile( path.string( ) ), std::runtime_error );
        std::filesystem::remove( path );
    }
    // Missing CZ row.
    {
        std::string headerNoCz =
                "IMGREQ\n2015 JUN 05 07:24:42.053\n1024 1024 500 65535 NPX, NLN, THRSH\n"
                "100.0 512.0 512.0 MMFL, CTR\n1.0 0.0 0.0 CX\n0.0 1.0 0.0 CY\n";
        const std::filesystem::path path = makeTemporaryPath( ".sum" );
        writeTextFile( path, headerNoCz + kMatrix + landmarks );
        BOOST_CHECK_THROW( input_output::sum_lmk::readSumFile( path.string( ) ), std::runtime_error );
        std::filesystem::remove( path );
    }
    // Semantically invalid calibration and attitude records must be rejected as well.
    checkInvalidSum( "IMGREQ\n2015 JUN 05 07:24:42.053\n0 1024 500 65535 NPX, NLN, THRSH\n" + header.substr( header.find( "100.0" ) ) +
                     kMatrix + landmarks );
    checkInvalidSum(
            "IMGREQ\n2015 JUN 05 07:24:42.053\n1024 1024 500 65535 NPX, NLN, THRSH\n"
            "0.0 512.0 512.0 MMFL, CTR\n"
            "1.0 0.0 0.0 CX\n0.0 1.0 0.0 CY\n0.0 0.0 1.0 CZ\n" +
            kMatrix + landmarks );
    checkInvalidSum(
            "IMGREQ\n2015 JUN 05 07:24:42.053\n1024 1024 500 65535 NPX, NLN, THRSH\n"
            "100.0 nan 512.0 MMFL, CTR\n"
            "1.0 0.0 0.0 CX\n0.0 1.0 0.0 CY\n0.0 0.0 1.0 CZ\n" +
            kMatrix + landmarks );
    checkInvalidSum(
            "IMGREQ\n2015 JUN 05 07:24:42.053\n1024 1024 500 65535 NPX, NLN, THRSH\n"
            "100.0 512.0 512.0 MMFL, CTR\n"
            "1.0 0.0 0.0 CX\n1.0 1.0 0.0 CY\n0.0 0.0 1.0 CZ\n" +
            kMatrix + landmarks );
    checkInvalidSum(
            "IMGREQ\n2015 JUN 05 07:24:42.053\n1024 1024 500 65535 NPX, NLN, THRSH\n"
            "100.0 512.0 512.0 MMFL, CTR\n"
            "1.0 0.0 0.0 CX\n0.0 -1.0 0.0 CY\n0.0 0.0 1.0 CZ\n" +
            kMatrix + landmarks );
    checkInvalidSum( header + "10.0 0.0 0.0 20.0 0.0 0.0 K-MATRIX\n" + landmarks );
    // Valid matrices in the full SPC archive have a Frobenius orthogonality error up to
    // approximately 7.9e-8 from stored numerical precision. This one has 7.8e-8.
    {
        const std::filesystem::path path = makeTemporaryPath( ".sum" );
        writeTextFile( path,
                       "IMGARCHIVEPRECISION\n2015 JUN 05 07:24:42.053\n1024 1024 500 65535 NPX, NLN, THRSH\n"
                       "100.0 512.0 512.0 MMFL, CTR\n1.0000000390 0.0 0.0 CX\n0.0 1.0 0.0 CY\n0.0 0.0 1.0 CZ\n" +
                               kMatrix + landmarks );
        BOOST_CHECK_NO_THROW( input_output::sum_lmk::readSumFile( path.string( ) ) );
        std::filesystem::remove( path );
    }
    // No landmark rows: a single SUM image parses, but the batch reader discards it.
    {
        const std::filesystem::path path = makeTemporaryPath( ".sum" );
        writeTextFile( path, header + kMatrix + "LANDMARKS\nEND FILE\n" );
        const input_output::sum_lmk::SumImageData image = input_output::sum_lmk::readSumFile( path.string( ) );
        BOOST_CHECK_EQUAL( image.landmarkObservations_.size( ), 0 );
        std::filesystem::remove( path );
    }
    // LMK missing VLM.
    {
        const std::filesystem::path path = makeTemporaryPath( ".lmk" );
        writeTextFile( path, "LMKNOVLM T\n0.1D-02 SIGMA_LMK\nEND FILE\n" );
        BOOST_CHECK_THROW( input_output::sum_lmk::readLmkFile( path.string( ) ), std::runtime_error );
        std::filesystem::remove( path );
    }
}

//! SUM images without landmark observations are valid parser inputs but are discarded in batches.
BOOST_AUTO_TEST_CASE( testSumFilesDiscardEmptyImages )
{
    const std::string baseSum =
            "IMGEMPTY\n2015 JUN 05 07:24:42.053\n1024 1024 500 65535 NPX, NLN, THRSH\n"
            "100.0 512.0 512.0 MMFL, CTR\n1.0 0.0 0.0 CX\n0.0 1.0 0.0 CY\n0.0 0.0 1.0 CZ\n"
            "10.0 0.0 0.0 0.0 10.0 0.0 K-MATRIX\n";

    {
        const std::filesystem::path emptyPath = makeTemporaryPath( ".sum" );
        const std::filesystem::path observedPath = makeTemporaryPath( ".sum" );
        writeTextFile( emptyPath, baseSum + "LANDMARKS\nEND FILE\n" );
        writeTextFile( observedPath,
                       "IMGOBS\n2015 JUN 05 07:24:42.053\n1024 1024 500 65535 NPX, NLN, THRSH\n"
                       "100.0 512.0 512.0 MMFL, CTR\n1.0 0.0 0.0 CX\n0.0 1.0 0.0 CY\n0.0 0.0 1.0 CZ\n"
                       "10.0 0.0 0.0 0.0 10.0 0.0 K-MATRIX\nLANDMARKS\nLMK0001 512.0 512.0\nEND FILE\n" );

        const std::vector< input_output::sum_lmk::SumImageData > images =
                input_output::sum_lmk::readSumFiles( { emptyPath.string( ), observedPath.string( ) } );
        BOOST_REQUIRE_EQUAL( images.size( ), 1 );
        BOOST_CHECK_EQUAL( images.at( 0 ).imageId_, "IMGOBS" );

        std::filesystem::remove( emptyPath );
        std::filesystem::remove( observedPath );
    }

    {
        const std::filesystem::path emptyPathA = makeTemporaryPath( ".sum" );
        const std::filesystem::path emptyPathB = makeTemporaryPath( ".sum" );
        writeTextFile( emptyPathA, baseSum + "END FILE\n" );
        writeTextFile( emptyPathB,
                       "IMGEMPTY2\n2015 JUN 05 07:24:42.053\n1024 1024 500 65535 NPX, NLN, THRSH\n"
                       "100.0 512.0 512.0 MMFL, CTR\n1.0 0.0 0.0 CX\n0.0 1.0 0.0 CY\n0.0 0.0 1.0 CZ\n"
                       "10.0 0.0 0.0 0.0 10.0 0.0 K-MATRIX\nLANDMARKS\nEND FILE\n" );

        const std::vector< input_output::sum_lmk::SumImageData > images =
                input_output::sum_lmk::readSumFiles( { emptyPathA.string( ), emptyPathB.string( ) } );
        BOOST_CHECK_EQUAL( images.size( ), 0 );

        std::filesystem::remove( emptyPathA );
        std::filesystem::remove( emptyPathB );
    }
}

//! Production SUM landmark rows may carry an optional trailing '-' flag.
BOOST_AUTO_TEST_CASE( testSumLandmarkRowsWithOptionalFlag )
{
    const std::filesystem::path path = makeTemporaryPath( ".sum" );
    writeTextFile( path,
                   "IMGFLAG\n"
                   "2015 JUN 05 07:24:42.053\n"
                   "1024 1024 500 65535 NPX, NLN, THRSH\n"
                   "100.0 512.0 512.0 MMFL, CTR\n"
                   "0.0 0.0 -10.0 SCOBJ\n"
                   "1.0 0.0 0.0 CX\n"
                   "0.0 1.0 0.0 CY\n"
                   "0.0 0.0 1.0 CZ\n"
                   "10.0 0.0 0.0 0.0 10.0 0.0 K-MATRIX\n"
                   "LANDMARKS\n"
                   "LMK0001 512.0 512.0 -\n"
                   "END FILE\n" );

    const input_output::sum_lmk::SumImageData image = input_output::sum_lmk::readSumFile( path.string( ) );
    BOOST_REQUIRE_EQUAL( image.landmarkObservations_.size( ), 1 );
    BOOST_CHECK_EQUAL( image.landmarkObservations_.at( 0 ).landmarkId_, "LMK0001" );
    checkClose( image.landmarkObservations_.at( 0 ).pixelCoordinates_( 0 ), 512.0 );
    checkClose( image.landmarkObservations_.at( 0 ).pixelCoordinates_( 1 ), 512.0 );
    std::filesystem::remove( path );
}

//! Fortran D-exponent notation including negative mantissas must parse.
BOOST_AUTO_TEST_CASE( testFortranExponentParsing )
{
    const std::filesystem::path path = makeTemporaryPath( ".lmk" );
    writeTextFile( path,
                   "LMKEXP T\n"
                   "-0.1545700894D+01 0.8640125542D+00 0.1009105347D+01 VLM\n"
                   "END FILE\n" );
    const input_output::sum_lmk::LmkLandmarkData landmark = input_output::sum_lmk::readLmkFile( path.string( ) );
    // km -> m, negative mantissa preserved.
    checkClose( landmark.bodyFixedPosition_( 0 ), -1545.700894, 1.0E-6 );
    checkClose( landmark.bodyFixedPosition_( 1 ), 864.0125542, 1.0E-6 );
    checkClose( landmark.bodyFixedPosition_( 2 ), 1009.105347, 1.0E-6 );
    std::filesystem::remove( path );
}

//! Cross-file duplicate handling: identical duplicates accepted; conflicting (image, landmark) throws.
BOOST_AUTO_TEST_CASE( testSumCrossFileDuplicateHandling )
{
    const std::string baseMetadata =
            "IMGDUPX\n2015 JUN 05 07:24:42.053\n1024 1024 500 65535 NPX, NLN, THRSH\n"
            "100.0 512.0 512.0 MMFL, CTR\n1.0 0.0 0.0 CX\n0.0 1.0 0.0 CY\n0.0 0.0 1.0 CZ\n"
            "10.0 0.0 0.0 0.0 10.0 0.0 K-MATRIX\n";
    const std::string baseSum = baseMetadata + "LANDMARKS\n";

    // Identical (image, landmark) across two files: accepted (no double counting).
    {
        const std::filesystem::path a = makeTemporaryPath( ".sum" );
        const std::filesystem::path b = makeTemporaryPath( ".sum" );
        writeTextFile( a, baseSum + "LMK0001 512.0 512.0\nEND FILE\n" );
        writeTextFile( b, baseSum + "LMK0001 512.0 512.0\nEND FILE\n" );
        const auto images = input_output::sum_lmk::readSumFiles( { a.string( ), b.string( ), a.string( ) } );
        BOOST_REQUIRE_EQUAL( images.size( ), 1 );
        BOOST_CHECK_EQUAL( images.front( ).landmarkObservations_.size( ), 1 );
        std::filesystem::remove( a );
        std::filesystem::remove( b );
    }
    // Conflicting pixel coordinates for the same (image, landmark) across files: throws.
    {
        const std::filesystem::path a = makeTemporaryPath( ".sum" );
        const std::filesystem::path b = makeTemporaryPath( ".sum" );
        writeTextFile( a, baseSum + "LMK0001 512.0 512.0\nEND FILE\n" );
        writeTextFile( b, baseSum + "LMK0001 600.0 512.0\nEND FILE\n" );
        BOOST_CHECK_THROW( input_output::sum_lmk::readSumFiles( { a.string( ), b.string( ) } ), std::runtime_error );
        std::filesystem::remove( a );
        std::filesystem::remove( b );
    }

    // All optional image metadata that affects conversion or its a-priori must also agree.
    const std::vector< std::pair< std::string, std::string > > conflictingMetadataRows = {
        { "0.0 0.0 -10.0 SCOBJ\n", "0.0 0.0 -11.0 SCOBJ\n" },
        { "1.0 0.0 0.0 SZ\n", "0.0 1.0 0.0 SZ\n" },
        { "0.0 0.0 0.0 0.0 DISTORTION\n", "1.0 0.0 0.0 0.0 DISTORTION\n" },
        { "0.1 SIGMA_VSO\n", "0.2 SIGMA_VSO\n" },
        { "1.0E-4 SIGMA_PTG\n", "2.0E-4 SIGMA_PTG\n" }
    };
    for( const auto& conflictingRows : conflictingMetadataRows )
    {
        const std::filesystem::path a = makeTemporaryPath( ".sum" );
        const std::filesystem::path b = makeTemporaryPath( ".sum" );
        writeTextFile( a, baseMetadata + conflictingRows.first + "LANDMARKS\nLMK0001 512.0 512.0\nEND FILE\n" );
        writeTextFile( b, baseMetadata + conflictingRows.second + "LANDMARKS\nLMK0002 513.0 512.0\nEND FILE\n" );
        BOOST_CHECK_THROW( input_output::sum_lmk::readSumFiles( { a.string( ), b.string( ) } ), std::runtime_error );
        std::filesystem::remove( a );
        std::filesystem::remove( b );
    }
}

//! Tracking data and camera/landmark metadata are available without constructing a body system.
BOOST_AUTO_TEST_CASE( testSumLmkTrackingDataWithoutEnvironment )
{
    using namespace input_output::sum_lmk;
    const std::string path = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const std::vector< std::string > sumFiles = { path + "/sample.sum" };
    const std::vector< std::string > lmkFiles = { path + "/LMK0001.lmk", path + "/LMK0002.lmk" };
    const auto images = readSumFiles( sumFiles );
    const auto landmarks = readLmkFiles( lmkFiles );
    const auto fromFiles = readSumLmkData< double, double >( sumFiles, lmkFiles, "Target", "Spacecraft" );
    const auto fromRaw = convertSumLmkData< double, double >( images, landmarks, "Target", "Spacecraft" );
    BOOST_REQUIRE_EQUAL( fromFiles.first.size( ), 2 );
    BOOST_REQUIRE_EQUAL( fromRaw.first.size( ), fromFiles.first.size( ) );
    for( unsigned int i = 0; i < fromFiles.first.size( ); ++i )
    {
        const auto entry = fromFiles.first.at( i );
        BOOST_CHECK_EQUAL( entry->getObservableType( ), "PixelCoordinates" );
        BOOST_CHECK_EQUAL( entry->getTimeScale( ), "UTC" );
        BOOST_CHECK_EQUAL( entry->getReferenceLinkEnd( ), "receiver" );
        BOOST_CHECK_EQUAL( entry->getNumberOfObservations( ), 1 );
        BOOST_CHECK_SMALL( ( entry->getObservationsVector( ) - fromRaw.first.at( i )->getObservationsVector( ) ).norm( ), 1.0E-14 );
        BOOST_CHECK_EQUAL( entry->getObservationEpochs( ).front( ), fromRaw.first.at( i )->getObservationEpochs( ).front( ) );
        BOOST_CHECK( entry->getLinkEnds( ) == fromRaw.first.at( i )->getLinkEnds( ) );
    }
    BOOST_REQUIRE_EQUAL( fromFiles.second.size( ), 3 );
    BOOST_REQUIRE_EQUAL( fromRaw.second.size( ), fromFiles.second.size( ) );
    unsigned int cameraCount = 0, landmarkCount = 0;
    for( const auto& entry : fromFiles.second )
    {
        if( entry->getReferencePointSupplementaryData( ) != nullptr )
        {
            ++landmarkCount;
            BOOST_CHECK_EQUAL( entry->getBodyName( ), "Target" );
            BOOST_CHECK_SMALL( ( entry->getReferencePointSupplementaryData( )->getBodyFixedPosition( ) -
                                 landmarks.at( entry->getReferencePointName( ) ).bodyFixedPosition_ )
                                       .norm( ),
                               1.0E-12 );
        }
        for( const auto& instrument : entry->getInstrumentSupplementaryData( ) )
        {
            const auto camera = std::dynamic_pointer_cast< data::CameraInstrumentSupplementaryData >( instrument );
            BOOST_REQUIRE( camera != nullptr );
            ++cameraCount;
            BOOST_CHECK_EQUAL( entry->getBodyName( ), "Spacecraft" );
            BOOST_CHECK_EQUAL( camera->getCameraId( ), "Camera_IMG001" );
            BOOST_CHECK_EQUAL( camera->getPointingReferenceBodyName( ), "Target" );
            BOOST_CHECK_SMALL( ( camera->getPointingUncertainty( ) - Eigen::Vector3d::Constant( 1.0E-4 ) ).norm( ), 1.0E-15 );
            BOOST_CHECK_SMALL( ( camera->getKMatrix( ) - images.front( ).kMatrix_ ).norm( ), 1.0E-12 );
        }
    }
    BOOST_CHECK_EQUAL( cameraCount, 1 );
    BOOST_CHECK_EQUAL( landmarkCount, 2 );
}

BOOST_AUTO_TEST_CASE( testIdenticalSumLmkDataIsLoadedOnce )
{
    using namespace input_output::sum_lmk;
    const std::string path = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const auto images = readSumFiles( { path + "/sample.sum", path + "/sample.sum" } );
    BOOST_REQUIRE_EQUAL( images.size( ), 1 );
    BOOST_CHECK_EQUAL( images.front( ).landmarkObservations_.size( ), 2 );
    BOOST_CHECK_EQUAL( images.front( ).limbFitObservations_.size( ), 1 );
    const auto landmarks = readLmkFiles( { path + "/LMK0001.lmk", path + "/LMK0002.lmk", path + "/LMK0001.lmk" } );
    BOOST_CHECK_EQUAL( landmarks.size( ), 2 );
    auto repeated = images.front( );
    repeated.landmarkObservations_.push_back( repeated.landmarkObservations_.front( ) );
    const auto loaded = convertSumLmkData< double, double >( { repeated, repeated }, landmarks, "Target", "Spacecraft" );
    BOOST_CHECK_EQUAL( loaded.first.size( ), 2 );
    BOOST_CHECK_EQUAL( loaded.second.size( ), 3 );
    repeated.landmarkObservations_.back( ).pixelCoordinates_( 0 ) += 1.0;
    BOOST_CHECK_THROW( ( convertSumLmkData< double, double >( { repeated }, landmarks, "Target", "Spacecraft" ) ), std::runtime_error );

    const auto minimalLmk = makeTemporaryPath( ".lmk" );
    writeTextFile( minimalLmk, "LMK0001 T\n0.0 0.0 0.0 VLM\nEND FILE\n" );
    const auto minimalLandmarks = readLmkFiles( { minimalLmk.string( ), minimalLmk.string( ) } );
    BOOST_REQUIRE_EQUAL( minimalLandmarks.size( ), 1 );
    BOOST_CHECK( minimalLandmarks.begin( )->second.localXAxis_.array( ).isNaN( ).all( ) );
    std::filesystem::remove( minimalLmk );

    const auto repeatedRows = makeTemporaryPath( ".sum" );
    writeTextFile( repeatedRows,
                   "IMGROWS\n2015 JUN 05 07:24:42.053\n1024 1024 500 65535 NPX, NLN, THRSH\n"
                   "100.0 512.0 512.0 MMFL, CTR\n1.0 0.0 0.0 CX\n0.0 1.0 0.0 CY\n0.0 0.0 1.0 CZ\n"
                   "10.0 0.0 0.0 0.0 10.0 0.0 K-MATRIX\nLANDMARKS\n"
                   "LMK0001 512.0 512.0\nLMK0001 512.0 512.0\nEND FILE\n" );
    BOOST_CHECK_EQUAL( readSumFile( repeatedRows.string( ) ).landmarkObservations_.size( ), 1 );
    std::filesystem::remove( repeatedRows );
}

BOOST_AUTO_TEST_CASE( testNonfiniteSumMeasurementsAndPointingUncertainties )
{
    using namespace input_output::sum_lmk;
    const std::string path = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const auto validImage = readSumFile( path + "/sample.sum" );
    const auto landmarks = readLmkFiles( { path + "/LMK0001.lmk", path + "/LMK0002.lmk" } );
    for( const double invalid : { std::numeric_limits< double >::infinity( ), -std::numeric_limits< double >::infinity( ), TUDAT_NAN } )
    {
        auto image = validImage;
        image.landmarkObservations_.front( ).pixelCoordinates_( 0 ) = invalid;
        BOOST_CHECK_THROW( ( convertSumLmkData< double, double >( { image }, landmarks, "Target", "Spacecraft" ) ), std::runtime_error );
    }
    for( const Eigen::Vector3d& sigma :
         std::vector< Eigen::Vector3d >{ Eigen::Vector3d::Constant( std::numeric_limits< double >::infinity( ) ),
                                         Eigen::Vector3d( 1.0E-4, TUDAT_NAN, TUDAT_NAN ),
                                         Eigen::Vector3d::Zero( ),
                                         Eigen::Vector3d::Constant( -1.0E-4 ) } )
    {
        auto image = validImage;
        image.pointingSigma_ = sigma;
        BOOST_CHECK_THROW( ( convertSumLmkData< double, double >( { image }, landmarks, "Target", "Spacecraft" ) ), std::runtime_error );
    }
    auto image = validImage;
    image.pointingSigma_ = Eigen::Vector3d::Constant( TUDAT_NAN );
    BOOST_CHECK_NO_THROW( ( convertSumLmkData< double, double >( { image }, landmarks, "Target", "Spacecraft" ) ) );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
