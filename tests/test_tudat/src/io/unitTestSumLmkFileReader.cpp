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
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/readSumLmkFiles.h"

namespace tudat
{
namespace unit_tests
{

using namespace input_output::sum_lmk;

BOOST_AUTO_TEST_SUITE( test_sum_lmk_file_reader )

//! Read the synthetic SUM/LMK files, including the LIMB FITS section and a single-valued SIGMA_PTG.
BOOST_AUTO_TEST_CASE( testReadSumAndLmkFiles )
{
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const SumImageData image = readSumFile( testDataPath + "/sample.sum" );

    BOOST_CHECK_EQUAL( image.imageId_, "IMG001" );
    BOOST_CHECK_EQUAL( image.utcEpochString_, "2015 JUN 05 07:24:42.053" );
    BOOST_CHECK_EQUAL( image.imageSize_( 0 ), 1024 );
    BOOST_CHECK_EQUAL( image.imageSize_( 1 ), 1024 );
    BOOST_CHECK_EQUAL( image.threshold_, 500 );
    BOOST_CHECK_EQUAL( image.maxDn_, 65535 );
    BOOST_CHECK_SMALL( image.focalLengthMm_ - 100.0, 1.0E-12 );
    BOOST_CHECK_SMALL( ( image.opticalCenter_ - Eigen::Vector2d( 512.0, 512.0 ) ).norm( ), 1.0E-12 );
    BOOST_CHECK_SMALL( ( image.spacecraftObjectVector_ - Eigen::Vector3d( 0.0, 0.0, 10000.0 ) ).norm( ), 1.0E-9 );
    BOOST_CHECK_SMALL( ( image.cameraAxes_ - Eigen::Matrix3d::Identity( ) ).norm( ), 1.0E-12 );
    BOOST_CHECK_SMALL( image.kMatrix_( 0, 0 ) - 10.0, 1.0E-12 );
    BOOST_CHECK_SMALL( image.kMatrix_( 1, 1 ) - 10.0, 1.0E-12 );
    BOOST_CHECK_SMALL( ( image.pointingSigma_ - Eigen::Vector3d::Constant( 1.0E-4 ) ).norm( ), 1.0E-18 );
    BOOST_REQUIRE_EQUAL( image.landmarkObservations_.size( ), 2 );
    BOOST_CHECK_EQUAL( image.landmarkObservations_.at( 1 ).landmarkId_, "LMK0002" );
    BOOST_CHECK_SMALL( ( image.landmarkObservations_.at( 1 ).pixelCoordinates_ - Eigen::Vector2d( 612.0, 512.0 ) ).norm( ), 1.0E-12 );
    BOOST_REQUIRE_EQUAL( image.limbFitObservations_.size( ), 1 );
    BOOST_CHECK_EQUAL( image.limbFitObservations_.at( 0 ).featureId_, "LIMB01" );
    BOOST_CHECK_SMALL( image.limbFitObservations_.at( 0 ).sigma_ - 2.0, 1.0E-12 );

    const std::map< std::string, LmkLandmarkData > landmarks =
            readLmkFiles( { testDataPath + "/LMK0001.lmk", testDataPath + "/LMK0002.lmk" } );
    BOOST_REQUIRE_EQUAL( landmarks.size( ), 2 );
    BOOST_CHECK_SMALL( landmarks.at( "LMK0002" ).bodyFixedPosition_( 0 ) - 1000.0, 1.0E-9 );
    BOOST_CHECK_SMALL( landmarks.at( "LMK0001" ).landmarkPositionSigma_( 0 ) - 1.0, 1.0E-12 );
    BOOST_REQUIRE_EQUAL( landmarks.at( "LMK0001" ).pictures_.size( ), 1 );
    BOOST_CHECK( landmarks.at( "LMK0001" ).pictures_.at( 0 ).flagged_ );
}

//! Read a real SPC SUM file (Rosetta/OSIRIS), with Fortran D exponents, a three-valued SIGMA_PTG and km -> m conversion.
BOOST_AUTO_TEST_CASE( testReadRealSumFile )
{
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const SumImageData image = readSumFile( testDataPath + "/W48230079013.SUM" );

    BOOST_CHECK_EQUAL( image.imageId_, "W48230079013" );
    BOOST_CHECK_EQUAL( image.utcEpochString_, "2015 APR 14 16:26:27.974" );
    BOOST_CHECK_EQUAL( image.imageSize_( 0 ), 1024 );
    BOOST_CHECK_EQUAL( image.imageSize_( 1 ), 1024 );
    BOOST_CHECK_SMALL( image.focalLengthMm_ - 135.68, 1.0E-9 );
    BOOST_CHECK_SMALL( image.opticalCenter_( 0 ) - 511.7, 1.0E-9 );
    BOOST_CHECK_SMALL( image.opticalCenter_( 1 ) - 510.6, 1.0E-9 );

    // SCOBJ, converted from km to m.
    BOOST_CHECK_SMALL( image.spacecraftObjectVector_( 0 ) - -143058.9240, 1.0E-4 );
    BOOST_CHECK_SMALL( image.spacecraftObjectVector_( 1 ) - -38931.15914, 1.0E-4 );
    BOOST_CHECK_SMALL( image.spacecraftObjectVector_( 2 ) - 84588.63431, 1.0E-4 );

    // CX/CY/CZ are stored as rows.
    BOOST_CHECK_SMALL( image.cameraAxes_( 0, 0 ) - 0.5399383659, 1.0E-12 );
    BOOST_CHECK_SMALL( image.cameraAxes_( 0, 2 ) - 0.8126088065, 1.0E-12 );
    BOOST_CHECK_SMALL( image.cameraAxes_( 1, 1 ) - 0.9503508287, 1.0E-12 );
    BOOST_CHECK_SMALL( image.cameraAxes_( 2, 0 ) - -0.8387783408, 1.0E-12 );
    BOOST_CHECK_SMALL( image.cameraAxes_( 2, 2 ) - 0.4977460006, 1.0E-12 );
    BOOST_CHECK_SMALL( image.sunDirectionBodyFixed_( 1 ) - 0.9750955554, 1.0E-12 );

    // K-MATRIX, row-major.
    BOOST_CHECK_SMALL( image.kMatrix_( 0, 0 ) - 37.0371, 1.0E-12 );
    BOOST_CHECK_SMALL( image.kMatrix_( 0, 1 ) - 0.1240, 1.0E-12 );
    BOOST_CHECK_SMALL( image.kMatrix_( 0, 2 ), 1.0E-12 );
    BOOST_CHECK_SMALL( image.kMatrix_( 1, 0 ) - -0.1240, 1.0E-12 );
    BOOST_CHECK_SMALL( image.kMatrix_( 1, 1 ) - 37.0371, 1.0E-12 );

    BOOST_CHECK_SMALL( image.pointingSigma_( 0 ) - 0.2231524849E-3, 1.0E-15 );
    BOOST_CHECK_SMALL( image.pointingSigma_( 2 ) - 0.5515340338E-3, 1.0E-15 );
    BOOST_CHECK_SMALL( image.spacecraftObjectSigma_( 0 ) - 39.07025850, 1.0E-7 );
    BOOST_CHECK_SMALL( image.distortionCoefficients_.norm( ), 1.0E-15 );

    BOOST_REQUIRE_EQUAL( image.landmarkObservations_.size( ), 3 );
    BOOST_CHECK_EQUAL( image.landmarkObservations_.at( 0 ).landmarkId_, "FA0115" );
    BOOST_CHECK_SMALL( ( image.landmarkObservations_.at( 0 ).pixelCoordinates_ - Eigen::Vector2d( 533.53, 475.22 ) ).norm( ), 1.0E-12 );
    BOOST_CHECK_EQUAL( image.landmarkObservations_.at( 2 ).landmarkId_, "GD0036" );
    BOOST_CHECK_EQUAL( image.limbFitObservations_.size( ), 0 );
}

//! Read real SPC LMK files, including the combined labels ("SIZE, SCALE(KM)" and "SIGKM, RMSLMK") of production files.
BOOST_AUTO_TEST_CASE( testReadRealLmkFile )
{
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const std::map< std::string, LmkLandmarkData > landmarks =
            readLmkFiles( { testDataPath + "/FA0115.LMK", testDataPath + "/FC0098.LMK", testDataPath + "/GD0036.LMK" } );
    BOOST_REQUIRE_EQUAL( landmarks.size( ), 3 );

    const LmkLandmarkData& landmark = landmarks.at( "FA0115" );
    BOOST_CHECK_EQUAL( landmark.landmarkId_, "FA0115" );
    BOOST_CHECK_EQUAL( landmark.typeFlag_, 'T' );

    // VLM and SIGMA_LMK, converted from km to m.
    BOOST_CHECK_SMALL( ( landmark.bodyFixedPosition_ - Eigen::Vector3d( 2182.154009, 479.8117875, -424.4737920 ) ).norm( ), 1.0E-9 );
    BOOST_CHECK_SMALL( landmark.landmarkPositionSigma_( 0 ) - 0.5522198847, 1.0E-10 );
    BOOST_CHECK_SMALL( landmark.landmarkPositionSigma_( 2 ) - 0.6491137634, 1.0E-10 );

    BOOST_CHECK_SMALL( landmark.localXAxis_( 0 ) - -0.5593237074, 1.0E-12 );
    BOOST_CHECK_SMALL( landmark.localZAxis_( 2 ) - -0.6904040917, 1.0E-12 );
    BOOST_CHECK_EQUAL( landmark.horizonFlags_( 0 ), -1 );
    BOOST_CHECK_EQUAL( landmark.patchSize_, 49 );
    BOOST_CHECK_SMALL( landmark.patchScaleKm_ - 0.01, 1.0E-12 );
    BOOST_CHECK_SMALL( landmark.sigKm_ - 0.5E-2, 1.0E-15 );
    BOOST_CHECK_SMALL( landmark.rmsLmk_ - 1.458559361, 1.0E-12 );

    BOOST_REQUIRE_EQUAL( landmark.pictures_.size( ), 3 );
    BOOST_CHECK_EQUAL( landmark.pictures_.at( 0 ).imageId_, "N46299874022" );
    BOOST_CHECK_SMALL( landmark.pictures_.at( 1 ).pixelCoordinates_( 1 ) - -481.78, 1.0E-12 );
    BOOST_CHECK( !landmark.pictures_.at( 0 ).flagged_ );
}

//! Conversion of SUM/LMK contents to tracking data and supplementary data: one observation per (image, landmark), with the
//! measured pixels, the UTC image epoch and the landmark/camera link ends, and a camera carrying the SUM calibration, the
//! CX/CY/CZ attitude (relative to the target body-fixed frame) and SIGMA_PTG.
BOOST_AUTO_TEST_CASE( testSumLmkConversion )
{
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const std::vector< std::string > sumFiles = { testDataPath + "/W48230079013.SUM" };
    const std::vector< std::string > lmkFiles = { testDataPath + "/FA0115.LMK",
                                                  testDataPath + "/FC0098.LMK",
                                                  testDataPath + "/GD0036.LMK" };
    const SumImageData image = readSumFile( sumFiles.at( 0 ) );
    const std::map< std::string, LmkLandmarkData > landmarks = readLmkFiles( lmkFiles );

    const auto convertedData = readSumLmkData< double, double >( sumFiles, lmkFiles, "Comet", "Spacecraft", "OSIRIS_WAC" );

    // Calendar epoch of the image, computed independently of the epoch string parser.
    const double expectedEpoch = basic_astrodynamics::DateTime( 2015, 4, 14, 16, 26, 27.974 ).epoch< double >( );
    const std::string cameraName = "OSIRIS_WAC";

    BOOST_REQUIRE_EQUAL( convertedData.first.size( ), 3 );
    for( unsigned int i = 0; i < convertedData.first.size( ); ++i )
    {
        const auto trackingData = convertedData.first.at( i );
        const SumLandmarkObservation& observation = image.landmarkObservations_.at( i );
        BOOST_CHECK_EQUAL( trackingData->getObservableType( ), "PixelCoordinates" );
        BOOST_CHECK_EQUAL( trackingData->getTimeScale( ), "UTC" );
        BOOST_CHECK_EQUAL( trackingData->getReferenceLinkEnd( ), "receiver" );
        BOOST_REQUIRE_EQUAL( trackingData->getNumberOfObservations( ), 1 );
        BOOST_CHECK_SMALL( trackingData->getObservationEpochs( ).at( 0 ) - expectedEpoch, 1.0E-6 );
        BOOST_CHECK_SMALL( ( trackingData->getObservations( ).at( 0 ) - observation.pixelCoordinates_ ).norm( ), 1.0E-12 );

        const data::PlainLinkDefinition& linkEnds = trackingData->getLinkEnds( );
        BOOST_REQUIRE_EQUAL( linkEnds.size( ), 2 );
        BOOST_CHECK( linkEnds.at( 0 ) ==
                     std::make_pair( std::make_pair( std::string( "Comet" ), observation.landmarkId_ ), std::string( "transmitter" ) ) );
        BOOST_CHECK( linkEnds.at( 1 ) ==
                     std::make_pair( std::make_pair( std::string( "Spacecraft" ), cameraName ), std::string( "receiver" ) ) );
    }

    unsigned int numberOfLandmarks = 0, numberOfCameras = 0;
    for( const auto& supplementaryData : convertedData.second )
    {
        if( supplementaryData->getReferencePointSupplementaryData( ) != nullptr )
        {
            ++numberOfLandmarks;
            BOOST_CHECK_EQUAL( supplementaryData->getBodyName( ), "Comet" );
            BOOST_CHECK_SMALL( ( supplementaryData->getReferencePointSupplementaryData( )->getBodyFixedPosition( ) -
                                 landmarks.at( supplementaryData->getReferencePointName( ) ).bodyFixedPosition_ )
                                       .norm( ),
                               1.0E-12 );
        }
        for( const auto& instrument : supplementaryData->getInstrumentSupplementaryData( ) )
        {
            const auto camera = std::dynamic_pointer_cast< data::CameraInstrumentSupplementaryData >( instrument );
            BOOST_REQUIRE( camera != nullptr );
            ++numberOfCameras;
            BOOST_CHECK_EQUAL( supplementaryData->getBodyName( ), "Spacecraft" );
            BOOST_CHECK_EQUAL( camera->getCameraId( ), cameraName );
            BOOST_CHECK_EQUAL( camera->getPointingReferenceBodyName( ), "Comet" );
            BOOST_CHECK_SMALL( camera->getFocalLength( ) - image.focalLengthMm_, 1.0E-12 );
            BOOST_CHECK_SMALL( ( camera->getPrincipalPoint( ) - image.opticalCenter_ ).norm( ), 1.0E-12 );
            BOOST_CHECK_SMALL( ( camera->getKMatrix( ) - image.kMatrix_ ).norm( ), 1.0E-12 );
            BOOST_REQUIRE_EQUAL( camera->getPointingUncertaintyHistory( ).size( ), 1 );
            BOOST_CHECK_SMALL( ( camera->getPointingUncertaintyHistory( ).begin( )->second - image.pointingSigma_ ).norm( ), 1.0E-15 );

            const std::map< double, Eigen::Quaterniond >& pointingHistory = camera->getRotationFromInertialToCameraFrameHistory( );
            BOOST_REQUIRE_EQUAL( pointingHistory.size( ), 1 );
            BOOST_CHECK_SMALL( pointingHistory.begin( )->first - expectedEpoch, 1.0E-6 );
            BOOST_CHECK_SMALL( ( pointingHistory.begin( )->second.toRotationMatrix( ) - image.cameraAxes_ ).norm( ), 1.0E-7 );
            BOOST_CHECK( camera->getPictureIdHistory( ) ==
                         ( std::map< double, std::string >{ { pointingHistory.begin( )->first, image.imageId_ } } ) );
        }
    }
    BOOST_CHECK_EQUAL( numberOfLandmarks, 3 );
    BOOST_CHECK_EQUAL( numberOfCameras, 1 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
