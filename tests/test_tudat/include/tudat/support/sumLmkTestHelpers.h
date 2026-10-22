/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SUM_LMK_TEST_HELPERS_H
#define TUDAT_SUM_LMK_TEST_HELPERS_H

#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/io/readSumLmkFiles.h"
#include "tudat/simulation/estimation_setup/createObservationDataset.h"
#include "tudat/simulation/estimation_setup/createInverseAprioriCovariance.h"
#include "tudat/simulation/estimation_setup/estimatableParameterSettings.h"
#include "tudat/simulation/estimation_setup/createObservationModelSettings.h"

namespace tudat
{
namespace unit_tests
{

//! Name of the camera that takes the SUM images in the tests.
const std::string sumLmkTestCameraName = "Camera";

//! Shared fixture setup; all conversion and environment changes use the public tracking-data interface.
struct LoadedPixelTestData {
    std::shared_ptr< observation_models::ObservationDataset< double, double > > observationDataset_;
    std::vector< std::shared_ptr< data::TrackingSupplementaryData > > supplementaryData_;
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettings_;
    std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >
            inverseAprioriCovarianceDiagonalEntries_;
};

//! Target at the origin and spacecraft on its -Z axis, both with identity rotation.
inline simulation_setup::SystemOfBodies createSyntheticSumLmkBodies( const double spacecraftZ = -10000.0 )
{
    simulation_setup::SystemOfBodies bodies( "SSB", "J2000" );
    bodies.createEmptyBody< double, double >( "Target", false );
    bodies.createEmptyBody< double, double >( "Spacecraft", false );
    bodies.at( "Target" )->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "J2000" ) );
    Eigen::Vector6d spacecraftState = Eigen::Vector6d::Zero( );
    spacecraftState( 2 ) = spacecraftZ;
    bodies.at( "Spacecraft" )->setEphemeris( std::make_shared< ephemerides::ConstantEphemeris >( spacecraftState, "SSB", "J2000" ) );
    bodies.at( "Target" )
            ->setRotationalEphemeris( std::make_shared< ephemerides::ConstantRotationalEphemeris >(
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "Target_Fixed" ) );
    bodies.at( "Spacecraft" )
            ->setRotationalEphemeris( std::make_shared< ephemerides::ConstantRotationalEphemeris >(
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "Spacecraft_Fixed" ) );
    bodies.processBodyFrameDefinitions< double, double >( );
    return bodies;
}

//! Minimal SUM image for the bodies of createSyntheticSumLmkBodies: identity camera axes (+Z boresight towards the target),
//! observing the given landmarks, with the measured pixels at the optical centre.
inline input_output::sum_lmk::SumImageData makeSyntheticSumImage( const std::string& imageId,
                                                                  const std::vector< std::string >& landmarkIds = { "LMK0001" },
                                                                  const double spacecraftZ = -10000.0 )
{
    input_output::sum_lmk::SumImageData image;
    image.imageId_ = imageId;
    image.utcEpochString_ = "2015 JUN 05 07:24:42.053";
    image.imageSize_ = Eigen::Vector2i( 1024, 1024 );
    image.focalLengthMm_ = 100.0;
    image.opticalCenter_ = Eigen::Vector2d( 512.0, 512.0 );
    // SCOBJ = spacecraft-to-object (target at origin), i.e. the negation of the spacecraft position.
    image.spacecraftObjectVector_ = Eigen::Vector3d( 0.0, 0.0, -spacecraftZ );
    image.cameraAxes_ = Eigen::Matrix3d::Identity( );
    image.kMatrix_ << 10.0, 0.0, 0.0, 0.0, 10.0, 0.0;
    image.pointingSigma_ = Eigen::Vector3d::Constant( 1.0E-4 );
    for( const std::string& landmarkId : landmarkIds )
    {
        input_output::sum_lmk::SumLandmarkObservation observation;
        observation.landmarkId_ = landmarkId;
        observation.pixelCoordinates_ = Eigen::Vector2d( 512.0, 512.0 );
        image.landmarkObservations_.push_back( observation );
    }
    return image;
}

inline std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > pixelModelSettingsForTest(
        const std::shared_ptr< observation_models::ObservationDataset< double, double > >& dataset )
{
    const std::map< observation_models::ObservableType, std::vector< observation_models::LinkEnds > > linkEndsPerObservable =
            dataset->getLinkEndsPerObservableType( );
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > settings;
    for( const auto& link : linkEndsPerObservable.at( observation_models::pixel_coordinates ) )
    {
        settings.push_back( observation_models::pixelCoordinatesSettings( link ) );
    }
    return settings;
}

inline LoadedPixelTestData loadPixelTestData( const std::vector< input_output::sum_lmk::SumImageData >& images,
                                              const std::map< std::string, input_output::sum_lmk::LmkLandmarkData >& landmarks,
                                              simulation_setup::SystemOfBodies& bodies,
                                              const std::string& targetBody = "Target",
                                              const std::string& receiverBody = "Spacecraft",
                                              const bool skipMissingLandmarks = false )
{
    const auto loaded = input_output::sum_lmk::convertSumLmkData< double, double >(
            { { sumLmkTestCameraName, images } }, landmarks, targetBody, receiverBody, skipMissingLandmarks );
    observation_models::setTrackingSupplementaryDataInBodies( bodies, loaded.second );
    LoadedPixelTestData result;
    result.supplementaryData_ = loaded.second;
    result.observationDataset_ = observation_models::createObservationDatasetFromTrackingData< double, double >( loaded.first, bodies );
    result.observationModelSettings_ = pixelModelSettingsForTest( result.observationDataset_ );
    result.inverseAprioriCovarianceDiagonalEntries_ = simulation_setup::getCameraPointingInverseAprioriCovarianceEntries( loaded.second );
    return result;
}

inline double sumImageTdbTimeForTest( const input_output::sum_lmk::SumImageData& image )
{
    return earth_orientation::createDefaultTimeConverter( )->getCurrentTime< double >(
            basic_astrodynamics::utc_scale,
            basic_astrodynamics::tdb_scale,
            input_output::sum_lmk::getSumImageObservationTime( image ),
            Eigen::Vector3d::Zero( ) );
}

inline std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > pointingSettingsForTest(
        const LoadedPixelTestData& loaded )
{
    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > settings;
    for( const auto& entry : loaded.supplementaryData_ )
    {
        for( const auto& instrument : entry->getInstrumentSupplementaryData( ) )
        {
            const auto camera = std::dynamic_pointer_cast< data::CameraInstrumentSupplementaryData >( instrument );
            if( camera != nullptr )
            {
                settings.push_back( estimatable_parameters::cameraPointingCorrection( entry->getBodyName( ), camera->getCameraId( ) ) );
            }
        }
    }
    return settings;
}

}  // namespace unit_tests
}  // namespace tudat

#endif  // TUDAT_SUM_LMK_TEST_HELPERS_H
