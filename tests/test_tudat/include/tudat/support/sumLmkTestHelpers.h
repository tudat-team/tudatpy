/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved. See http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SUM_LMK_TEST_HELPERS_H
#define TUDAT_SUM_LMK_TEST_HELPERS_H

#include "tudat/io/readSumLmkFiles.h"
#include "tudat/simulation/estimation_setup/createObservationCollection.h"
#include "tudat/simulation/estimation_setup/createInverseAprioriCovariance.h"
#include "tudat/simulation/estimation_setup/estimatableParameterSettings.h"
#include "tudat/simulation/estimation_setup/createObservationModelSettings.h"

namespace tudat
{
namespace unit_tests
{

//! Shared fixture setup; all conversion and environment changes use the public tracking-data interface.
struct LoadedPixelTestData
{
    std::shared_ptr< observation_models::ObservationCollection< double, double > > observationCollection_;
    std::vector< std::shared_ptr< data::TrackingSupplementaryData > > supplementaryData_;
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettings_;
    std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >
            inverseAprioriCovarianceDiagonalEntries_;
};

inline std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > pixelModelSettingsForTest(
        const std::shared_ptr< observation_models::ObservationCollection< double, double > >& collection )
{
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > settings;
    for( const auto& link : collection->getLinkDefinitionsForSingleObservable( observation_models::pixel_coordinates ) )
    {
        settings.push_back( observation_models::pixelCoordinatesSettings( link ) );
    }
    return settings;
}

inline LoadedPixelTestData loadPixelTestData(
        const std::vector< input_output::sum_lmk::SumImageData >& images,
        const std::map< std::string, input_output::sum_lmk::LmkLandmarkData >& landmarks,
        simulation_setup::SystemOfBodies& bodies,
        const std::string& targetBody = "Target",
        const std::string& receiverBody = "Spacecraft",
        const bool skipMissingLandmarks = false )
{
    const auto loaded = input_output::sum_lmk::convertSumLmkData< double, double >(
            images, landmarks, targetBody, receiverBody, skipMissingLandmarks );
    observation_models::setTrackingSupplementaryDataInBodies( bodies, loaded.second );
    LoadedPixelTestData result;
    result.supplementaryData_ = loaded.second;
    result.observationCollection_ = observation_models::createObservationCollection< double, double >( loaded.first, bodies );
    result.observationModelSettings_ = pixelModelSettingsForTest( result.observationCollection_ );
    result.inverseAprioriCovarianceDiagonalEntries_ =
            simulation_setup::getCameraPointingInverseAprioriCovarianceEntries( loaded.second );
    return result;
}

inline double sumImageTdbTimeForTest( const input_output::sum_lmk::SumImageData& image )
{
    return earth_orientation::createDefaultTimeConverter( )->getCurrentTime< double >(
            basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale,
            input_output::sum_lmk::getSumImageObservationTime( image ), Eigen::Vector3d::Zero( ) );
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
