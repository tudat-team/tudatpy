/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROCESS_SUM_LMK_FILES_H
#define TUDAT_PROCESS_SUM_LMK_FILES_H

#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/io/readSumLmkFiles.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameterSet.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createObservationCollection.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/estimation_setup/createObservationModelSettings.h"
#include "tudat/simulation/estimation_setup/estimatableParameterSettings.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace tudat
{

namespace observation_models
{

class SumLmkObservationConversionSettings
{
public:
    SumLmkObservationConversionSettings( const std::string& targetBodyName, const std::string& receiverBodyName ):
        targetBodyName_( targetBodyName ), receiverBodyName_( receiverBodyName )
    {}

    std::string targetBodyName_;
    std::string receiverBodyName_;
    bool validateSpacecraftObjectGeometry_ = true;
    bool skipObservationsWithMissingLandmarks_ = false;
};

template< typename ObservationScalarType = double, typename TimeType = double >
struct SumLmkObservationConversionResult {
    std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > observationCollection_;
    std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >
            inverseAprioriCovarianceDiagonalEntries_;
    std::map< std::string, std::string > imageIdToCameraName_;
    // Body carrying the per-image cameras, i.e. the body the pointing-correction parameters belong to.
    std::string receiverBodyName_;
    // Observation model settings matching observationCollection_ (one pixel_coordinates setting per
    // (image, landmark) link end), ready for createObservationSimulators / residual computation.
    std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettings_;
};

//! Observation epoch of a SUM image, in seconds since J2000 in the TDB time scale.
template< typename TimeType = double >
TimeType getSumImageTdbObservationTime( const input_output::sum_lmk::SumImageData& image )
{
    return earth_orientation::createDefaultTimeConverter( )->getCurrentTime< TimeType >(
            basic_astrodynamics::utc_scale,
            basic_astrodynamics::tdb_scale,
            input_output::sum_lmk::getSumImageObservationTime( image ),
            Eigen::Vector3d::Zero( ) );
}

namespace detail
{

inline void validateReceiverAndTargetBodies( const simulation_setup::SystemOfBodies& bodies,
                                             const SumLmkObservationConversionSettings& conversionSettings )
{
    for( const std::string& bodyName : { conversionSettings.targetBodyName_, conversionSettings.receiverBodyName_ } )
    {
        if( bodies.count( bodyName ) == 0 )
        {
            throw std::runtime_error( "Error when converting SUM/LMK observations: body '" + bodyName + "' not found." );
        }
        if( bodies.at( bodyName )->getRotationalEphemeris( ) == nullptr )
        {
            throw std::runtime_error( "Error when converting SUM/LMK observations: body '" + bodyName +
                                      "' does not have a rotational ephemeris." );
        }
    }
}

//! Inverse a-priori variance of each camera pointing correction, from the SUM SIGMA_PTG values.
inline std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > > getSumLmkPointingAprioriEntries(
        const std::vector< std::shared_ptr< data::TrackingSupplementaryData > >& supplementaryData )
{
    std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > > entries;
    for( const std::shared_ptr< data::TrackingSupplementaryData >& currentData : supplementaryData )
    {
        for( const std::shared_ptr< data::InstrumentSupplementaryData >& instrument : currentData->getInstrumentSupplementaryData( ) )
        {
            auto camera = std::dynamic_pointer_cast< data::CameraInstrumentSupplementaryData >( instrument );
            if( camera == nullptr || !camera->getPointingUncertainty( ).array( ).isFinite( ).all( ) )
            {
                continue;
            }
            entries.push_back( std::make_pair( estimatable_parameters::EstimatebleParameterIdentifier(
                                                       estimatable_parameters::camera_pointing_correction,
                                                       std::make_pair( currentData->getBodyName( ), camera->getCameraId( ) ) ),
                                               Eigen::VectorXd( camera->getPointingUncertainty( ).array( ).square( ).inverse( ) ) ) );
        }
    }
    return entries;
}

}  // namespace detail

//! Build the pixel-coordinate observation model settings matching every (image, landmark) link end
//! in a SUM/LMK observation collection. Light-time is geometric and single-leg (an empty correction
//! list still triggers the model's light-time iteration); stellar aberration is off, per the SUM/LMK
//! conventions. No bodies are needed - the camera is resolved from VehicleSystems when simulators are
//! created. The resulting settings plug directly into createObservationSimulators and the generic
//! residual/estimation machinery, exactly like any other observable type.
template< typename ObservationScalarType = double, typename TimeType = double >
std::vector< std::shared_ptr< ObservationModelSettings > > createSumLmkObservationModelSettings(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >& observationCollection,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ) )
{
    if( observationCollection == nullptr )
    {
        throw std::runtime_error( "Error when creating SUM/LMK observation model settings: observation collection is null." );
    }
    std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettings;
    for( const LinkDefinition& linkDefinition : observationCollection->getLinkDefinitionsForSingleObservable( pixel_coordinates ) )
    {
        observationModelSettings.push_back( pixelCoordinatesSettings( linkDefinition, lightTimeCorrections ) );
    }
    return observationModelSettings;
}

//! Convert parsed SUM/LMK data and add the landmarks and per-image cameras to the bodies.
/*!
 *  This is a convenience wrapper around the general tracking-data path: convertSumLmkData, then
 *  setTrackingSupplementaryDataInBodies and createObservationCollection.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
SumLmkObservationConversionResult< ObservationScalarType, TimeType > createSumLmkObservationCollection(
        const std::vector< input_output::sum_lmk::SumImageData >& sumImages,
        const std::map< std::string, input_output::sum_lmk::LmkLandmarkData >& landmarks,
        simulation_setup::SystemOfBodies& bodies,
        const SumLmkObservationConversionSettings& conversionSettings )
{
    detail::validateReceiverAndTargetBodies( bodies, conversionSettings );

    // Conversion checks all input before the bodies are changed.
    const auto trackingData = input_output::sum_lmk::convertSumLmkData< ObservationScalarType, TimeType >(
            sumImages,
            landmarks,
            conversionSettings.targetBodyName_,
            conversionSettings.receiverBodyName_,
            conversionSettings.skipObservationsWithMissingLandmarks_,
            conversionSettings.validateSpacecraftObjectGeometry_ );
    setTrackingSupplementaryDataInBodies( bodies, trackingData.second );

    SumLmkObservationConversionResult< ObservationScalarType, TimeType > result;
    result.receiverBodyName_ = conversionSettings.receiverBodyName_;
    result.observationCollection_ = createObservationCollection< ObservationScalarType, TimeType >( trackingData.first, bodies );
    result.inverseAprioriCovarianceDiagonalEntries_ = detail::getSumLmkPointingAprioriEntries( trackingData.second );
    for( const input_output::sum_lmk::SumImageData& image : input_output::sum_lmk::getSumImagesToConvert(
                 sumImages, landmarks, conversionSettings.skipObservationsWithMissingLandmarks_, false ) )
    {
        result.imageIdToCameraName_[ image.imageId_ ] = input_output::sum_lmk::getSumImageCameraName( image.imageId_ );
    }
    result.observationModelSettings_ =
            createSumLmkObservationModelSettings< ObservationScalarType, TimeType >( result.observationCollection_ );
    return result;
}

template< typename ObservationScalarType = double, typename TimeType = double >
SumLmkObservationConversionResult< ObservationScalarType, TimeType > createSumLmkObservationCollection(
        const std::vector< std::string >& sumFiles,
        const std::vector< std::string >& lmkFiles,
        simulation_setup::SystemOfBodies& bodies,
        const SumLmkObservationConversionSettings& conversionSettings )
{
    return createSumLmkObservationCollection< ObservationScalarType, TimeType >(
            input_output::sum_lmk::readSumFiles( sumFiles ), input_output::sum_lmk::readLmkFiles( lmkFiles ), bodies, conversionSettings );
}

//! Build the camera_pointing_correction parameter settings for the per-image cameras registered by a
//! SUM/LMK conversion: one 3-vector pointing parameter per image, on the receiver body that carries the
//! cameras. The settings are returned ordered by camera name, so that the resulting parameter vector has a
//! reproducible layout. Pass a non-empty imageIds to restrict the settings to a subset of the images (for
//! instance to estimate pointing only for the images that have a SIGMA_PTG a-priori); an unknown image ID
//! is an error rather than being silently ignored.
template< typename ObservationScalarType = double, typename TimeType = double >
std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > createSumLmkPointingParameterSettings(
        const SumLmkObservationConversionResult< ObservationScalarType, TimeType >& conversionResult,
        const std::vector< std::string >& imageIds = std::vector< std::string >( ) )
{
    if( conversionResult.receiverBodyName_.empty( ) )
    {
        throw std::runtime_error( "Error when creating SUM/LMK pointing parameter settings: conversion result has no receiver body name." );
    }

    std::set< std::string > cameraNames;
    if( imageIds.empty( ) )
    {
        for( const auto& imageEntry : conversionResult.imageIdToCameraName_ )
        {
            cameraNames.insert( imageEntry.second );
        }
    }
    else
    {
        for( const std::string& imageId : imageIds )
        {
            if( conversionResult.imageIdToCameraName_.count( imageId ) == 0 )
            {
                throw std::runtime_error( "Error when creating SUM/LMK pointing parameter settings: image '" + imageId +
                                          "' is not part of the converted observations." );
            }
            cameraNames.insert( conversionResult.imageIdToCameraName_.at( imageId ) );
        }
    }

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameterSettings > > parameterSettings;
    for( const std::string& cameraName : cameraNames )
    {
        parameterSettings.push_back( estimatable_parameters::cameraPointingCorrection( conversionResult.receiverBodyName_, cameraName ) );
    }
    return parameterSettings;
}

//! Assemble an inverse a-priori covariance matrix for an EstimationInput from per-parameter inverse-variance
//! diagonal entries, such as the per-image SIGMA_PTG entries produced by a SUM/LMK conversion. Each entry
//! identifies a parameter by its (type, body, reference point) identifier; the entry's values are written onto
//! the diagonal of that parameter's block in the full parameter vector.
//!
//! Entries whose parameter is not part of parametersToEstimate are skipped: a conversion typically produces an
//! a-priori for every image, while only a subset of the images may have their pointing estimated. Pass
//! baseInverseAprioriCovariance to add these entries on top of an existing a-priori (for instance one already
//! constraining the initial state); it must then be square with the estimated parameter-set size.
template< typename InitialStateParameterType = double >
Eigen::MatrixXd createInverseAprioriCovarianceFromDiagonalEntries(
        const std::vector< std::pair< estimatable_parameters::EstimatebleParameterIdentifier, Eigen::VectorXd > >&
                inverseAprioriCovarianceDiagonalEntries,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parametersToEstimate,
        const Eigen::MatrixXd& baseInverseAprioriCovariance = Eigen::MatrixXd::Zero( 0, 0 ) )
{
    if( parametersToEstimate == nullptr )
    {
        throw std::runtime_error( "Error when creating inverse a-priori covariance: parameter set is null." );
    }

    const int numberOfParameters = parametersToEstimate->getEstimatedParameterSetSize( );
    Eigen::MatrixXd inverseAprioriCovariance;
    if( baseInverseAprioriCovariance.rows( ) == 0 && baseInverseAprioriCovariance.cols( ) == 0 )
    {
        inverseAprioriCovariance = Eigen::MatrixXd::Zero( numberOfParameters, numberOfParameters );
    }
    else if( baseInverseAprioriCovariance.rows( ) != numberOfParameters || baseInverseAprioriCovariance.cols( ) != numberOfParameters )
    {
        throw std::runtime_error( "Error when creating inverse a-priori covariance: base matrix is " +
                                  std::to_string( baseInverseAprioriCovariance.rows( ) ) + "x" +
                                  std::to_string( baseInverseAprioriCovariance.cols( ) ) + ", but the estimated parameter set has size " +
                                  std::to_string( numberOfParameters ) + "." );
    }
    else
    {
        inverseAprioriCovariance = baseInverseAprioriCovariance;
    }

    std::set< int > alreadySetStartIndices;
    for( const auto& entry : inverseAprioriCovarianceDiagonalEntries )
    {
        const std::vector< std::pair< int, int > > parameterIndices = parametersToEstimate->getIndicesForParameterType( entry.first );
        if( parameterIndices.empty( ) )
        {
            // Parameter is not estimated; its a-priori simply does not apply to this parameter set.
            continue;
        }
        if( parameterIndices.size( ) > 1 )
        {
            throw std::runtime_error( "Error when creating inverse a-priori covariance: parameter " +
                                      estimatable_parameters::getParameterTypeString( entry.first.first ) + "of (" +
                                      entry.first.second.first + ", " + entry.first.second.second +
                                      ") matches more than one parameter block." );
        }

        const int startIndex = parameterIndices.at( 0 ).first;
        const int parameterSize = parameterIndices.at( 0 ).second;
        if( entry.second.size( ) != parameterSize )
        {
            throw std::runtime_error( "Error when creating inverse a-priori covariance: a-priori for parameter " +
                                      estimatable_parameters::getParameterTypeString( entry.first.first ) + "of (" +
                                      entry.first.second.first + ", " + entry.first.second.second + ") has size " +
                                      std::to_string( entry.second.size( ) ) + ", but the parameter has size " +
                                      std::to_string( parameterSize ) + "." );
        }
        if( !alreadySetStartIndices.insert( startIndex ).second )
        {
            throw std::runtime_error( "Error when creating inverse a-priori covariance: parameter " +
                                      estimatable_parameters::getParameterTypeString( entry.first.first ) + "of (" +
                                      entry.first.second.first + ", " + entry.first.second.second +
                                      ") is given an a-priori more than once." );
        }

        for( int i = 0; i < parameterSize; ++i )
        {
            inverseAprioriCovariance( startIndex + i, startIndex + i ) += entry.second( i );
        }
    }

    return inverseAprioriCovariance;
}

//! Convenience overload assembling the inverse a-priori covariance directly from a SUM/LMK conversion result,
//! i.e. from the per-image SIGMA_PTG pointing a-priori it carries.
template< typename ObservationScalarType = double, typename TimeType = double, typename InitialStateParameterType = double >
Eigen::MatrixXd createSumLmkInverseAprioriCovariance(
        const SumLmkObservationConversionResult< ObservationScalarType, TimeType >& conversionResult,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< InitialStateParameterType > >& parametersToEstimate,
        const Eigen::MatrixXd& baseInverseAprioriCovariance = Eigen::MatrixXd::Zero( 0, 0 ) )
{
    return createInverseAprioriCovarianceFromDiagonalEntries< InitialStateParameterType >(
            conversionResult.inverseAprioriCovarianceDiagonalEntries_, parametersToEstimate, baseInverseAprioriCovariance );
}

//! Compute observed-minus-computed (O-C) pixel residuals for a SUM/LMK observation collection given a
//! fixed environment (e.g. a spacecraft SPICE trajectory). Builds the matching pixel-coordinate
//! observation simulators, evaluates the model at each observation's link end/time, and stores the
//! residuals into the collection (readable via getResiduals / getResidualStatistics). Returns the
//! concatenated residual vector for convenience. This is the pixel-landmark analogue of computing
//! Doppler residuals against a reference trajectory.
template< typename ObservationScalarType = double, typename TimeType = double >
Eigen::VectorXd computeSumLmkResiduals(
        const std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > >& observationCollection,
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections =
                std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ) )
{
    if( observationCollection == nullptr )
    {
        throw std::runtime_error( "Error when computing SUM/LMK residuals: observation collection is null." );
    }
    const std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettings =
            createSumLmkObservationModelSettings< ObservationScalarType, TimeType >( observationCollection, lightTimeCorrections );
    const std::vector< std::shared_ptr< ObservationSimulatorBase< ObservationScalarType, TimeType > > > observationSimulators =
            createObservationSimulators< ObservationScalarType, TimeType >( observationModelSettings, bodies );
    simulation_setup::computeResidualsAndDependentVariables< ObservationScalarType, TimeType >(
            observationCollection, observationSimulators, bodies );
    return observationCollection->getConcatenatedResiduals( ).template cast< double >( );
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_PROCESS_SUM_LMK_FILES_H
