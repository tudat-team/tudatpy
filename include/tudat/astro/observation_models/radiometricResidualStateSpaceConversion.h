/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RADIOMETRIC_RESIDUAL_STATE_SPACE_CONVERSION_H
#define TUDAT_RADIOMETRIC_RESIDUAL_STATE_SPACE_CONVERSION_H

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/observationAncillarySettings.h"
#include "tudat/astro/observation_models/observationFrequencies.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace observation_models
{

inline bool isFinitePositiveConversionValue( const double value )
{
    return std::isfinite( value ) && value > 0.0;
}

inline bool isValidFrequencyBandDouble( const double frequencyBandDouble )
{
    if( !std::isfinite( frequencyBandDouble ) )
    {
        return false;
    }
    const int bandIndex = static_cast< int >( std::lround( frequencyBandDouble ) );
    return std::fabs( frequencyBandDouble - static_cast< double >( bandIndex ) ) < 1.0e-12 && bandIndex >= static_cast< int >( s_band ) &&
            bandIndex <= static_cast< int >( ku_band );
}

inline std::string getUnsupportedResidualStateSpaceConversionMessage( const ObservableType observableType )
{
    return "Error when computing residual state-space conversion factor: observable type " + getObservableName( observableType ) +
            " is not supported. Supported types are dsn_n_way_range, dsn_n_way_averaged_doppler, "
            "doppler_measured_frequency, and one_way_doppler_measured_frequency.";
}

inline double computeDopplerResidualStateSpaceConversionFactorFromReferenceCarrierFrequency(
        const double downlinkReferenceCarrierFrequency,
        const unsigned int numberOfWayFactors )
{
    if( numberOfWayFactors != 1 && numberOfWayFactors != 2 )
    {
        throw std::runtime_error(
                "Error when computing Doppler residual state-space conversion factor: the round-trip factor must be 1 "
                "(one-way) or 2 (two-/three-way)." );
    }
    if( !isFinitePositiveConversionValue( downlinkReferenceCarrierFrequency ) )
    {
        throw std::runtime_error(
                "Error when computing Doppler residual state-space conversion factor: the downlink reference carrier "
                "frequency f_C must be finite and positive." );
    }
    return physical_constants::SPEED_OF_LIGHT / ( static_cast< double >( numberOfWayFactors ) * downlinkReferenceCarrierFrequency );
}

inline double reconstructDsnDownlinkReferenceCarrierFrequency( ObservationAncillarySimulationSettings* ancillarySettings )
{
    if( ancillarySettings == nullptr )
    {
        throw std::runtime_error(
                "Error when reconstructing DSN downlink reference carrier frequency: no ancillary settings were provided." );
    }

    const double referenceFrequency = ancillarySettings->getAncillaryDoubleData( doppler_reference_frequency, false );
    const double receptionBandDouble = ancillarySettings->getAncillaryDoubleData( reception_reference_frequency_band, false );
    const std::vector< double > frequencyBandsDouble = ancillarySettings->getAncillaryDoubleVectorData( frequency_bands, false );

    if( !isFinitePositiveConversionValue( referenceFrequency ) )
    {
        throw std::runtime_error(
                "Error when reconstructing DSN downlink reference carrier frequency: doppler_reference_frequency is "
                "missing or invalid." );
    }
    if( !isValidFrequencyBandDouble( receptionBandDouble ) )
    {
        throw std::runtime_error(
                "Error when reconstructing DSN downlink reference carrier frequency: reception_reference_frequency_band "
                "is missing or invalid." );
    }
    if( frequencyBandsDouble.empty( ) || !isValidFrequencyBandDouble( frequencyBandsDouble.back( ) ) )
    {
        throw std::runtime_error(
                "Error when reconstructing DSN downlink reference carrier frequency: frequency_bands is missing the "
                "downlink band." );
    }

    const FrequencyBands receptionReferenceBand = convertDoubleToFrequencyBand( receptionBandDouble );
    const FrequencyBands downlinkBand = convertDoubleToFrequencyBand( frequencyBandsDouble.back( ) );
    const double turnaroundRatio = getDsnDefaultTurnaroundRatios( receptionReferenceBand, downlinkBand );
    const double downlinkReferenceCarrierFrequency = turnaroundRatio * referenceFrequency;
    if( !isFinitePositiveConversionValue( downlinkReferenceCarrierFrequency ) )
    {
        throw std::runtime_error(
                "Error when reconstructing DSN downlink reference carrier frequency: M_2R * f_ref is not a finite "
                "positive frequency." );
    }
    return downlinkReferenceCarrierFrequency;
}

inline double computeDopplerResidualStateSpaceConversionFactor( const ObservableType observableType,
                                                                ObservationAncillarySimulationSettings* ancillarySettings )
{
    if( ancillarySettings != nullptr )
    {
        const double storedFactor = ancillarySettings->getAncillaryDoubleData( doppler_conversion_factor, false );
        if( isFinitePositiveConversionValue( storedFactor ) )
        {
            return storedFactor;
        }
    }

    switch( observableType )
    {
        case dsn_n_way_averaged_doppler:
        case doppler_measured_frequency: {
            // Set-level conversion only: f_C ≈ M_2R * f_ref, then c / (2 f_C). This is a reference conversion for
            // grouped ODF/DSN sets, not a per-observation ramp-aware conversion.
            const double downlinkReferenceCarrierFrequency = reconstructDsnDownlinkReferenceCarrierFrequency( ancillarySettings );
            return computeDopplerResidualStateSpaceConversionFactorFromReferenceCarrierFrequency( downlinkReferenceCarrierFrequency, 2 );
        }
        case one_way_doppler_measured_frequency: {
            if( ancillarySettings == nullptr )
            {
                throw std::runtime_error(
                        "Error when computing residual state-space conversion factor for "
                        "one_way_doppler_measured_frequency: no ancillary settings were provided. Conversion requires a "
                        "stored doppler_conversion_factor or a set-level doppler_reference_frequency used as f_C. This "
                        "is a set-level reference conversion, not a per-observation ramp-aware conversion." );
            }
            const double referenceCarrierFrequency = ancillarySettings->getAncillaryDoubleData( doppler_reference_frequency, false );
            if( !isFinitePositiveConversionValue( referenceCarrierFrequency ) )
            {
                throw std::runtime_error(
                        "Error when computing residual state-space conversion factor for "
                        "one_way_doppler_measured_frequency: neither doppler_conversion_factor nor a set-level "
                        "doppler_reference_frequency (used as f_C) is available. This is a set-level reference "
                        "conversion, not a per-observation ramp-aware conversion." );
            }
            return computeDopplerResidualStateSpaceConversionFactorFromReferenceCarrierFrequency( referenceCarrierFrequency, 1 );
        }
        default:
            throw std::runtime_error( getUnsupportedResidualStateSpaceConversionMessage( observableType ) );
    }
}

inline double computeDopplerResidualStateSpaceConversionFactor(
        const ObservableType observableType,
        const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings )
{
    return computeDopplerResidualStateSpaceConversionFactor( observableType, ancillarySettings.get( ) );
}

inline double computeResidualStateSpaceConversionFactor( const ObservableType observableType,
                                                         ObservationAncillarySimulationSettings* ancillarySettings )
{
    switch( observableType )
    {
        case dsn_n_way_range: {
            if( ancillarySettings == nullptr )
            {
                throw std::runtime_error(
                        "Error when computing residual state-space conversion factor for dsn_n_way_range: no ancillary "
                        "settings were provided. range_conversion_factor depends on the uplink frequency and is typically "
                        "available only after the observation has been computed." );
            }
            const double rangeConversionFactor = ancillarySettings->getAncillaryDoubleData( range_conversion_factor, false );
            if( !isFinitePositiveConversionValue( rangeConversionFactor ) )
            {
                throw std::runtime_error(
                        "Error when computing residual state-space conversion factor for dsn_n_way_range: "
                        "range_conversion_factor is missing or invalid. This factor depends on the uplink frequency and "
                        "is typically available only after the observation has been computed." );
            }
            return rangeConversionFactor;
        }
        case dsn_n_way_averaged_doppler:
        case doppler_measured_frequency:
        case one_way_doppler_measured_frequency:
            return computeDopplerResidualStateSpaceConversionFactor( observableType, ancillarySettings );
        default:
            throw std::runtime_error( getUnsupportedResidualStateSpaceConversionMessage( observableType ) );
    }
}

inline double computeResidualStateSpaceConversionFactor(
        const ObservableType observableType,
        const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings )
{
    return computeResidualStateSpaceConversionFactor( observableType, ancillarySettings.get( ) );
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_RADIOMETRIC_RESIDUAL_STATE_SPACE_CONVERSION_H
