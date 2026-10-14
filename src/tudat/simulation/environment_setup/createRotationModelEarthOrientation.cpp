/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/createRotationModel.h"

#include <iostream>

#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"
#include "tudat/astro/earth_orientation/shortPeriodEarthOrientationCorrectionCalculator.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/io/readEarthOrientationData.h"

namespace tudat
{

namespace simulation_setup
{

// Construct the selected Earth rotation directly from neutral EOP observations and its model policies.
std::shared_ptr< ephemerides::RotationalEphemeris > createGcrsToItrsRotationModel(
        const std::shared_ptr< GcrsToItrsRotationModelSettings >& gcrsToItrsRotationSettings )
{
    // Use supplied observations, or read the requested file without building legacy UTC maps.
    std::shared_ptr< const input_output::EarthOrientationData > data = gcrsToItrsRotationSettings->getEopData( );

    if( !data )
    {
        // Preserve the file-format names accepted by the existing rotation settings.
        const std::string formatName = gcrsToItrsRotationSettings->getEopFileFormat( );

        if( formatName != "C04" && formatName != "finals2000A" )
        {
            throw std::runtime_error( "Unsupported Earth orientation format: " + formatName );
        }

        const input_output::EarthOrientationFormat format =
                formatName == "C04" ? input_output::EarthOrientationFormat::c04 : input_output::EarthOrientationFormat::finals2000a;

        data = std::make_shared< input_output::EarthOrientationData >(
                input_output::readEarthOrientationData( gcrsToItrsRotationSettings->getEopFile( ), format ) );

        // Retain the legacy file-reader warning for unsupported measured-pole conventions.
        if( gcrsToItrsRotationSettings->getNutationTheory( ) != basic_astrodynamics::iau_2000_a &&
            gcrsToItrsRotationSettings->getNutationTheory( ) != basic_astrodynamics::iau_2006 )
        {
            std::cerr << "Warning, only IAU2000 nutation theory format currently supported by reader." << std::endl;
        }
    }

    // One evaluator supplies the selected daily observations to clocks and orientation models.
    const auto evaluator = std::make_shared< earth_orientation::EarthOrientationDataEvaluator >(
            *data, gcrsToItrsRotationSettings->getEopEvaluationSettings( ) );

    // Selected sub-daily corrections accompany the daily observations.
    const auto polarSettings = gcrsToItrsRotationSettings->getPolarMotionCorrectionSettings( );

    const auto ut1Settings = gcrsToItrsRotationSettings->getUt1CorrectionSettings( );

    if( !polarSettings || !ut1Settings )
    {
        throw std::runtime_error( "Missing short-period EOP correction settings" );
    }

    // Store the selected polar-motion and UT1 correction calculators.
    earth_orientation::HighFrequencyEopCorrectionCalculators corrections;

    if( gcrsToItrsRotationSettings->getOceanTideEopModel( ) == earth_orientation::OceanTideEopModel::desai_sibois_2016 )
    {
        // Custom files have no reliable component identity. Reject an ambiguous
        // combination instead of dropping custom terms or counting ocean terms twice.
        GcrsToItrsRotationModelSettings conventional;
        const auto matches = []( const std::shared_ptr< EopCorrectionSettings >& actual,
                                 const std::shared_ptr< EopCorrectionSettings >& expected ) {
            return actual->conversionFactor_ == expected->conversionFactor_ && actual->amplitudesFiles_ == expected->amplitudesFiles_ &&
                    actual->argumentMultipliersFile_ == expected->argumentMultipliersFile_;
        };
        if( !matches( polarSettings, conventional.getPolarMotionCorrectionSettings( ) ) ||
            !matches( ut1Settings, conventional.getUt1CorrectionSettings( ) ) )
            throw std::runtime_error(
                    "Desai--Sibois replaces the conventional ocean-tide EOP terms; "
                    "it cannot be combined automatically with custom short-period correction files or units" );
        corrections = earth_orientation::createHighFrequencyEopCorrectionCalculators(
                earth_orientation::OceanTideEopModel::desai_sibois_2016, polarSettings->minimumAmplitude_, ut1Settings->minimumAmplitude_ );
    }
    else
    {
        // Preserve the existing conventional/custom path and summation order.
        corrections.polarMotion = std::make_shared< earth_orientation::ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > >(
                gcrsToItrsRotationSettings->getPolarMotionCorrectionSettings( )->conversionFactor_,
                gcrsToItrsRotationSettings->getPolarMotionCorrectionSettings( )->minimumAmplitude_,
                gcrsToItrsRotationSettings->getPolarMotionCorrectionSettings( )->amplitudesFiles_,
                gcrsToItrsRotationSettings->getPolarMotionCorrectionSettings( )->argumentMultipliersFile_,
                std::bind( &sofa_interface::calculateApproximateDelaunayFundamentalArgumentsWithGmst, std::placeholders::_1 ),
                gcrsToItrsRotationSettings->getShortTermInterpolatorSettings( ) );
        corrections.ut1 = std::make_shared< earth_orientation::ShortPeriodEarthOrientationCorrectionCalculator< double > >(
                gcrsToItrsRotationSettings->getUt1CorrectionSettings( )->conversionFactor_,
                gcrsToItrsRotationSettings->getUt1CorrectionSettings( )->minimumAmplitude_,
                gcrsToItrsRotationSettings->getUt1CorrectionSettings( )->amplitudesFiles_,
                gcrsToItrsRotationSettings->getUt1CorrectionSettings( )->argumentMultipliersFile_,
                std::bind( &sofa_interface::calculateApproximateDelaunayFundamentalArgumentsWithGmst, std::placeholders::_1 ),
                gcrsToItrsRotationSettings->getShortTermInterpolatorSettings( ) );
    }

    // Create IAU 2006 precession/nutation calculator
    std::shared_ptr< earth_orientation::PrecessionNutationCalculator > precessionNutationCalculator =
            std::make_shared< earth_orientation::PrecessionNutationCalculator >( gcrsToItrsRotationSettings->getNutationTheory( ),
                                                                                 nullptr,
                                                                                 gcrsToItrsRotationSettings->getCioInterpolatorSettings( ),
                                                                                 evaluator );

    // Create default time scale converter
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > terrestrialTimeScaleConverter =
            std::make_shared< earth_orientation::TerrestrialTimeScaleConverter >(
                    nullptr, corrections.ut1, gcrsToItrsRotationSettings->getTdbToTtInterpolatorSettings( ), evaluator );

    const auto polarMotionCalculator = std::make_shared< earth_orientation::PolarMotionCalculator >(
            nullptr, corrections.polarMotion, evaluator, [ terrestrialTimeScaleConverter ]( const Time& tai ) {
                return terrestrialTimeScaleConverter->getUt1FromTai( tai );
            } );

    // Create rotation model
    std::shared_ptr< earth_orientation::EarthOrientationAnglesCalculator > earthOrientationCalculator =
            std::make_shared< earth_orientation::EarthOrientationAnglesCalculator >(
                    polarMotionCalculator, precessionNutationCalculator, terrestrialTimeScaleConverter );
    return std::make_shared< ephemerides::GcrsToItrsRotationModel >(
            earthOrientationCalculator, gcrsToItrsRotationSettings->getInputTimeScale( ), gcrsToItrsRotationSettings->getOriginalFrame( ) );
}

}  // namespace simulation_setup

}  // namespace tudat
