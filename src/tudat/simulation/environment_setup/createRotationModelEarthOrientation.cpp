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
    std::shared_ptr< EopCorrectionSettings > polarSettings = gcrsToItrsRotationSettings->getPolarMotionCorrectionSettings( );

    std::shared_ptr< EopCorrectionSettings > ut1Settings = gcrsToItrsRotationSettings->getUt1CorrectionSettings( );

    if( !polarSettings || !ut1Settings )
    {
        throw std::runtime_error( "Missing short-period EOP correction settings" );
    }

    if( gcrsToItrsRotationSettings->getOceanTideEopModel( ) == earth_orientation::OceanTideEopModel::desai_sibois_2016 )
    {
        // The Desai and Sibois (2016) tables replace the conventional ocean-tide files; the libration files are retained.
        // Custom correction files cannot be identified as ocean or libration terms, so they are rejected.
        GcrsToItrsRotationModelSettings conventional;
        const auto isConventional = []( const std::shared_ptr< EopCorrectionSettings >& actual,
                                        const std::shared_ptr< EopCorrectionSettings >& expected ) {
            return actual->conversionFactor_ == expected->conversionFactor_ && actual->amplitudesFiles_ == expected->amplitudesFiles_ &&
                    actual->argumentMultipliersFile_ == expected->argumentMultipliersFile_;
        };
        if( !isConventional( polarSettings, conventional.getPolarMotionCorrectionSettings( ) ) ||
            !isConventional( ut1Settings, conventional.getUt1CorrectionSettings( ) ) )
        {
            throw std::runtime_error(
                    "Desai--Sibois replaces the conventional ocean-tide EOP terms; "
                    "it cannot be combined automatically with custom short-period correction files or units" );
        }

        // The second file of each conventional list contains the ocean-tide terms.
        const std::string path = paths::getEarthOrientationDataFilesPath( ) + "/";
        polarSettings = std::make_shared< EopCorrectionSettings >( *polarSettings );
        polarSettings->amplitudesFiles_.at( 1 ) = path + "desaiSibois2016PolarMotionOceanTidesAmplitudes.txt";
        polarSettings->argumentMultipliersFile_.at( 1 ) = path + "desaiSibois2016PolarMotionOceanTidesFundamentalArgumentMultipliers.txt";
        ut1Settings = std::make_shared< EopCorrectionSettings >( *ut1Settings );
        ut1Settings->amplitudesFiles_.at( 1 ) = path + "desaiSibois2016Ut1OceanTidesAmplitudes.txt";
        ut1Settings->argumentMultipliersFile_.at( 1 ) = path + "desaiSibois2016Ut1OceanTidesFundamentalArgumentMultipliers.txt";
    }

    // Short-period polar motion and UT1 corrections, summed in the order of the correction files.
    const auto polarMotionCorrections =
            std::make_shared< earth_orientation::ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > >(
                    polarSettings->conversionFactor_,
                    polarSettings->minimumAmplitude_,
                    polarSettings->amplitudesFiles_,
                    polarSettings->argumentMultipliersFile_,
                    std::bind( &sofa_interface::calculateApproximateDelaunayFundamentalArgumentsWithGmst, std::placeholders::_1 ),
                    gcrsToItrsRotationSettings->getShortTermInterpolatorSettings( ) );

    const auto ut1Corrections = std::make_shared< earth_orientation::ShortPeriodEarthOrientationCorrectionCalculator< double > >(
            ut1Settings->conversionFactor_,
            ut1Settings->minimumAmplitude_,
            ut1Settings->amplitudesFiles_,
            ut1Settings->argumentMultipliersFile_,
            std::bind( &sofa_interface::calculateApproximateDelaunayFundamentalArgumentsWithGmst, std::placeholders::_1 ),
            gcrsToItrsRotationSettings->getShortTermInterpolatorSettings( ) );

    // Create IAU 2006 precession/nutation calculator
    std::shared_ptr< earth_orientation::PrecessionNutationCalculator > precessionNutationCalculator =
            std::make_shared< earth_orientation::PrecessionNutationCalculator >( gcrsToItrsRotationSettings->getNutationTheory( ),
                                                                                 nullptr,
                                                                                 gcrsToItrsRotationSettings->getCioInterpolatorSettings( ),
                                                                                 evaluator );

    // Create default time scale converter
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > terrestrialTimeScaleConverter =
            std::make_shared< earth_orientation::TerrestrialTimeScaleConverter >(
                    nullptr, ut1Corrections, gcrsToItrsRotationSettings->getTdbToTtInterpolatorSettings( ), evaluator );

    const auto polarMotionCalculator = std::make_shared< earth_orientation::PolarMotionCalculator >(
            nullptr, polarMotionCorrections, evaluator, [ terrestrialTimeScaleConverter ]( const Time& tai ) {
                return terrestrialTimeScaleConverter->getUt1FromTai( tai );
            } );

    // Create rotation model
    std::shared_ptr< earth_orientation::EarthOrientationAnglesCalculator > earthOrientationCalculator =
            std::make_shared< earth_orientation::EarthOrientationAnglesCalculator >(
                    polarMotionCalculator, precessionNutationCalculator, terrestrialTimeScaleConverter );

    // Inside an environment, every rotation model (and every station and deformation model that uses it) is
    // evaluated at TDB epochs. A different input time scale is kept for compatibility, but is then misinterpreted.
    if( gcrsToItrsRotationSettings->getInputTimeScale( ) != basic_astrodynamics::tdb_scale )
    {
        std::cerr << "Warning when creating GCRS-ITRS rotation model: input time scale is not TDB. Environment models evaluate "
                     "body rotations at TDB epochs, so these epochs will be interpreted in the selected time scale instead."
                  << std::endl;
    }
    return std::make_shared< ephemerides::GcrsToItrsRotationModel >(
            earthOrientationCalculator, gcrsToItrsRotationSettings->getInputTimeScale( ), gcrsToItrsRotationSettings->getOriginalFrame( ) );
}

}  // namespace simulation_setup

}  // namespace tudat
