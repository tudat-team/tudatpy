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

#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"
#include "tudat/astro/earth_orientation/shortPeriodEarthOrientationCorrectionCalculator.h"
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/math/interpolators/jumpDataLinearInterpolator.h"

namespace tudat
{

namespace simulation_setup
{

std::shared_ptr< ephemerides::RotationalEphemeris > createGcrsToItrsRotationModel(
        const std::shared_ptr< GcrsToItrsRotationModelSettings >& gcrsToItrsRotationSettings )
{
    if( gcrsToItrsRotationSettings->getShortTermInterpolatorSettings( ) )
        throw std::runtime_error(
                "Strict GCRS/ITRS models do not support short-term EOP interpolation settings: "
                "short-period corrections must use continuous UT1/TT arguments" );
    std::shared_ptr< const input_output::EarthOrientationData > data = gcrsToItrsRotationSettings->getEopData( );
    if( !data )
    {
        const earth_orientation::EOPReader reader( gcrsToItrsRotationSettings->getEopFile( ),
                                                   gcrsToItrsRotationSettings->getEopFileFormat( ),
                                                   gcrsToItrsRotationSettings->getNutationTheory( ) );
        data = std::make_shared< input_output::EarthOrientationData >( reader.getRawEarthOrientationData( ) );
    }
    const auto evaluator = std::make_shared< earth_orientation::EarthOrientationDataEvaluator >(
            *data, gcrsToItrsRotationSettings->getEopEvaluationSettings( ) );
    const auto cipInItrsInterpolator = std::make_shared< earth_orientation::EarthOrientationDataInterpolator< Eigen::Vector2d > >(
            evaluator, earth_orientation::EarthOrientationDataField::polar_motion );
    const auto cipInGcrsCorrectionInterpolator = std::make_shared< earth_orientation::EarthOrientationDataInterpolator< Eigen::Vector2d > >(
            evaluator, earth_orientation::EarthOrientationDataField::celestial_pole );
    // Create polar motion correction (sub-diural frequencies) object
    std::shared_ptr< earth_orientation::ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > >
            shortPeriodPolarMotionCalculator =
                    std::make_shared< earth_orientation::ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > >(
                            gcrsToItrsRotationSettings->getPolarMotionCorrectionSettings( )->conversionFactor_,
                            gcrsToItrsRotationSettings->getPolarMotionCorrectionSettings( )->minimumAmplitude_,
                            gcrsToItrsRotationSettings->getPolarMotionCorrectionSettings( )->amplitudesFiles_,
                            gcrsToItrsRotationSettings->getPolarMotionCorrectionSettings( )->argumentMultipliersFile_,
                            std::bind( &sofa_interface::calculateApproximateDelaunayFundamentalArgumentsWithGmst, std::placeholders::_1 ),
                            gcrsToItrsRotationSettings->getShortTermInterpolatorSettings( ) );

    // Create IAU 2006 precession/nutation calculator
    std::shared_ptr< earth_orientation::PrecessionNutationCalculator > precessionNutationCalculator =
            std::make_shared< earth_orientation::PrecessionNutationCalculator >( gcrsToItrsRotationSettings->getNutationTheory( ),
                                                                                 cipInGcrsCorrectionInterpolator,
                                                                                 gcrsToItrsRotationSettings->getCioInterpolatorSettings( ),
                                                                                 evaluator );

    // Create UT1 correction (sub-diural frequencies) object
    std::shared_ptr< earth_orientation::ShortPeriodEarthOrientationCorrectionCalculator< double > > ut1CorrectionSettings =
            std::make_shared< earth_orientation::ShortPeriodEarthOrientationCorrectionCalculator< double > >(
                    gcrsToItrsRotationSettings->getUt1CorrectionSettings( )->conversionFactor_,
                    gcrsToItrsRotationSettings->getUt1CorrectionSettings( )->minimumAmplitude_,
                    gcrsToItrsRotationSettings->getUt1CorrectionSettings( )->amplitudesFiles_,
                    gcrsToItrsRotationSettings->getUt1CorrectionSettings( )->argumentMultipliersFile_,
                    std::bind( &sofa_interface::calculateApproximateDelaunayFundamentalArgumentsWithGmst, std::placeholders::_1 ),
                    gcrsToItrsRotationSettings->getShortTermInterpolatorSettings( ) );

    const auto dailyUtcUt1CorrectionInterpolator = std::make_shared< earth_orientation::EarthOrientationDataInterpolator< double > >(
            evaluator, earth_orientation::EarthOrientationDataField::ut1_minus_utc );

    // Create default time scale converter
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > terrestrialTimeScaleConverter =
            std::make_shared< earth_orientation::TerrestrialTimeScaleConverter >(
                    dailyUtcUt1CorrectionInterpolator,
                    ut1CorrectionSettings,
                    gcrsToItrsRotationSettings->getTdbToTtInterpolatorSettings( ),
                    evaluator );

    const auto polarMotionCalculator = std::make_shared< earth_orientation::PolarMotionCalculator >(
            cipInItrsInterpolator, shortPeriodPolarMotionCalculator, evaluator, [ terrestrialTimeScaleConverter ]( const Time& tai ) {
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
