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
    std::shared_ptr< earth_orientation::EOPReader > eopReader =
            std::make_shared< earth_orientation::EOPReader >( gcrsToItrsRotationSettings->getEopFile( ),
                                                              gcrsToItrsRotationSettings->getEopFileFormat( ),
                                                              gcrsToItrsRotationSettings->getNutationTheory( ) );

    // Load polar motion corrections
    std::shared_ptr< interpolators::LinearInterpolator< double, Eigen::Vector2d > > cipInItrsInterpolator =
            std::make_shared< interpolators::LinearInterpolator< double, Eigen::Vector2d > >(
                    eopReader->getCipInItrsMapInSecondsSinceJ2000( ),
                    interpolators::huntingAlgorithm,
                    interpolators::use_default_value,
                    std::make_pair( Eigen::Vector2d::Zero( ), Eigen::Vector2d::Zero( ) ) );

    // Load nutation corrections
    std::shared_ptr< interpolators::LinearInterpolator< double, Eigen::Vector2d > > cipInGcrsCorrectionInterpolator =
            std::make_shared< interpolators::LinearInterpolator< double, Eigen::Vector2d > >(
                    eopReader->getCipInGcrsCorrectionMapInSecondsSinceJ2000( ),
                    interpolators::huntingAlgorithm,
                    interpolators::use_default_value,
                    std::make_pair( Eigen::Vector2d::Zero( ), Eigen::Vector2d::Zero( ) ) );

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

    // Create full polar motion calculator
    std::shared_ptr< earth_orientation::PolarMotionCalculator > polarMotionCalculator =
            std::make_shared< earth_orientation::PolarMotionCalculator >( cipInItrsInterpolator, shortPeriodPolarMotionCalculator );

    // Create IAU 2006 precession/nutation calculator
    std::shared_ptr< earth_orientation::PrecessionNutationCalculator > precessionNutationCalculator =
            std::make_shared< earth_orientation::PrecessionNutationCalculator >(
                    gcrsToItrsRotationSettings->getNutationTheory( ),
                    cipInGcrsCorrectionInterpolator,
                    gcrsToItrsRotationSettings->getCioInterpolatorSettings( ) );

    // Create UT1 correction (sub-diural frequencies) object
    std::shared_ptr< earth_orientation::ShortPeriodEarthOrientationCorrectionCalculator< double > > ut1CorrectionSettings =
            std::make_shared< earth_orientation::ShortPeriodEarthOrientationCorrectionCalculator< double > >(
                    gcrsToItrsRotationSettings->getUt1CorrectionSettings( )->conversionFactor_,
                    gcrsToItrsRotationSettings->getUt1CorrectionSettings( )->minimumAmplitude_,
                    gcrsToItrsRotationSettings->getUt1CorrectionSettings( )->amplitudesFiles_,
                    gcrsToItrsRotationSettings->getUt1CorrectionSettings( )->argumentMultipliersFile_,
                    std::bind( &sofa_interface::calculateApproximateDelaunayFundamentalArgumentsWithGmst, std::placeholders::_1 ),
                    gcrsToItrsRotationSettings->getShortTermInterpolatorSettings( ) );

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > dailyUtcUt1CorrectionInterpolator =
            std::make_shared< interpolators::JumpDataLinearInterpolator< double, double > >(
                    eopReader->getUt1MinusUtcMapInSecondsSinceJ2000( ),
                    0.5,
                    1.0,
                    interpolators::huntingAlgorithm,
                    interpolators::use_default_value,
                    0.0 );

    // Create default time scale converter
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > terrestrialTimeScaleConverter =
            std::make_shared< earth_orientation::TerrestrialTimeScaleConverter >(
                    dailyUtcUt1CorrectionInterpolator,
                    ut1CorrectionSettings,
                    gcrsToItrsRotationSettings->getTdbToTtInterpolatorSettings( ) );

    // Create rotation model
    std::shared_ptr< earth_orientation::EarthOrientationAnglesCalculator > earthOrientationCalculator =
            std::make_shared< earth_orientation::EarthOrientationAnglesCalculator >(
                    polarMotionCalculator, precessionNutationCalculator, terrestrialTimeScaleConverter );
    return std::make_shared< ephemerides::GcrsToItrsRotationModel >(
            earthOrientationCalculator, gcrsToItrsRotationSettings->getInputTimeScale( ), gcrsToItrsRotationSettings->getOriginalFrame( ) );
}

}  // namespace simulation_setup

}  // namespace tudat
