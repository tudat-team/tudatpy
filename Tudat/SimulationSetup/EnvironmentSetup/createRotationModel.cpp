/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>
#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"

#if USE_CSPICE
#include "Tudat/External/SpiceInterface/spiceRotationalEphemeris.h"
#endif

#include "Tudat/SimulationSetup/EnvironmentSetup/createRotationModel.h"

#if USE_SOFA
#include "Tudat/Astrodynamics/Ephemerides/itrsToGcrsRotationModel.h"
#include "Tudat/Astrodynamics/EarthOrientation/earthOrientationCalculator.h"
#include "Tudat/Astrodynamics/EarthOrientation/shortPeriodEarthOrientationCorrectionCalculator.h"
#include "Tudat/Mathematics/Interpolators/jumpDataLinearInterpolator.h"
#endif

namespace tudat
{

namespace simulation_setup
{

//! Function to create a rotation model.
std::shared_ptr< ephemerides::RotationalEphemeris > createRotationModel(
        const std::shared_ptr< RotationModelSettings > rotationModelSettings,
        const std::string& body )
{
    using namespace tudat::ephemerides;

    // Declare return object.
    std::shared_ptr< RotationalEphemeris > rotationalEphemeris;

    // Check which type of rotation model is to be created.
    switch( rotationModelSettings->getRotationType( ) )
    {
    case simple_rotation_model:
    {
        // Check whether settings for simple rotation model are consistent with its type.
        std::shared_ptr< SimpleRotationModelSettings > simpleRotationSettings =
                std::dynamic_pointer_cast< SimpleRotationModelSettings >( rotationModelSettings );
        if( simpleRotationSettings == nullptr )
        {
            throw std::runtime_error(
                        "Error, expected simple rotation model settings for " + body );
        }
        else
        {
            // Create and initialize simple rotation model.
            rotationalEphemeris = std::make_shared< SimpleRotationalEphemeris >(
                        simpleRotationSettings->getInitialOrientation( ),
                        simpleRotationSettings->getRotationRate( ),
                        simpleRotationSettings->getInitialTime( ),
                        simpleRotationSettings->getOriginalFrame( ),
                        simpleRotationSettings->getTargetFrame( ) );
        }
        break;
    }
#if USE_SOFA
    case gcrs_to_itrs_rotation_model:
    {

        // Check whether settings for simple rotation model are consistent with its type.
        std::shared_ptr< GcrsToItrsRotationModelSettings > gcrsToItrsRotationSettings =
                std::dynamic_pointer_cast< GcrsToItrsRotationModelSettings >( rotationModelSettings );
        if( gcrsToItrsRotationSettings == nullptr )
        {
            throw std::runtime_error(
                        "Error, expected GCRS to ITRS rotation model settings for " + body );
        }
        else
        {
            std::shared_ptr< earth_orientation::EOPReader > eopReader = std::make_shared< earth_orientation::EOPReader >(
                        gcrsToItrsRotationSettings->getEopFile( ),
                        gcrsToItrsRotationSettings->getEopFileFormat( ),
                        gcrsToItrsRotationSettings->getNutationTheory( ) );

            // Load polar motion corrections
            std::shared_ptr< interpolators::LinearInterpolator< double, Eigen::Vector2d > > cipInItrsInterpolator =
                    std::make_shared< interpolators::LinearInterpolator< double, Eigen::Vector2d > >(
                        eopReader->getCipInItrsMapInSecondsSinceJ2000( ) );

            // Load nutation corrections
            std::shared_ptr< interpolators::LinearInterpolator< double, Eigen::Vector2d > > cipInGcrsCorrectionInterpolator =
                    std::make_shared< interpolators::LinearInterpolator< double, Eigen::Vector2d > >(
                        eopReader->getCipInGcrsCorrectionMapInSecondsSinceJ2000( ) );

            // Create polar motion correction (sub-diural frequencies) object
            std::shared_ptr< earth_orientation::ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > >
                    shortPeriodPolarMotionCalculator =
                    std::make_shared< earth_orientation::ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > >(
                        gcrsToItrsRotationSettings->getPolarMotionCorrectionSettings( )->conversionFactor_,
                        gcrsToItrsRotationSettings->getPolarMotionCorrectionSettings( )->minimumAmplitude_,
                        gcrsToItrsRotationSettings->getPolarMotionCorrectionSettings( )->amplitudesFiles_,
                        gcrsToItrsRotationSettings->getPolarMotionCorrectionSettings( )->argumentMultipliersFile_ );

            // Create full polar motion calculator
            std::shared_ptr< earth_orientation::PolarMotionCalculator > polarMotionCalculator =
                    std::make_shared< earth_orientation::PolarMotionCalculator >
                    ( cipInItrsInterpolator, shortPeriodPolarMotionCalculator );

            // Create IAU 2006 precession/nutation calculator
            std::shared_ptr< earth_orientation::PrecessionNutationCalculator > precessionNutationCalculator =
                    std::make_shared< earth_orientation::PrecessionNutationCalculator >(
                        gcrsToItrsRotationSettings->getNutationTheory( ), cipInGcrsCorrectionInterpolator );

            // Create UT1 correction (sub-diural frequencies) object
            std::shared_ptr< earth_orientation::ShortPeriodEarthOrientationCorrectionCalculator< double > >
                    ut1CorrectionSettings =
                    std::make_shared< earth_orientation::ShortPeriodEarthOrientationCorrectionCalculator< double > >(
                        gcrsToItrsRotationSettings->getUt1CorrectionSettings( )->conversionFactor_,
                        gcrsToItrsRotationSettings->getUt1CorrectionSettings( )->minimumAmplitude_,
                        gcrsToItrsRotationSettings->getUt1CorrectionSettings( )->amplitudesFiles_,
                        gcrsToItrsRotationSettings->getUt1CorrectionSettings( )->argumentMultipliersFile_ );

            std::shared_ptr< interpolators::OneDimensionalInterpolator < double, double > > dailyUtcUt1CorrectionInterpolator =
                    std::make_shared< interpolators::JumpDataLinearInterpolator< double, double > >(
                        eopReader->getUt1MinusUtcMapInSecondsSinceJ2000( ), 0.5, 1.0 );

            // Create default time scale converter
            std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > terrestrialTimeScaleConverter =
                    std::make_shared< earth_orientation::TerrestrialTimeScaleConverter >
                    (  dailyUtcUt1CorrectionInterpolator, ut1CorrectionSettings );

            // Create rotation model
            std::shared_ptr< earth_orientation::EarthOrientationAnglesCalculator > earthOrientationCalculator =
                    std::make_shared< earth_orientation::EarthOrientationAnglesCalculator >(
                        polarMotionCalculator, precessionNutationCalculator, terrestrialTimeScaleConverter );
            rotationalEphemeris = std::make_shared< ephemerides::GcrsToItrsRotationModel >(
                        earthOrientationCalculator, gcrsToItrsRotationSettings->getInputTimeScale( ),
                        gcrsToItrsRotationSettings->getOriginalFrame( ) );

            break;
        }

    }
#endif

#if USE_CSPICE
    case spice_rotation_model:
    {
        // Create rotational ephemeris directly from Spice.
        rotationalEphemeris = std::make_shared< SpiceRotationalEphemeris >(
                    rotationModelSettings->getOriginalFrame( ),
                    rotationModelSettings->getTargetFrame( ) );
        break;
    }
#endif
    default:
        throw std::runtime_error(
                    "Error, did not recognize rotation model settings type " +
                    std::to_string( rotationModelSettings->getRotationType( ) ) );
    }

    return rotationalEphemeris;
}

} // namespace simulation_setup

} // namespace tudat
