/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/fullPlanetaryRotationModel.h"
#include "tudat/astro/ephemerides/tabulatedRotationalEphemeris.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"

#include "tudat/simulation/environment_setup/createRotationModel.h"

#if TUDAT_BUILD_WITH_SOFA_INTERFACE
#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"
#include "tudat/astro/earth_orientation/shortPeriodEarthOrientationCorrectionCalculator.h"
#include "tudat/math/interpolators/jumpDataLinearInterpolator.h"
#endif

#include "tudat/astro/ephemerides/synchronousRotationalEphemeris.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to retrieve a state from one of two functions
Eigen::Vector6d getStateFromSelectedStateFunction(
        const double currentTime,
        const bool useFirstFunction,
        const std::function< Eigen::Vector6d( const double ) > stateFunction1,
        const std::function< Eigen::Vector6d( const double ) > stateFunction2 )
{
    return ( useFirstFunction ) ? ( stateFunction1( currentTime ) ) : ( stateFunction2( currentTime ) );
}


//! Function to create a state function for a body, valid both during propagation, and outside propagation
std::function< Eigen::Vector6d( const double, bool ) > createRelativeStateFunction(
        const SystemOfBodies& bodies,
        const std::string orbitingBody,
        const std::string centralBody )
{
    // Retrieve state functions for relevant bodies (obtained from current state of body objects)
    std::function< Eigen::Vector6d( const double ) > bodyInertialStateFunction =
            std::bind( &Body::getState, bodies.at( orbitingBody ) );
    std::function< Eigen::Vector6d( const double ) > centralBodyInertialStateFunction =
            std::bind( &Body::getState, bodies.at(  centralBody ) );

    // Define relative state function from body object
    std::function< Eigen::Vector6d( const double ) > fromBodyStateFunction =
            std::bind(
                &ephemerides::getDifferenceBetweenStates, bodyInertialStateFunction,
                centralBodyInertialStateFunction, std::placeholders::_1 );

    // Define state function from ephemeris
    std::function< Eigen::Vector6d( const double ) > fromEphemerisStateFunction;

    if( bodies.at( orbitingBody )->getEphemeris( )->getReferenceFrameOrigin( ) == centralBody )
    {
        fromEphemerisStateFunction = std::bind( &ephemerides::Ephemeris::getCartesianState,
                                                bodies.at( orbitingBody )->getEphemeris( ), std::placeholders::_1 );

    }
    else
    {
        std::function< Eigen::Vector6d( const double ) > ephemerisInertialStateFunction =
                std::bind( &Body::getStateInBaseFrameFromEphemeris< double, double >, bodies.at( orbitingBody ),
                           std::placeholders::_1 );
        std::function< Eigen::Vector6d( const double ) > ephemerisCentralBodyInertialStateFunction =
                std::bind( &Body::getStateInBaseFrameFromEphemeris< double, double >, bodies.at( centralBody ),
                           std::placeholders::_1 );
        fromEphemerisStateFunction = std::bind(
                    &ephemerides::getDifferenceBetweenStates,
                    ephemerisInertialStateFunction,
                    ephemerisCentralBodyInertialStateFunction, std::placeholders::_1 );
    }

    return std::bind( &getStateFromSelectedStateFunction, std::placeholders::_1, std::placeholders::_2,
                      fromBodyStateFunction, fromEphemerisStateFunction );
}

//! Function to create a rotation model.
std::shared_ptr< ephemerides::RotationalEphemeris > createRotationModel(
        const std::shared_ptr< RotationModelSettings > rotationModelSettings,
        const std::string& body,
        const SystemOfBodies& bodies )
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
#if TUDAT_BUILD_WITH_SOFA_INTERFACE
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

    case spice_rotation_model:
    {
        // Create rotational ephemeris directly from Spice.
        rotationalEphemeris = std::make_shared< SpiceRotationalEphemeris >(
                    rotationModelSettings->getOriginalFrame( ),
                    rotationModelSettings->getTargetFrame( ) );
        break;
    }
    case planetary_rotation_model:
    {
        std::shared_ptr< PlanetaryRotationModelSettings > planetaryRotationModelSettings =
                std::dynamic_pointer_cast< PlanetaryRotationModelSettings >( rotationModelSettings );
        if( planetaryRotationModelSettings == nullptr )
        {
            std::cerr<<"Error, expected planetary rotation model settings for "<<body<<std::endl;
        }
        else
        {

            Eigen::VectorXd bodyKeplerElements = orbital_element_conversions::convertCartesianToKeplerianElements(
                        spice_interface::getBodyCartesianStateAtEpoch(
                            body, planetaryRotationModelSettings->getCentralBody( ), "ECLIPJ2000", "NONE", 0.0 ),
                        spice_interface::getBodyGravitationalParameter( planetaryRotationModelSettings->getCentralBody( ) ) );

            double meanAnomalyAtJ2000 = orbital_element_conversions::convertEccentricAnomalyToMeanAnomaly(
                        orbital_element_conversions::convertTrueAnomalyToEccentricAnomaly(
                            bodyKeplerElements( 5 ), bodyKeplerElements( 1 ) ), bodyKeplerElements( 1 ) );

            double meanMotion = orbital_element_conversions::convertSemiMajorAxisToEllipticalMeanMotion(
                        bodyKeplerElements( 0 ), spice_interface::getBodyGravitationalParameter(
                            planetaryRotationModelSettings->getCentralBody( ) ) + spice_interface::getBodyGravitationalParameter( body ) );

            std::shared_ptr< PlanetaryOrientationAngleCalculator > planetaryOrientationAnglesCalculator
                    = std::make_shared< PlanetaryOrientationAngleCalculator >(
                        planetaryRotationModelSettings->getAnglePsiAtEpoch( ),
                        planetaryRotationModelSettings->getAnglePsiRateAtEpoch( ),
                        planetaryRotationModelSettings->getAngleIAtEpoch( ),
                        planetaryRotationModelSettings->getAngleIRateAtEpoch( ),
                        planetaryRotationModelSettings->getAnglePhiAtEpoch( ),
                        planetaryRotationModelSettings->getAnglePhiRateAtEpoch( ),
                        planetaryRotationModelSettings->getCoreFactor( ),
                        planetaryRotationModelSettings->getFreeCoreNutationRate( ),
                        meanMotion, meanAnomalyAtJ2000,
                        planetaryRotationModelSettings->getOriginalFrame( ),
                        planetaryRotationModelSettings->getMeanMotionDirectNutationCorrections( ),
                        planetaryRotationModelSettings->getMeanMotionTimeDependentPhaseNutationCorrections( ),
                        planetaryRotationModelSettings->getTimeDependentPhaseCorrectionFunctions( ),
                        planetaryRotationModelSettings->getRotationRateCorrections( ),
                        planetaryRotationModelSettings->getxPolarMotionCoefficients( ),
                        planetaryRotationModelSettings->getyPolarMotionCoefficients( ) );


            rotationalEphemeris = std::make_shared< PlanetaryRotationModel >(
                        planetaryRotationModelSettings->getAngleN( ),
                        planetaryRotationModelSettings->getAngleJ( ),
                        planetaryOrientationAnglesCalculator,
                        planetaryRotationModelSettings->getOriginalFrame( ), planetaryRotationModelSettings->getTargetFrame( ) );
        }

        break;
    }
    case synchronous_rotation_model:
    {
        std::shared_ptr< SynchronousRotationModelSettings > synchronousRotationSettings =
                std::dynamic_pointer_cast< SynchronousRotationModelSettings >( rotationModelSettings );
        if( synchronousRotationSettings == NULL )
        {
            throw std::runtime_error( "Error, expected synchronous rotation model settings for " + body );
        }
        else
        {
            if( bodies.count( synchronousRotationSettings->getCentralBodyName( ) ) == 0 )
            {
                throw std::runtime_error( "Error, when making synchronous rotation model for " + body + ", central body " +
                                          synchronousRotationSettings->getCentralBodyName( ) + " not found." );
            }
            if( bodies.at( body )->getEphemeris( )->getReferenceFrameOrigin( ) == synchronousRotationSettings->getCentralBodyName( ) )
            {
                if( bodies.at( body )->getEphemeris( )->getReferenceFrameOrientation( ) !=
                        synchronousRotationSettings->getOriginalFrame( ) )
                {
                    throw std::runtime_error( "Error, ephemeris of body " + body + " is in " +
                                              bodies.at( body )->getEphemeris( )->getReferenceFrameOrientation( ) +
                                              " frame when making synchronous rotation model, expected " +
                                              synchronousRotationSettings->getOriginalFrame( ) + " frame." );
                }
            }
            std::shared_ptr< SynchronousRotationalEphemeris > synchronousRotationalEphemeris = std::make_shared< SynchronousRotationalEphemeris >(
                        createRelativeStateFunction( bodies, body, synchronousRotationSettings->getCentralBodyName( ) ),
                        synchronousRotationSettings->getCentralBodyName( ),
                        synchronousRotationSettings->getOriginalFrame( ),
                        synchronousRotationSettings->getTargetFrame( ) );

            rotationalEphemeris = synchronousRotationalEphemeris;
        }
        break;
    }
    case tabulated_rotation_model:
    {
        // Check whether settings for simple rotation model are consistent with its type.
        std::shared_ptr< TabulatedRotationSettings > tabulatedRotationSettings =
                std::dynamic_pointer_cast< TabulatedRotationSettings >( rotationModelSettings );
        if( tabulatedRotationSettings == nullptr )
        {
            throw std::runtime_error(
                        "Error, expected tabulated rotation model settings for " + body );
        }
        else
        {
            // Create and initialize simple rotation model.
            rotationalEphemeris = std::make_shared< TabulatedRotationalEphemeris< double, double > >(
                        interpolators::createOneDimensionalInterpolator(
                            tabulatedRotationSettings->getBodyStateHistory( ),
                            tabulatedRotationSettings->getInterpolatorSettings( ) ),
                        tabulatedRotationSettings->getOriginalFrame( ),
                        tabulatedRotationSettings->getTargetFrame( ) );
        }
        break;
    }
    default:
        throw std::runtime_error(
                    "Error, did not recognize rotation model settings type " +
                    std::to_string( rotationModelSettings->getRotationType( ) ) );
    }

    return rotationalEphemeris;
}

} // namespace simulation_setup

} // namespace tudat
