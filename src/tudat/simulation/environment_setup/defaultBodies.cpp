/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define DEFAULT_MERCURY_GRAVITY_FIELD_SETTINGS std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( jgmess160a )
#define DEFAULT_VENUS_GRAVITY_FIELD_SETTINGS std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( shgj180u )
#define DEFAULT_EARTH_GRAVITY_FIELD_SETTINGS std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( goco05c, 200 )
#define DEFAULT_MOON_GRAVITY_FIELD_SETTINGS std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( gggrx1200, 200 )
#define DEFAULT_MARS_GRAVITY_FIELD_SETTINGS std::make_shared< FromFileSphericalHarmonicsGravityFieldSettings >( jgmro120d )

#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createRotationModel.h"
#include "tudat/simulation/environment_setup/createGravityField.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"

namespace tudat
{
namespace simulation_setup
{

//! Function to create default settings for a body's atmosphere model.
std::shared_ptr< AtmosphereSettings > getDefaultAtmosphereModelSettings( const std::string& bodyName,
                                                                         const double initialTime,
                                                                         const double finalTime )
{
    std::shared_ptr< AtmosphereSettings > atmosphereSettings;

    // A default atmosphere is only implemented for Earth.
    if( bodyName == "Earth" )
    {
        std::string atmosphereTableFile = paths::getAtmosphereTablesPath( ) + "/USSA1976Until100kmPer100mUntil1000kmPer1000m.dat";
        atmosphereSettings = std::make_shared< TabulatedAtmosphereSettings >( atmosphereTableFile );
    }

    return atmosphereSettings;
}

std::shared_ptr< RadiationSourceModelSettings > getKnockeEarthRadiationPressureSettings( )
{
    // Model from Knocke (1988)
    return extendedRadiationSourceModelSettings(
            { albedoPanelRadiosityModelSettings( KnockeTypeSurfacePropertyDistributionModel::albedo_knocke, "Sun" ),
              delayedThermalPanelRadiosityModelSettings( KnockeTypeSurfacePropertyDistributionModel::emissivity_knocke, "Sun" ) },
            { 6, 12 } );
}
//! Function to create default settings for a body's radiation source model.
std::shared_ptr< RadiationSourceModelSettings > getDefaultRadiationSourceModelSettings( const std::string& bodyName,
                                                                                        const double initialTime,
                                                                                        const double finalTime )
{
    std::shared_ptr< RadiationSourceModelSettings > radiationSourceModelSettings;

    if( bodyName == "Sun" )
    {
        radiationSourceModelSettings =
                isotropicPointRadiationSourceModelSettings( constantLuminosityModelSettings( celestial_body_constants::SUN_LUMINOSITY ) );
    }
    // Do not use Earth default source model since they require Sun to be present
    //    else if( bodyName == "Earth" )
    //    {
    //        // Model from Knocke (1988)
    //        radiationSourceModelSettings =
    //                extendedRadiationSourceModelSettings({
    //                    albedoPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::albedo_knocke, "Sun"),
    //                    delayedThermalPanelRadiosityModelSettings(KnockeTypeSurfacePropertyDistributionModel::emissivity_knocke, "Sun")
    //                }, {6, 12});
    //    }

    return radiationSourceModelSettings;
}

//! Function to create default settings for a body's ephemeris.
std::shared_ptr< EphemerisSettings > getDefaultEphemerisSettings( const std::string& bodyName,
                                                                  const std::string& baseFrameOrientation,
                                                                  const std::string& originatingNameBodyName )
{
    std::string bodyNameToUse = ( originatingNameBodyName == "" ) ? bodyName : originatingNameBodyName;
    if( originatingNameBodyName == "Uranus" || originatingNameBodyName == "Neptune" || originatingNameBodyName == "Pluto" )
    {
        bodyNameToUse += "_BARYCENTER";
    }
    // Create settings for an interpolated Spice ephemeris.
    return std::make_shared< DirectSpiceEphemerisSettings >( "SSB", baseFrameOrientation, bodyNameToUse );
}

//! Function to create default settings for a body's ephemeris.
std::shared_ptr< EphemerisSettings > getDefaultEphemerisSettings( const std::string& bodyName,
                                                                  const double initialTime,
                                                                  const double finalTime,
                                                                  const std::string& baseFrameOrientation,
                                                                  const std::string& originatingNameBodyName,
                                                                  const double timeStep )
{
    std::string bodyNameToUse = ( originatingNameBodyName == "" ) ? bodyName : originatingNameBodyName;
    if( originatingNameBodyName == "Uranus" || originatingNameBodyName == "Neptune" || originatingNameBodyName == "Pluto" )
    {
        bodyNameToUse += "_BARYCENTER";
    }
    // Create settings for an interpolated Spice ephemeris.
    return std::make_shared< InterpolatedSpiceEphemerisSettings >( initialTime,
                                                                   finalTime,
                                                                   timeStep,
                                                                   "SSB",
                                                                   baseFrameOrientation,
                                                                   std::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 ),
                                                                   bodyNameToUse );
}

//! Function to create default settings for a body's gravity field model.
std::shared_ptr< GravityFieldSettings > getDefaultGravityFieldSettings( const std::string& bodyName,
                                                                        const double initialTime,
                                                                        const double finalTime )
{
    if( bodyName == "Earth" )
    {
        return DEFAULT_EARTH_GRAVITY_FIELD_SETTINGS;
    }
    else if( bodyName == "Moon" )
    {
        return DEFAULT_MOON_GRAVITY_FIELD_SETTINGS;
    }
    else if( bodyName == "Mars" )
    {
        return DEFAULT_MARS_GRAVITY_FIELD_SETTINGS;
    }
    else if( bodyName == "Venus" )
    {
        return DEFAULT_VENUS_GRAVITY_FIELD_SETTINGS;
    }
    else if( bodyName == "Mercury" )
    {
        return DEFAULT_MERCURY_GRAVITY_FIELD_SETTINGS;
    }
    else if( bodyName == "Jupiter" )
    {
        double jupiterJ2 = 14.696572E-3;
        double jupiterJ3 = -0.042E-6;
        double jupiterJ4 = -586.609E-6;
        double jupiterJ5 = -0.069E-6;
        double jupiterJ6 = 34.198E-6;
        double jupiterJ7 = 0.124E-6;
        double jupiterJ8 = -2.426E-6;

        Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 21, 21 );
        Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 21, 21 );

        cosineCoefficients( 0, 0 ) = 1.0;
        cosineCoefficients( 2, 0 ) = -jupiterJ2 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
        cosineCoefficients( 3, 0 ) = -jupiterJ3 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 3, 0 );
        cosineCoefficients( 4, 0 ) = -jupiterJ4 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 4, 0 );
        cosineCoefficients( 5, 0 ) = -jupiterJ5 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 5, 0 );
        cosineCoefficients( 6, 0 ) = -jupiterJ6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 6, 0 );
        cosineCoefficients( 7, 0 ) = -jupiterJ7 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 7, 0 );
        cosineCoefficients( 8, 0 ) = -jupiterJ8 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 8, 0 );

        return std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                1.266865341960128E17, 71492.0E3, cosineCoefficients, sineCoefficients, "IAU_Jupiter" );  // Mass from jup329.cmt
    }
    else if( bodyName == "Io" )
    {
        Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
        Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );

        cosineCoefficients( 0, 0 ) = 1.0;
        cosineCoefficients( 2, 0 ) = -1845.9E-6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
        cosineCoefficients( 2, 2 ) = 553.7E-6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );

        return std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                5.959924010272514E+12, 1821.6E3, cosineCoefficients, sineCoefficients, "IAU_Io" );
    }
    else if( bodyName == "Europa" )
    {
        Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
        Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );

        cosineCoefficients( 0, 0 ) = 1.0;
        cosineCoefficients( 2, 0 ) = -435.5E-6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
        cosineCoefficients( 2, 2 ) = 131.0E-6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );

        return std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                3.202739815114734E+12, 1565.0E3, cosineCoefficients, sineCoefficients, "IAU_Europa" );
    }
    else if( bodyName == "Ganymede" )
    {
        Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 33, 33 );
        Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 33, 33 );

        cosineCoefficients( 0, 0 ) = 1.0;
        cosineCoefficients( 2, 0 ) = -127.8E-6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
        cosineCoefficients( 2, 2 ) = 38.3E-6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );

        return std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                9.887819980080976E+12, 2634.0E3, cosineCoefficients, sineCoefficients, "IAU_Ganymede" );
    }
    else if( bodyName == "Callisto" )
    {
        Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );
        Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 13, 13 );

        cosineCoefficients( 0, 0 ) = 1.0;
        cosineCoefficients( 2, 0 ) = -32.7E-6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 );
        cosineCoefficients( 2, 2 ) = 10.2E-6 / basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );

        return std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                7.179304867611079E+12, 2410.3E3, cosineCoefficients, sineCoefficients, "IAU_Callisto" );
    }
    else
    {
        // Create settings for a point mass gravity with data from Spice
        return centralGravityFromSpiceSettings( );
    }
}

//! Function to create default settings from which to create a single body object.
std::shared_ptr< RotationModelSettings > getDefaultRotationModelSettings( const std::string& bodyName,
                                                                          const double initialTime,
                                                                          const double finalTime,
                                                                          const std::string& baseFrameOrientation,
                                                                          const std::string& spiceBodyName )
{
    TUDAT_UNUSED_PARAMETER( initialTime );
    TUDAT_UNUSED_PARAMETER( finalTime );

    std::string spiceFrameName = "IAU_" + bodyName;
    if( spiceBodyName != "" )
    {
        spiceFrameName = "IAU_" + spiceBodyName;
    }
    // Create settings for a rotation model taken directly from Spice.
    return std::make_shared< SpiceRotationModelSettings >( baseFrameOrientation, "IAU_" + bodyName, spiceFrameName );
}

double marsTimeDependentPhaseAngleCorrectionFunction( const double secondsSinceJ2000 )
{
    double centuriesSinceJ2000 = secondsSinceJ2000 / ( 100.0 * physical_constants::JULIAN_YEAR );
    return ( 142.0 + 1.3 * centuriesSinceJ2000 ) * mathematical_constants::PI / 180.0;
}

// Mars orientatio// Implementation of the fully customizable version
std::shared_ptr< RotationModelSettings > getHighAccuracyMarsRotationModel(
        const std::string& baseFrameOrientation,
        const std::string& targetFrameOrientation,
        const double angleN,
        const double angleJ,
        const double anglePsiAtEpoch,
        const double anglePsiRateAtEpoch,
        const std::map< double, std::pair< double, double > >& nutationCorrectionSettings,
        const std::vector< std::map< double, std::pair< double, double > > >& meanMotionTimeDependentPhaseNutationCorrections,
        const std::map< double, std::pair< double, double > >& rotationRateCorrections,
        const std::map< double, std::pair< double, double > >& xPolarMotionCoefficients,
        const std::map< double, std::pair< double, double > >& yPolarMotionCoefficients )
{
    std::shared_ptr< RotationModelSettings > rotationModelSettings;

    using namespace tudat::unit_conversions;

    double milliArcSecondToRadian = mathematical_constants::PI / ( 180.0 * 1000.0 * 3600.0 );

    // Set up phase correction functions (not customizable via parameters)
    std::vector< std::function< double( const double ) > > timeDependentPhaseCorrectionFunctions;
    timeDependentPhaseCorrectionFunctions.push_back(
            std::bind( &tudat::simulation_setup::marsTimeDependentPhaseAngleCorrectionFunction, std::placeholders::_1 ) );

    // Create the rotation model settings with all provided parameters
    rotationModelSettings =
            std::make_shared< PlanetaryRotationModelSettings >( angleN,
                                                                angleJ,
                                                                anglePsiAtEpoch,
                                                                anglePsiRateAtEpoch,
                                                                convertDegreesToRadians( 25.1893823 ),
                                                                ( -2.0 * milliArcSecondToRadian ) / physical_constants::JULIAN_YEAR,
                                                                convertDegreesToRadians( 133.386277 ),
                                                                convertDegreesToRadians( 350.891985307 ) / physical_constants::JULIAN_DAY,
                                                                0.07,
                                                                convertDegreesToRadians( -1.5 ) / physical_constants::JULIAN_DAY,
                                                                baseFrameOrientation,
                                                                targetFrameOrientation,
                                                                "Sun",
                                                                nutationCorrectionSettings,
                                                                meanMotionTimeDependentPhaseNutationCorrections,
                                                                timeDependentPhaseCorrectionFunctions,
                                                                rotationRateCorrections,
                                                                xPolarMotionCoefficients,
                                                                yPolarMotionCoefficients );

    return rotationModelSettings;
}

// Implementation of the version that allows customizing only the basic angles
std::shared_ptr< RotationModelSettings > getHighAccuracyMarsRotationModel( const std::string& baseFrameOrientation,
                                                                           const std::string& targetFrameOrientation,
                                                                           const double angleN,
                                                                           const double angleJ,
                                                                           const double anglePsiAtEpoch,
                                                                           const double anglePsiRateAtEpoch )
{
    using namespace tudat::unit_conversions;

    double milliArcSecondToRadian = mathematical_constants::PI / ( 180.0 * 1000.0 * 3600.0 );

    // Create default nutation correction settings
    std::map< double, std::pair< double, double > > nutationCorrectionSettings;
    nutationCorrectionSettings[ 0.0 ] = std::make_pair( -1.4 * milliArcSecondToRadian, 0.0 );
    nutationCorrectionSettings[ 1.0 ] = std::make_pair( -0.4 * milliArcSecondToRadian, -632.6 * milliArcSecondToRadian );
    nutationCorrectionSettings[ 2.0 ] = std::make_pair( 0.0, -44.2 * milliArcSecondToRadian );
    nutationCorrectionSettings[ 3.0 ] = std::make_pair( 0.0, -4.0 * milliArcSecondToRadian );

    // Create default mean motion time dependent phase nutation corrections
    std::vector< std::map< double, std::pair< double, double > > > meanMotionTimeDependentPhaseNutationCorrections;
    std::map< double, std::pair< double, double > > meanMotionTimeDependentPhaseNutationCorrection;
    meanMotionTimeDependentPhaseNutationCorrection[ 1.0 ] =
            std::make_pair( -49.1 * milliArcSecondToRadian, -104.5 * milliArcSecondToRadian );
    meanMotionTimeDependentPhaseNutationCorrection[ 2.0 ] =
            std::make_pair( 515.7 * milliArcSecondToRadian, 1097.0 * milliArcSecondToRadian );
    meanMotionTimeDependentPhaseNutationCorrection[ 3.0 ] =
            std::make_pair( 112.8 * milliArcSecondToRadian, 240.1 * milliArcSecondToRadian );
    meanMotionTimeDependentPhaseNutationCorrection[ 4.0 ] = std::make_pair( 19.2 * milliArcSecondToRadian, 40.9 * milliArcSecondToRadian );
    meanMotionTimeDependentPhaseNutationCorrection[ 5.0 ] = std::make_pair( 3.0 * milliArcSecondToRadian, 6.5 * milliArcSecondToRadian );
    meanMotionTimeDependentPhaseNutationCorrection[ 6.0 ] = std::make_pair( 0.4 * milliArcSecondToRadian, 1.0 * milliArcSecondToRadian );
    meanMotionTimeDependentPhaseNutationCorrections.push_back( meanMotionTimeDependentPhaseNutationCorrection );

    // Create default rotation rate corrections
    std::map< double, std::pair< double, double > > rotationRateCorrections;
    rotationRateCorrections[ 1.0 ] =
            std::make_pair( 481.0 * milliArcSecondToRadian, -155.0 * milliArcSecondToRadian - 176 * milliArcSecondToRadian );
    rotationRateCorrections[ 2.0 ] =
            std::make_pair( -103.0 * milliArcSecondToRadian, -93.0 * milliArcSecondToRadian - 8 * milliArcSecondToRadian );
    rotationRateCorrections[ 3.0 ] =
            std::make_pair( -35.0 * milliArcSecondToRadian, -3.0 * milliArcSecondToRadian - 1 * milliArcSecondToRadian );
    rotationRateCorrections[ 4.0 ] = std::make_pair( -10.0 * milliArcSecondToRadian, -8.0 * milliArcSecondToRadian );

    // Create default polar motion coefficients
    std::map< double, std::pair< double, double > > xPolarMotionCoefficients;
    xPolarMotionCoefficients[ 1.0 ] = std::make_pair( 2.8 * milliArcSecondToRadian * std::sin( convertDegreesToRadians( 46.5 ) ),
                                                      2.8 * milliArcSecondToRadian * std::cos( convertDegreesToRadians( 46.5 ) ) );
    xPolarMotionCoefficients[ 2.0 ] = std::make_pair( 8.9 * milliArcSecondToRadian * std::sin( convertDegreesToRadians( -150.1 ) ),
                                                      8.9 * milliArcSecondToRadian * std::cos( convertDegreesToRadians( -150.1 ) ) );
    xPolarMotionCoefficients[ 3.0 ] = std::make_pair( 0.0, 0.0 );
    xPolarMotionCoefficients[ 4.0 ] = std::make_pair( 0.0, 0.0 );
    xPolarMotionCoefficients[ 3.34 ] = std::make_pair( 0.0, 50.0 * milliArcSecondToRadian );  // Mars's Chandler wobble T=205 dd

    std::map< double, std::pair< double, double > > yPolarMotionCoefficients;
    yPolarMotionCoefficients[ 1.0 ] = std::make_pair( 11.7 * milliArcSecondToRadian * std::sin( convertDegreesToRadians( 118.7 ) ),
                                                      11.7 * milliArcSecondToRadian * std::cos( convertDegreesToRadians( 118.7 ) ) );
    yPolarMotionCoefficients[ 2.0 ] = std::make_pair( 3.9 * milliArcSecondToRadian * std::sin( convertDegreesToRadians( 172.5 ) ),
                                                      3.9 * milliArcSecondToRadian * std::cos( convertDegreesToRadians( 118.7 ) ) );
    yPolarMotionCoefficients[ 3.0 ] = std::make_pair( 0.0, 0.0 );
    yPolarMotionCoefficients[ 4.0 ] = std::make_pair( 0.0, 0.0 );
    yPolarMotionCoefficients[ 3.34 ] = std::make_pair( 0.0, 50.0 * milliArcSecondToRadian );  // Mars's Chandler wobble T=205 dd

    // Call the fully customizable version with the default correction parameters
    return getHighAccuracyMarsRotationModel( baseFrameOrientation,
                                             targetFrameOrientation,
                                             angleN,
                                             angleJ,
                                             anglePsiAtEpoch,
                                             anglePsiRateAtEpoch,
                                             nutationCorrectionSettings,
                                             meanMotionTimeDependentPhaseNutationCorrections,
                                             rotationRateCorrections,
                                             xPolarMotionCoefficients,
                                             yPolarMotionCoefficients );
}

// Implementation of the original function - uses all default values
std::shared_ptr< RotationModelSettings > getHighAccuracyMarsRotationModel( const std::string& baseFrameOrientation,
                                                                           const std::string& targetFrameOrientation )
{
    using namespace tudat::unit_conversions;

    double milliArcSecondToRadian = mathematical_constants::PI / ( 180.0 * 1000.0 * 3600.0 );

    // Default values from Konopliv et al. 2016
    const double defaultAngleN = convertDegreesToRadians( 3.37919183 );
    const double defaultAngleJ = convertDegreesToRadians( 24.67682669 );
    const double defaultAnglePsiAtEpoch = convertDegreesToRadians( 81.9683988 );
    const double defaultAnglePsiRateAtEpoch = ( -7608.3 * milliArcSecondToRadian ) / physical_constants::JULIAN_YEAR;

    // Call the overloaded function with default values
    return getHighAccuracyMarsRotationModel( baseFrameOrientation,
                                             targetFrameOrientation,
                                             defaultAngleN,
                                             defaultAngleJ,
                                             defaultAnglePsiAtEpoch,
                                             defaultAnglePsiRateAtEpoch );
}

//! Function to create default settings for a body's shape model.
std::shared_ptr< BodyShapeSettings > getDefaultBodyShapeSettings( const std::string& bodyName,
                                                                  const double initialTime,
                                                                  const double finalTime )
{
    TUDAT_UNUSED_PARAMETER( initialTime );
    TUDAT_UNUSED_PARAMETER( finalTime );

    return std::make_shared< SphericalBodyShapeSettings >( spice_interface::getAverageRadius( bodyName ) );
}

//! Function to create default settings for a body's rotation model.
std::shared_ptr< BodySettings > getDefaultSingleBodySettings( const std::string& bodyName,
                                                              const double initialTime,
                                                              const double finalTime,
                                                              const std::string& baseFrameOrientation,
                                                              const double timeStep )
{
    return getDefaultSingleAlternateNameBodySettings( bodyName, bodyName, initialTime, finalTime, baseFrameOrientation, timeStep );
}

std::shared_ptr< BodySettings > getDefaultSingleAlternateNameBodySettings( const std::string& bodyName,
                                                                           const std::string& originatingName,
                                                                           const double initialTime,
                                                                           const double finalTime,
                                                                           const std::string& baseFrameOrientation,
                                                                           const double timeStep )
{
    std::shared_ptr< BodySettings > singleBodySettings = std::make_shared< BodySettings >( );

    // Get default settings for each of the environment models in the body.
    singleBodySettings->atmosphereSettings = getDefaultAtmosphereModelSettings( originatingName, initialTime, finalTime );
    singleBodySettings->radiationSourceModelSettings = getDefaultRadiationSourceModelSettings( bodyName, initialTime, finalTime );
    singleBodySettings->rotationModelSettings =
            getDefaultRotationModelSettings( bodyName, initialTime, finalTime, baseFrameOrientation, originatingName );

    if( !( std::isnan( initialTime ) == std::isnan( finalTime ) ) )
    {
        throw std::runtime_error( "Error when getting default body settings, some but not all input times are NaN" );
    }
    else if( std::isnan( initialTime ) )
    {
        singleBodySettings->ephemerisSettings = getDefaultEphemerisSettings( originatingName, baseFrameOrientation, originatingName );
    }
    else
    {
        if( std::isnan( timeStep ) )
        {
            throw std::runtime_error( "Error when getting default body settings, time step is NaN" );
        }

        singleBodySettings->ephemerisSettings =
                getDefaultEphemerisSettings( originatingName, initialTime, finalTime, baseFrameOrientation, originatingName, timeStep );
    }

    singleBodySettings->gravityFieldSettings = getDefaultGravityFieldSettings( originatingName, initialTime, finalTime );
    if( std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >( singleBodySettings->gravityFieldSettings ) != nullptr )
    {
        std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >( singleBodySettings->gravityFieldSettings )
                ->resetAssociatedReferenceFrame( singleBodySettings->rotationModelSettings->getTargetFrame( ) );
    }
    singleBodySettings->shapeModelSettings = getDefaultBodyShapeSettings( originatingName, initialTime, finalTime );

    return singleBodySettings;
}

std::shared_ptr< BodySettings > getDefaultSingleBodySettings( const std::string& bodyName, const std::string& baseFrameOrientation )
{
    return getDefaultSingleBodySettings( bodyName, TUDAT_NAN, TUDAT_NAN, baseFrameOrientation );
}

std::shared_ptr< BodySettings > getDefaultSingleAlternateNameBodySettings( const std::string& bodyName,
                                                                           const std::string& originatingName,
                                                                           const std::string& baseFrameOrientation )
{
    return getDefaultSingleAlternateNameBodySettings( bodyName, originatingName, TUDAT_NAN, TUDAT_NAN, baseFrameOrientation );
}

//! Function to create default settings from which to create a set of body objects.
BodyListSettings getDefaultBodySettings( const std::vector< std::string >& bodies,
                                         const double initialTime,
                                         const double finalTime,
                                         const std::string baseFrameOrigin,
                                         const std::string baseFrameOrientation,
                                         const double timeStep )
{
    std::map< std::string, std::shared_ptr< BodySettings > > settingsMap;

    // Iterative over all bodies and get default settings.
    for( unsigned int i = 0; i < bodies.size( ); i++ )
    {
        settingsMap[ bodies.at( i ) ] = getDefaultSingleBodySettings(
                bodies.at( i ), initialTime - 10.0 * timeStep, finalTime + 10.0 * timeStep, baseFrameOrientation, timeStep );
    }
    return BodyListSettings( settingsMap, baseFrameOrigin, baseFrameOrientation );
}

//! Function to create default settings from which to create a set of body objects, without stringent limitations on
//! time-interval of validity of environment.
BodyListSettings getDefaultBodySettings( const std::vector< std::string >& bodies,
                                         const std::string baseFrameOrigin,
                                         const std::string baseFrameOrientation )
{
    std::map< std::string, std::shared_ptr< BodySettings > > settingsMap;

    // Iterative over all bodies and get default settings.
    for( unsigned int i = 0; i < bodies.size( ); i++ )
    {
        settingsMap[ bodies.at( i ) ] = getDefaultSingleBodySettings( bodies.at( i ), TUDAT_NAN, TUDAT_NAN, baseFrameOrientation );
    }
    return BodyListSettings( settingsMap, baseFrameOrigin, baseFrameOrientation );
}

}  // namespace simulation_setup

}  // namespace tudat
