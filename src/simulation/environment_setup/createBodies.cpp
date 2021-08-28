/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <iostream>
#include <cmath>

#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"

#include "tudat/simulation/environment_setup/createBodies.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
//#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/interface/spice/spiceRotationalEphemeris.h"

namespace tudat
{

namespace simulation_setup
{

using namespace ephemerides;
using namespace gravitation;
using namespace basic_astrodynamics;

void addAerodynamicCoefficientInterface(
        const SystemOfBodies& bodies, const std::string bodyName,
        const std::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings )
{
    if( bodies.count( bodyName ) == 0 )
    {
        throw std::runtime_error( "Error when setting aerodynamic coefficients for body "+ bodyName + ", body is not found in system of bodies" );
    }
    bodies.at( bodyName )->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, bodyName) );
}

void addRadiationPressureInterface(
        const SystemOfBodies& bodies, const std::string bodyName,
        const std::shared_ptr< RadiationPressureInterfaceSettings > radiationPressureSettings )
{
    if( bodies.count( bodyName ) == 0 )
    {
        throw std::runtime_error( "Error when setting radiation pressure interface for body "+ bodyName + ", body is not found in system of bodies" );
    }
    bodies.at( bodyName )->setRadiationPressureInterface(
                radiationPressureSettings->getSourceBody( ), createRadiationPressureInterface(
                    radiationPressureSettings, bodyName, bodies ) );
}

void setSimpleRotationSettingsFromSpice(
        const BodyListSettings& bodySettings, const std::string& bodyName, const double spiceEvaluationTime )
{
    if( bodySettings.count( bodyName ) == 0 )
    {
        throw std::runtime_error( "Error when setting simple rotation model settings for body " +
                                  bodyName + ", no settings found for this body." );
    }

    Eigen::Quaterniond rotationAtReferenceTime =
            spice_interface::computeRotationQuaternionBetweenFrames(
                bodySettings.getFrameOrientation( ), "IAU_" + bodyName, spiceEvaluationTime );
    double rotationRateAtReferenceTime =
            spice_interface::getAngularVelocityVectorOfFrameInOriginalFrame(
                bodySettings.getFrameOrientation( ), "IAU_" + bodyName, spiceEvaluationTime ).norm( );

    bodySettings.at( bodyName )->rotationModelSettings =
            std::make_shared< SimpleRotationModelSettings >(
                bodySettings.getFrameOrientation( ), "IAU_" + bodyName, rotationAtReferenceTime,
                spiceEvaluationTime, rotationRateAtReferenceTime );
}

void addEmptyTabulateEphemeris(
        const SystemOfBodies& bodies, const std::string& bodyName, const std::string& ephemerisOrigin )
{
    if( bodies.count( bodyName ) ==  0 )
    {
        throw std::runtime_error( "Error when setting empty tabulated ephemeris for body " + bodyName + ", no such body found" );
    }
    std::string ephemerisOriginToUse = ( ephemerisOrigin == "" ) ? bodies.getFrameOrigin( ) : ephemerisOrigin;
    bodies.at( bodyName )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                            std::shared_ptr< interpolators::OneDimensionalInterpolator
                                            < double, Eigen::Vector6d > >( ), ephemerisOriginToUse, bodies.getFrameOrientation( ) ) );
}


//! Function that determines the order in which bodies are to be created
std::vector< std::pair< std::string, std::shared_ptr< BodySettings > > > determineBodyCreationOrder(
        const std::map< std::string, std::shared_ptr< BodySettings > >& bodySettings )
{
    std::vector< std::pair< std::string, std::shared_ptr< BodySettings > > > outputVector;

    // Create vector of pairs (body name and body settings) that is to be created.
    for( std::map< std::string, std::shared_ptr< BodySettings > >::const_iterator bodyIterator
         = bodySettings.begin( );
         bodyIterator != bodySettings.end( ); bodyIterator ++ )
    {
        outputVector.push_back( std::make_pair( bodyIterator->first, bodyIterator->second ) );
    }

    return outputVector;
}


//! Function to create a map of bodies objects.
SystemOfBodies createSystemOfBodies(
        const BodyListSettings& bodySettings )
{
    std::vector< std::pair< std::string, std::shared_ptr< BodySettings > > > orderedBodySettings
            = determineBodyCreationOrder( bodySettings.getMap( ) );

    // Declare map of bodies that is to be returned.
    SystemOfBodies bodyList = SystemOfBodies(
                bodySettings.getFrameOrigin( ), bodySettings.getFrameOrientation( ) );

    // Create empty body objects.
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        bodyList.createEmptyBody( orderedBodySettings.at( i ).first, false );
    }

    // Define constant mass for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        const double constantMass = orderedBodySettings.at( i ).second->constantMass;
        if ( constantMass == constantMass )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setConstantBodyMass( constantMass );
        }
    }

    // Create ephemeris objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->ephemerisSettings != nullptr )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setEphemeris(
                        createBodyEphemeris( orderedBodySettings.at( i ).second->ephemerisSettings,
                                             orderedBodySettings.at( i ).first ) );
        }
    }

    // Create atmosphere model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->atmosphereSettings != nullptr )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setAtmosphereModel(
                        createAtmosphereModel( orderedBodySettings.at( i ).second->atmosphereSettings,
                                               orderedBodySettings.at( i ).first ) );
        }
    }

    // Create body shape model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->shapeModelSettings != nullptr )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setShapeModel(
                        createBodyShapeModel( orderedBodySettings.at( i ).second->shapeModelSettings,
                                              orderedBodySettings.at( i ).first ) );
        }
    }

    // Create rotation model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->rotationModelSettings != nullptr )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setRotationalEphemeris(
                        createRotationModel( orderedBodySettings.at( i ).second->rotationModelSettings,
                                             orderedBodySettings.at( i ).first, bodyList ) );
        }
    }

    // Create gravity field model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->gravityFieldSettings != nullptr )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setGravityFieldModel(
                        createGravityFieldModel( orderedBodySettings.at( i ).second->gravityFieldSettings,
                                                 orderedBodySettings.at( i ).first, bodyList,
                                                 orderedBodySettings.at( i ).second->gravityFieldVariationSettings ) );
        }
    }

    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->gravityFieldVariationSettings.size( ) > 0 )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setGravityFieldVariationSet(
                        createGravityFieldModelVariationsSet(
                            orderedBodySettings.at( i ).first, bodyList,
                            orderedBodySettings.at( i ).second->gravityFieldVariationSettings ) );
        }
    }

    // Create aerodynamic coefficient interface objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        if( orderedBodySettings.at( i ).second->aerodynamicCoefficientSettings != nullptr )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface(
                            orderedBodySettings.at( i ).second->aerodynamicCoefficientSettings,
                            orderedBodySettings.at( i ).first ) );
        }
    }


    // Create radiation pressure coefficient objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        std::map< std::string, std::shared_ptr< RadiationPressureInterfaceSettings > >
                radiationPressureSettings
                = orderedBodySettings.at( i ).second->radiationPressureSettings;
        for( std::map< std::string, std::shared_ptr< RadiationPressureInterfaceSettings > >::iterator
             radiationPressureSettingsIterator = radiationPressureSettings.begin( );
             radiationPressureSettingsIterator != radiationPressureSettings.end( );
             radiationPressureSettingsIterator++ )
        {
            bodyList.at( orderedBodySettings.at( i ).first )->setRadiationPressureInterface(
                        radiationPressureSettingsIterator->first,
                        createRadiationPressureInterface(
                            radiationPressureSettingsIterator->second,
                            orderedBodySettings.at( i ).first, bodyList ) );
        }

    }

    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        for( unsigned int j = 0; j < orderedBodySettings.at( i ).second->groundStationSettings.size( ); j++ )
        {
            createGroundStation( bodyList.at( orderedBodySettings.at( i ).first ), orderedBodySettings.at( i ).first,
                     orderedBodySettings.at( i ).second->groundStationSettings.at( j ) );
        }
    }

    bodyList.processBodyFrameDefinitions( );

    return bodyList;

}


simulation_setup::SystemOfBodies createSimplifiedSystemOfBodies(const double initialTime)
{
    using namespace ephemerides;
    using namespace gravitation;

    SystemOfBodies bodies("SSB","ECLIPJ2000");
    bodies.createEmptyBody( "Sun" );
    bodies.createEmptyBody( "Mercury" );
    bodies.createEmptyBody( "Venus" );
    bodies.createEmptyBody( "Earth" );
    bodies.createEmptyBody( "Mars" );
    bodies.createEmptyBody( "Jupiter" );
    bodies.createEmptyBody( "Saturn" );
    bodies.createEmptyBody( "Uranus" );
    bodies.createEmptyBody( "Neptune" );
    bodies.createEmptyBody( "Pluto" );

    bodies.getBody( "Sun" )->setEphemeris( std::make_shared< ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "ECLIPJ2000" ) );
    bodies.getBody( "Mercury" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Mercury") );
    bodies.getBody( "Venus" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Venus" ) );
    bodies.getBody( "Earth" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Earth" ) );
    bodies.getBody( "Mars" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Mars" ) );
    bodies.getBody( "Jupiter" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Jupiter" ) );
    bodies.getBody( "Saturn" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Saturn" ) );
    bodies.getBody( "Uranus" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Uranus" ) );
    bodies.getBody( "Neptune" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Neptune" ) );
    bodies.getBody( "Pluto" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Pluto" ) );

    bodies.getBody( "Sun" )->setGravityFieldModel( std::make_shared< GravityFieldModel >(celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER ) );
    bodies.getBody( "Mercury" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( celestial_body_constants::MERCURY_GRAVITATIONAL_PARAMETER ) );
    bodies.getBody( "Venus" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( celestial_body_constants::VENUS_GRAVITATIONAL_PARAMETER ) );
    bodies.getBody( "Earth" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( celestial_body_constants::EARTH_GRAVITATIONAL_PARAMETER ) );
    bodies.getBody( "Mars" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( celestial_body_constants::MARS_GRAVITATIONAL_PARAMETER ) );
    bodies.getBody( "Jupiter" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( celestial_body_constants::JUPITER_GRAVITATIONAL_PARAMETER ) );
    bodies.getBody( "Saturn" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( celestial_body_constants::SATURN_GRAVITATIONAL_PARAMETER ) );
    bodies.getBody( "Uranus" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( celestial_body_constants::URANUS_GRAVITATIONAL_PARAMETER ) );
    bodies.getBody( "Neptune" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( celestial_body_constants::NEPTUNE_GRAVITATIONAL_PARAMETER ) );
    bodies.getBody( "Pluto" )->setGravityFieldModel( std::make_shared< GravityFieldModel >( celestial_body_constants::PLUTO_GRAVITATIONAL_PARAMETER ) );

    bodies.getBody( "Earth" )->setShapeModel( std::make_shared< SphericalBodyShapeModel >(celestial_body_constants::EARTH_EQUATORIAL_RADIUS ) );

    //bodies.getBody( "Earth" )->setRotationalEphemeris( std::make_shared< SimpleRotationalEphemeris >(Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ),
    //                                                                                                 7.2921159E-5, 0.0, "ECLIPJ2000", "Earth_Fixed" ) );

    // Calculate position of rotation axis at initialTime, with respect to J2000 frame. Values from
    // "Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009", B.A. Archinal et al.(2011)
    /*
    const double initialTimeDays = initialTime / 86400;
    const double initialTimeCenturies = initialTimeDays / 36525;
    const double poleRightAscension = (0 - 0.641 * initialTimeCenturies) * mathematical_constants::PI/180;
    const double poleDeclination = (90 - 0.557 * initialTimeCenturies) * mathematical_constants::PI/180;
    const double initialRotationAngle = 190.147 * mathematical_constants::PI/180;
    const double rotationRate = 360.985235 * mathematical_constants::PI/180 / physical_constants::JULIAN_DAY;

    bodies.getBody( "Earth" )->setRotationalEphemeris(
            std::make_shared< SimpleRotationalEphemeris >(
                    reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion (
                            poleDeclination, poleRightAscension, initialRotationAngle),
                            rotationRate, initialTime, "J2000", "Earth_Fixed" ) );
    */

    std::string originalFrame = "ECLIPJ2000";
    std::string targetFrame = "IAU_Earth";
    bodies.getBody( "Earth" )->setRotationalEphemeris(std::make_shared< SpiceRotationalEphemeris >( originalFrame, targetFrame,
                                                                                                    initialTime/86400 + 2451545.0));

    return bodies;
}

/*
Eigen::Matrix3d rotationMatrixAroundX (const double angle)
    {
        return (Eigen::Matrix3d() <<
            1, 0, 0,
            0, std::cos(angle), std::sin(angle),
            0, -std::sin(angle), std::cos(angle)).finished() ;
    }

Eigen::Matrix3d rotationMatrixAroundZ (const double angle)
    {
        return (Eigen::Matrix3d() <<
            std::cos(angle), std::sin(anglea), 0,
            -std::sin(angle), std::cos(angle), 0,
            0, 0, 1).finished() ;
    }
*/


} // namespace simulation_setup

} // namespace tudat
