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

void addRotationModel(
        const SystemOfBodies& bodies, const std::string bodyName,
        const std::shared_ptr< RotationModelSettings > rotationModelSettings )
{
    if( bodies.count( bodyName ) == 0 )
    {
        throw std::runtime_error( "Error when setting rotation model for body "+ bodyName + ", body is not found in system of bodies" );
    }
    bodies.at( bodyName )->setRotationalEphemeris( createRotationModel(
                    rotationModelSettings, bodyName, bodies ) );
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

void addEmptyTabulatedEphemeris(
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

    bodies.processBodyFrameDefinitions( );
}

void addEmptyTabulatedRotationalEphemeris(
        const SystemOfBodies& bodies, const std::string& bodyName, const std::string& bodyFixedFrameName )
{
    if( bodies.count( bodyName ) ==  0 )
    {
        throw std::runtime_error( "Error when setting empty tabulated rotational ephemeris for body " + bodyName + ", no such body found" );
    }
    std::string bodyFixedFrameNameToUse = ( bodyFixedFrameName == "" ) ? ( bodyName + "_fixed" ) : bodyFixedFrameName;

    bodies.at( bodyName )->setRotationalEphemeris( std::make_shared< ephemerides::TabulatedRotationalEphemeris< > >(
                                            std::shared_ptr< interpolators::OneDimensionalInterpolator
                                            < double, Eigen::Vector7d > >( ), bodies.getFrameOrientation( ), bodyFixedFrameNameToUse ) );
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

//! Function to create a simplified system of bodies
simulation_setup::SystemOfBodies createSimplifiedSystemOfBodies(const double secondsSinceJ2000)
{
    using namespace ephemerides;
    using namespace gravitation;

    // Creation of bodies
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

    // Ephemerides
    bodies.getBody( "Sun" )->setEphemeris( std::make_shared< ConstantEphemeris >( Eigen::Vector6d::Zero( )) );
    bodies.getBody( "Mercury" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Mercury") );
    bodies.getBody( "Venus" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Venus" ) );
    bodies.getBody( "Earth" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Earth" ) );
    bodies.getBody( "Mars" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Mars" ) );
    bodies.getBody( "Jupiter" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Jupiter" ) );
    bodies.getBody( "Saturn" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Saturn" ) );
    bodies.getBody( "Uranus" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Uranus" ) );
    bodies.getBody( "Neptune" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Neptune" ) );
    bodies.getBody( "Pluto" )->setEphemeris( std::make_shared< ApproximateGtopEphemeris >("Pluto" ) );

    // Gravity field
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

    // Earth's shape model
    bodies.getBody( "Earth" )->setShapeModel( std::make_shared< SphericalBodyShapeModel >(celestial_body_constants::EARTH_EQUATORIAL_RADIUS ) );

    // Calculate position of rotation axis at initialTime, with respect to J2000 frame. Values from:
    // "Report of the IAU Working Group on Cartographic Coordinates and Rotational Elements: 2009", B.A. Archinal et al.(2011)
    const double daysSinceJ2000 = secondsSinceJ2000 / 86400;
    const double centuriesSinceJ2000 = daysSinceJ2000 / 36525;
    const double poleRightAscension = (0 - 0.641 * centuriesSinceJ2000) * mathematical_constants::PI/180;
    const double poleDeclination = (90 - 0.557 * centuriesSinceJ2000) * mathematical_constants::PI/180;
    const double initialRotationAngle = 190.147 * mathematical_constants::PI/180;
    const double rotationRate = 360.985235 * mathematical_constants::PI/180 / physical_constants::JULIAN_DAY;

    Eigen::Matrix3d J2000toPlanetocentricMatrix = reference_frames::getInertialToPlanetocentricFrameTransformationMatrix (
            poleDeclination, poleRightAscension, initialRotationAngle);
    Eigen::Matrix3d ECLIPJ2000toJ2000Matrix = reference_frames::getECLIPJ2000toJ2000TransformationMatrix();
    Eigen::Matrix3d ECLIPJ2000toPlanetocentricMatrix =  J2000toPlanetocentricMatrix * ECLIPJ2000toJ2000Matrix;

    bodies.getBody( "Earth" )->setRotationalEphemeris(std::make_shared< SimpleRotationalEphemeris >(
                    Eigen::Quaterniond(ECLIPJ2000toPlanetocentricMatrix),
                    rotationRate, secondsSinceJ2000, "ECLIPJ2000", "Earth_Fixed" ) );

    return bodies;
}

} // namespace simulation_setup

} // namespace tudat
