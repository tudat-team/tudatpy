/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/assign/list_of.hpp>

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"

namespace tudat
{

namespace simulation_setup
{

using namespace ephemerides;
using namespace gravitation;

std::vector< std::pair< std::string, boost::shared_ptr< BodySettings > > > determineBodyCreationOrder(
        const std::map< std::string, boost::shared_ptr< BodySettings > >& bodySettings )
{
    std::vector< std::pair< std::string, boost::shared_ptr< BodySettings > > > outputVector;

    std::map< std::string, std::vector< boost::shared_ptr< GravityFieldVariationSettings > > >
            gravityFieldVariationSettingsList;

    for( std::map< std::string, boost::shared_ptr< BodySettings > >::const_iterator bodyIterator
         = bodySettings.begin( );
         bodyIterator != bodySettings.end( ); bodyIterator ++ )
    {
        if( bodyIterator->second->gravityFieldVariationSettings.size( ) != 0 )
        {
            gravityFieldVariationSettingsList[ bodyIterator->first ] =
                    bodyIterator->second->gravityFieldVariationSettings;
        }
        outputVector.push_back( std::make_pair( bodyIterator->first, bodyIterator->second ) );
    }

    if( gravityFieldVariationSettingsList.size( ) != 0 )
    {
        for( std::map< std::string,
             std::vector< boost::shared_ptr< GravityFieldVariationSettings > > >::iterator
             variationIterator = gravityFieldVariationSettingsList.begin( );
             variationIterator != gravityFieldVariationSettingsList.end( ); variationIterator++ )
        {
            unsigned int currentBodyIndex = -1;

            for( unsigned int k = 0; k < variationIterator->second.size( ); k++ )
            {
                boost::shared_ptr< BasicSolidBodyGravityFieldVariationSettings > variationSettings1
                        = boost::dynamic_pointer_cast<
                        BasicSolidBodyGravityFieldVariationSettings >( variationIterator->second[ k ] );

                if( variationSettings1 != NULL )
                {

                    for( unsigned int i = 0; i < outputVector.size( ); i++ )
                    {
                        if( outputVector[ i ].first == variationIterator->first )
                        {
                            currentBodyIndex = i;
                        }
                    }

                    std::vector< std::string > deformingBodies
                            = variationSettings1->getDeformingBodies( );
                    for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
                    {
                        for( unsigned int j = 0; j < outputVector.size( ); j++ )
                        {
                            if( deformingBodies[ i ] == outputVector[ j ].first
                                    && j > currentBodyIndex )
                            {
                                std::pair< std::string, boost::shared_ptr< BodySettings > >
                                        entryToMove = outputVector[ j ];
                                outputVector.erase( outputVector.begin( ) + j );
                                outputVector.insert( outputVector.begin( ) + i, entryToMove );
                            }
                        }
                    }
                }
            }
        }
    }
    return outputVector;
}

//! Function to create a single Body object from BodySettings and add it to a NamedBodyMap.
void addBody( NamedBodyMap& bodyMap, const std::string& bodyName,
              const boost::shared_ptr< BodySettings >& bodySettings )
{
    // Create empty body object.
    boost::shared_ptr< Body > body = boost::make_shared< Body >( );
    bodyMap[ bodyName ] = body;

    // Create ephemeris object (if required).
    if( bodySettings->ephemerisSettings != NULL )
    {
        body->setEphemeris( createBodyEphemeris( bodySettings->ephemerisSettings, bodyName ) );
    }

    // Create atmosphere model object (if required).
    if( bodySettings->atmosphereSettings != NULL )
    {
        body->setAtmosphereModel( createAtmosphereModel( bodySettings->atmosphereSettings, bodyName ) );
    }

    // Create body shape model object (if required).
    if( bodySettings->shapeModelSettings != NULL )
    {
        body->setShapeModel( createBodyShapeModel( bodySettings->shapeModelSettings, bodyName ) );
    }

    // Create rotation model object (if required).
    if( bodySettings->rotationModelSettings != NULL )
    {
        body->setRotationalEphemeris( createRotationModel( bodySettings->rotationModelSettings, bodyName ) );
    }

    // Create gravity field model object (if required).
    if( bodySettings->gravityFieldSettings != NULL )
    {
        body->setGravityFieldModel( createGravityFieldModel( bodySettings->gravityFieldSettings,  bodyName, bodyMap,
                                                             bodySettings->gravityFieldVariationSettings ) );
    }

    if( bodySettings->gravityFieldVariationSettings.size( ) > 0 )
    {
        body->setGravityFieldVariationSet( createGravityFieldModelVariationsSet(
                                               bodyName, bodyMap, bodySettings->gravityFieldVariationSettings ) );
    }

    // Create aerodynamic coefficient interface object (if required).
    if( bodySettings->aerodynamicCoefficientSettings != NULL )
    {
        body->setAerodynamicCoefficientInterface( createAerodynamicCoefficientInterface(
                                                      bodySettings->aerodynamicCoefficientSettings, bodyName ) );
    }

    // Create radiation pressure coefficient object (if required).
    std::map< std::string, boost::shared_ptr< RadiationPressureInterfaceSettings > >
            radiationPressureSettings
            = bodySettings->radiationPressureSettings;
    for( std::map< std::string, boost::shared_ptr< RadiationPressureInterfaceSettings > >::iterator
         radiationPressureSettingsIterator = radiationPressureSettings.begin( );
         radiationPressureSettingsIterator != radiationPressureSettings.end( );
         radiationPressureSettingsIterator++ )
    {
        body->setRadiationPressureInterface(
                    radiationPressureSettingsIterator->first, createRadiationPressureInterface(
                        radiationPressureSettingsIterator->second, bodyName, bodyMap ) );
    }

    // Constant mass (if required).
    const double constantMass = bodySettings->constantMass;
    if ( constantMass == constantMass )
    {
        body->setConstantBodyMass( constantMass );
    }
}

//! Function to create a map of bodies objects.
NamedBodyMap createBodies(
        const std::map< std::string, boost::shared_ptr< BodySettings > >& bodySettings )
{
    std::vector< std::pair< std::string, boost::shared_ptr< BodySettings > > > orderedBodySettings
            = determineBodyCreationOrder( bodySettings );

    // Declare map of bodies that is to be returned.
    NamedBodyMap bodyMap;

    // Create bodies
    for( std::pair< std::string, boost::shared_ptr< BodySettings > > entry : orderedBodySettings )
    {
        addBody( bodyMap, entry.first, entry.second );
    }

    return bodyMap;
}

} // namespace simulation_setup

} // namespace tudat
