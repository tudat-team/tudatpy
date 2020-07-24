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

#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/basic/mathematicalConstants.h"

#include "tudat/simulation/environment/createBodies.h"

namespace tudat
{

namespace simulation_setup
{

using namespace ephemerides;
using namespace gravitation;

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
NamedBodyMap createBodies(
        const BodyListSettings& bodySettings )
{
    std::vector< std::pair< std::string, std::shared_ptr< BodySettings > > > orderedBodySettings
            = determineBodyCreationOrder( bodySettings.get( ) );

    // Declare map of bodies that is to be returned.
    NamedBodyMap bodyList = NamedBodyMap(
                bodySettings.getFrameOrigin( ), bodySettings.getFrameOrientation( ) );

    // Create empty body objects.
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )
    {
        bodyList.addNewBody( orderedBodySettings.at( i ).first, false );
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

} // namespace simulation_setup

} // namespace tudat
