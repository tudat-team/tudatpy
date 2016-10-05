/*    Copyright (c) 2010-2016, Delft University of Technology
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


//! Function to create a map of bodies objects.
NamedBodyMap createBodies(
        const std::map< std::string, boost::shared_ptr< BodySettings > >& bodySettings )
{
    std::vector< std::pair< std::string, boost::shared_ptr< BodySettings > > > orderedBodySettings
           = determineBodyCreationOrder( bodySettings );

    // Declare map of bodies that is to be returned.
    NamedBodyMap bodyMap;

    // Create empty body objects.
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        bodyMap[ orderedBodySettings.at( i ).first ] = boost::make_shared< Body >( );
    }

    // Create ephemeris objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        if( orderedBodySettings.at( i ).second->ephemerisSettings != NULL )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setEphemeris(
                        createBodyEphemeris( orderedBodySettings.at( i ).second->ephemerisSettings,
                                             orderedBodySettings.at( i ).first ) );
        }
    }

    // Create atmosphere model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        if( orderedBodySettings.at( i ).second->atmosphereSettings != NULL )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setAtmosphereModel(
                        createAtmosphereModel( orderedBodySettings.at( i ).second->atmosphereSettings,
                                               orderedBodySettings.at( i ).first ) );
        }
    }

    // Create body shape model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        if( orderedBodySettings.at( i ).second->shapeModelSettings != NULL )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setShapeModel(
                        createBodyShapeModel( orderedBodySettings.at( i ).second->shapeModelSettings,
                                              orderedBodySettings.at( i ).first ) );
        }
    }

    // Create rotation model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        if( orderedBodySettings.at( i ).second->rotationModelSettings != NULL )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setRotationalEphemeris(
                        createRotationModel( orderedBodySettings.at( i ).second->rotationModelSettings,
                                             orderedBodySettings.at( i ).first ) );
        }
    }

    // Create gravity field model objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        if( orderedBodySettings.at( i ).second->gravityFieldSettings != NULL )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setGravityFieldModel(
                        createGravityFieldModel( orderedBodySettings.at( i ).second->gravityFieldSettings,
                                                 orderedBodySettings.at( i ).first, bodyMap,
                                                 orderedBodySettings.at( i ).second->gravityFieldVariationSettings ) );
        }
    }

    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        if( orderedBodySettings.at( i ).second->gravityFieldVariationSettings.size( ) > 0 )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setGravityFieldVariationSet(
                        createGravityFieldModelVariationsSet(
                            orderedBodySettings.at( i ).first, bodyMap,
                            orderedBodySettings.at( i ).second->gravityFieldVariationSettings ) );
        }
    }

    // Create aerodynamic coefficient interface objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        if( orderedBodySettings.at( i ).second->aerodynamicCoefficientSettings != NULL )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setAerodynamicCoefficientInterface(
                        createAerodynamicCoefficientInterface(
                            orderedBodySettings.at( i ).second->aerodynamicCoefficientSettings,
                            orderedBodySettings.at( i ).first ) );
        }
    }


    // Create radiation pressure coefficient objects for each body (if required).
    for( unsigned int i = 0; i < orderedBodySettings.size( ); i++ )         
    {
        std::map< std::string, boost::shared_ptr< RadiationPressureInterfaceSettings > >
                radiationPressureSettings
                = orderedBodySettings.at( i ).second->radiationPressureSettings;
        for( std::map< std::string, boost::shared_ptr< RadiationPressureInterfaceSettings > >::iterator
             radiationPressureSettingsIterator = radiationPressureSettings.begin( );
             radiationPressureSettingsIterator != radiationPressureSettings.end( );
             radiationPressureSettingsIterator++ )
        {
            bodyMap[ orderedBodySettings.at( i ).first ]->setRadiationPressureInterface(
                        radiationPressureSettingsIterator->first,
                        createRadiationPressureInterface(
                            radiationPressureSettingsIterator->second,
                            orderedBodySettings.at( i ).first, bodyMap ) );
        }

    }
    return bodyMap;

}

} // namespace simulation_setup

} // namespace tudat
