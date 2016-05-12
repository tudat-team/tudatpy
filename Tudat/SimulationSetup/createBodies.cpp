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

#include "Tudat/SimulationSetup/createBodies.h"

namespace tudat
{

namespace simulation_setup
{

using namespace ephemerides;
using namespace gravitation;

//! Function to create a map of bodies objects.
NamedBodyMap createBodies(
        const std::map< std::string, boost::shared_ptr< BodySettings > >& bodySettings )
{
    // Declare map of bodies that is to be returned.
    NamedBodyMap bodyMap;

    // Create empty body objects.
    for( std::map< std::string, boost::shared_ptr< BodySettings > >::const_iterator settingIterator
         = bodySettings.begin( ); settingIterator != bodySettings.end( ); settingIterator++ )
    {
        bodyMap[ settingIterator->first ] = boost::make_shared< Body >( );
    }

    // Create ephemeris objects for each body (if required).
    for( std::map< std::string, boost::shared_ptr< BodySettings > >::const_iterator settingIterator
         = bodySettings.begin( ); settingIterator != bodySettings.end( ); settingIterator++ )
    {
        if( settingIterator->second->ephemerisSettings != NULL )
        {
            bodyMap[ settingIterator->first ]->setEphemeris(
                        createBodyEphemeris( settingIterator->second->ephemerisSettings,
                                                 settingIterator->first ) );
        }
        else
        {
            std::cerr<<"Warning, no ephemeris data found for body "<<
                       settingIterator->first<<std::endl;
        }
    }

    // Create atmosphere model objects for each body (if required).
    for( std::map< std::string, boost::shared_ptr< BodySettings > >::const_iterator settingIterator
         = bodySettings.begin( ); settingIterator != bodySettings.end( ); settingIterator++ )
    {
        if( settingIterator->second->atmosphereSettings != NULL )
        {
            bodyMap[ settingIterator->first ]->setAtmosphereModel(
                        createAtmosphereModel( settingIterator->second->atmosphereSettings,
                                               settingIterator->first ) );
        }
    }

    // Create rotation model objects for each body (if required).
    for( std::map< std::string, boost::shared_ptr< BodySettings > >::const_iterator settingIterator
         = bodySettings.begin( ); settingIterator != bodySettings.end( ); settingIterator++ )
    {
        if( settingIterator->second->rotationModelSettings != NULL )
        {
            bodyMap[ settingIterator->first ]->setRotationalEphemeris(
                        createRotationModel( settingIterator->second->rotationModelSettings,
                                             settingIterator->first ) );
        }
    }

    // Create gravity field model objects for each body (if required).
    for( std::map< std::string, boost::shared_ptr< BodySettings > >::const_iterator settingIterator
         = bodySettings.begin( ); settingIterator != bodySettings.end( ); settingIterator++ )
    {
        if( settingIterator->second->gravityFieldSettings != NULL )
        {
            bodyMap[ settingIterator->first ]->setGravityFieldModel(
                        createGravityFieldModel( settingIterator->second->gravityFieldSettings,
                                                 settingIterator->first ) );
        }
    }

    return bodyMap;

}

} // namespace simulation_setup

} // namespace tudat
