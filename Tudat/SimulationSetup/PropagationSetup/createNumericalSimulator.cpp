/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include "Tudat/SimulationSetup/PropagationSetup/createNumericalSimulator.h"
#include "Tudat/Astrodynamics/Ephemerides/frameManager.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a set of acceleration models from a map of bodies and acceleration model types.
basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
        const NamedBodyMap& bodyMap,
        const SelectedAccelerationMap& selectedAccelerationPerBody,
        const std::map< std::string, std::string >& centralBodies )
{
    // Declare return map.
    basic_astrodynamics::AccelerationMap accelerationModelMap;

    // Put selectedAccelerationPerBody in correct order
    SelectedAccelerationMap orderedAccelerationPerBody =
            orderSelectedAccelerationMap( selectedAccelerationPerBody );

    // Iterate over all bodies which are undergoing acceleration
    for( SelectedAccelerationMap::const_iterator bodyIterator =
         orderedAccelerationPerBody.begin( ); bodyIterator != orderedAccelerationPerBody.end( );
         bodyIterator++ )
    {
        boost::shared_ptr< Body > currentCentralBody;

        // Retrieve name of body undergoing acceleration.
        std::string bodyUndergoingAcceleration = bodyIterator->first;

        // Retrieve name of current central body.
        std::string currentCentralBodyName = centralBodies.at( bodyUndergoingAcceleration );

        if( !ephemerides::isFrameInertial( currentCentralBodyName ) )
        {
            if( bodyMap.count( currentCentralBodyName ) == 0 )
            {
                throw std::runtime_error(
                            std::string( "Error, could not find non-inertial central body ") +
                            currentCentralBodyName + " of " + bodyUndergoingAcceleration +
                            " when making acceleration model." );
            }
            else
            {
                currentCentralBody = bodyMap.at( currentCentralBodyName );
            }
        }

        // Check if body undergoing acceleration is included in bodyMap
        if( bodyMap.count( bodyUndergoingAcceleration ) ==  0 )
        {
            throw std::runtime_error(
                        std::string( "Error when making acceleration models, requested forces" ) +
                        "acting on body " + bodyUndergoingAcceleration  +
                        ", but no such body found in map of bodies" );
        }

        // Declare map of acceleration models acting on current body.
        basic_astrodynamics::SingleBodyAccelerationMap mapOfAccelerationsForBody;

        // Retrieve list of required acceleration model types and bodies exerting accelerationd on
        // current body.
        std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > >
                accelerationsForBody = bodyIterator->second;

        // Iterate over all bodies exerting an acceleration
        for( std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > >::
             iterator body2Iterator = accelerationsForBody.begin( );
             body2Iterator != accelerationsForBody.end( ); body2Iterator++ )
        {
            // Retrieve name of body exerting acceleration.
            std::string bodyExertingAcceleration = body2Iterator->first;

            // Check if body exerting acceleration is included in bodyMap
            if( bodyMap.count( bodyExertingAcceleration ) ==  0 )
            {
                throw std::runtime_error(
                            std::string( "Error when making acceleration models, requested forces ")
                            + "acting on body " + bodyUndergoingAcceleration  + " due to body " +
                            bodyExertingAcceleration +
                            ", but no such body found in map of bodies" );
            }

            // Retrieve list of accelerations due to current body.
            std::vector< boost::shared_ptr< AccelerationSettings > > accelerationList =
                    body2Iterator->second;

            for( unsigned int i = 0; i < accelerationList.size( ); i++ )
            {
                // Create acceleration model.
                mapOfAccelerationsForBody[ bodyExertingAcceleration ].push_back(
                            createAccelerationModel( bodyMap.at( bodyUndergoingAcceleration ),
                                                     bodyMap.at( bodyExertingAcceleration ),
                                                     accelerationList.at( i ),
                                                     bodyUndergoingAcceleration,
                                                     bodyExertingAcceleration,
                                                     currentCentralBody,
                                                     currentCentralBodyName,
                                                     bodyMap ) );
            }
        }

        // Put acceleration models on current body in return map.
        accelerationModelMap[ bodyUndergoingAcceleration ] = mapOfAccelerationsForBody;
    }

    return accelerationModelMap;
}

//! Function to create acceleration models from a map of bodies and acceleration model types.
basic_astrodynamics::AccelerationMap createAccelerationModelsMap(
        const NamedBodyMap& bodyMap,
        const SelectedAccelerationMap& selectedAccelerationPerBody,
        const std::vector< std::string >& propagatedBodies,
        const std::vector< std::string >& centralBodies )
{
    if( centralBodies.size( ) != propagatedBodies.size( ) )
    {
        throw std::runtime_error( "Error, number of propagated bodies must equal number of central bodies" );
    }

    std::map< std::string, std::string > centralBodyMap;
    for( unsigned int i = 0; i < propagatedBodies.size( ); i++ )
    {
        centralBodyMap[ propagatedBodies.at( i ) ] = centralBodies.at( i );
    }

    return createAccelerationModelsMap( bodyMap, selectedAccelerationPerBody, centralBodyMap );
}

}

}
