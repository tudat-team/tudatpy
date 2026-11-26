/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GETPROCESSEDSTATESTRINGS_H
#define TUDAT_GETPROCESSEDSTATESTRINGS_H

#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "tudat/astro/propagators/nBodyStateDerivative.h"
#include "tudat/astro/propagators/rotationalMotionStateDerivative.h"
#include "tudat/astro/propagators/singleStateTypeDerivative.h"
#include "tudat/simulation/propagation_setup/propagatorType.h"

namespace tudat
{

namespace propagators
{

inline std::map< std::pair< int, int >, std::string > getProcessedStateStrings(
        const std::map< IntegratedStateType, std::vector< std::tuple< std::string, std::string, PropagatorType > > >
                integratedTypeAndBodyList )
{
    unsigned int stateVectorIndex = 0;
    std::map< std::pair< int, int >, std::string > stringPerIndex;

    for( auto integratedTypeAndBody : integratedTypeAndBodyList )
    {
        // Extract state type and list of body names
        IntegratedStateType stateType = integratedTypeAndBody.first;
        std::vector< std::tuple< std::string, std::string, PropagatorType > > bodyList = integratedTypeAndBody.second;

        int stateSize = -1;

        // Loop trough list of body names
        for( unsigned int i = 0; i < bodyList.size( ); i++ )
        {
            std::string currentString = getIntegratedStateTypString( stateType );
            switch( stateType )
            {
                case translational_state:
                    stateSize = getSingleIntegrationSize( stateType );
                    currentString += ", in Cartesian elements, of body " + std::get< 0 >( bodyList.at( i ) ) + " w.r.t. " +
                            std::get< 1 >( bodyList.at( i ) );
                    break;
                case rotational_state:
                    stateSize = getSingleIntegrationSize( stateType );
                    currentString += ", in quaternions and body-fixed angular velocity, of body " + std::get< 0 >( bodyList.at( i ) );
                    break;
                case body_mass_state:
                    stateSize = getSingleIntegrationSize( stateType );
                    currentString += " of body " + std::get< 0 >( bodyList.at( i ) );
                    break;
                case proper_time:
                    stateSize = getSingleIntegrationSize( stateType );
                    currentString += " of body " + std::get< 0 >( bodyList.at( i ) );
                    break;
                case custom_state: {
                    stateSize = std::get< 2 >( bodyList.at( i ) ).customStateSize_;
                    if( stateSize <= 0 )
                    {
                        throw std::runtime_error( "Error when getting custom state size; size is <= 0" );
                    }
                    if( std::get< 0 >( bodyList.at( i ) ) != "" )
                    {
                        currentString += " of body " + std::get< 0 >( bodyList.at( i ) );
                    }
                    break;
                }
                default:
                    throw std::runtime_error( "Error when getting processed state strings, type not recognized" );
            }
            stringPerIndex[ std::make_pair( stateVectorIndex, stateSize ) ] = currentString;
            // Remember where we are at trough the state vector
            stateVectorIndex += stateSize;
        }
    }
    return stringPerIndex;
}

inline std::map< std::pair< int, int >, std::string > getPropagatedStateStrings(
        const std::map< IntegratedStateType, std::vector< std::tuple< std::string, std::string, PropagatorType > > >
                integratedTypeAndBodyList )
{
    unsigned int stateVectorIndex = 0;
    std::map< std::pair< int, int >, std::string > stringPerIndex;

    for( auto integratedTypeAndBody : integratedTypeAndBodyList )
    {
        // Extract state type and list of body names
        IntegratedStateType stateType = integratedTypeAndBody.first;
        std::vector< std::tuple< std::string, std::string, PropagatorType > > bodyList = integratedTypeAndBody.second;

        int stateSize = -1;

        // Loop trough list of body names
        for( unsigned int i = 0; i < bodyList.size( ); i++ )
        {
            std::string currentString = getIntegratedStateTypString( stateType );
            switch( stateType )
            {
                case translational_state: {
                    TranslationalPropagatorType propagatorType = std::get< 2 >( bodyList.at( i ) ).translationalPropagatorType_;
                    stateSize = getTranslationalStateSize( propagatorType );
                    currentString += " (" + getTranslationalPropagatorName( propagatorType ) + " formulation)";
                    currentString += " of body " + std::get< 0 >( bodyList.at( i ) ) + " w.r.t. " + std::get< 1 >( bodyList.at( i ) );
                    break;
                }
                case rotational_state: {
                    RotationalPropagatorType propagatorType = std::get< 2 >( bodyList.at( i ) ).rotationalPropagatorType_;
                    stateSize = getRotationalStateSize( propagatorType );
                    currentString += " (" + getRotationalPropagatorName( propagatorType ) + " formulation)";
                    currentString += " of body " + std::get< 0 >( bodyList.at( i ) );
                    break;
                }
                case body_mass_state:
                    stateSize = getSingleIntegrationSize( stateType );
                    currentString += " of body " + std::get< 0 >( bodyList.at( i ) );
                    break;
                case proper_time:
                    stateSize = getSingleIntegrationSize( stateType );
                    currentString += " of body " + std::get< 0 >( bodyList.at( i ) );
                    break;
                case custom_state: {
                    stateSize = std::get< 2 >( bodyList.at( i ) ).customStateSize_;
                    if( stateSize <= 0 )
                    {
                        throw std::runtime_error( "Error when getting custom state size; size is <= 0" );
                    }
                    if( std::get< 0 >( bodyList.at( i ) ) != "" )
                    {
                        currentString += " of body " + std::get< 0 >( bodyList.at( i ) );
                    }
                    break;
                }
                default:
                    throw std::runtime_error( "Error when getting processed state strings, type not recognized" );
            }
            stringPerIndex[ std::make_pair( stateVectorIndex, stateSize ) ] = currentString;
            // Remember where we are at trough the state vector
            stateVectorIndex += stateSize;
        }
    }
    return stringPerIndex;
}

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_GETPROCESSEDSTATESTRINGS_H
