/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"

namespace tudat
{

namespace observation_models
{

//! Function to get a string identifier for a link end type
std::string getLinkEndTypeString( const LinkEndType linkEndType )
{
    std::string linkEndString = "";
    switch( linkEndType )
    {
    case transmitter:
        linkEndString = "transmitter";
        break;
    case reflector1:
        linkEndString = "reflector_1";
        break;
    case reflector2:
        linkEndString = "reflector_2";
        break;
    case reflector3:
        linkEndString = "reflector_3";
        break;
    case reflector4:
        linkEndString = "reflector_4";
        break;
    case receiver:
        linkEndString = "receiver";
        break;
    case observed_body:
        linkEndString = "observed body";
        break;
    default:
        std::string errorMessage = "Error when getting link end string for type " +
                std::to_string( linkEndType ) + ", type not found.";
        throw std::runtime_error( errorMessage );
    }
    return linkEndString;
}

//! Function to get a string identifier for a set of link ends
std::string getLinkEndsString( const LinkEnds linkEnds )
{
    std::string linkEndsString = "";
    for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( );
         linkEndIterator++ )
    {
        linkEndsString += getLinkEndTypeString( linkEndIterator->first ) + ": (" + linkEndIterator->second.first;
        if( linkEndIterator->second.second == "" )
        {
            linkEndsString +=  ")";
        }
        else
        {
            linkEndsString += ", " + linkEndIterator->second.second +  ")";
        }

        if( linkEndIterator != (--linkEnds.end( ) ) )
        {
            linkEndsString += "; ";
        }
    }
    return linkEndsString;
}

//! Function to get the link end index (0=transmitter, numberOfLinkEnds-1=receiver) of a link end in n-way observable
int getNWayLinkIndexFromLinkEndType( const LinkEndType linkEndType, const int numberOfLinkEnds )
{
    int linkEndIndex;
    // If index is first or last, set 0 or numberOfLinkEnds - 1, respectively
    if( linkEndType == transmitter )
    {
        linkEndIndex = 0;
    }
    else if( linkEndType == receiver )
    {
        linkEndIndex = numberOfLinkEnds - 1;
    }
    else
    {
        linkEndIndex = static_cast< int >( linkEndType ) - static_cast< int >( reflector1 ) + 1;
        if( linkEndIndex > numberOfLinkEnds - 2 )
        {
            throw std::runtime_error( "Error when getting n-way link end index; value too large." );
        }
    }
    return linkEndIndex;
}

//! Function to get the link end type enum of a link end in n-way observable from link index
LinkEndType getNWayLinkEnumFromIndex( const int linkEndIndex, const int numberOfLinkEnds )
{
    LinkEndType linkEndType;

    // If index is first or last, set transmitter or receiver, respectively
    if( linkEndIndex == 0 )
    {
        linkEndType = transmitter;
    }
    else if( linkEndIndex == numberOfLinkEnds - 1 )
    {
        linkEndType = receiver;
    }
    // Check feasibility of inner link end
    else if( linkEndIndex >= numberOfLinkEnds )
    {
        throw std::runtime_error(
                    "Error, found link end index " + std::to_string( linkEndIndex ) +
                    " when getting n-way link end index for " + std::to_string(
                   numberOfLinkEnds ) + " link end total." );
    }
    else
    {
        linkEndType = static_cast< LinkEndType >( static_cast< int >( reflector1 ) + ( linkEndIndex - 1 ) );
    }

    return linkEndType;
}

//! Function to get the list of indices in link-end list for n-way observables that matches a given link end id.
std::vector< int > getNWayLinkEndIndicesFromLinkEndId( const LinkEndId& linkEndid, const LinkEnds& linkEnds )
{
    std::vector< LinkEndType >  matchingLinkEndTypes = getNWayLinkIndicesFromLinkEndId(
                linkEndid, linkEnds );
    return getNWayLinkEndIndicesFromLinkEndId( matchingLinkEndTypes, linkEnds );
}

//! Function to get the list of indices in link-end list for n-way observables that matches a list of link-end types.
std::vector< int > getNWayLinkEndIndicesFromLinkEndId( const std::vector< LinkEndType >& linkEndTypes, const LinkEnds& linkEnds )
{
    std::vector< int > linkEndIndices;
    for( unsigned int i = 0; i < linkEndTypes.size( ); i++ )
    {
        linkEndIndices.push_back( getNWayLinkIndexFromLinkEndType( linkEndTypes.at( i ), linkEnds.size( ) ) );
    }
    return linkEndIndices;
}

//! Function to get the list of link end types in link-end list for n-way observables that match a given link end id.
std::vector< LinkEndType > getNWayLinkIndicesFromLinkEndId( const LinkEndId& linkEndid, const LinkEnds& linkEnds )
{
    std::vector< LinkEndType > matchingLinkEndTypes;

    for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( ); linkEndIterator++ )
    {
        if( linkEndIterator->second == linkEndid || ( ( linkEndIterator->second.first == linkEndid.first ) &&
                                                                        linkEndid.second == "" ) )
        {
            matchingLinkEndTypes.push_back( linkEndIterator->first );
        }
    }
    return matchingLinkEndTypes;
}

} // namespace observation_models

} // namespace tudat
