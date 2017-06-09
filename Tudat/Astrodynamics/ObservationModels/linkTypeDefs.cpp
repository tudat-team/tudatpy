/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/lexical_cast.hpp>
#include <iostream>

#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"

namespace tudat
{

namespace observation_models
{

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
                boost::lexical_cast< std::string >( linkEndType ) + ", type not found.";
        throw std::runtime_error( errorMessage );
    }
    return linkEndString;
}

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

int getNWayLinkIndexFromLinkEndType( const LinkEndType linkEndType, const int numberOfLinkEnds )
{
    int linkEndIndex;
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

LinkEndType getNWayLinkEnumFromIndex( const int linkEndIndex, const int numberOfLinkEnds )
{
    LinkEndType linkEndType;
    if( linkEndIndex == 0 )
    {
        linkEndType = transmitter;
    }
    else if( linkEndIndex == numberOfLinkEnds - 1 )
    {
        linkEndType = receiver;
    }
    else if( linkEndIndex >= numberOfLinkEnds )
    {
        std::cerr<<"Error, found link end index "<<linkEndIndex<<" when getting n-way link end index for "<<
                   numberOfLinkEnds<<" link end total."<<std::endl;
    }
    else
    {
        linkEndType = static_cast< LinkEndType >( static_cast< int >( reflector1 ) + ( linkEndIndex - 1 ) );
    }

    return linkEndType;
}


} // namespace observation_models

} // namespace tudat
