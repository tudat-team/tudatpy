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
    case reflector:
        linkEndString = "reflector";
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


} // namespace observation_models

} // namespace tudat
