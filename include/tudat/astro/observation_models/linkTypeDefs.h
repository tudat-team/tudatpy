/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_LINKTYPEDEFS_H
#define TUDAT_LINKTYPEDEFS_H

#include <stdexcept>
#include <map>
#include <string>
#include <vector>

namespace tudat
{

namespace observation_models
{

//! Enum defining different link end types.
/*!
 *  Enum defining different roles that a given link end can play in an observation model.
 */
enum LinkEndType
{
    unidentified_link_end = -1,
    default_link_end = unidentified_link_end,
    transmitter = 0,
    reflector1 = 1,
    reflector = reflector1,
    retransmitter = reflector,
    retransmitter1 = retransmitter,
    reflector2 = 2,
    retransmitter2 = reflector2,
    reflector3 = 3,
    retransmitter3 = reflector3,
    reflector4 = 4,
    retransmitter4 = reflector4,
    receiver = 5,
    observed_body = 6,
    transmitter2 = 7,
    observer = 8
};

////! Typedef for the identifier of a given link-end (body and reference points)
//typedef std::pair< std::string, std::string > LinkEndId;


struct LinkEndId
{
    LinkEndId( ){ }

    LinkEndId( const std::pair< std::string, std::string >& linkEnd ):
        bodyName_( linkEnd.first ),
        stationName_( linkEnd.second ){ }

    LinkEndId( const std::string& bodyName,
               const std::string& stationName ):
        bodyName_( bodyName ),
        stationName_( stationName ){ }

    LinkEndId( const std::string& bodyName ):
        bodyName_( bodyName ),
        stationName_( "" ){ }

    std::pair< std::string, std::string > getDualStringLinkEnd( ) const
    {
        return std::make_pair( bodyName_, stationName_ );
    }

    friend bool operator==( const LinkEndId& linkEnd1, const LinkEndId& linkEnd2 )
    {
        return ( ( linkEnd1.bodyName_ == linkEnd2.bodyName_ ) && ( linkEnd1.stationName_ == linkEnd2.stationName_ ) );
    }

    friend bool operator!=( const LinkEndId& linkEnd1, const LinkEndId& linkEnd2 )
    {
        return !operator==( linkEnd1, linkEnd2 );
    }

    friend bool operator< ( const LinkEndId& linkEnd1, const LinkEndId& linkEnd2 )
    {
        if( linkEnd1.bodyName_ < linkEnd2.bodyName_ )
        {
            return true;
        }
        else if( linkEnd1.bodyName_ > linkEnd2.bodyName_ )
        {
            return false;
        }
        else if( linkEnd1.stationName_ < linkEnd2.stationName_ )
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    std::string bodyName_;

    std::string stationName_;

    std::string getBodyName( )
    {
        return bodyName_;
    }

    std::string getStationName( )
    {
        return stationName_;
    }

};

inline LinkEndId linkEndId( const std::string& bodyName,
                            const std::string& stationName )
{
    return  LinkEndId( bodyName, stationName );
}

inline LinkEndId linkEndId( const std::string& bodyName )
{
    return LinkEndId( bodyName );
}


//! Function to get a string identifier for a link end type
/*!
 * Function to get a string identifier for a link end type
 * \param linkEndType Enum identifier for a link end type
 * \return String identifier for a link end type
 */
std::string getLinkEndTypeString( const LinkEndType linkEndType );

//! Typedef for list of link ends, with associated role, used for a single observation (model).
typedef std::map< LinkEndType, LinkEndId > LinkEnds;

struct LinkDefinition
{
    LinkDefinition( ){ }

    LinkDefinition( const std::map< LinkEndType, LinkEndId >& linkEnds ):
        linkEnds_( linkEnds ){ }

    LinkDefinition( const std::map< LinkEndType, std::pair< std::string, std::string > >& linkEnds )
    {
        for( auto it : linkEnds )
        {
            linkEnds_[ it.first ] = LinkEndId( it.second );
        }
    }

    LinkEndId at( const LinkEndType linkEndType ) const
    {
        if( linkEnds_.count( linkEndType ) == 0 )
        {
            throw std::runtime_error( "Error when extracing link end, requested link end type " +
                                      getLinkEndTypeString( linkEndType ) +
                                      " not found" );
        }
        return linkEnds_.at( linkEndType );
    }

    LinkEndId &operator[](LinkEndType linkEndType)
    {
        return linkEnds_[ linkEndType ];
    }

    std::map< LinkEndType, LinkEndId > linkEnds_;

    unsigned int size( ) const { return static_cast< unsigned int >( linkEnds_.size( ) ); }

    friend bool operator==( const LinkDefinition& linkEnds1, const LinkDefinition& linkEnds2 )
    {
        bool isEqual = true;
        std::map< LinkEndType, LinkEndId > firstLinkEnds = linkEnds1.linkEnds_;
        std::map< LinkEndType, LinkEndId > secondLinkEnds = linkEnds2.linkEnds_;

        std::map< LinkEndType, LinkEndId >::iterator firstLinkEndIterator = firstLinkEnds.begin( );
        std::map< LinkEndType, LinkEndId >::iterator secondLinkEndIterator = secondLinkEnds.begin( );

        if( firstLinkEnds.size( ) == secondLinkEnds.size( ) )
        {
            for( unsigned int i = 0; i < firstLinkEnds.size( ); i++ )
            {
                if( ( firstLinkEndIterator->first != secondLinkEndIterator->first ) ||
                        ( firstLinkEndIterator->second != secondLinkEndIterator->second ) )
                {
                    isEqual = false;
                    break;
                }

                firstLinkEndIterator++;
                secondLinkEndIterator++;
            }
        }
        else
        {
            isEqual = false;
        }

        return isEqual;
    }

//    friend bool operator< ( const LinkEnds& linkEnds1, const LinkEnds& linkEnds2 )
//    {
//        return linkEnds1.linkEnds_ < linkEnds2.linkEnds_;
//    }

    friend bool operator!=( const LinkEnds& linkEnds1, const LinkEnds& linkEnds2 )
    {
        return !operator==( linkEnds1, linkEnds2 );
    }

};

inline LinkDefinition linkDefinition( const std::map< LinkEndType, LinkEndId >& linkEnds )
{
    return LinkDefinition( linkEnds );
}

//! Function to get a string identifier for a set of link ends
/*!
 * Function to get a string identifier for a set of link ends
 * \param linkEnds Link ends for which string identifier is to be created.
 * \return String identifier for a link ends
 */
std::string getLinkEndsString( const LinkEnds linkEnds );

//! Function to get the link end index (0=transmitter, numberOfLinkEnds-1=receiver) of a link end in n-way observable
/*!
 * Function to get the link end index (0=transmitter, numberOfLinkEnds-1=receiver) of a link end in n-way observable
 * \param linkEndType Link end type for which index is to be retrieved
 * \param numberOfLinkEnds Total number of link ends in n-way observable
 * \return Requested link end index
 */
int getNWayLinkIndexFromLinkEndType( const LinkEndType linkEndType, const int numberOfLinkEnds );

//! Function to get the link end type enum of a link end in n-way observable from link index
/*!
 *  Function to get the link end type enum of a link end in n-way observable from link index
 *  (0=transmitter, numberOfLinkEnds-1=receiver)
 *  \param linkEndIndex Link end index for which type is to be retrieved
 *  \param numberOfLinkEnds Total number of link ends in n-way observable
 *  \return Requested link end type
 */
LinkEndType getNWayLinkEnumFromIndex( const int linkEndIndex, const int numberOfLinkEnds );

//! Function to get the list of indices in link-end list for n-way observables that match a given link end id.
/*!
 * Function to get the list of indices in link-end list for n-way observables that match a given link end id.
 * \param linkEndId Link end id for which the link-end indices in linkEnds are to be determined
 * \param linkEnds N-way link ends
 * \return List of indices (transmitter = 0, receiver = linkEnds.size( ) - 1 ) at which linkEndId occurs in linkEnds.
 */
std::vector< int > getNWayLinkEndIndicesFromLinkEndId( const LinkEndId& linkEndId, const LinkEnds& linkEnds );

//! Function to get the list of indices in link-end list for n-way observables that match a list of link-end types.
/*!
 * Function to get the list of indices in link-end list for n-way observables that match a list of link-end types.
 * \param linkEndTypes List of link end types for which the link-end indices in linkEnds are to be determined
 * \param linkEnds N-way link ends
 * \return List of indices (transmitter = 0, receiver = linkEnds.size( ) - 1 ) at which linkEndTypes occurs in linkEnds.
 */
std::vector< int > getNWayLinkEndIndicesFromLinkEndId( const std::vector< LinkEndType >& linkEndTypes, const LinkEnds& linkEnds );

//! Function to get the list of link end types in link-end list for n-way observables that match a given link end id.
/*!
 * Function to get the list of link end types in link-end list for n-way observables that match a given link end id.
 * \param linkEndId Link end id for which the link-end indices in linkEnds are to be determined
 * \param linkEnds N-way link ends
 * \return List of link end type at which linkEndId occurs in linkEnds.
 */
std::vector< LinkEndType > getNWayLinkIndicesFromLinkEndId( const LinkEndId& linkEndId, const LinkEnds& linkEnds );

LinkEnds mergeUpDownLink( const LinkEnds& uplink, const LinkEnds& downlink );

LinkDefinition mergeUpDownLink( const LinkDefinition& uplink, const LinkDefinition& downlink );

LinkEnds mergeOneWayLinkEnds( const std::vector< LinkEnds >& linkEnds );

LinkDefinition mergeOneWayLinkEnds( const std::vector< LinkDefinition >& linkEnds );

LinkEnds getUplinkFromTwoWayLinkEnds(
        const LinkEnds& twoWayLinkEnds );

LinkDefinition getUplinkFromTwoWayLinkEnds(
        const LinkDefinition& twoWayLinkEnds );

LinkEnds getDownlinkFromTwoWayLinkEnds(
        const LinkEnds& twoWayLinkEnds );

LinkDefinition getDownlinkFromTwoWayLinkEnds(
        const LinkDefinition& twoWayLinkEnds );

LinkEnds getSingleLegLinkEnds(
        const LinkEnds& nWayLinkEnds, const unsigned int legIndex );

LinkEnds mergeOneWayLinkEnds(
        const std::vector< LinkEnds >& linkEnds );

std::vector< LinkEnds > getOneWayDownlinkLinkEndsList(
        const LinkEndId singleTransmitter,
        const std::vector< LinkEndId >& listOfReceivers );

std::vector< LinkEnds > getOneWayUplinkLinkEndsList(
        const std::vector< LinkEndId > listOfTransmitters,
        const LinkEndId singleReceivers );

std::vector< LinkEnds > getSameStationTwoWayLinkEndsList(
        const std::vector< LinkEndId > listOfStations,
        const LinkEndId spacecraft );

std::vector< LinkEnds > getTwoWayLinkEndsList(
        const std::vector< LinkEndId > listOfStations,
        const LinkEndId spacecraft );


bool isLinkEndPresent(
        const LinkEnds linkEnds,
        const LinkEndId linkEndToSearch );



} // namespace observation_models

} // namespace tudat

#endif // TUDAT_LINKTYPEDEFS_H
