/*    Copyright (c) 2010-2017, Delft University of Technology
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
    transmitter = 0,
    reflector1 = 1,
    reflector2 = 2,
    reflector3 = 3,
    reflector4 = 4,
    receiver = 5,
    observed_body = 6
};

//! Typedef for the identifier of a given link-end (body and reference points)
typedef std::pair< std::string, std::string > LinkEndId;

//! Typedef for list of link ends, with associated role, used for a single observation (model).
typedef std::map< LinkEndType, LinkEndId > LinkEnds;

//! Function to get a string identifier for a link end type
/*!
 * Function to get a string identifier for a link end type
 * \param linkEndType Enum identifier for a link end type
 * \return String identifier for a link end type
 */
std::string getLinkEndTypeString( const LinkEndType linkEndType );

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


std::vector< int > getNWayLinkEndIndicesFromLinkEndId( const LinkEndId& linkEndid, const LinkEnds& linkEnds );

std::vector< int > getNWayLinkEndIndicesFromLinkEndId( const std::vector< LinkEndType >& linkEndTypes, const LinkEnds& linkEnds );

std::vector< LinkEndType > getNWayLinkIndicesFromLinkEndId( const LinkEndId& linkEndid, const LinkEnds& linkEnds );

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_LINKTYPEDEFS_H
