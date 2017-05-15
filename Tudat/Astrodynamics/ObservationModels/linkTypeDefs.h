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
    reflector = 1,
    receiver = 2,
    observed_body = 3
};

//! Typedef for the identifier of a given link-end (body and reference points)
typedef std::pair< std::string, std::string > LinkEndId;

//! Typedef for list of link ends, with associated role, used for a single observation (model).
typedef std::map< LinkEndType, LinkEndId > LinkEnds;

std::string getLinkEndTypeString( const LinkEndType linkEndType );

std::string getLinkEndsString( const LinkEnds linkEnds );

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_LINKTYPEDEFS_H
