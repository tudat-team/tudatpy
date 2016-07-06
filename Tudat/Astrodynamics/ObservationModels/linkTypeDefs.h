/*    Copyright (c) 2010-2016, Delft University of Technology
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
 *  Enum defining different link end types. Enum is used for book-keeping purposes when setting link characteristics
 *  so that one does not need to remember the numerical order of types of link ends in input vector.
 */
enum LinkEndType
{
    transmitter = 0,
    reflector = 1,
    receiver = 2,
    observed_body = 3
};

typedef std::pair< std::string, std::string > LinkEndId;

typedef std::map< LinkEndType, LinkEndId > LinkEnds;

}

}

#endif // TUDAT_LINKTYPEDEFS_H
