/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DSNSTATIONDATA_H
#define TUDAT_DSNSTATIONDATA_H

#include <map>
#include <string>
#include <vector>

namespace tudat
{

namespace input_output
{

/*!
 * Returns the default DSN station names per DSN station complex id. Stations are named as "DSS-i".
 */
std::map< int, std::vector< std::string > > getDefaultDsnStationNamesPerComplex( );

}  // namespace input_output

}  // namespace tudat

#endif  // TUDAT_DSNSTATIONDATA_H
