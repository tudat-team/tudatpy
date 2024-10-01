/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References: 820-013, TRK-2-18 Tracking System Interfaces Orbit Data File Interface, Revision E, 2008, JPL/DSN
 *
 */

#ifndef TUDAT_READ_VARIOUS_PDS_FILES_H
#define TUDAT_READ_VARIOUS_PDS_FILES_H

#include <iostream>
#include <fstream>
#include <bitset>
#include <cmath>
#include <vector>
#include <map>

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

namespace tudat
{
namespace input_output
{


std::pair< std::vector< double >, std::vector< double > > grailAntennaFileReader( const std::string antennaFile );

std::map< double, double > grailMassLevel0FileReader( const std::string massFile );
std::map< double, double > grailMassLevel1FileReader( const std::string massFile, const std::string dataLevel = "1b" );


} // namespace input_output

} // namespace tudat

#endif // TUDAT_READ_VARIOUS_PDS_FILES_H
