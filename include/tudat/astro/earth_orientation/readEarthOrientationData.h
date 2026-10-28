/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_READ_EARTH_ORIENTATION_DATA_H
#define TUDAT_READ_EARTH_ORIENTATION_DATA_H

#include <array>
#include <istream>
#include <limits>
#include <map>
#include <optional>
#include <string>

namespace tudat
{
namespace earth_orientation
{

//! File formats supported by the Earth orientation data reader.
enum class EarthOrientationFormat { c04, finals2000a };

//! Earth orientation parameters at a single epoch (one line of an EOP file), in SI units.
struct EarthOrientationRecord {
    //! Epoch of the record, as UTC modified Julian date.
    double utcModifiedJulianDate = std::numeric_limits< double >::quiet_NaN( );

    //! Polar motion x_p and y_p (radians).
    std::array< double, 2 > polarMotion = { std::numeric_limits< double >::quiet_NaN( ), std::numeric_limits< double >::quiet_NaN( ) };

    //! UT1 - UTC (seconds).
    double ut1MinusUtc = std::numeric_limits< double >::quiet_NaN( );

    //! Excess length of day (seconds); empty if not in the file.
    std::optional< double > lengthOfDay;

    //! Celestial pole offsets dX and dY (radians); empty if not in the file.
    std::optional< std::array< double, 2 > > celestialPoleOffsets;
};

//! Read Earth orientation records from legacy 16-column (14 C04) or modern 21-column (20 C04) C04 data, or from USNO
//! Bulletin A finals2000A data (celestial pole offsets as dX, dY). Returns the records with their UTC modified Julian date
//! as key. Throws with source and line number on malformed, empty, unsorted or duplicate data.
std::map< double, EarthOrientationRecord > readEarthOrientationData( std::istream& stream,
                                                                     const EarthOrientationFormat format,
                                                                     const std::string& source = "<stream>" );

//! Read Earth orientation records from a file (see stream version for details).
std::map< double, EarthOrientationRecord > readEarthOrientationData( const std::string& fileName,
                                                                     const EarthOrientationFormat format = EarthOrientationFormat::c04 );

}  // namespace earth_orientation
}  // namespace tudat

#endif  // TUDAT_READ_EARTH_ORIENTATION_DATA_H
