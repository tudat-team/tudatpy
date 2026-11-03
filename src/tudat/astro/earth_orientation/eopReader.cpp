/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/earth_orientation/eopReader.h"
#include "tudat/astro/earth_orientation/readEarthOrientationData.h"

namespace tudat
{

namespace earth_orientation
{

//! Constructor
EOPReader::EOPReader( const std::string& eopFile, const std::string& format, const basic_astrodynamics::IAUConventions nutationTheory )
{
    EarthOrientationFormat fileFormat;
    if( format == "C04" )
    {
        fileFormat = EarthOrientationFormat::c04;
    }
    else if( format == "finals2000A" )
    {
        fileFormat = EarthOrientationFormat::finals2000a;
    }
    else
    {
        throw std::runtime_error( "Error, only C04 and finals2000A EOP file formats currently supported by reader, found " + format );
    }
    if( !( nutationTheory == basic_astrodynamics::iau_2000_a || nutationTheory == basic_astrodynamics::iau_2006 ) )
    {
        std::cerr << ( "Warning, only IAU2000 nutation theory format currently supported by reader." ) << std::endl;
    }

    for( const auto& [ modifiedJulianDate, record ] : readEarthOrientationData( eopFile, fileFormat ) )
    {
        cipInItrs[ modifiedJulianDate ] = Eigen::Vector2d( record.polarMotion[ 0 ], record.polarMotion[ 1 ] );
        ut1MinusUtc[ modifiedJulianDate ] = record.ut1MinusUtc;
        if( record.lengthOfDay )
        {
            lengthOfDayOffset[ modifiedJulianDate ] = *record.lengthOfDay;
        }
        if( record.celestialPoleOffsets )
        {
            cipInGcrsCorrection[ modifiedJulianDate ] =
                    Eigen::Vector2d( ( *record.celestialPoleOffsets )[ 0 ], ( *record.celestialPoleOffsets )[ 1 ] );
        }
    }

    // Add one hour buffer time
    cipInItrs[ cipInItrs.begin( )->first - 1.0 / 24.0 ] = cipInItrs.begin( )->second;
    ut1MinusUtc[ ut1MinusUtc.begin( )->first - 1.0 / 24.0 ] = ut1MinusUtc.begin( )->second;
    if( !cipInGcrsCorrection.empty( ) )
    {
        cipInGcrsCorrection[ cipInGcrsCorrection.begin( )->first - 1.0 / 24.0 ] = cipInGcrsCorrection.begin( )->second;
    }
    if( !lengthOfDayOffset.empty( ) )
    {
        lengthOfDayOffset[ lengthOfDayOffset.begin( )->first - 1.0 / 24.0 ] = lengthOfDayOffset.begin( )->second;
    }
}

}  // namespace earth_orientation

}  // namespace tudat
