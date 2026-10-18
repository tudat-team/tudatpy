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
#include "tudat/io/readEarthOrientationData.h"

namespace tudat
{
namespace earth_orientation
{

namespace
{
input_output::EarthOrientationData loadEarthOrientationFile( const std::string& fileName, const std::string& format )
{
    input_output::EarthOrientationFormat inputFormat;
    if( format == "C04" )
    {
        inputFormat = input_output::EarthOrientationFormat::c04;
    }
    else if( format == "finals2000A" )
    {
        inputFormat = input_output::EarthOrientationFormat::finals2000a;
    }
    else
    {
        throw std::runtime_error( "Unsupported Earth orientation format: " + format );
    }
    return input_output::readEarthOrientationData( fileName, inputFormat );
}
}  // namespace

EOPReader::EOPReader( const std::string& eopFile, const std::string& format, const basic_astrodynamics::IAUConventions nutationTheory ):
    EOPReader( loadEarthOrientationFile( eopFile, format ) )
{
    if( !( nutationTheory == basic_astrodynamics::iau_2000_a || nutationTheory == basic_astrodynamics::iau_2006 ) )
    {
        std::cerr << "Warning, only IAU2000 nutation theory format currently supported by reader." << std::endl;
    }
}

EOPReader::EOPReader( const input_output::EarthOrientationData& data ): rawData_( data )
{
    for( const auto& entry : data.records )
    {
        const double mjd = entry.first;
        const auto& record = entry.second;
        cipInItrs[ mjd ] = Eigen::Vector2d( record.polarMotion[ 0 ], record.polarMotion[ 1 ] );
        ut1MinusUtc[ mjd ] = record.ut1MinusUtc;
        if( record.lengthOfDay )
        {
            lengthOfDayOffset[ mjd ] = *record.lengthOfDay;
        }
        if( record.celestialPoleOffsets )
        {
            cipInGcrsCorrection[ mjd ] = Eigen::Vector2d( ( *record.celestialPoleOffsets )[ 0 ], ( *record.celestialPoleOffsets )[ 1 ] );
        }
    }
}

}  // namespace earth_orientation
}  // namespace tudat
