/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved.
 *    This file is part of Tudat and is distributed under the Modified BSD license.
 */

#include "tudat/io/readDsnAmcFile.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "tudat/astro/basic_astro/dateTime.h"

namespace tudat
{
namespace input_output
{

std::shared_ptr< DsnAmcData > readDsnAmcFile( const std::string& fileName )
{
    std::ifstream file( fileName );
    if( !file.is_open( ) )
    {
        throw std::runtime_error( "Could not open DSN AMC file: " + fileName );
    }

    const double unixEpoch = basic_astrodynamics::DateTime( 1970, 1, 1, 0, 0, 0.0 ).epoch< double >( );
    std::vector< DsnAmcRecord > records;
    std::string line;
    unsigned int lineNumber = 0;
    while( std::getline( file, line ) )
    {
        ++lineNumber;
        const auto firstCharacter = line.find_first_not_of( " \t\r" );
        if( firstCharacter == std::string::npos || line[ firstCharacter ] == '#' )
        {
            continue;
        }

        DsnAmcRecord record;
        std::string identifier, trailingField;
        std::istringstream fields( line );
        if( !( fields >> identifier >> record.utcTime >> record.azimuth >> record.elevation >> record.zenithWetDelay >>
               record.zenithDryDelay >> record.cloudLiquid >> record.algorithmId >> record.cloudFlag >> record.wvrFlag >> record.mtpFlag >>
               record.surfaceFlag >> record.windFlag >> record.pathDelayFlag >> record.windSpeed >> record.windDirection ) ||
            identifier != "P" || fields >> trailingField )
        {
            throw std::runtime_error( "Malformed DSN AMC record in " + fileName + " at line " + std::to_string( lineNumber ) );
        }

        for( const double value : { record.utcTime,
                                    record.azimuth,
                                    record.elevation,
                                    record.zenithWetDelay,
                                    record.zenithDryDelay,
                                    record.cloudLiquid,
                                    record.windSpeed,
                                    record.windDirection } )
        {
            if( !std::isfinite( value ) )
            {
                throw std::runtime_error( "Non-finite DSN AMC field in " + fileName + " at line " + std::to_string( lineNumber ) );
            }
        }
        record.utcTime += unixEpoch;
        record.zenithWetDelay *= 0.01;
        record.zenithDryDelay *= 0.01;
        if( !records.empty( ) && record.utcTime <= records.back( ).utcTime )
        {
            throw std::runtime_error( "Non-increasing DSN AMC epoch in " + fileName + " at line " + std::to_string( lineNumber ) );
        }
        records.push_back( record );
    }
    if( file.bad( ) || records.empty( ) )
    {
        throw std::runtime_error( "Could not read DSN AMC records from: " + fileName );
    }
    return std::make_shared< DsnAmcData >( records, fileName );
}

}  // namespace input_output
}  // namespace tudat
