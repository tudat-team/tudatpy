/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/io/readDsnAmcFile.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace input_output
{

std::vector< DsnAmcRecord > readDsnAmcFile( const std::string& fileName )
{
    std::ifstream file( fileName );
    if( !file.is_open( ) )
    {
        throw std::runtime_error( "Error when reading DSN AMC file, could not open " + fileName );
    }

    // AMC times are given in seconds since 1970
    const double unixEpoch = basic_astrodynamics::DateTime( 1970, 1, 1, 0, 0, 0.0 ).epoch< double >( );

    std::vector< DsnAmcRecord > records;
    std::string line;
    while( std::getline( file, line ) )
    {
        // Skip empty and comment lines
        std::size_t firstCharacter = line.find_first_not_of( " \t\r" );
        if( firstCharacter == std::string::npos || line.at( firstCharacter ) == '#' )
        {
            continue;
        }

        // Wind speed and direction (last two entries) are only given in later (e.g. Juno-era) files, not in earlier (e.g.
        // Cassini-era) files
        std::istringstream entryStream( line );
        const long numberOfEntries =
                std::distance( std::istream_iterator< std::string >( entryStream ), std::istream_iterator< std::string >( ) );

        DsnAmcRecord record;
        record.windSpeed = TUDAT_NAN;
        record.windDirection = TUDAT_NAN;
        std::string recordType;
        std::istringstream lineStream( line );
        if( ( numberOfEntries != 14 && numberOfEntries != 16 ) ||
            !( lineStream >> recordType >> record.utcTime >> record.azimuth >> record.elevation >> record.zenithWetDelay >>
               record.zenithDryDelay >> record.cloudLiquid >> record.algorithmId >> record.cloudFlag >> record.wvrFlag >> record.mtpFlag >>
               record.surfaceFlag >> record.windFlag >> record.pathDelayFlag ) ||
            recordType != "P" || ( numberOfEntries == 16 && !( lineStream >> record.windSpeed >> record.windDirection ) ) )
        {
            throw std::runtime_error( "Error when reading DSN AMC file " + fileName + ", could not parse line: " + line );
        }

        record.utcTime += unixEpoch;
        record.zenithWetDelay /= 100.0;
        record.zenithDryDelay /= 100.0;
        records.push_back( record );
    }
    return records;
}

std::map< double, Eigen::Vector2d > getValidDsnAmcZenithDelays( const std::vector< DsnAmcRecord >& records,
                                                                const int minimumRejectedCloudFlag,
                                                                const bool rejectOnElevationFlag,
                                                                const bool rejectOnMtpFlag,
                                                                const bool rejectOnSurfaceFlag,
                                                                const bool rejectOnWindFlag )
{
    std::map< double, Eigen::Vector2d > zenithDelays;
    for( const DsnAmcRecord& record : records )
    {
        // Path delay flag 2 (elevation outside prescribed range) may be accepted for statistical retrievals only (no retrieval is
        // performed by the Bayesian algorithm 200)
        const bool isPathDelayFlagAccepted =
                record.pathDelayFlag == 0 || ( record.pathDelayFlag == 2 && !rejectOnElevationFlag && record.algorithmId != "200" );

        // Retain records with acceptable quality flags and valid wet path delay (9.99 m is the 999 cm placeholder for a missing
        // retrieval)
        if( record.wvrFlag == 0 && isPathDelayFlagAccepted && record.cloudFlag < minimumRejectedCloudFlag &&
            !( rejectOnMtpFlag && record.mtpFlag != 0 ) && !( rejectOnSurfaceFlag && record.surfaceFlag != 0 ) &&
            !( rejectOnWindFlag && record.windFlag != 0 ) && record.zenithWetDelay >= 0.0 && record.zenithWetDelay < 9.99 )
        {
            // Store invalid dry path delay as NaN, so that the wet path delay can still be used
            const double zenithDryDelay =
                    ( record.zenithDryDelay >= 0.0 && record.zenithDryDelay < 9.99 ) ? record.zenithDryDelay : TUDAT_NAN;
            zenithDelays[ record.utcTime ] = ( Eigen::Vector2d( ) << zenithDryDelay, record.zenithWetDelay ).finished( );
        }
    }
    return zenithDelays;
}

}  // namespace input_output
}  // namespace tudat
