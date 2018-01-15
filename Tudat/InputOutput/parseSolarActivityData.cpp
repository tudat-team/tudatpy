/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Data file:
 *                        http://celestrak.com/SpaceData/sw19571001.txt
 *      Data format explanation:
 *                        http://celestrak.com/SpaceData/SpaceWx-format.asp
 *
 */

#include "Tudat/InputOutput/parseSolarActivityData.h"
#include "Tudat/InputOutput/fixedWidthParser.h"

#include <vector>
#include <string>

namespace tudat
{
namespace input_output
{
namespace solar_activity
{

//! Parses the stream of text.
void ParseSolarActivityData::parseStream( std::istream& fileContent)
{

    using namespace tudat::input_output::field_types::solar_activity;
    using namespace tudat::input_output::field_types::time;

    // Construct FixedWidthParser containing the fieldtypes occuring in the solar activity file
    tudat::input_output::FixedWidthParser solarParser(
            34, year, month, day, bartelsSolarRotationNumber, dayOfBartelsCycle,
            planetaryRangeIndex0to3, planetaryRangeIndex3to6, planetaryRangeIndex6to9,
            planetaryRangeIndex9to12, planetaryRangeIndex12to15, planetaryRangeIndex15to18,
            planetaryRangeIndex18to21, planetaryRangeIndex21to24, planetaryRangeIndexSum,
            planetaryEquivalentAmplitude0to3, planetaryEquivalentAmplitude3to6,
            planetaryEquivalentAmplitude6to9, planetaryEquivalentAmplitude9to12,
            planetaryEquivalentAmplitude12to15, planetaryEquivalentAmplitude15to18,
            planetaryEquivalentAmplitude18to21, planetaryEquivalentAmplitude21to24,
            planetaryEquivalentAmplitudeAverage, planetaryDailyCharacterFigure,
            planetaryDailyCharacterFigureConverted, internationalSunspotNumber,
            solarRadioFlux107Adjusted, fluxQualifier, centered81DaySolarRadioFlux107Adjusted,
            last81DaySolarRadioFlux107Adjusted, solarRadioFlux107Observed,
            centered81DaySolarRadioFlux107Observed, last81DaySolarRadioFlux107Observed, dataType,
            4, 3, 3, 5, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 4, 6, 2, 6,
            6, 6, 6, 6, 2 );

    // String containing line to be parsed
    std::string line;

    int dataType = 0;
    bool validdata = false;

    while ( !fileContent.fail( ) && !fileContent.eof( ) )   // Read stream line by line
    {
        std::getline( fileContent,line );

        // Determine dataType of line (observed/daily predicted/monthly predicted/monthly fit)
        if ( line.substr( 0, 14 ).compare( "BEGIN OBSERVED" ) == 0)
        {
            dataType = 1;
            validdata = true;
            continue;
        }

        if ( line.substr( 0, 21 ).compare( "BEGIN DAILY_PREDICTED" ) == 0 )
        {
            dataType = 2;
            validdata = true;
            continue;
        }

        if ( line.substr( 0, 23 ).compare( "BEGIN MONTHLY_PREDICTED" ) == 0 )
        {
            dataType = 3;
            validdata = true;
            continue;
        }

        if ( line.substr( 0, 17 ).compare( "BEGIN MONTHLY_FIT" ) == 0)
        {
            dataType = 4;
            validdata = true;
            continue;
        }

        if ( line.substr( 0, 12 ).compare( "END OBSERVED" )== 0  )
        {
            validdata = false;
            continue;
        }

        if ( line.substr( 0, 19 ).compare( "END DAILY_PREDICTED" ) == 0 )
        {
            validdata = false;
            continue;
        }

        if ( line.substr( 0, 21 ).compare( "END MONTHLY_PREDICTED" ) == 0 )
        {
            validdata = false;
            continue;
        }

        if ( line.substr( 0, 15 ).compare( "END MONTHLY_FIT" ) == 0 )
        {
            validdata = false;
            continue;
        }

        if ( validdata == true )
        {
            // add datatype at the end of the parsed line
            line = line + " " + std::to_string( dataType );

            parsedData->push_back( solarParser.parse( line )->at( 0 ) );
        }
    }
}

}   // namespace solar_activity
}   // namespace input_output
}   // namespace tudat
