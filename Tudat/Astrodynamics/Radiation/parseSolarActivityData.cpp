/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120607    A. Ronse          Creation of code.
 *
 *    References
 *      Data file:
 *                        http://celestrak.com/SpaceData/sw19571001.txt
 *      Data format explanation:
 *                        http://celestrak.com/SpaceData/SpaceWx-format.asp
 *
 */

#include "Tudat/Astrodynamics/Radiation/parseSolarActivityData.h"
#include "Tudat/InputOutput/fixedWidthParser.h"

#include <vector>

#include <boost/lexical_cast.hpp>

namespace tudat
{
namespace radiation
{
namespace solar_activity
{

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
        std::getline (fileContent,line);

        // Determine dataType of line (observed/daily predicted/monthly predicted/monthly fit)
        if (line.substr( 0, 14 ).compare("BEGIN OBSERVED") == 0)
        {
            dataType = 1;
            validdata = true;
            continue;
        }

        if (line.substr( 0, 21 ).compare("BEGIN DAILY_PREDICTED")== 0)
        {
            dataType = 2;
            validdata = true;
            continue;
        }

        if (line.substr( 0, 23 ).compare("BEGIN MONTHLY_PREDICTED")== 0)
        {
            dataType = 3;
            validdata = true;
            continue;
        }

        if (line.substr( 0, 17 ).compare("BEGIN MONTHLY_FIT") == 0)
        {
            dataType = 4;
            validdata = true;
            continue;
        }

        if (line.substr( 0, 12 ).compare("END OBSERVED")== 0  )
        {
            validdata = false;
            continue;
        }
        if (line.substr( 0, 19 ).compare("END DAILY_PREDICTED")== 0)
        {
            validdata = false;
            continue;
        }

        if (line.substr( 0, 21 ).compare("END MONTHLY_PREDICTED")== 0)
        {
            validdata = false;
            continue;
        }

        if (line.substr( 0, 15 ).compare("END MONTHLY_FIT")== 0)
        {
            validdata = false;
            continue;
        }

        if (validdata == true)
        {
            // add datatype at the end of the parsed line
            line = line + " " + boost::lexical_cast<string>( dataType );

            parsedData->push_back( solarParser.parse( line )->at(0) );
        }
    }
}

}   // namespace solar_activity
}   // namespace radiation
}   // namespace tudat
