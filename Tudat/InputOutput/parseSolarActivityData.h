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
 *                        http://celestrak.com/SpaceData/sw20110101.txt
 *      Data format explanation:
 *                        http://celestrak.com/SpaceData/SpaceWx-format.asp
 *
 */

#ifndef TUDAT_PARSESOLARACTIVITY_H
#define TUDAT_PARSESOLARACTIVITY_H

#include <iostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "Tudat/InputOutput/parsedDataVectorUtilities.h"
#include "Tudat/InputOutput/textParser.h"

namespace tudat
{
namespace input_output
{
namespace field_types
{
namespace solar_activity
{
    static const FieldType bartelsSolarRotationNumber = hash_constructor(
            "Solar: Bartels_solar_rotation number" );
    static const FieldType dayOfBartelsCycle = hash_constructor(
            "Solar: Day_of_Bartels_cycle" );
    static const FieldType planetaryRangeIndex0to3 = hash_constructor(
            "Solar: Planetary_range_index_0to3" );
    static const FieldType planetaryRangeIndex3to6 = hash_constructor(
            "Solar: Planetary_range_index_3to6" );
    static const FieldType planetaryRangeIndex6to9 = hash_constructor(
            "Solar: Planetary_range_index_6to9" );
    static const FieldType planetaryRangeIndex9to12 = hash_constructor(
            "Solar: Planetary_range_index_9to12" );
    static const FieldType planetaryRangeIndex12to15 = hash_constructor(
            "Solar: Planetary_range_index_12to15" );
    static const FieldType planetaryRangeIndex15to18 = hash_constructor(
            "Solar: Planetary_range_index_15to18" );
    static const FieldType planetaryRangeIndex18to21 = hash_constructor(
            "Solar: Planetary_range_index_18to21" );
    static const FieldType planetaryRangeIndex21to24 = hash_constructor(
            "Solar: Planetary_range_index_21to24" );
    static const FieldType planetaryRangeIndexSum = hash_constructor(
            "Solar: Planetary_range_index_Sum" );
    static const FieldType planetaryEquivalentAmplitude0to3 = hash_constructor(
            "Solar: Planetary_equivalent_amplitude_0to3" );
    static const FieldType planetaryEquivalentAmplitude3to6 = hash_constructor(
            "Solar: Planetary_equivalent_amplitude_3to6" );
    static const FieldType planetaryEquivalentAmplitude6to9 = hash_constructor(
            "Solar: Planetary_equivalent_amplitude_6to9" );
    static const FieldType planetaryEquivalentAmplitude9to12 = hash_constructor(
            "Solar: Planetary_equivalent_amplitude_9to12" );
    static const FieldType planetaryEquivalentAmplitude12to15 = hash_constructor(
            "Solar: Planetary_equivalent_amplitude_12to15" );
    static const FieldType planetaryEquivalentAmplitude15to18 = hash_constructor(
            "Solar: Planetary_equivalent_amplitude_15to18" );
    static const FieldType planetaryEquivalentAmplitude18to21 = hash_constructor(
            "Solar: Planetary_equivalent_amplitude_18to21" );
    static const FieldType planetaryEquivalentAmplitude21to24 = hash_constructor(
            "Solar: Planetary_equivalent_amplitude_21to24" );
    static const FieldType planetaryEquivalentAmplitudeAverage = hash_constructor(
            "Solar: Planetary_equivalent_amplitude_average" );
    static const FieldType planetaryDailyCharacterFigure = hash_constructor(
            "Solar: Planetary_daily_character_figure" );
    static const FieldType planetaryDailyCharacterFigureConverted =  hash_constructor(
            "Solar: Planetary_daily_character_figure_converted" );
    static const FieldType internationalSunspotNumber = hash_constructor(
            "Solar: International_Sunspot_Number" );
    static const FieldType solarRadioFlux107Adjusted = hash_constructor(
            "Solar: Solar_radio_Flux_10.7_Adjusted" );
    static const FieldType fluxQualifier = hash_constructor( "Solar: Flux_qualifier" );
    static const FieldType centered81DaySolarRadioFlux107Adjusted = hash_constructor(
            "Solar: Centered_81_Day_Solar_Radio_Flux_10.7_Adjusted" );
    static const FieldType last81DaySolarRadioFlux107Adjusted = hash_constructor(
            "Solar: Last_81_Day_Solar_Radio_Flux_10.7_Adjusted" );
    static const FieldType solarRadioFlux107Observed = hash_constructor(
            "Solar: Solar_Radio_Flux_10.7_Observed" );
    static const FieldType centered81DaySolarRadioFlux107Observed = hash_constructor(
            "Solar: Centered_81_Day_Solar_Radio_Flux_10.7_Observed" );
    static const FieldType last81DaySolarRadioFlux107Observed = hash_constructor(
            "Solar: Last_81_Day_Solar_Radio_Flux_10.7_Observed" );
    static const FieldType dataType = hash_constructor( "Data: type");

}   // solar_activity

namespace time
{
    static const FieldType year = hash_constructor( "Time: Year" );
    static const FieldType month = hash_constructor( "Time: Month" );
    static const FieldType day = hash_constructor( "Time: Day" );

}   // time
}   // field_types
}   // input_output

namespace input_output
{
namespace solar_activity
{

using std::string;

//! Solar activity parser class.
/*!
 * This class implements a fixed width parser specifically designed for parsing solar activity
 * files (space weather data) in the format of http://celestrak.com/SpaceData/sw19571001.txt.
 */
class ParseSolarActivityData : public input_output::TextParser
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ParseSolarActivityData(): TextParser( true ) { }

protected:

    //! Parses the stream of text.
    /*!
     * Parses the stream of (solarActivity) data by detecting datatype
     * and subsequently applying a fixed width parse process.
     */
    void parseStream( std::istream& fileContent );

private:
};

} // namespace solar_activity
} // namespace input_output
} // namerspace tudat

#endif // TUDAT_PARSESOLARACTIVITY_H
