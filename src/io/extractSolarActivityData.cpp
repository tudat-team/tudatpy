/*    Copyright (c) 2010-2019, Delft University of Technology
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
 *                        http://celestrak.com/SpaceData/sw20110101.txt
 *      Data format explanation:
 *                        http://celestrak.com/SpaceData/SpaceWx-format.asp
 *
 */

#include "tudat/io/extractSolarActivityData.h"
#include "tudat/io/parseSolarActivityData.h"  // For the FieldTypes

namespace tudat
{

namespace input_output
{

namespace solar_activity
{

//! Extracts the solar activity data to a SolarActivityData container.
std::shared_ptr< SolarActivityData > ExtractSolarActivityData::extract(
        tudat::input_output::parsed_data_vector_utilities::ParsedDataLineMapPtr data )
{

    using namespace tudat::input_output::parsed_data_vector_utilities;
    using namespace tudat::input_output::field_types::solar_activity;
    using namespace tudat::input_output::field_types::time;

    // Check if the required fields are available
    checkRequiredFieldType( data, 33, year, month, day, bartelsSolarRotationNumber,
                            dayOfBartelsCycle, planetaryRangeIndex0to3, planetaryRangeIndex3to6,
                            planetaryRangeIndex6to9, planetaryRangeIndex9to12,
                            planetaryRangeIndex12to15, planetaryRangeIndex15to18,
                            planetaryRangeIndex18to21, planetaryRangeIndex21to24,
                            planetaryRangeIndexSum, planetaryEquivalentAmplitude0to3,
                            planetaryEquivalentAmplitude3to6, planetaryEquivalentAmplitude6to9,
                            planetaryEquivalentAmplitude9to12, planetaryEquivalentAmplitude12to15,
                            planetaryEquivalentAmplitude15to18, planetaryEquivalentAmplitude18to21,
                            planetaryEquivalentAmplitude21to24, planetaryEquivalentAmplitudeAverage,
                            planetaryDailyCharacterFigure, planetaryDailyCharacterFigureConverted,
                            internationalSunspotNumber, solarRadioFlux107Adjusted, fluxQualifier,
                            centered81DaySolarRadioFlux107Adjusted,
                            last81DaySolarRadioFlux107Adjusted, solarRadioFlux107Observed,
                            centered81DaySolarRadioFlux107Observed,
                            last81DaySolarRadioFlux107Observed );


    // Convert string data and append to solar activity data object
    int currentYear = getField< unsigned int >( data, year );
    int currentMonth = getField< unsigned int >( data, month );
    int currentDay = getField< unsigned int >( data, day );

    // Create the resulting solar activity data object (will be returned at the end)
    std::shared_ptr< SolarActivityData > solarActivityContainer =
        std::make_shared< SolarActivityData >( currentYear, currentMonth, currentDay );

    solarActivityContainer->bartelsSolarRotationNumber = getField< unsigned int >(
                data, bartelsSolarRotationNumber );
    solarActivityContainer->dayOfBartelsCycle = getField< unsigned int >( data, dayOfBartelsCycle );
    solarActivityContainer->solarRadioFlux107Adjusted = getField< double >(
                data, solarRadioFlux107Adjusted );
    solarActivityContainer->centered81DaySolarRadioFlux107Adjusted = getField< double >(
                data, centered81DaySolarRadioFlux107Adjusted );
    solarActivityContainer->last81DaySolarRadioFlux107Adjusted = getField< double >(
                data, last81DaySolarRadioFlux107Adjusted );
    solarActivityContainer->solarRadioFlux107Observed = getField< double >(
                data, solarRadioFlux107Observed );
    solarActivityContainer->centered81DaySolarRadioFlux107Observed = getField< double >(
                data, centered81DaySolarRadioFlux107Observed );
    solarActivityContainer->last81DaySolarRadioFlux107Observed = getField< double >(
                data, last81DaySolarRadioFlux107Observed );
    try
    {
        solarActivityContainer->dataType = getField< unsigned int >( data, dataType );
    }
    catch( std::bad_cast& )
    {
        solarActivityContainer->dataType = std::numeric_limits< unsigned int >::max( );
    }


    // Make sure only non-empty fields are extracted
    if ( !data->find( planetaryRangeIndex0to3 )->second->getRaw( ).empty( ) ) // check if string is empty
    {
        solarActivityContainer->planetaryRangeIndexSum = getField< unsigned int >(
                    data, planetaryRangeIndexSum);
        solarActivityContainer->planetaryEquivalentAmplitudeAverage = getField< unsigned int >(
                    data, planetaryEquivalentAmplitudeAverage);
        solarActivityContainer->planetaryRangeIndexVector = Eigen::VectorXd::Zero( 8 );
        solarActivityContainer->planetaryEquivalentAmplitudeVector = Eigen::VectorXd::Zero( 8 );
        solarActivityContainer->planetaryRangeIndexVector
                << getField< unsigned int >( data, planetaryRangeIndex0to3 ),
                getField< unsigned int >( data, planetaryRangeIndex3to6 ),
                getField< unsigned int >( data, planetaryRangeIndex6to9 ),
                getField< unsigned int >( data, planetaryRangeIndex9to12 ),
                getField< unsigned int >( data, planetaryRangeIndex12to15 ),
                getField< unsigned int >( data, planetaryRangeIndex15to18 ),
                getField< unsigned int >( data, planetaryRangeIndex18to21 ),
                getField< unsigned int >( data, planetaryRangeIndex21to24 );

        solarActivityContainer->planetaryEquivalentAmplitudeVector
                << getField< unsigned int >( data, planetaryEquivalentAmplitude0to3 ),
                getField< unsigned int >( data, planetaryEquivalentAmplitude3to6 ),
                getField< unsigned int >( data, planetaryEquivalentAmplitude6to9 ),
                getField< unsigned int >( data, planetaryEquivalentAmplitude9to12 ),
                getField< unsigned int >( data, planetaryEquivalentAmplitude12to15 ),
                getField< unsigned int >( data, planetaryEquivalentAmplitude15to18 ),
                getField< unsigned int >( data, planetaryEquivalentAmplitude18to21 ),
                getField< unsigned int >( data, planetaryEquivalentAmplitude21to24 );
    }

    if ( !data->find( planetaryDailyCharacterFigure )->second->getRaw( ).empty( ) )
    {
        solarActivityContainer->planetaryDailyCharacterFigure = getField< double >(
                    data, planetaryDailyCharacterFigure);
        solarActivityContainer->planetaryDailyCharacterFigureConverted = getField< unsigned int >(
                    data, planetaryDailyCharacterFigureConverted);
    }

    if  ( !data->find( internationalSunspotNumber )->second->getRaw( ).empty( ) )
    {
        solarActivityContainer->internationalSunspotNumber = getField< unsigned int >(
                    data, internationalSunspotNumber);
    }

    if ( !data->find( fluxQualifier )->second->getRaw( ).empty( ) )
    {
        solarActivityContainer->fluxQualifier = getField< unsigned int >( data, fluxQualifier );
    }

    return solarActivityContainer;
}

} // namespace solar_activity

} // namespace input_output

} // namespace tudat
