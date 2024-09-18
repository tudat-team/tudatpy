/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/estimation_setup/createAtmosphericLightTimeCorrection.h"

namespace tudat
{

namespace observation_models
{

std::shared_ptr< TabulatedMediaReferenceCorrection > createReferenceCorrection(
        const double startTime,
        const double endTime,
        const std::vector< double >& coefficients,
        const std::string& computationSpecifier )
{

    std::shared_ptr< TabulatedMediaReferenceCorrection > mediaCorrection;

    if ( computationSpecifier == "BY CONST" || computationSpecifier == "BY DCONST" )
    {
        mediaCorrection =std::make_shared< ConstantReferenceCorrection >(
                        startTime, endTime, coefficients.front( ) );
    }
    else if ( computationSpecifier == "BY NRMPOW" || computationSpecifier == "BY DNRMPOW" )
    {
        mediaCorrection = std::make_shared< PowerSeriesReferenceCorrection >(
                        startTime, endTime, coefficients );
    }
    else if ( computationSpecifier == "BY TRIG" || computationSpecifier == "BY DTRIG" )
    {
        mediaCorrection = std::make_shared< FourierSeriesReferenceCorrection >(
                        startTime, endTime, coefficients );
    }
    else
    {
        throw std::runtime_error( "Error when creating media reference correction object: computation specifier " +
            computationSpecifier + " is not valid." );
    }

    return mediaCorrection;
}

AtmosphericCorrectionPerStationAndSpacecraftType extractAtmosphericCorrection(
        std::vector< std::shared_ptr< input_output::CspRawFile > > rawCspFiles,
        const std::string& modelIdentifier,
        const std::map< int, std::string >& spacecraftNamePerSpacecraftId,
        const std::map< int, std::string >& quasarNamePerQuasarId )
{
    // Sort CSP files by start date
    std::stable_sort( rawCspFiles.begin( ), rawCspFiles.end( ), &input_output::compareAtmosphericCspFileStartDate );

    observation_models::AtmosphericCorrectionPerStationAndSpacecraftType troposphericCorrection;

    // Keys: std::pair< station, source (i.e. spacecraft or quasar) >, data type identifier
    std::map< std::pair< std::string, std::string >, std::map< std::string,
        std::shared_ptr< observation_models::TabulatedMediaReferenceCorrectionManager > > >
        correctionsPerStationsAndSourcePerType;

    for ( unsigned int i = 0; i < rawCspFiles.size( ); ++i )
    {
        for ( unsigned int j = 0; j < rawCspFiles.at( i )->getCspCommands( ).size( ); ++j )
        {
            std::shared_ptr< input_output::AtmosphericCorrectionCspCommand > cspCommand = std::dynamic_pointer_cast< input_output::AtmosphericCorrectionCspCommand >(
                    rawCspFiles.at( i )->getCspCommands( ).at( j ) );
            if ( cspCommand == nullptr )
            {
                throw std::runtime_error( "Error when creating tabulated atmospheric corrections: inconsistent CSP command type." );
            }

            if ( modelIdentifier != cspCommand->modelIdentifier_ )
            {
                continue;
            }

            std::string sourceName = input_output::getSourceName(  cspCommand->sourceSpecifier_, cspCommand->sourceId_,
                                                     spacecraftNamePerSpacecraftId, quasarNamePerQuasarId );
            if ( ( modelIdentifier == "DRY NUPART" || modelIdentifier == "WET NUPART" ) && !sourceName.empty( ) )
            {
                throw std::runtime_error( "Error when creating tabulated atmospheric corrections: invalid source name was"
                                          "created for tropospheric corrections." );
            }
            std::pair< std::string, std::string > stationsAndSource = std::make_pair(
                    cspCommand->groundStationsId_, sourceName );

            // Check whether new media correction manager object should be created
            bool createNewObject = false;
            if ( correctionsPerStationsAndSourcePerType.count( stationsAndSource ) == 0 ||
                correctionsPerStationsAndSourcePerType.at( stationsAndSource ).count( cspCommand->dataTypesIdentifier_ ) == 0 )
            {
                createNewObject = true;
            }

            if ( createNewObject )
            {
                correctionsPerStationsAndSourcePerType[ stationsAndSource ][ cspCommand->dataTypesIdentifier_ ] =
                        std::make_shared< observation_models::TabulatedMediaReferenceCorrectionManager >( );
            }

            // Create correction calculator and save it to correction manager
            correctionsPerStationsAndSourcePerType.at( stationsAndSource ).at(
                    cspCommand->dataTypesIdentifier_ )->pushReferenceCorrectionCalculator( createReferenceCorrection(
                            cspCommand->startTime_, cspCommand->endTime_,
                            cspCommand->computationCoefficients_, cspCommand->computationSpecifier_ ) );
        }
    }

    // Loop over created calculator managers and assign them to map with ground station and observation type as keys
    for ( auto stationsAndSourceIt = correctionsPerStationsAndSourcePerType.begin( );
            stationsAndSourceIt != correctionsPerStationsAndSourcePerType.end( ); ++stationsAndSourceIt )
    {
        std::vector< std::string > groundStations = input_output::getGroundStationsNames( stationsAndSourceIt->first.first );
        std::string source = stationsAndSourceIt->first.second;

        for ( auto observableIt = stationsAndSourceIt->second.begin( ); observableIt != stationsAndSourceIt->second.end( ); ++observableIt )
        {
            std::vector< observation_models::ObservableType > observableTypes = input_output::getBaseObservableTypes(
                observableIt->first );

            for ( const std::string& groundStation : groundStations )
            {
                for ( const observation_models::ObservableType& observableType : observableTypes )
                {
                    troposphericCorrection[ std::make_pair( groundStation, source ) ][ observableType ] =
                            correctionsPerStationsAndSourcePerType.at( stationsAndSourceIt->first ).at( observableIt->first );
                }
            }
        }
    }

    return troposphericCorrection;
}

AtmosphericCorrectionPerStationAndSpacecraftType extractTroposphericDryCorrectionAdjustment(
        const std::vector< std::shared_ptr< input_output::CspRawFile > >& rawCspFiles )
{
    std::string modelIdentifier = "DRY NUPART";
    return extractAtmosphericCorrection( rawCspFiles, modelIdentifier );
}

AtmosphericCorrectionPerStationAndSpacecraftType extractTroposphericWetCorrectionAdjustment(
        const std::vector< std::shared_ptr< input_output::CspRawFile > >& rawCspFiles )
{
    std::string modelIdentifier = "WET NUPART";
    return extractAtmosphericCorrection( rawCspFiles, modelIdentifier );
}

AtmosphericCorrectionPerStationAndSpacecraftType extractIonosphericCorrection(
        const std::vector< std::shared_ptr< input_output::CspRawFile > >& rawCspFiles,
        const std::map< int, std::string >& spacecraftNamePerSpacecraftId,
        const std::map< int, std::string >& quasarNamePerQuasarId )
{
    std::string modelIdentifier = "CHPART";
    return extractAtmosphericCorrection( rawCspFiles, modelIdentifier, spacecraftNamePerSpacecraftId,
                                         quasarNamePerQuasarId );
}

AtmosphericCorrectionPerStationAndSpacecraftType extractDefaultTroposphericDryCorrection( )
{
    std::string modelIdentifier = "DRY NUPART";
    return extractAtmosphericCorrection(
            { input_output::getDsnDefaultTroposphericSeasonalModelCspFile( ) }, modelIdentifier );
}

AtmosphericCorrectionPerStationAndSpacecraftType extractDefaultTroposphericWetCorrection( )
{
    std::string modelIdentifier = "WET NUPART";
    return extractAtmosphericCorrection(
            { input_output::getDsnDefaultTroposphericSeasonalModelCspFile( ) }, modelIdentifier );
}

} // namespace observation_models

} // namespace tudat