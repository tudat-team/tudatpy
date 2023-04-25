/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readTabulatedMediaCorrections.h"

#include <map>
#include <fstream>
#include <string>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include "tudat/astro/basic_astro/timeConversions.h"

namespace tudat
{

namespace input_output
{

AtmosphericCorrectionCspCommand::AtmosphericCorrectionCspCommand( std::vector< std::string >& cspCommand )
{
    for ( unsigned int i = 0; i < cspCommand.size(); ++i )
    {
        boost::algorithm::trim( cspCommand.at( i ) );

        if ( cspCommand.at( i ) == "MODEL" )
        {
            modelIdentifier_ = cspCommand.at( i + 1 );
            ++i;
        }
        else if ( cspCommand.at( i ) == "FROM" )
        {
            startTime_ = convertTime( cspCommand.at( i + 1 ), cspCommand.at( i + 2 ) );
            i = i + 2;
        }
        else if ( cspCommand.at( i ) == "TO" )
        {
            endTime_ = convertTime( cspCommand.at( i + 1 ), cspCommand.at( i + 2 ) );
            i = i + 2;
        }
        else if ( !boost::algorithm::find_first( cspCommand.at( i ), "BY " ).empty( )  )
        {
            computationSpecifier_ = cspCommand.at( i );
            for ( ++i; i < cspCommand.size(); ++i )
            {
                try
                {
                    computationCoefficients_.push_back( std::stod( cspCommand.at( i ) ) );
                }
                catch( ... )
                {
                    --i;
                    break;
                }
            }
        }
        else if ( cspCommand.at( i ) == "DSN" )
        {
            groundStationsId_ = cspCommand.at( i + 1 );
            ++i;
        }
        else if ( cspCommand.at( i ) == "ADJUST" )
        {
            dataTypesIdentifier_ = cspCommand.at( i + 1 );
            ++i;
        }
        else
        {
            throw std::runtime_error( "Invalid CSP command for atmospheric correction file." );
        }

    }

}

double AtmosphericCorrectionCspCommand::convertTime( std::string yearMonthDay, std::string hoursMinutesSeconds )
{
    std::vector< std::string > vectorOfIndividualStrings;

    // Extract years, months, days
    boost::algorithm::split( vectorOfIndividualStrings, yearMonthDay, boost::algorithm::is_any_of( "/" ),
                             boost::algorithm::token_compress_on );
    if ( vectorOfIndividualStrings.size( ) != 3 )
    {
        throw std::runtime_error("Error when extracting year/month/day from CSP command: number of extracted elements (" +
            std::to_string( vectorOfIndividualStrings.size( ) ) + ") is invalid.");
    }
    int year = std::stoi( vectorOfIndividualStrings.at( 0 ) );
    if ( year >= 00 && year <= 68 )
    {
        year += 2000;
    }
    else if ( year >= 69 && year <= 99 )
    {
        year += 1900;
    }
    else
    {
        throw std::runtime_error("Error when extracting year from CSP command: " + std::to_string( year ) +
            " is not a valid year (valid years are 00 through 99).");
    }
    int month = std::stoi( vectorOfIndividualStrings.at( 1 ) );
    int day = std::stoi( vectorOfIndividualStrings.at( 2 ) );

    // Extract hours, minutes and seconds
    boost::algorithm::split( vectorOfIndividualStrings, hoursMinutesSeconds, boost::algorithm::is_any_of( ":" ),
                             boost::algorithm::token_compress_on );
    if ( vectorOfIndividualStrings.size( ) != 2 && vectorOfIndividualStrings.size( ) != 3 )
    {
        throw std::runtime_error("Error when extracting hours:minutes:seconds from CSP command: number of extracted elements (" +
            std::to_string( vectorOfIndividualStrings.size( ) ) + ") is invalid.");
    }
    int hours = std::stoi( vectorOfIndividualStrings.at( 0 ) );
    int minutes = std::stoi( vectorOfIndividualStrings.at( 1 ) );
    double seconds = 0.0;
    if ( vectorOfIndividualStrings.size( ) == 3 )
    {
        seconds = std::stod( vectorOfIndividualStrings.at( 2 ) );
    }

    return basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
            year, month, day, hours, minutes, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * physical_constants::JULIAN_DAY;
}

IonosphericCorrectionCspCommand::IonosphericCorrectionCspCommand( std::vector< std::string >& cspCommand ):
    AtmosphericCorrectionCspCommand( cspCommand )
{
    for ( unsigned int i = 0; i < cspCommand.size(); ++i )
    {
        boost::algorithm::trim( cspCommand.at( i ) );

        if ( cspCommand.at( i ) == "SCID" || cspCommand.at( i ) == "QUASAR" )
        {
            sourceSpecifier_ = cspCommand.at( i );
            sourceId_ = std::stoi( cspCommand.at( i + 1) );
            break;
        }
    }
}

CspRawFile::CspRawFile( const std::string& cspFile ):
    fileName_( cspFile )
{
    // Read CSP commands. Each command saved as one string
    std::vector< std::string > cspCommandsVector = readCspCommandsFile( fileName_ );

    for ( unsigned int i = 0; i < cspCommandsVector.size( ); ++i )
    {
        std::vector< std::string > cspVectorOfIndividualStrings;
        boost::algorithm::split( cspVectorOfIndividualStrings, cspCommandsVector.at( i ),
                                 boost::algorithm::is_any_of( ",()" ), boost::algorithm::token_compress_on );
        // Check if last string is empty and remove it if so
        if ( cspVectorOfIndividualStrings.back( ).empty( ) )
        {
            cspVectorOfIndividualStrings.pop_back( );
        }

        std::shared_ptr< CspCommand > cspCommand;

        // Loop over strings and find type of model
        for ( unsigned int j = 0; j < cspVectorOfIndividualStrings.size( ); ++ j )
        {
            if ( cspVectorOfIndividualStrings.at( j ) == "MODEL" )
            {
                if ( cspVectorOfIndividualStrings.at( j + 1 ) == "DRY NUPART" ||
                    cspVectorOfIndividualStrings.at( j + 1 ) == "WET NUPART" )
                {
                    cspCommand = std::make_shared< AtmosphericCorrectionCspCommand >( cspVectorOfIndividualStrings );
                }
                else if ( cspVectorOfIndividualStrings.at( j + 1 ) == "CHPART" )
                {
                    cspCommand = std::make_shared< IonosphericCorrectionCspCommand >( cspVectorOfIndividualStrings );
                }
                else
                {
                    throw std::runtime_error( "Error when reading CSP file: invalid model type." );
                }

                cspCommands_.push_back( cspCommand );
                break;
            }
        }
    }
}

std::vector< std::string > CspRawFile::readCspCommandsFile( std::string file )
{

    // Open file
    std::ifstream stream( file, std::ios_base::binary );
    if ( !stream.good( ) )
    {
        throw std::runtime_error( "Error when opening file: " + file );
    }

    // Line based parsing
    std::string line;
    std::vector< std::string > commandVector;
    bool newCommand = true;
    while ( stream.good( ) && !stream.eof( ) )
    {
        // Get line from stream
        std::getline( stream, line );

        // Trim input string (removes all leading and trailing whitespaces).
        boost::algorithm::trim( line );

        // Skip empty and comment lines
        if ( line.size( ) <= 0 || ( line.at( 0 ) == '#' ) )
        {
            continue;
        }

        // Search and remove in-line comments
        auto iteratorRange = boost::algorithm::find_first( line, "#" );
        if( !iteratorRange.empty( ) )
        {
            line.erase( iteratorRange.begin( ), line.end( ) );
            // Re-trim string
            boost::algorithm::trim( line );
        }

        // If continuation of previous command, merge it with previous portions
        if ( !newCommand )
        {
            // Check if it's the last part of the command
            auto iteratorRange = boost::algorithm::find_first( line, ")." );
            if ( !iteratorRange.empty( ) )
            {
                // Remove command terminator (i.e. the period)
                line.erase( ++iteratorRange.begin( ), line.end( ) );
                newCommand = true;
            }

            // Merge with previous part of command
            commandVector.back( ) += line;
        }
        // If new command, push it to vector
        else
        {
            commandVector.push_back( line );
            newCommand = false;
        }
    }

    // Close file
    stream.close();

    return commandVector;
}

std::shared_ptr< observation_models::TabulatedMediaReferenceCorrection > createReferenceCorrection(
        double startTime,
        double endTime,
        std::vector< double > coefficients,
        std::string computationSpecifier )
{

    std::shared_ptr< observation_models::TabulatedMediaReferenceCorrection > mediaCorrection;

    if ( computationSpecifier == "BY CONST" || computationSpecifier == "BY DCONST" )
    {
        mediaCorrection =std::make_shared< observation_models::ConstantReferenceCorrection >(
                        startTime, endTime, coefficients.front( ) );
    }
    else if ( computationSpecifier == "BY NRMPOW" || computationSpecifier == "BY DNRMPOW" )
    {
        mediaCorrection = std::make_shared< observation_models::PowerSeriesReferenceCorrection >(
                        startTime, endTime, coefficients );
    }
    else if ( computationSpecifier == "BY TRIG" || computationSpecifier == "BY DTRIG" )
    {
        mediaCorrection = std::make_shared< observation_models::FourierSeriesReferenceCorrection >(
                        startTime, endTime, coefficients );
    }
    else
    {
        throw std::runtime_error( "Error when creating media reference correction object: computation specifier " +
            computationSpecifier + " is not valid." );
    }

    return mediaCorrection;
}

std::vector< std::string > getGroundStationsNames( std::string groundStationIdentifier )
{
    std::vector< std::string > groundStations;

    if ( groundStationIdentifier == "C10" )
    {
        groundStations = { "DSS-13", "DSS-14", "DSS-15", "DSS-24", "DSS-25", "DSS-26" };
    }
    else if ( groundStationIdentifier == "C40" )
    {
        groundStations = { "DSS-34", "DSS-35", "DSS-36", "DSS-43", "DSS-45" };
    }
    else if ( groundStationIdentifier == "C60" )
    {
        groundStations = { "DSS-54", "DSS-55", "DSS-63", "DSS-65" };
    }
    else
    {
        groundStations = { "DSS-" + groundStationIdentifier };
    }

    return groundStations;
}

std::vector< observation_models::ObservableType > getBaseObservableTypes( std::string observableTypeIdentifier )
{
    std::vector< observation_models::ObservableType > baseObservableTypes;

    if ( observableTypeIdentifier == "DOPRNG" )
    {
        baseObservableTypes = { observation_models::one_way_range, observation_models::one_way_doppler };
    }
    else if ( observableTypeIdentifier == "DOPPLER" )
    {
        baseObservableTypes = { observation_models::one_way_doppler };
    }
    else if ( observableTypeIdentifier == "RANGE" )
    {
        baseObservableTypes = { observation_models::one_way_range };
    }
    else if ( observableTypeIdentifier == "VLBI" )
    {

    }
    else if ( observableTypeIdentifier == "ALL" )
    {
        baseObservableTypes = { observation_models::one_way_range, observation_models::one_way_doppler };
    }
    else
    {
        throw std::runtime_error( "Invalid CSP command observable type identifier." );
    }

    return baseObservableTypes;
}

bool compareAtmosphericCspFileStartDate( std::shared_ptr< CspRawFile > rawCspData1,
                                         std::shared_ptr< CspRawFile > rawCspData2 )
{

    std::shared_ptr< AtmosphericCorrectionCspCommand > cspCommand1 = std::dynamic_pointer_cast< AtmosphericCorrectionCspCommand >(
                rawCspData1->getCspCommands( ).at( 0 ) );
    std::shared_ptr< AtmosphericCorrectionCspCommand > cspCommand2 = std::dynamic_pointer_cast< AtmosphericCorrectionCspCommand >(
                rawCspData2->getCspCommands( ).at( 0 ) );

    if ( cspCommand1 == nullptr || cspCommand2 == nullptr )
    {
        throw std::runtime_error( "Error when comparing atmospheric CSP files: inconsistent CSP command type." );
    }

    if ( cspCommand1->startTime_ < cspCommand2->startTime_ )
    {
        return true;
    }
    else
    {
        return false;
    }
}

observation_models::AtmosphericCorrectionPerStationAndSpacecraftType createTroposphericCorrection(
        std::vector< std::shared_ptr< CspRawFile > > rawCspFiles,
        const std::string& modelIdentifier )
{
    // Sort CSP files by start date
    std::stable_sort( rawCspFiles.begin( ), rawCspFiles.end( ), &compareAtmosphericCspFileStartDate );

    observation_models::AtmosphericCorrectionPerStationAndSpacecraftType troposphericCorrection;

    std::map< std::string, std::map< std::string, std::shared_ptr< observation_models::TabulatedMediaReferenceCorrectionManager > > >
            troposphericCorrectionPerStationPerType;

    for ( unsigned int i = 0; i < rawCspFiles.size( ); ++i )
    {
        for ( unsigned int j = 0; j < rawCspFiles.at( i )->getCspCommands( ).size( ); ++j )
        {
            std::shared_ptr< AtmosphericCorrectionCspCommand > cspCommand = std::dynamic_pointer_cast< AtmosphericCorrectionCspCommand >(
                    rawCspFiles.at( i )->getCspCommands( ).at( j ) );
            if ( cspCommand == nullptr )
            {
                throw std::runtime_error( "Error when creating tabulated atmospheric corrections: inconsistent CSP command type." );
            }

            if ( modelIdentifier != cspCommand->modelIdentifier_ )
            {
                continue;
            }

            // Check whether new media correction manager object should be created
            bool createNewObject = false;
            if ( troposphericCorrectionPerStationPerType.count( cspCommand->groundStationsId_ ) == 0 )
            {
                createNewObject = true;
            }
            else if ( troposphericCorrectionPerStationPerType.at( cspCommand->groundStationsId_ ).count(
                    cspCommand->dataTypesIdentifier_ ) == 0 )
            {
                createNewObject = true;
            }

            if ( createNewObject )
            {
                troposphericCorrectionPerStationPerType[ cspCommand->groundStationsId_ ][ cspCommand->dataTypesIdentifier_ ] =
                        std::make_shared< observation_models::TabulatedMediaReferenceCorrectionManager >( );
            }

            // Create correction calculator and save it to correction manager
            troposphericCorrectionPerStationPerType.at( cspCommand->groundStationsId_ ).at(
                    cspCommand->dataTypesIdentifier_ )->pushReferenceCorrectionCalculator( createReferenceCorrection(
                            cspCommand->startTime_, cspCommand->endTime_,
                            cspCommand->computationCoefficients_, cspCommand->computationSpecifier_ ) );
        }
    }

    // Loop over created calculator managers and assign them to map with ground station and observation type as keys
    for ( auto stationIt = troposphericCorrectionPerStationPerType.begin( );
            stationIt != troposphericCorrectionPerStationPerType.end( ); ++stationIt )
    {
        std::vector< std::string > groundStations = getGroundStationsNames( stationIt->first );

        for ( auto observableIt = stationIt->second.begin( ); observableIt != stationIt->second.end( ); ++observableIt )
        {
            std::vector< observation_models::ObservableType > observableTypes = getBaseObservableTypes(
                observableIt->first );

            for ( const std::string& groundStation : groundStations )
            {
                for ( const observation_models::ObservableType& observableType : observableTypes )
                {
                    troposphericCorrection[ std::make_pair( groundStation, "" ) ][ observableType ] =
                            troposphericCorrectionPerStationPerType.at( stationIt->first ).at( observableIt->first );
                }
            }
        }
    }

    return troposphericCorrection;
}

observation_models::AtmosphericCorrectionPerStationAndSpacecraftType createTroposphericDryCorrection(
        const std::vector< std::shared_ptr< CspRawFile > >& rawCspFiles )
{
    std::string modelIdentifier = "DRY NUPART";
    return createTroposphericCorrection( rawCspFiles, modelIdentifier );
}

observation_models::AtmosphericCorrectionPerStationAndSpacecraftType createTroposphericWetCorrection(
        const std::vector< std::shared_ptr< CspRawFile > >& rawCspFiles )
{
    std::string modelIdentifier = "WET NUPART";
    return createTroposphericCorrection( rawCspFiles, modelIdentifier );
}

} // namespace input_output

} // namespace tudat