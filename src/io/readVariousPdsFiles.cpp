/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readVariousPdsFiles.h"
#include "boost/algorithm/string.hpp"

namespace tudat
{
namespace input_output
{

std::pair< std::vector< double >, std::vector< double > > grailAntennaFileReader( const std::string antennaFile )
{
    std::vector< double > switchTimes;
    std::vector< double > antennaPositions;

    // Open file
    std::ifstream dataFile( antennaFile );
    if ( !dataFile.good( ) )
    {
        throw std::runtime_error( "Error when opening GRAIL antenna file, file " + antennaFile +  " could not be opened." );
    }

    unsigned int numberHeaderLines = 20;
    unsigned int numberLinesParsed = 0;
    bool isHeaderParsed = false;

    // Line based parsing
    std::string line;
    while ( dataFile.good( ) && !dataFile.eof( ) )
    {
        // Get line from stream
        std::getline( dataFile, line );

        // Skip empty lines
        std::string trimmedLine = line;
        boost::algorithm::trim( trimmedLine );
        if ( trimmedLine.empty( ) )
        {
            continue;
        }

        if ( isHeaderParsed )
        {
            std::vector< std::string > vectorOfIndividualStrings;
            boost::algorithm::split( vectorOfIndividualStrings, line,
                                     boost::algorithm::is_any_of( " " ), boost::algorithm::token_compress_on );
            // Check if last string is empty and remove it if so
            if ( vectorOfIndividualStrings.back( ).empty( ) )
            {
                vectorOfIndividualStrings.pop_back( );
            }

            if ( vectorOfIndividualStrings.size( ) != 7 )
            {
                throw std::runtime_error( "Error when reading antenna switch file for GRAIL, inconsistent format. Number of entries is "
                    + std::to_string( vectorOfIndividualStrings.size( ) ) + ", should be 7." );
            }

            std::string emptyFlags = vectorOfIndividualStrings.at( 6 );
            if ( emptyFlags.compare( "00000000" ) != 1 )
            {
                throw std::runtime_error( "Error when reading antenna switch file for GRAIL, inconsistent entry. The last eights digits should be set to zero." );
            }

            double switchTimeTbd = std::stod( vectorOfIndividualStrings.at( 0 ) );
            double antennaPositionNorm = std::stod( vectorOfIndividualStrings.at( 2 ) );
            double xAxisCosinus = std::stod( vectorOfIndividualStrings.at( 3 ) );
            double yAxisCosinus = std::stod( vectorOfIndividualStrings.at( 4 ) );
            double zAxisCosinus = std::stod( vectorOfIndividualStrings.at( 5 ) );
            Eigen::Vector3d antennaPosition = antennaPositionNorm * ( Eigen::Vector3d( ) << xAxisCosinus, yAxisCosinus, zAxisCosinus  ).finished( );

            switchTimes.push_back( switchTimeTbd );
            antennaPositions.push_back( antennaPosition[ 0 ] );
            antennaPositions.push_back( antennaPosition[ 1 ] );
            antennaPositions.push_back( antennaPosition[ 2 ] );
        }

        numberLinesParsed++;
        if ( numberLinesParsed >= numberHeaderLines  )
        {
            isHeaderParsed = true;
        }
    }
    dataFile.close( );

    return std::make_pair( switchTimes, antennaPositions );
}

std::map< double, double > grailMassLevel0FileReader( const std::string massFile )
{
    std::map< double, double > massHistory;

    // Open file
    std::ifstream dataFile( massFile );
    if ( !dataFile.good( ) )
    {
        throw std::runtime_error( "Error when opening GRAIL mass file (level 0), file " + massFile +  " could not be opened." );
    }

    unsigned int numberHeaderLines = 4;
    unsigned int numberLinesParsed = 0;

    bool isHeaderParsed = false;

    // Line based parsing
    std::string line;
    while ( dataFile.good( ) && !dataFile.eof( ) )
    {
        // Get line from stream
        std::getline( dataFile, line );

        // Skip empty lines
        std::string trimmedLine = line;
        boost::algorithm::trim( trimmedLine );
        if ( trimmedLine.empty( ) )
        {
            continue;
        }

        if ( isHeaderParsed )
        {
            std::vector< std::string > vectorOfIndividualStrings;

            std::string delimiter = "   ";
            unsigned int positionInLine = 0;
            while( positionInLine <= line.size( ) )
            {
                positionInLine = line.find( delimiter );
                vectorOfIndividualStrings.push_back( line.substr( 0, positionInLine ) );
                line.erase( 0, positionInLine + delimiter.size( ) );
            }

//            // Check if last string is empty and remove it if so
//            if ( vectorOfIndividualStrings.back( ).empty( ) )
//            {
//                vectorOfIndividualStrings.pop_back( );
//            }

            if ( vectorOfIndividualStrings.size( ) != 13 )
            {
                throw std::runtime_error( "Error when reading mass level 0 file for GRAIL, inconsistent format. Number of entries is "
                                          + std::to_string( vectorOfIndividualStrings.size( ) ) + ", should be 13." );
            }

            // Retrieve manoeuvre date
            std::string dateUtc = vectorOfIndividualStrings.at( 1 );
            std::vector< std::string > individualStrDate;
            boost::algorithm::split( individualStrDate, dateUtc, boost::algorithm::is_any_of( "/" ), boost::algorithm::token_compress_on );

            // Retrieve manoeuvre end time (UTC)
            std::string manoeuvreEndTimeUtc = vectorOfIndividualStrings.at( 3 );
            std::vector< std::string > individualStrTime;
            boost::algorithm::split( individualStrTime, manoeuvreEndTimeUtc, boost::algorithm::is_any_of( ":" ), boost::algorithm::token_compress_on );

            // Compute time UTC in seconds
            int year = std::stoi( individualStrDate.at( 2 ) );
            int month = std::stoi( individualStrDate.at( 0 ) );
            int day = std::stoi( individualStrDate.at( 1 ) );
            int hour = std::stoi( individualStrTime.at( 0 ) );
            int minutes = std::stoi( individualStrTime.at( 1 ) );
            double seconds = stod( individualStrTime.at( 2 ) );
            double timeInSecondsUtc = basic_astrodynamics::timeFromDecomposedDateTime< double >( year, month, day, hour, minutes, seconds );


            double timeInSecondsTdb = earth_orientation::TerrestrialTimeScaleConverter( ).getCurrentTime< double >(
                    basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, timeInSecondsUtc, Eigen::Vector3d::Zero( ) );
            double mass = std::stod( vectorOfIndividualStrings.at( 6 ) );
            massHistory[ timeInSecondsTdb ] = mass;
        }

        numberLinesParsed++;
        if ( numberLinesParsed >= numberHeaderLines  )
        {
            isHeaderParsed = true;
        }
    }

    dataFile.close( );

    return massHistory;
}

std::map< double, double > grailMassLevel1FileReader( const std::string massFile,
                                                      const std::string dataLevel )
{
    if ( dataLevel != "1a" && dataLevel != "1b" )
    {
        throw std::runtime_error( "Error when reading level 1 mass file for GRAIL, the input dataLevel should be set to either 1a or 1b." );
    }

    std::map< double, double > massHistory;

    // Open file
    std::ifstream dataFile( massFile );
    if ( !dataFile.good( ) )
    {
        throw std::runtime_error( "Error when opening GRAIL mass file (level 1), file " + massFile +  " could not be opened." );
    }

    unsigned int numberHeaderLines = 20;
    unsigned int numberLinesParsed = 0;

    bool isHeaderParsed = false;

    // Line based parsing
    std::string line;
    while ( dataFile.good( ) && !dataFile.eof( ) )
    {
        // Get line from stream
        std::getline( dataFile, line );

        // Skip empty lines
        std::string trimmedLine = line;
        boost::algorithm::trim( trimmedLine );
        if ( trimmedLine.empty( ) )
        {
            continue;
        }

        if ( isHeaderParsed )
        {
            std::vector< std::string > vectorOfIndividualStrings;
            boost::algorithm::split( vectorOfIndividualStrings, line,
                                     boost::algorithm::is_any_of( " " ), boost::algorithm::token_compress_on );
            // Check if last string is empty and remove it if so
            if ( vectorOfIndividualStrings.back( ).empty( ) )
            {
                vectorOfIndividualStrings.pop_back( );
            }

            if ( vectorOfIndividualStrings.size( ) != 7 )
            {
                throw std::runtime_error( "Error when reading mass level 1 file for GRAIL, inconsistent format. Number of entries is "
                                          + std::to_string( vectorOfIndividualStrings.size( ) ) + ", should be 7." );
            }

            std::string firstFlagSet = vectorOfIndividualStrings.at( 4 );
            std::string secondFlagSet = vectorOfIndividualStrings.at( 5 );
            if ( firstFlagSet != "00000000" )
            {
                throw std::runtime_error( "Error when reading mass level 1 file for GRAIL, inconsistent entry. "
                                          "In the first set of flags, all digits should be set to zero." );
            }
            if ( secondFlagSet != "00000001" )
            {
                throw std::runtime_error( "Error when reading mass level 1 file for GRAIL, inconsistent entry. "
                                          "In the second set of flags, only digit 0 should be set to 1." );
            }

            // Retrieve mass and time
            double timeTdb;
            if ( dataLevel == "1a" )
            {
               double timeUtcScet = std::stod( vectorOfIndividualStrings.at( 0 ) ) + std::stod( vectorOfIndividualStrings.at( 1 ) );
               timeTdb = earth_orientation::TerrestrialTimeScaleConverter( ).getCurrentTime< double >(
                       basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, timeUtcScet, Eigen::Vector3d::Zero( ) );
            }
            if ( dataLevel == "1b" )
            {
                timeTdb = std::stod( vectorOfIndividualStrings.at( 0 ) ) + std::stod( vectorOfIndividualStrings.at( 1 ) );
            }
            double mass = std::stod( vectorOfIndividualStrings.at( 6 ) );

            massHistory[ timeTdb ] = mass;
        }

        numberLinesParsed++;
        if ( numberLinesParsed >= numberHeaderLines  )
        {
            isHeaderParsed = true;
        }
    }

    dataFile.close( );

    return massHistory;
}



} // namespace input_output

} // namespace tudat
