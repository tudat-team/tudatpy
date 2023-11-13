/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readJplRangeFile.h"


namespace tudat
{
namespace input_output
{

unsigned int getMonthIndex( std::string monthName )
{
    std::map< std::string, unsigned int > months {{ "JAN", 1 },
                                                  { "FEB", 2 },
                                                  { "MAR", 3 },
                                                  { "APR", 4 },
                                                  { "MAY", 5 },
                                                  { "JUN", 6 },
                                                  { "JUL", 7 },
                                                  { "AUG", 8 },
                                                  { "SEP", 9 },
                                                  { "OCT", 10 },
                                                  { "NOV", 11 },
                                                  { "DEC", 12 }};

    const auto iter = months.find( monthName );

    if ( iter != months.cend( )) {
        return iter->second;
    }
    throw std::runtime_error( "Invalid month name found in file" );
}

JplRangeDataBlock::JplRangeDataBlock( std::string rawBlock )
{
    boost::algorithm::trim( rawBlock );
    boost::algorithm::split( splitRawBlock_,
                             rawBlock,
                             boost::is_any_of( ",: \t" ),
                             boost::algorithm::token_compress_on );

    spacecraftId_ = std::stoi( splitRawBlock_.at( 0 ));
    transmittingStationId_ = std::stoi( splitRawBlock_.at( 1 ));
    receivingStationId_ = std::stoi( splitRawBlock_.at( 2 ));
    utcYear_ = std::stoi( splitRawBlock_.at( 3 ));
    utcMonth_ = getMonthIndex( splitRawBlock_.at( 4 ));
    utcDay_ = std::stoi( splitRawBlock_.at( 5 ));
    utcHour_ = std::stoi( splitRawBlock_.at( 6 ));
    utcMinute_ = std::stoi( splitRawBlock_.at( 7 ));
    utcSecond_ = std::stoi( splitRawBlock_.at( 8 ));
    roundTripLightTime_ = std::stod( splitRawBlock_.at( 9 ));
    transmissionDelay_ += std::stod( splitRawBlock_.at( 10 ));
}

JplRangeFileContents::JplRangeFileContents( const std::string& rangeFileName )
{
    fileName_ = rangeFileName;

    std::ifstream dataFile( rangeFileName );
    if ( !dataFile.good( )) {
        throw std::runtime_error(
                "Error when opening Jpl Range file, file " + rangeFileName + " could not be opened." );
    }

    readHeader( dataFile );
    checkFileFormat( );
    readDataBlocks( dataFile );
    removeTransponderDelay( );

}

void JplRangeFileContents::readHeader( std::ifstream& dataFile )
{
    fileHeaderInfo_.clear( );
    std::string currentLine;

    // Read the file title
    std::getline( dataFile, currentLine );
    fileTitle_ = currentLine.substr( currentLine.find( commentSymbol_ ) + 1, currentLine.find( '\n' ));
    boost::algorithm::trim( fileTitle_ );

    // Read the header
    while ( std::getline( dataFile, currentLine ) && !currentLine.empty( ) && currentLine.at( 0 ) == commentSymbol_ ) {
        fileHeaderInfo_ += ( currentLine + '\n' );
    }



    // Put the file pointer back to the beginning to avoid first line of data not being read.
    dataFile.seekg( 0, std::ios::beg );
}

void JplRangeFileContents::readDataBlocks( std::ifstream& dataFile )
{
    dataBlocks_.clear( );
    std::string currentLine;
    while ( std::getline( dataFile, currentLine )) {
        if ( !currentLine.empty( ) && currentLine.at( 0 ) != commentSymbol_ ) {
            std::shared_ptr< JplRangeDataBlock > newDataBlock = std::make_shared< JplRangeDataBlock >( currentLine );
            dataBlocks_.push_back( newDataBlock );
        }
    }

}

void JplRangeFileContents::checkFileFormat( )
{
    if ( fileHeaderInfo_.empty( )) {
        throw std::runtime_error( "Error when reading JPL file: file header information not found." );
        formatSupported_ = false;
    }

    if ( supportedFileTitles.find( fileTitle_ ) == supportedFileTitles.end( )) {
        throw std::runtime_error(
                "Error when reading JPL file: File Title not known. You might make slight adjustments to the "
                "format and add `# Manual Ranging Data` to the top of the file." );
    }

    if ( !formatSupported_ ) {
        throw std::runtime_error( "The file format of " + fileName_ + "is not currently supported." );
    }
}

void JplRangeFileContents::removeTransponderDelay( )
{
    if ( fileHeaderInfo_.find( "transponder delay of ~420 has NOT been removed" ) != std::string::npos ) {
        for ( auto dataBlock: dataBlocks_ ) {
            dataBlock->transmissionDelay_ += 420.0;
        }
    }
}

}
}