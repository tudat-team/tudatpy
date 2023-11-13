//
// Created by simon on 11/7/23.
//

#ifndef TUDAT_READJPLRANGEFILE_H
#define TUDAT_READJPLRANGEFILE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

#include "tudat/io/basicInputOutput.h"
#include "tudat/basics/timeType.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include <boost/algorithm/string.hpp>


namespace tudat
{
namespace input_output
{

unsigned int getMonthIndex( std::string monthName );


class JplRangeDataBlock
{
public:
    explicit JplRangeDataBlock( std::string line );

    unsigned int spacecraftId_;
    unsigned int transmittingStationId_;
    unsigned int receivingStationId_;
    unsigned int utcYear_;
    unsigned int utcMonth_;
    unsigned int utcDay_;
    unsigned int utcHour_;
    unsigned int utcMinute_;
    unsigned int utcSecond_;
    double roundTripLightTime_;
    double transmissionDelay_ = 0.0;
    std::vector< std::string > splitRawBlock_;

private:

};


class JplRangeFileContents
{
public:
    explicit JplRangeFileContents( const std::string& rangeFileName );

    void readHeader( std::ifstream& dataFile );

    void readDataBlocks( std::ifstream& dataFile );

    void checkFileFormat( );

    void removeTransponderDelay();


public:
    std::string fileName_;
    std::string fileTitle_;
    std::string fileHeaderInfo_;
    bool formatSupported_ = true;
    char commentSymbol_ = '#';
    std::vector< std::shared_ptr< JplRangeDataBlock>> dataBlocks_;

private:
    std::set< std::string > supportedFileTitles = { "Mars Pathfinder Ranging Data", "Viking lander range data","MESSENGER range data", "Mars Global Surveyor range data", "Mars Odyssey range data", "Mars Reconnaissance Orbiter range data" "Manual Ranging Data" };

};


inline std::shared_ptr< JplRangeFileContents > readJplRangeFile( std::string& fileName )
{
    return std::make_shared< JplRangeFileContents >( fileName );
}

}
}

#endif //TUDAT_READJPLRANGEFILE_H
