//
// Created by simon on 11/7/23.
//

#ifndef TUDAT_READJPLRANGEFILE_H
#define TUDAT_READJPLRANGEFILE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "tudat/io/basicInputOutput.h"
#include "tudat/basics/timeType.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include <boost/algorithm/string.hpp>


namespace tudat
{
namespace input_output
{

int getMonthIndex( std::string monthName );

class JplRangeDataBlock
{
public:
    JplRangeDataBlock( std::string line );

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
    double measurementAccuracy_;
    std::vector< std::string > splitBlock_;

private:

};


class JplRangeFileContents
{
public:
    JplRangeFileContents( const std::string& rangeFileName );

    std::string getFileName( )
    {
        return fileName_;
    }

private:
    std::string fileName_;
    std::vector< std::shared_ptr< JplRangeDataBlock>> dataBlocks_;

};

}
}

#endif //TUDAT_READJPLRANGEFILE_H
