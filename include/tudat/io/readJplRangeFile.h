//
// Created by simon on 11/7/23.
//

#ifndef TUDAT_READJPLRANGEFILE_H
#define TUDAT_READJPLRANGEFILE_H

#include <iostream>
#include <fstream>


namespace tudat
{
namespace input_output
{

class JPLRangeFileContents
{
public:
    JPLRangeFileContents( const std::string& rangeFileName );
private:
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




};

}
}

#endif //TUDAT_READJPLRANGEFILE_H
