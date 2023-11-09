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

int getMonthIndex( std::string monthName )
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
    boost::algorithm::split( splitBlock_, rawBlock, boost::is_any_of( ",: \t" ), boost::algorithm::token_compress_on );

    spacecraftId_ = std::stoi( splitBlock_.at( 0 ));
    transmittingStationId_ = std::stoi( splitBlock_.at( 1 ));
    receivingStationId_ = std::stoi( splitBlock_.at( 2 ));
    utcYear_ = std::stoi( splitBlock_.at( 3 ));
    utcMonth_ = getMonthIndex( splitBlock_.at( 4 ));
    utcDay_ = std::stoi( splitBlock_.at( 5 ));
    utcHour_ = std::stoi( splitBlock_.at( 6 ));
    utcMinute_ = std::stoi( splitBlock_.at( 7 ));
    utcSecond_ = std::stoi( splitBlock_.at( 8 ));
    roundTripLightTime_ = std::stod( splitBlock_.at( 9 ));
    measurementAccuracy_ = std::stod( splitBlock_.at( 10 ));
}

JplRangeFileContents::JplRangeFileContents( const std::string& rangeFileName )
{
    fileName_ = rangeFileName;
}

}
}