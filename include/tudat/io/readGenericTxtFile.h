/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */


#ifndef TUDAT_READ_GENERIC_TXT_FILE_H
#define TUDAT_READ_GENERIC_TXT_FILE_H

#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdarg>
#include <memory>
#include <boost/algorithm/string.hpp>

#include "tudat/io/fieldType.h"


namespace tudat
{
namespace input_output
{

class TxtFieldType
{
public:
    explicit TxtFieldType( std::string fieldName ): fieldName_( std::move( fieldName )) { }

    virtual ~TxtFieldType( ) = default;

    virtual double toDouble( const std::string& stringValue )
    {
        return std::stod( stringValue );
    }

    bool operator<( const TxtFieldType& other ) const
    {
        return ( this->fieldName_ < other.fieldName_ );
    }

private:
    std::string fieldName_;
};


std::shared_ptr< TxtFieldType > newTxtFieldType( std::string fieldName )
{
    return std::make_shared< TxtFieldType >( std::move( fieldName ));
}


class TxtMonthFieldType: public TxtFieldType
{
public:
    explicit TxtMonthFieldType( std::string fieldName ): TxtFieldType( std::move( fieldName )) { }

    double toDouble( const std::string& stringValue ) override
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
        const auto iter = months.find( stringValue );
        if ( iter != months.cend( )) {
            return iter->second;
        }
        throw std::runtime_error( "Invalid month name found in file" );
    }
};


std::shared_ptr< TxtMonthFieldType > newTxtMonthFieldType( std::string fieldName )
{
    return std::make_shared< TxtMonthFieldType >( std::move( fieldName ));
}

namespace txt_field_types
{

std::shared_ptr< TxtFieldType > spacecraftIdentifier = newTxtFieldType( "Spacecraft ID" );
std::shared_ptr< TxtFieldType > dsnTransmittingStation = newTxtFieldType( "DSN Transmitting Station" );
std::shared_ptr< TxtFieldType > dsnReceivingStation = newTxtFieldType( "DSN Receiving Station" );
std::shared_ptr< TxtFieldType > utcYear = newTxtFieldType( "UTC Year" );
std::shared_ptr< TxtFieldType > utcMonth = newTxtMonthFieldType( "UTC Month" );
std::shared_ptr< TxtFieldType > utcDay = newTxtFieldType( "UTC Day" );
std::shared_ptr< TxtFieldType > utcHour = newTxtFieldType( "UTC Hour" );
std::shared_ptr< TxtFieldType > utcMinute = newTxtFieldType( "UTC Minute" );
std::shared_ptr< TxtFieldType > utcSecond = newTxtFieldType( "UTC Second" );
std::shared_ptr< TxtFieldType > roundTripLightTimeMicroSec = newTxtFieldType( "Round Trip Light Time [micro s]" );
std::shared_ptr< TxtFieldType > roundTripLightTimeSec = newTxtFieldType( "Round Trip Light Time [s]" );
std::shared_ptr< TxtFieldType > lightTimeMeasurementAccuracyMicroSec = newTxtFieldType(
        "Measurement Accuracy of the Light Time [micro s]" );
std::shared_ptr< TxtFieldType > lightTimeMeasurementDelayMicroSec = newTxtFieldType(
        "Measurement Delay of the Light Time [micro s]" );
std::shared_ptr< TxtFieldType > tbdTimeJ2000 = newTxtFieldType( "TBD Time since J2000.0 [s]" );
std::shared_ptr< TxtFieldType > planetNumber = newTxtFieldType( "Number of the planet" );
std::shared_ptr< TxtFieldType > xPlanetFrame = newTxtFieldType( "x coordinate in planet reference frame [km]" );
std::shared_ptr< TxtFieldType > yPlanetFrame = newTxtFieldType( "y coordinate in planet reference frame [km]" );
std::shared_ptr< TxtFieldType > zPlanetFrame = newTxtFieldType( "z coordinate in planet reference frame [km]" );
std::shared_ptr< TxtFieldType > vXPlanetFrame = newTxtFieldType( "x velocity in planet reference frame [km/s]" );
std::shared_ptr< TxtFieldType > vYPlanetFrame = newTxtFieldType( "y velocity in planet reference frame [km/s]" );
std::shared_ptr< TxtFieldType > vZPlanetFrame = newTxtFieldType( "z velocity in planet reference frame [km/s]" );
}

class TxtFileContents
{
public:

    TxtFileContents( ) = default;

    TxtFileContents( std::string fileName,
                     std::vector< std::shared_ptr< TxtFieldType > > columnTypes,
                     char commentSymbol = '#',
                     std::string valueSeparators = ",: \t" )
            : fileName_( std::move( fileName )),
              columnTypes_( std::move( columnTypes )),
              commentSymbol_( commentSymbol ),
              valueSeparators_( std::move( valueSeparators ))
    {
        parseData( );
    }


    void parseData( )
    {
        std::ifstream dataFile( fileName_ );
        if ( !dataFile.good( )) {
            throw std::runtime_error(
                    "Error when opening Jpl Range file, file " + fileName_ + " could not be opened." );
        }
        initialiseDataMap( );
        fillDataMap( dataFile );
    }

    void initialiseDataMap( )
    {
        dataMap_.clear( );

        for ( std::shared_ptr< TxtFieldType >& columnType: columnTypes_ ) {
            dataMap_[columnType] = std::vector< double >( );
        }

    }

    void fillDataMap( std::ifstream& dataFile )
    {
        std::string currentLine;
        while ( std::getline( dataFile, currentLine )) {
            if ( !currentLine.empty( ) && currentLine.at( 0 ) != commentSymbol_ ) {
                addLineToDataMap( currentLine );
            }
        }
    }

    void addLineToDataMap( std::string& rawLine )
    {
        size_t numColumns = getNumColumns( );

        // Trim the line and split based on the separators
        boost::algorithm::trim( rawLine );
        boost::algorithm::split( splitRawLine_,
                                 rawLine,
                                 boost::is_any_of( valueSeparators_ ),
                                 boost::algorithm::token_compress_on );

        // Check if the expected number of columns is present in this line
        if ( splitRawLine_.size( ) != numColumns ) {
            unsigned int columnsFound = splitRawLine_.size( );
            throw std::runtime_error(

                    "The current line in file " + fileName_ + " has " + std::to_string( columnsFound ) +
                    " columns but " + std::to_string( numColumns ) + " columns were expected.\nRaw line:" + rawLine );
        }

        std::map< std::shared_ptr< TxtFieldType >, double > currentLineDataMap;

        // Populate the dataMap_ with a new row on each of the vectors
        for ( std::size_t i = 0; i < numColumns; ++i ) {
            std::shared_ptr< TxtFieldType > currentField = columnTypes_.at( i );
            double currentValue = currentField->toDouble( splitRawLine_.at( i ));
            dataMap_[currentField].push_back( currentValue );
            currentLineDataMap[currentField] = currentValue;
        }

        // Also add to a vector that stores a map per data point.
        dataVector_.push_back( currentLineDataMap );
    }

    size_t getNumColumns( ) const
    {
        return columnTypes_.size( );
    }

public:
    std::string fileName_ = "None";
    std::string separators_ = ":, \t";
    std::vector< std::shared_ptr< TxtFieldType > > columnTypes_;

    // FIXME: Maybe just one of the two is enough (dataMap_, dataVector_)
    std::map< std::shared_ptr< TxtFieldType >, std::vector< double >> dataMap_;
    std::vector< std::map< std::shared_ptr< TxtFieldType >, double>> dataVector_;
    char commentSymbol_;
    std::string valueSeparators_;

    // TODO: Maybe something to add metadata (delay, ...)

private:
    std::vector< std::string > splitRawLine_;
};


static inline std::unique_ptr< TxtFileContents > createTxtFileContents( const std::string& fileName,
                                                                        std::vector< std::shared_ptr< TxtFieldType > >& columnTypes,
                                                                        char commentSymbol = '#',
                                                                        const std::string& valueSeparators = ",: \t" )
{
    return std::make_unique< TxtFileContents >( fileName, columnTypes, commentSymbol, valueSeparators );
}

} // input_output
} // tudat

#endif // TUDAT_READ_GENERIC_TXT_FILE_H