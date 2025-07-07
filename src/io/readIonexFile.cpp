/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readIonexFile.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <boost/algorithm/string.hpp>

namespace tudat
{
namespace input_output
{

double parseIonexEpochToSecondsSinceJ2000( const std::string& line )
{
    int year, month, day, hour, minute;
    double second;
    std::stringstream ss( line.substr( 0, 60 ) );
    ss >> year >> month >> day >> hour >> minute >> second;

    double julianDate = basic_astrodynamics::convertCalendarDateToJulianDay< double >(
        year, month, day, hour, minute, second );
    return basic_astrodynamics::convertJulianDayToSecondsSinceEpoch(
        julianDate, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
}

void readIonexFile( const std::string& filePath, IonexTecMap& data )
{
    std::ifstream file( filePath );
    if ( !file )
    {
        throw std::runtime_error( "Cannot open IONEX file: " + filePath );
    }

    std::string line;
    std::vector< double > epochQueue;
    std::vector< Eigen::MatrixXd > mapQueue;

    double currentEpoch = -1.0;

    while ( std::getline( file, line ) )
    {
        boost::algorithm::trim( line );

        if ( line.find( "LAT1 / LAT2 / DLAT" ) != std::string::npos )
        {
            std::stringstream ss( line );
            ss >> data.latMax >> data.latMin >> data.dLat;

            data.latitudes.clear( );
            for ( double lat = data.latMax; lat >= data.latMin - 1e-3; lat += data.dLat )
            {
                data.latitudes.push_back( lat );
            }
        }
        else if ( line.find( "LON1 / LON2 / DLON" ) != std::string::npos )
        {
            std::stringstream ss( line );
            ss >> data.lonMin >> data.lonMax >> data.dLon;

            data.longitudes.clear( );
            for ( double lon = data.lonMin; lon <= data.lonMax + 1e-3; lon += data.dLon )
            {
                data.longitudes.push_back( lon );
            }
        }
        else if ( line.find( "HGT1 / HGT2 / DHGT" ) != std::string::npos )
        {
            std::stringstream ss( line );
            ss >> data.hgtMin >> data.hgtMax >> data.dHgt;
            data.referenceIonosphereHeight_ = data.hgtMin;  // assume single-shell height = HGT1
        }
        else if ( line.find( "EPOCH OF FIRST MAP" ) != std::string::npos )
        {
            data.epochStart = parseIonexEpochToSecondsSinceJ2000( line );
        }
        else if ( line.find( "EPOCH OF LAST MAP" ) != std::string::npos )
        {
            data.epochEnd = parseIonexEpochToSecondsSinceJ2000( line );
        }
        else if ( line.find( "EPOCH OF CURRENT MAP" ) != std::string::npos )
        {
            currentEpoch = parseIonexEpochToSecondsSinceJ2000( line );
            epochQueue.push_back( currentEpoch );
        }
        else if ( line.find( "START OF TEC MAP" ) != std::string::npos )
        {
            Eigen::MatrixXd tecMap = Eigen::MatrixXd::Zero( data.latitudes.size( ), data.longitudes.size( ) );
            std::vector< std::vector< bool > > fillMask(
                data.latitudes.size( ), std::vector< bool >( data.longitudes.size( ), false ) );

            while ( std::getline( file, line ) )
            {
                if ( line.find( "END OF TEC MAP" ) != std::string::npos )
                {
                    break;
                }

                if ( line.find( "LAT/LON1/LON2/DLON/H" ) != std::string::npos )
                {
                    double lat, lonStart, lonEnd, lonStep, height;
                    std::stringstream ss( line.substr( 0, 60 ) );
                    ss >> lat >> lonStart >> lonEnd >> lonStep >> height;

                    std::size_t latIndex = std::distance( data.latitudes.begin( ),
                        std::find( data.latitudes.begin( ), data.latitudes.end( ), lat ) );

                    std::vector< double > tecValues;
                    while ( tecValues.size( ) < data.longitudes.size( ) )
                    {
                        std::getline( file, line );
                        for ( std::size_t i = 0; i + 5 <= line.size( ); i += 5 )
                        {
                            std::string val = line.substr( i, 5 );
                            tecValues.push_back( std::stod( val ) * 0.1 );
                        }
                    }

                    for ( std::size_t j = 0; j < tecValues.size( ); ++j )
                    {
                        double lon = lonStart + lonStep * j;
                        std::size_t lonIndex = std::distance( data.longitudes.begin( ),
                            std::find( data.longitudes.begin( ), data.longitudes.end( ), lon ) );

                        if ( latIndex >= data.latitudes.size( ) || lonIndex >= data.longitudes.size( ) )
                        {
                            throw std::runtime_error( "LAT/LON index out of bounds in TEC map." );
                        }

                        tecMap( latIndex, lonIndex ) = tecValues.at( j );
                        fillMask[ latIndex ][ lonIndex ] = true;
                    }
                }
            }

            // Confirm completeness
            for ( std::size_t i = 0; i < data.latitudes.size( ); ++i )
            {
                for ( std::size_t j = 0; j < data.longitudes.size( ); ++j )
                {
                    if ( !fillMask[ i ][ j ] )
                    {
                        std::cerr << "Warning: TEC map at epoch " << currentEpoch
                                  << " missing value at lat=" << data.latitudes[ i ]
                                  << ", lon=" << data.longitudes[ j ] << std::endl;
                    }
                }
            }

            mapQueue.push_back( tecMap );
        }
    }

    file.close( );

    // Reverse latitudes to ascending order (required by Tudat interpolators)
    std::reverse( data.latitudes.begin(), data.latitudes.end( ) );
    for ( Eigen::MatrixXd& mat : mapQueue )
    {
        mat = mat.colwise().reverse().eval();  // flip rows to match new lat order
    }

    // Store to output structure
    if ( epochQueue.size( ) != mapQueue.size( ) )
    {
        throw std::runtime_error( "IONEX epoch list and TEC map count mismatch." );
    }

    for ( std::size_t i = 0; i < epochQueue.size( ); ++i )
    {
        data.epochs.push_back( epochQueue[ i ] );
        data.tecMaps[ epochQueue[ i ] ] = mapQueue[ i ];
    }

    data.validate( );
    data.printMetadata( );
}


void readIonexFiles( const std::vector< std::string >& filePaths, IonexTecMap& data )
{
    for ( const auto& path : filePaths )
    {
        readIonexFile( path, data );
    }
}

} // namespace input_output
} // namespace tudat
