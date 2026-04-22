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
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <set>

namespace tudat
{
namespace input_output
{

namespace
{

double parseIonexEpochToSecondsSinceJ2000( const std::string& line )
{
    int year, month, day, hour, minute;
    double second;
    std::stringstream ss( line.substr( 0, 60 ) );
    ss >> year >> month >> day >> hour >> minute >> second;

    double julianDate = basic_astrodynamics::convertCalendarDateToJulianDay< double >( year, month, day, hour, minute, second );
    return basic_astrodynamics::convertJulianDayToSecondsSinceEpoch( julianDate, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
}

//! Find the index of a value in a sorted grid using tolerance-based lookup
std::size_t findGridIndex( const std::vector< double >& grid, const double value, const double tolerance = 0.01 )
{
    for( std::size_t i = 0; i < grid.size( ); ++i )
    {
        if( std::fabs( grid[ i ] - value ) < tolerance )
        {
            return i;
        }
    }
    throw std::runtime_error( "IONEX: Value " + std::to_string( value ) + " not found in grid." );
}

//! Parse a single map block (TEC or RMS) from the current file position.
//! Expects the file stream to be positioned just after the START OF ... MAP line.
//! Also extracts the epoch from the EPOCH OF CURRENT MAP line within the block.
Eigen::MatrixXd parseMapBlock( std::ifstream& file, const IonexTecMap& data,
                               const std::string& endMarker, const int exponent,
                               double& epoch )
{
    Eigen::MatrixXd mapData = Eigen::MatrixXd::Zero( data.latitudes.size( ), data.longitudes.size( ) );
    std::string line;
    const double scaleFactor = std::pow( 10.0, exponent );

    while( std::getline( file, line ) )
    {
        if( line.find( endMarker ) != std::string::npos )
        {
            break;
        }

        if( line.find( "EPOCH OF CURRENT MAP" ) != std::string::npos )
        {
            epoch = parseIonexEpochToSecondsSinceJ2000( line );
        }
        else if( line.find( "LAT/LON1/LON2/DLON/H" ) != std::string::npos )
        {
            double lat, lonStart, lonEnd, lonStep, height;
            std::stringstream ss( line.substr( 0, 60 ) );
            ss >> lat >> lonStart >> lonEnd >> lonStep >> height;

            std::size_t latIndex = findGridIndex( data.latitudes, lat );

            std::vector< double > values;
            while( values.size( ) < data.longitudes.size( ) )
            {
                std::getline( file, line );
                for( std::size_t i = 0; i + 5 <= line.size( ) && values.size( ) < data.longitudes.size( ); i += 5 )
                {
                    std::string val = line.substr( i, 5 );
                    boost::algorithm::trim( val );
                    if( !val.empty( ) )
                    {
                        values.push_back( std::stod( val ) * scaleFactor );
                    }
                }
            }

            for( std::size_t j = 0; j < values.size( ) && j < data.longitudes.size( ); ++j )
            {
                mapData( latIndex, j ) = values[ j ];
            }
        }
    }

    return mapData;
}

}  // anonymous namespace

void readIonexFile( const std::string& filePath, IonexTecMap& data, const bool loadRmsMaps )
{
    std::ifstream file( filePath );
    if( !file )
    {
        throw std::runtime_error( "Cannot open IONEX file: " + filePath );
    }

    std::string line;
    std::vector< double > epochQueue;
    std::vector< Eigen::MatrixXd > tecMapQueue;
    std::vector< Eigen::MatrixXd > rmsMapQueue;

    while( std::getline( file, line ) )
    {
        boost::algorithm::trim( line );

        if( line.find( "LAT1 / LAT2 / DLAT" ) != std::string::npos )
        {
            std::stringstream ss( line );
            ss >> data.latMax >> data.latMin >> data.dLat;

            data.latitudes.clear( );
            for( double lat = data.latMax; lat >= data.latMin - 1e-3; lat += data.dLat )
            {
                data.latitudes.push_back( lat );
            }
        }
        else if( line.find( "LON1 / LON2 / DLON" ) != std::string::npos )
        {
            std::stringstream ss( line );
            ss >> data.lonMin >> data.lonMax >> data.dLon;

            data.longitudes.clear( );
            for( double lon = data.lonMin; lon <= data.lonMax + 1e-3; lon += data.dLon )
            {
                data.longitudes.push_back( lon );
            }
        }
        else if( line.find( "HGT1 / HGT2 / DHGT" ) != std::string::npos )
        {
            std::stringstream ss( line );
            ss >> data.hgtMin >> data.hgtMax >> data.dHgt;
            data.referenceIonosphereHeight_ = data.hgtMin * 1.0e3;  // convert km to m
        }
        else if( line.find( "EXPONENT" ) != std::string::npos )
        {
            std::stringstream ss( line );
            ss >> data.exponent;
        }
        else if( line.find( "EPOCH OF FIRST MAP" ) != std::string::npos )
        {
            data.epochStart = parseIonexEpochToSecondsSinceJ2000( line );
        }
        else if( line.find( "EPOCH OF LAST MAP" ) != std::string::npos )
        {
            data.epochEnd = parseIonexEpochToSecondsSinceJ2000( line );
        }
        else if( line.find( "START OF TEC MAP" ) != std::string::npos )
        {
            double epoch = -1.0;
            tecMapQueue.push_back( parseMapBlock( file, data, "END OF TEC MAP", data.exponent, epoch ) );
            epochQueue.push_back( epoch );
        }
        else if( line.find( "START OF RMS MAP" ) != std::string::npos )
        {
            if( loadRmsMaps )
            {
                double rmsEpoch = -1.0;
                rmsMapQueue.push_back( parseMapBlock( file, data, "END OF RMS MAP", data.exponent, rmsEpoch ) );
            }
            else
            {
                // Skip RMS block
                while( std::getline( file, line ) )
                {
                    if( line.find( "END OF RMS MAP" ) != std::string::npos )
                    {
                        break;
                    }
                }
            }
        }
    }

    file.close( );

    // Reverse latitudes to ascending order (required by Tudat interpolators)
    std::reverse( data.latitudes.begin( ), data.latitudes.end( ) );
    for( Eigen::MatrixXd& mat: tecMapQueue )
    {
        mat = mat.colwise( ).reverse( ).eval( );
    }
    for( Eigen::MatrixXd& mat: rmsMapQueue )
    {
        mat = mat.colwise( ).reverse( ).eval( );
    }

    // Store to output structure
    if( epochQueue.size( ) != tecMapQueue.size( ) )
    {
        throw std::runtime_error( "IONEX epoch list and TEC map count mismatch." );
    }
    if( loadRmsMaps && !rmsMapQueue.empty( ) && rmsMapQueue.size( ) != tecMapQueue.size( ) )
    {
        throw std::runtime_error( "IONEX TEC and RMS map count mismatch." );
    }

    // Apply microsecond adjustments to avoid boundary overlaps
    if( epochQueue.size( ) >= 1 )
    {
        epochQueue.front( ) += 1.0e-6;
        epochQueue.back( ) -= 1.0e-6;
    }

    for( std::size_t i = 0; i < epochQueue.size( ); ++i )
    {
        data.epochs.push_back( epochQueue[ i ] );
        data.tecMaps[ epochQueue[ i ] ] = tecMapQueue[ i ];
        if( loadRmsMaps && i < rmsMapQueue.size( ) )
        {
            data.rmsMaps[ epochQueue[ i ] ] = rmsMapQueue[ i ];
        }
    }
}

void readIonexFiles( const std::vector< std::string >& filePaths, IonexTecMap& data, const bool loadRmsMaps )
{
    for( const auto& path: filePaths )
    {
        readIonexFile( path, data, loadRmsMaps );
    }

    std::set< double > uniqueHeights;
    for( const auto& it: data.tecMaps )
    {
        uniqueHeights.insert( data.referenceIonosphereHeight_ );
    }

    if( uniqueHeights.size( ) > 1 )
    {
        std::cerr << "Warning: IONEX files use multiple reference heights for the ionospheric shell:\n";
        for( const auto& h: uniqueHeights )
        {
            std::cerr << "    - " << h << " m\n";
        }
    }

    std::sort( data.epochs.begin( ), data.epochs.end( ) );
    for( std::size_t i = 1; i < data.epochs.size( ); ++i )
    {
        double delta = data.epochs[ i ] - data.epochs[ i - 1 ];
        if( delta > 24.0 * 3600.0 )
        {
            std::cerr << "Warning: Gap of " << delta / 3600.0 << " hours between TEC maps at epochs: " << data.epochs[ i - 1 ] << " and "
                      << data.epochs[ i ] << "\n";
        }
    }
}

}  // namespace input_output
}  // namespace tudat
