/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readViennaMappingFunctionData.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/basics/utilities.h"

namespace tudat
{

namespace input_output
{

void VMFData::validate( )
{
    if( ( delayData.size( ) != meteoData.size( ) ) && ( meteoData.size( ) > 0 ) )
    {
        throw std::runtime_error( "Error, meteo data is of incompatible size in VMF data" );
    }

    if( ( delayData.size( ) != gradientData.size( ) ) && ( gradientData.size( ) > 0 ) )
    {
        throw std::runtime_error( "Error, gradient data is of incompatible size in VMF data" );
    }

    for( auto it: delayData )
    {
        if( meteoData.size( ) > 0 )
        {
            if( meteoData.count( it.first ) == 0 )
            {
                throw std::runtime_error( "Error, meteo data is at incompatible times in VMF data" );
            }
        }

        if( gradientData.size( ) > 0 )
        {
            if( gradientData.count( it.first ) == 0 )
            {
                throw std::runtime_error( "Error, gradient data is at incompatible times in VMF data" );
            }
        }
    }
}

void VMFData::getFullDataSet( std::map< double, Eigen::VectorXd >& processedTroposphereData,
                              std::map< double, Eigen::VectorXd >& processedMeteoData )
{
    validate( );

    bool hasMeteoData = ( meteoData.size( ) > 0 );
    bool hasGradientData = ( gradientData.size( ) > 0 );
    int singleEntrySize = 4 + ( hasGradientData ? 4 : 0 );

    Eigen::VectorXd currentData = Eigen::VectorXd::Zero( singleEntrySize );
    double mjdAtJ2000 = basic_astrodynamics::getModifiedJulianDayOnJ2000< double >( );

    for( auto it: delayData )
    {
        double currentMjd = it.first;
        double currentTt = ( it.first - mjdAtJ2000 ) * physical_constants::JULIAN_DAY;
        double currentUtc = sofa_interface::convertTTtoUTC< double >( currentTt );

        currentData.setZero( );
        currentData.segment( 0, 4 ) = utilities::convertStlVectorToEigenVector( it.second );
        if( hasGradientData )
        {
            currentData.segment( 4, 4 ) = utilities::convertStlVectorToEigenVector( gradientData.at( it.first ) );
        }
        processedTroposphereData[ currentUtc ] = currentData;

        if( hasMeteoData )
        {
            processedMeteoData[ currentUtc ] = utilities::convertStlVectorToEigenVector( meteoData.at( it.first ) );
        }
    }
}

void readVMFFile( const std::string& fileName,
                  std::map< std::string, VMFData >& vmfData,
                  const bool fileHasMeteo,
                  const bool fileHasGradient )
{
    // Open file
    std::ifstream stream( fileName, std::ios_base::in );
    if( !stream.good( ) )
    {
        throw std::runtime_error( "Error when opening VMF file: " + fileName );
    }

    // Line based parsing
    std::string currentLine;
    std::vector< std::string > currentSplitLine;

    std::vector< double > currentDelayData;
    currentDelayData.resize( 4 );
    std::vector< double > currentMeteoData;
    currentMeteoData.resize( 3 );
    std::vector< double > currentGradientData;
    currentGradientData.resize( 4 );

    while( stream.good( ) && !stream.eof( ) )
    {
        // Get line from stream
        std::getline( stream, currentLine );

        if( currentLine.size( ) > 0 )
        {
            if( currentLine.at( 0 ) != '#' )
            {
                // Trim input string (removes all leading and trailing whitespaces).
                boost::algorithm::trim( currentLine );
                boost::algorithm::split(
                        currentSplitLine, currentLine, boost::algorithm::is_any_of( " \t" ), boost::algorithm::token_compress_on );
                std::string stationName = currentSplitLine.at( 0 );
                if( vmfData.count( stationName ) == 0 )
                {
                    vmfData[ stationName ] = VMFData( );
                }

                double currentMjd = std::stod( currentSplitLine.at( 1 ) );
                int currentIndex = 2;
                int currentStartIndex = currentIndex;

                for( ; currentIndex < currentStartIndex + 4; currentIndex++ )
                {
                    currentDelayData[ currentIndex - currentStartIndex ] = std::stod( currentSplitLine.at( currentIndex ) );
                }
                vmfData[ stationName ].delayData[ currentMjd ] = currentDelayData;

                if( fileHasMeteo )
                {
                    currentStartIndex = currentIndex;
                    for( ; currentIndex < currentStartIndex + 3; currentIndex++ )
                    {
                        currentMeteoData[ currentIndex - currentStartIndex ] = std::stod( currentSplitLine.at( currentIndex ) );
                        switch( currentIndex - currentStartIndex )
                        {
                            case 0:
                                currentMeteoData[ currentIndex - currentStartIndex ] *= 100.0;
                                break;
                            case 1:
                                currentMeteoData[ currentIndex - currentStartIndex ] += 273.15;
                                break;
                            case 2:
                                currentMeteoData[ currentIndex - currentStartIndex ] *= 100.0;
                                break;
                            default:
                                throw std::runtime_error( "Error when reading meteo data from VMF file, index is incompatible" );
                        }
                    }
                    vmfData[ stationName ].meteoData[ currentMjd ] = currentMeteoData;
                }

                if( fileHasGradient )
                {
                    currentStartIndex = currentIndex;
                    for( ; currentIndex < currentStartIndex + 4; currentIndex++ )
                    {
                        currentGradientData[ currentIndex - currentStartIndex ] =
                                std::stod( currentSplitLine.at( currentIndex ) ) / 1000.0;  // convert mm to m
                    }
                    vmfData[ stationName ].gradientData[ currentMjd ] = currentGradientData;
                }
            }
        }
    }

    // Close file
    stream.close( );
}

void readVMFFiles( const std::vector< std::string >& fileName,
                   std::map< std::string, VMFData >& vmfData,
                   const bool fileHasMeteo,
                   const bool fileHasGradient )
{
    for( unsigned int i = 0; i < fileName.size( ); i++ )
    {
        readVMFFile( fileName.at( i ), vmfData, fileHasMeteo, fileHasGradient );
    }
}

}  // namespace input_output

}  // namespace tudat