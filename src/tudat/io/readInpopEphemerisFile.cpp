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

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <iomanip>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/math/interpolators/chebyshevApproximation.h"
#include "tudat/io/readInpopEphemerisFile.h"

namespace tudat
{

namespace input_output
{
std::map< double, std::vector< Eigen::Vector6d > > readInpopStateEphemeris(
        const std::string& positionFileName,
        const std::string& velocityFileName,
        const double referenceJulianDay )
{
    std::map< double, std::vector< Eigen::Vector6d > > chebyshevCoefficients;

    std::ifstream positionStream( positionFileName );
    std::ifstream velocityStream( velocityFileName );

    if ( !positionStream.is_open( ) )
    {
        throw std::runtime_error( "Could not open position file: " + positionFileName );
    }
    if ( !velocityStream.is_open( ) )
    {
        throw std::runtime_error( "Could not open velocity file: " + velocityFileName );
    }

    auto parseHeader = []( std::ifstream& stream ) -> std::vector<std::string>
    {
        std::string line;
        std::getline( stream, line );  // discard version line
        std::getline( stream, line );  // header line
        std::istringstream iss( line );
        return { std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>() };
    };

    auto splitLine = []( const std::string& line ) -> std::vector<std::string>
    {
        std::istringstream iss( line );
        return { std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>() };
    };

    std::vector<std::string> positionHeader = parseHeader( positionStream );
    int order = std::stoi( positionHeader.at(6) );
    int span = std::stoi( positionHeader.at(7) );
    int numberOfSteps = std::stoi( positionHeader.at(8) );

    double distanceMultiplier = ( positionHeader.at(4) == "km" ) ? 1000.0 :
        throw std::runtime_error("Unknown position unit: " + positionHeader.at(4));

    std::vector<std::string> velocityHeader = parseHeader( velocityStream );
    if ( std::stoi( velocityHeader.at(6) ) != order ||
         std::stoi( velocityHeader.at(7) ) != span ||
         std::stoi( velocityHeader.at(8) ) != numberOfSteps )
    {
        throw std::runtime_error( "INPOP state file headers inconsistent between position and velocity files." );
    }

    double velocityMultiplier = ( velocityHeader.at(4) == "km/day" ) ?
        1000.0 / physical_constants::JULIAN_DAY :
        throw std::runtime_error( "Unknown velocity unit: " + velocityHeader.at(4) );

    std::string line;
    while ( true )
    {
        std::vector< Eigen::Vector6d > currentCoefficients( order );
        double currentStartTime = 0.0;

        for ( int i = 0; i < 3; ++i )
        {
            if ( !std::getline( positionStream, line ) )
                return chebyshevCoefficients;
            std::replace( line.begin(), line.end(), 'D', 'E' );
            std::replace( line.begin(), line.end(), 'd', 'e' );
            std::vector<std::string> tokens = splitLine( line );

            if ( tokens.size( ) > 1 )
            {
                for ( int j = 0; j < order; ++j )
                {
                    currentCoefficients[ j ]( i ) = std::stod( tokens.at( j + 2 ) ) * distanceMultiplier;
                }
                currentStartTime = std::stod( tokens.at( 0 ) );
            }
        }

        for ( int i = 0; i < 3; ++i )
        {
            if ( !std::getline( velocityStream, line ) )
                return chebyshevCoefficients;
            std::replace( line.begin(), line.end(), 'D', 'E' );
            std::replace( line.begin(), line.end(), 'd', 'e' );
            std::vector<std::string> tokens = splitLine( line );

            if ( tokens.size( ) > 1 )
            {
                for ( int j = 0; j < order; ++j )
                {
                    currentCoefficients[ j ]( i + 3 ) = std::stod( tokens.at( j + 2 ) ) * velocityMultiplier;
                }
            }
        }

        // Store coefficients
        const double epoch = ( currentStartTime - referenceJulianDay ) * physical_constants::JULIAN_DAY;
        chebyshevCoefficients[ epoch ] = currentCoefficients;
    }
    return chebyshevCoefficients;
}


std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > createInpopTimeEphemerisInterpolator(
        const std::string& fileName,
        const double referenceJulianDay )
{
    auto timeDifferenceMap = readInpopTimeEphemeris< double >( fileName, referenceJulianDay );
    return std::make_shared< interpolators::ChebyshevInterpolator< double, double > >( timeDifferenceMap );
}

std::shared_ptr< interpolators::OneDimensionalInterpolator< double, long double > > createLongInpopTimeEphemerisInterpolator(
        const std::string& fileName,
        const double referenceJulianDay )
{
    auto timeDifferenceMap = readInpopTimeEphemeris< long double >( fileName, referenceJulianDay );
    return std::make_shared< interpolators::ChebyshevInterpolator< double, long double > >( timeDifferenceMap );
}

std::shared_ptr< ephemerides::Ephemeris > createInpopEphemerisFromFiles(
        const std::string& positionFileName,
        const std::string& velocityFileName,
        const double referenceJulianDay,
        const bool useGeocentricReference )
{
    auto stateMap = readInpopStateEphemeris( positionFileName, velocityFileName, referenceJulianDay );
    auto stateInterpolator = std::make_shared< interpolators::ChebyshevInterpolator< double, Eigen::Vector6d > >( stateMap );

    std::string referencePoint = useGeocentricReference ? "Earth" : "SSB";

    return std::make_shared< ephemerides::TabulatedCartesianEphemeris< double, double > >( stateInterpolator, referencePoint );
}

} // namespace input_output

} // namespace tudat
