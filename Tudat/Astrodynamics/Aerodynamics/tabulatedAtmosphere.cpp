/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <boost/make_shared.hpp>

#include "Tudat/InputOutput/matrixTextFileReader.h"

#include "Tudat/Astrodynamics/Aerodynamics/tabulatedAtmosphere.h"

namespace tudat
{
namespace aerodynamics
{

//! Initialize atmosphere table reader.
void TabulatedAtmosphere::initialize( const std::string& atmosphereTableFile )
{
    // Locally store the atmosphere table file name.
    atmosphereTableFile_ = atmosphereTableFile;

    Eigen::MatrixXd containerOfAtmosphereTableFileData
            = input_output::readMatrixFromFile( atmosphereTableFile_, " \t", "%" );

    // Check whether data is present in the file.
    if ( containerOfAtmosphereTableFileData.rows( ) < 1
         || containerOfAtmosphereTableFileData.cols( ) < 1 )
    {
        std::string errorMessage = "The atmosphere table file " + atmosphereTableFile_ + " is empty";
        throw std::runtime_error( errorMessage );
    }

    // Initialize vectors.
    altitudeData_.resize( containerOfAtmosphereTableFileData.rows( ) );
    densityData_.resize( containerOfAtmosphereTableFileData.rows( ) );
    pressureData_.resize( containerOfAtmosphereTableFileData.rows( ) );
    temperatureData_.resize( containerOfAtmosphereTableFileData.rows( ) );

    // Loop through all the strings stored in the container and store the data
    // in the right Eigen::VectorXd.
    for ( int i = 0; i < containerOfAtmosphereTableFileData.rows( ); i++  )
    {
        altitudeData_[ i ] = containerOfAtmosphereTableFileData( i, 0 );
        densityData_[ i ] = containerOfAtmosphereTableFileData( i, 1 );
        pressureData_[ i ] = containerOfAtmosphereTableFileData( i, 2 );
        temperatureData_[ i ] = containerOfAtmosphereTableFileData( i, 3 );
    }

    using namespace interpolators;

    cubicSplineInterpolationForDensity_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, densityData_ );
    cubicSplineInterpolationForPressure_
            = boost::make_shared< CubicSplineInterpolatorDouble >( altitudeData_, pressureData_ );
    cubicSplineInterpolationForTemperature_
            = boost::make_shared< CubicSplineInterpolatorDouble >(
                altitudeData_, temperatureData_ );
}

} // namespace aerodynamics
} // namespace tudat
