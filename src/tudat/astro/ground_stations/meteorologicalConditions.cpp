#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/math/basic/mathematicalConstants.h"

/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/ground_stations/meteorologicalConditions.h"
namespace tudat
{

namespace ground_stations
{

double calculateBeanAndDuttonWaterVaporPartialPressure( double relativeHumidity, double temperature )
{
    if( relativeHumidity < 0 || relativeHumidity > 1 )
    {
        throw std::runtime_error( "Error when computing partial vapor pressure: invalid relative humidity value (" +
                                  std::to_string( relativeHumidity ) + "), should be in [0,1]." );
    }

    // Estefan and Sovers (1994), eq. 16
    // 1e2 factor is conversion from mbar to Pa
    return 6.11 * relativeHumidity * std::pow( 10.0, 7.5 * ( ( temperature - 273.15 ) / ( temperature - 35.85 ) ) ) * 1e2;
}

std::function< double( const double ) > getBeanAndDuttonWaterVaporPartialPressureFunction(
        std::function< double( const double time ) > relativeHumidity,
        std::function< double( const double time ) > temperature )
{
    return [ = ]( double time ) {
        return calculateBeanAndDuttonWaterVaporPartialPressure( relativeHumidity( time ), temperature( time ) );
    };
}

// Function to compute saturation vapor pressure (Pa) using Tetens' formula
double computeSaturationWaterVaporPressure( const double temperature )
{
    double temperatureInCelsius = temperature - 273.15;  // Convert Kelvin to Celsius
    return 611.2 * exp( ( 17.62 * temperatureInCelsius ) / ( temperatureInCelsius + 243.12 ) );
}

// Magnus formula
double computeDewPoint( const double relativeHumidity, const double temperature )
{
    double temperatureInCelsius = temperature - 273.15;  // Convert Kelvin to Celsius
    double gamma = std::log( relativeHumidity ) + 17.625 * temperatureInCelsius / ( 243.04 + temperatureInCelsius );
    return 243.04 * gamma / ( 17.625 - gamma ) + 273.15;
}

void ContinuousInterpolatedMeteoData::updateData( const double currentUtc )
{
    if( !( currentUtc == currentUtc_ ) )
    {
        currentUtc_ = currentUtc;
        try
        {
            currentData_ = meteoDataInterpolator_->interpolate( currentUtc_ );
        }
        catch( std::runtime_error& caughtException )
        {
            throw std::runtime_error( "Error in continuous meteo data interpolator.\nOriginal error: " +
                                      std::string( caughtException.what( ) ) );
        }
    }
}

void PiecewiseInterpolatedMeteoData::updateData( const double currentUtc )
{
    if( !( currentUtc == currentUtc_ ) )
    {
        currentUtc_ = currentUtc;
        int nearestIndex = lookUpScheme_->findNearestLowerNeighbour( currentUtc_ );

        // If currentUtc falls in a gap between two segment's data, extrapolate from the segment
        // that is closer in time, rather than always extrapolating forward from the lower segment.
        if( currentUtc_ > endTimes_.at( nearestIndex ) && static_cast< unsigned int >( nearestIndex ) + 1 < startTimes_.size( ) )
        {
            double distanceToLowerSegment = currentUtc_ - endTimes_.at( nearestIndex );
            double distanceToUpperSegment = startTimes_.at( nearestIndex + 1 ) - currentUtc_;
            if( distanceToUpperSegment < distanceToLowerSegment )
            {
                nearestIndex += 1;
            }
        }

        currentInterpolator_ = nearestIndex;
        try
        {
            currentData_ = meteoDataInterpolators_.at( currentInterpolator_ )->interpolate( currentUtc_ );
        }
        catch( std::runtime_error& caughtException )
        {
            throw std::runtime_error( "Error in piecewise constant meteo data interpolator.\nOriginal error: " +
                                      std::string( caughtException.what( ) ) );
        }
    }
}
}  // namespace ground_stations

}  // namespace tudat
