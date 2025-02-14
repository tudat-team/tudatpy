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


// Function to compute saturation vapor pressure (Pa) using Tetens' formula
double computeSaturationWaterVaporPressure( const double temperature)
{
    double temperatureInCelsius = temperature - 273.15; // Convert Kelvin to Celsius
    return 611.2 * exp((17.62 * temperatureInCelsius) / (temperatureInCelsius + 243.12));
}

double computeDewPoint( const double waterVaporPartialPressure) {
    double dewTemperatureInCelsius = (243.12 * std::log(waterVaporPartialPressure / 6.112)) /
                                     (17.62 - std::log(waterVaporPartialPressure / 6.112));
    return dewTemperatureInCelsius + 273.15; // Convert Celsius to Kelvin
}


}  // namespace ground_stations

}  // namespace tudat
