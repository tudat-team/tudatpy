/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include "tudat/astro/aerodynamics/mcdAtmosphereModel.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/basics/utilities.h"
#include "tudat/io/basicInputOutput.h"
#include <cmath>
#include <stdexcept>

namespace tudat
{
namespace aerodynamics
{

// Get density
double McdAtmosphereModel::getDensity( double altitude, double longitude, double latitude, double time)
{
    return  marsClimateDatabaseClimateModel_->getCache( {longitude, latitude, time} )->density_;
}

// Get pressure
double McdAtmosphereModel::getPressure( double altitude, double longitude, double latitude, double time)
{
    return marsClimateDatabaseClimateModel_->getCache( {longitude, latitude, time} )->pressure_;
}

// Get temperature
double McdAtmosphereModel::getTemperature( double altitude, double longitude, double latitude, double time)
{
    return marsClimateDatabaseClimateModel_->getCache( {longitude, latitude, time} )->temperature_;
}

double McdAtmosphereModel::getZonalWind( double altitude, double longitude, double latitude, double time) const
{
    return marsClimateDatabaseClimateModel_->getCache( {longitude, latitude, time} )->zonalWind_;
}

double McdAtmosphereModel::getMeridionalWind( double altitude, double longitude, double latitude, double time) const
{
    return marsClimateDatabaseClimateModel_->getCache( {longitude, latitude, time} )->meridionalWind_;
}

// Get speed of sound
double McdAtmosphereModel::getSpeedOfSound( double altitude, double longitude, double latitude, double time)
{
    // Get gamma and R from extra variables
    double gamma = marsClimateDatabaseClimateModel_->getExtraVariable( 
        mcd_interface::ExtVar::ratio_of_specific_heats, {longitude, latitude, time} ); // extvar(60): gamma
    double R = marsClimateDatabaseClimateModel_->getExtraVariable( 
        mcd_interface::ExtVar::reduced_molecular_gas_constant, {longitude, latitude, time} ); // extvar(61): R (J/kg/K)

    if( gamma > 0.0 && R > 0.0 )
    {
        return std::sqrt( gamma * R * getTemperature( altitude, longitude, latitude, time ) );
    }
    else
    {
        // Fallback to Mars defaults
        const double defaultGamma = 1.3;
        const double defaultR = 192.0;  // J/kg/K for Mars CO2 atmosphere
        return std::sqrt( defaultGamma * defaultR * getTemperature( altitude, longitude, latitude, time ) );
    }
}

}  // namespace aerodynamics
}  // namespace tudat