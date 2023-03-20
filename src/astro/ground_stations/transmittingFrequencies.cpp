/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/ground_stations/transmittingFrequencies.h"

namespace tudat
{

namespace ground_stations
{

template< >
double StationFrequencyInterpolator::getTemplatedCurrentFrequency( const double& lookupTime )
{
    return getCurrentFrequency( lookupTime );
}

template< >
double StationFrequencyInterpolator::getTemplatedCurrentFrequency( const Time& lookupTime )
{
    return getCurrentFrequency( lookupTime );
}

template< >
long double StationFrequencyInterpolator::getTemplatedCurrentFrequency( const double& lookupTime )
{
    return getCurrentLongFrequency( lookupTime );
}

template< >
long double StationFrequencyInterpolator::getTemplatedCurrentFrequency( const Time& lookupTime )
{
    return getCurrentLongFrequency( lookupTime );
}

template< >
double StationFrequencyInterpolator::getTemplatedFrequencyIntegral(
        const double& quadratureStartTime, const double& quadratureEndTime )
{
    return getFrequencyIntegral( quadratureStartTime, quadratureEndTime );
}

template< >
double StationFrequencyInterpolator::getTemplatedFrequencyIntegral(
        const Time& quadratureStartTime, const Time& quadratureEndTime )
{
    return getFrequencyIntegral( quadratureStartTime, quadratureEndTime );
}

template< >
long double StationFrequencyInterpolator::getTemplatedFrequencyIntegral(
        const double& quadratureStartTime, const double& quadratureEndTime )
{
    return getLongFrequencyIntegral( quadratureStartTime, quadratureEndTime );
}

template< >
long double StationFrequencyInterpolator::getTemplatedFrequencyIntegral(
        const Time& quadratureStartTime, const Time& quadratureEndTime )
{
    return getLongFrequencyIntegral( quadratureStartTime, quadratureEndTime );
}



} // namespace ground_stations

} // namespace tudat
