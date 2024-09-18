/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/lambda/lambda.hpp>

#include "tudat/math/interpolators/lagrangeInterpolator.h"
#include "tudat/simulation/environment_setup/createEphemeris.h"

namespace tudat
{

namespace simulation_setup
{

using namespace ephemerides;

//! Function that retrieves the time interval at which an ephemeris can be safely interrogated
std::pair< double, double > getSafeInterpolationInterval( const std::shared_ptr< ephemerides::Ephemeris > ephemerisModel )
{
    // Make default output pair
    std::pair< double, double > safeInterval = std::make_pair(
                std::numeric_limits< double >::lowest( ),  std::numeric_limits< double >::max( ) );

    // Check if model is tabulated, and retrieve safe interval from model
    if( isTabulatedEphemeris( ephemerisModel ) )
    {
        safeInterval = getTabulatedEphemerisSafeInterval( ephemerisModel );
    }
    // Check if model is multi-arc, and retrieve safe intervals from first and last arc.
    else if( std::dynamic_pointer_cast< ephemerides::MultiArcEphemeris >( ephemerisModel ) != nullptr )
    {
        std::shared_ptr< ephemerides::MultiArcEphemeris > multiArcEphemerisModel  =
                std::dynamic_pointer_cast< ephemerides::MultiArcEphemeris >( ephemerisModel );
        safeInterval.first = getSafeInterpolationInterval( multiArcEphemerisModel->getSingleArcEphemerides( ).at( 0 ) ).first;
        safeInterval.second = getSafeInterpolationInterval(
                    multiArcEphemerisModel->getSingleArcEphemerides( ).at(
                        multiArcEphemerisModel->getSingleArcEphemerides( ).size( ) - 1 ) ).second;
    }
    return safeInterval;
}


} // namespace simulation_setup

} // namespace tudat
