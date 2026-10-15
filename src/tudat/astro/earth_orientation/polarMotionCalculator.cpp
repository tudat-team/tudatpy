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

#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/earth_orientation/polarMotionCalculator.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"

namespace tudat
{

namespace earth_orientation
{

Eigen::Vector2d getSecularPolePositionInMas( const double time )
{
    double timeInYears = time / physical_constants::JULIAN_YEAR;
    return ( Eigen::Vector2d( ) << 55.0 + 1.677 * timeInYears, 320.5 + 3.460 * timeInYears ).finished( );
}

// Daily values from the Earth orientation data, plus short-period variations.
Eigen::Vector2d PolarMotionCalculator::getPositionOfCipInItrsFromTai( const Time& tai, const Eigen::Vector6d& arguments )
{
    Eigen::Vector2d value = eopData_->evaluate( tai ).polarMotion;
    if( shortPeriodPolarMotionCalculator_ )
    {
        value += shortPeriodPolarMotionCalculator_->getCorrectionsFromFundamentalArgument( arguments );
    }
    return value;
}

}  // namespace earth_orientation

}  // namespace tudat
