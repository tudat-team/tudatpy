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


//! Calculate the position of the Celestial Intermediate Pole in the ITRS
Eigen::Vector2d PolarMotionCalculator::getPositionOfCipInItrs(
        const double ttSinceEpoch,
        const double utcSinceEpoch )
{
    // Initialize offset to zero.
    Eigen::Vector2d poleOffsetsInItrs = Eigen::Vector2d::Zero( );

    // Add interpolated measured offsets.
    poleOffsetsInItrs += dailyIersValueInterpolator_->interpolate( utcSinceEpoch );

    // Add short period motion
    poleOffsetsInItrs += shortPeriodPolarMotionCalculator_->getCorrections(
                ttSinceEpoch );

    return poleOffsetsInItrs;
}

//! Calculate the position of the Celestial Intermediate Pole in the ITRS
Eigen::Vector2d PolarMotionCalculator::getPositionOfCipInItrs(
        Eigen::Vector6d fundamentalArguments,
        const double utcSinceEpoch )
{
    // Initialize offset to zero.
    Eigen::Vector2d poleOffsetsInItrs = Eigen::Vector2d::Zero( );

    // Add interpolated measured offsets.
    poleOffsetsInItrs += dailyIersValueInterpolator_->interpolate( utcSinceEpoch );

    // Add short period motion
    poleOffsetsInItrs += shortPeriodPolarMotionCalculator_->getCorrectionsFromFundamentalArgument(
                fundamentalArguments );

    return poleOffsetsInItrs;
}

}

}
