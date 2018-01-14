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

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/EarthOrientation/polarMotionCalculator.h"
#include "Tudat/External/SofaInterface/sofaTimeConversions.h"

namespace tudat
{

namespace earth_orientation
{


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
    poleOffsetsInItrs += shortPeriodPolarMotionCalculator_->getCorrections(
                fundamentalArguments );

    return poleOffsetsInItrs;
}

}

}
