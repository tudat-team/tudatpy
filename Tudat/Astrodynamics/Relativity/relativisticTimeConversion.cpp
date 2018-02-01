/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/Relativity/relativisticTimeConversion.h"

namespace tudat
{

namespace relativity
{

//! Function to compute proper-time rate w.r.t. coordinate time, minus 1.0, from a speed and scalar potential
double calculateFirstCentralBodyProperTimeRateDifference(
        const double relativeSpeed, const double gravitationalScalarPotential,
        const double equivalencePrincipleLpiViolationParameter )
{
    return ( -( 0.5 * relativeSpeed * relativeSpeed + ( 1.0 + equivalencePrincipleLpiViolationParameter ) *
                gravitationalScalarPotential )  ) * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
}

//! Function to compute proper-time rate w.r.t. coordinate time, minus 1.0, from a speed and scalar potential
double calculateFirstCentralBodyProperTimeRateDifference(
        const Eigen::Vector6d relativeStateVector, const double centralBodyGravitationalParameter,
        const double equivalencePrincipleLpiViolationParameter )
{
    return calculateFirstCentralBodyProperTimeRateDifference(
                relativeStateVector.segment( 3, 3 ).norm( ), centralBodyGravitationalParameter /
                relativeStateVector.segment( 0, 3 ).norm( ), equivalencePrincipleLpiViolationParameter );
}

}

}

