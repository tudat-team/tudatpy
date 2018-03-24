/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <stdexcept>

#include "Tudat/Astrodynamics/Aerodynamics/exponentialAtmosphere.h"

namespace tudat
{
namespace aerodynamics
{

ExponentialAtmosphere::ExponentialAtmosphere(
        const BodiesWithPredefinedExponentialAtmospheres bodyWithPredefinedExponentialAtmosphere )
{
    switch( bodyWithPredefinedExponentialAtmosphere )
    {
    case earth:
        // Set local variables for Earth exponential atmosphere. Based on  lecture notes
        // Rocket Motion by Prof. Ir. B.A.C. Ambrosius, November 2009.

        // Set scale height.
        scaleHeight_ = 7.200e3;

        //Set density at zero altitude.
        densityAtZeroAltitude_ = 1.225;

        //Set atmosphere temperature.
        constantTemperature_ = 246.0;

        // Set specific gas constant.
        specificGasConstant_ = physical_constants::SPECIFIC_GAS_CONSTANT_AIR;

        ratioOfSpecificHeats_ = 1.4;

        break;

    default:
        throw std::runtime_error(
                    "Error when making exponential atmosphere, predefined atmosphere" +
                    std::to_string(
                        bodyWithPredefinedExponentialAtmosphere ) + "not recognized." );
    }
}

} // namespace aerodynamics

} // namespace tudat
