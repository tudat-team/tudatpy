/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"

namespace tudat
{

namespace observation_models
{

bool requiresMultiLegIterations( const LightTimeCorrectionType& lightTimeCorrectionType )
{
    bool requiresMultiLegIterations = false;
    switch( lightTimeCorrectionType )
    {
    case first_order_relativistic:
    case tabulated_tropospheric:
    case saastamoinen_tropospheric:
        requiresMultiLegIterations = false;
        break;
    case tabulated_ionospheric:
    case mapped_vtec_ionospheric:
        requiresMultiLegIterations = true;
        break;
    default:
        throw std::runtime_error(
                "Error when getting whether light time corrections require multi-leg iterations: could not find light"
                "time correction type " +  std::to_string( lightTimeCorrectionType ) + "." );
    }

    return requiresMultiLegIterations;
}

std::string getLightTimeCorrectionName( const LightTimeCorrectionType& lightTimeCorrectionType )
{
    std::string name;

    switch( lightTimeCorrectionType )
    {
    case first_order_relativistic:
        name = "first order relativistic";
        break;
    case tabulated_tropospheric:
        name = "tabulated tropospheric";
        break;
    case tabulated_ionospheric:
        name = "tabulated ionospheric";
        break;
    case saastamoinen_tropospheric:
        name = "Saastamoinen tropospheric";
        break;
    default:
        throw std::runtime_error(
                "Error when getting light time correction name: could not find light"
                "time correction type " +  std::to_string( lightTimeCorrectionType ) + "." );
    }

    return name;
}


} // namespace observation_models

} // namespace tudat
