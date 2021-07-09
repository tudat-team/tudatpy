/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/estimation_setup/observationOutputSettings.h"
#include "tudat/astro/observation_models/observationViabilityCalculator.h"

namespace tudat
{

namespace simulation_setup
{

std::string getObservationDependentVariableName(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings )
{
    return "Variable name not yet implemented.";
}


bool isObservationDependentVariableVectorial(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings )
{
    bool isVariableVectorial = false;
    switch( variableSettings->variableType_ )
    {
    case station_elevation_angle:
        isVariableVectorial = false;
        break;
    case station_azimuth_angle:
        isVariableVectorial = false;
        break;
    case target_range:
        isVariableVectorial = false;
        break;
    case body_avoidance_angle:
        isVariableVectorial = false;
        break;
    case integration_time:
        isVariableVectorial = false;
        break;
    case station_local_time:
        isVariableVectorial = false;
        break;
    case limb_separation_angle:
        isVariableVectorial = false;
        break;
    default:
        throw std::runtime_error( "Error when checking observation dependent variable. Type " +
                                  getObservationDependentVariableName( variableSettings ) +
                                  " not found when checking fif variable is vectorial." );

    }
    return isVariableVectorial;
}

bool isObservationDependentVariableGroundStationProperty(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings )
{
    bool isGroundStationProperty;
    switch( variableSettings->variableType_ )
    {
    case station_elevation_angle:
        isGroundStationProperty = true;
        break;
    case station_azimuth_angle:
        isGroundStationProperty = true;
        break;
    case target_range:
        isGroundStationProperty = false;
        break;
    case body_avoidance_angle:
        isGroundStationProperty = false;
        break;
    case integration_time:
        isGroundStationProperty = false;
        break;
    case station_local_time:
        isGroundStationProperty = true;
        break;
    case limb_separation_angle:
        isGroundStationProperty = false;
        break;
    default:
        throw std::runtime_error( "Error when checking observation dependent variable. Type " +
                                  getObservationDependentVariableName( variableSettings ) +
                                  " not found when checking for ground station dependency." );

    }
    return isGroundStationProperty;
}

int getObservationDependentVariableSize(
        const std::shared_ptr< ObservationDependentVariableSettings > variableSettings )
{
    int variableSize = 0;
    if( !isObservationDependentVariableVectorial( variableSettings ) )
    {
        variableSize = 1;
    }
    else
    {
        switch( variableSettings->variableType_ )
        {
        default:
            throw std::runtime_error( "Error when checking observation dependent variable. Type " +
                                      getObservationDependentVariableName( variableSettings ) +
                                      " not found when determining parameter size." );

        }
    }
    return variableSize;
}

}

}
