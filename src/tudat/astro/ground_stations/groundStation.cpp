/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/ground_stations/groundStation.h"
#include "tudat/astro/system_models/vehicleSystems.h"

namespace tudat
{

namespace ground_stations
{

std::shared_ptr< system_models::VehicleSystems > GroundStation::getVehicleSystems( )
{
    if( vehicleSystems_ == nullptr )
    {
        vehicleSystems_ = std::make_shared< system_models::VehicleSystems >( );
    }
    return vehicleSystems_;
}

}  // namespace ground_stations

}  // namespace tudat
