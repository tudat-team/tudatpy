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

#include "Tudat/JsonInterface/Propagation/torque.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `TorqueSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< TorqueSettings >& torqueSettings )
{
    if ( ! torqueSettings )
    {
        return;
    }
    using namespace basic_astrodynamics;
    using namespace json_interface;
    using K = Keys::Propagator::Torque;

    const AvailableTorque torqueType = torqueSettings->torqueType_;
    jsonObject[ K::type ] = torqueType;

    switch ( torqueType ) {
    case underfined_torque:
    case second_order_gravitational_torque:
    case aerodynamic_torque:
        return;
    default:
        handleUnimplementedEnumValue( torqueType, torqueTypes, unsupportedTorqueTypes );
    }
}

//! Create a shared pointer to a `TorqueSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< TorqueSettings >& torqueSettings )
{
    using namespace json_interface;
    using namespace basic_astrodynamics;
    using K = Keys::Propagator::Torque;

    // Get acceleration type
    const AvailableTorque torqueType = getValue< AvailableTorque >( jsonObject, K::type );

    switch ( torqueType ) {
    case underfined_torque:
    case second_order_gravitational_torque:
    case aerodynamic_torque:
    {
        torqueSettings = std::make_shared< TorqueSettings >( torqueType );
        return;
    }
    default:
        handleUnimplementedEnumValue( torqueType, torqueTypes, unsupportedTorqueTypes );
    }
}

} // namespace simulation_setup

} // namespace tudat
