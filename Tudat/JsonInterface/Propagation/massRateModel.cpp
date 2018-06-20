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

#include "Tudat/JsonInterface/Propagation/massRateModel.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `MassRateModelSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< MassRateModelSettings >& massRateModelSettings )
{
    if ( ! massRateModelSettings )
    {
        return;
    }
    using namespace basic_astrodynamics;
    using namespace json_interface;
    using K = Keys::Propagator::MassRateModel;

    const AvailableMassRateModels massRateType = massRateModelSettings->massRateType_;
    jsonObject[ K::type ] = massRateType;

    switch ( massRateType ) {
    case from_thrust_mass_rate_model:
    {
        std::shared_ptr< FromThrustMassModelSettings > fromThrustMassModelSettings =
                std::dynamic_pointer_cast< FromThrustMassModelSettings >( massRateModelSettings );
        assertNonnullptrPointer( fromThrustMassModelSettings );
        jsonObject[ K::useAllThrustModels ] = fromThrustMassModelSettings->useAllThrustModels_;
        assignIfNotEmpty( jsonObject, K::associatedThrustSource,
                          fromThrustMassModelSettings->associatedThrustSource_ );
        return;
    }
    default:
        handleUnimplementedEnumValue( massRateType, massRateTypes, unsupportedMassRateType );
    }
}

//! Create a shared pointer to a `MassRateModelSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< MassRateModelSettings >& massRateModelSettings )
{
    using namespace json_interface;
    using namespace basic_astrodynamics;
    using K = Keys::Propagator::MassRateModel;

    // Get mass-rate type
    const AvailableMassRateModels massRateType = getValue< AvailableMassRateModels >( jsonObject, K::type );

    switch ( massRateType ) {
    case from_thrust_mass_rate_model:
    {
        std::shared_ptr< FromThrustMassModelSettings > fromThrustMassModelSettings =
                std::make_shared< FromThrustMassModelSettings >( );
        updateFromJSONIfDefined( fromThrustMassModelSettings->useAllThrustModels_, jsonObject, K::useAllThrustModels );
        if ( ! fromThrustMassModelSettings->useAllThrustModels_ )
        {
            updateFromJSON( fromThrustMassModelSettings->associatedThrustSource_,
                            jsonObject, K::associatedThrustSource );
        }
        massRateModelSettings = fromThrustMassModelSettings;
        return;
    }
    default:
        handleUnimplementedEnumValue( massRateType, massRateTypes, unsupportedMassRateType );
    }
}

} // namespace simulation_setup

} // namespace tudat
