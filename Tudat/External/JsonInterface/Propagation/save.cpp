/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "save.h"

#include "variable.h"

namespace tudat
{

namespace propagators
{

//! Create a `json` object from a shared pointer to a `DependentVariableSaveSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< DependentVariableSaveSettings >& saveSettings )
{
    if ( ! saveSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Propagator::Save;

    jsonObject[ K::variables ] = saveSettings->dependentVariables_;
    jsonObject[ K::printVariablesNames ] = saveSettings->printDependentVariableTypes_;
}

//! Create a shared pointer to a `DependentVariableSaveSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< DependentVariableSaveSettings >& saveSettings )
{
    using namespace json_interface;
    using K = Keys::Propagator::Save;

    saveSettings = boost::make_shared< DependentVariableSaveSettings >(
        getValue< std::vector< boost::shared_ptr< SingleDependentVariableSaveSettings > > >(
            jsonObject, K::variables ) );
    updateFromJSONIfDefined( saveSettings->printDependentVariableTypes_, jsonObject, K::printVariablesNames );
}

} // namespace propagators

} // namespace tudat
