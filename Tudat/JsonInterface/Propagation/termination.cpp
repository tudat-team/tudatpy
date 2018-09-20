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

#include "Tudat/JsonInterface/Propagation/termination.h"

#include "Tudat/JsonInterface/Propagation/variable.h"

namespace tudat
{

namespace propagators
{

// PropagationHybridTerminationSettings

//! Create a `json` object from a shared pointer to a `PropagationHybridTerminationSettings` object.
void to_json( nlohmann::json& jsonObject,
              const std::shared_ptr< PropagationHybridTerminationSettings >& hybridTerminationSettings )
{
    if ( ! hybridTerminationSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Termination;

    jsonObject[ hybridTerminationSettings->fulfillSingleCondition_ ? K::anyOf : K::allOf ] =
            hybridTerminationSettings->terminationSettings_;
}

//! Create a shared pointer to a `PropagationHybridTerminationSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject,
                std::shared_ptr< PropagationHybridTerminationSettings >& hybridTerminationSettings )
{
    using namespace json_interface;
    using K = Keys::Termination;

    // Not-hybrid
    if ( ! isDefined( jsonObject, K::allOf ) && ! isDefined( jsonObject, K::anyOf ) )
    {
        hybridTerminationSettings = std::make_shared< PropagationHybridTerminationSettings >(
                    std::vector< std::shared_ptr< PropagationTerminationSettings > >(
        { getAs< std::shared_ptr< PropagationTerminationSettings > >( jsonObject ) } ), true );
    }
    else
    {
        const bool meetAnyCondition = isDefined( jsonObject, K::anyOf );
        hybridTerminationSettings = std::make_shared< PropagationHybridTerminationSettings >(
                    getValue< std::vector< std::shared_ptr< PropagationTerminationSettings > > >(
                        jsonObject, meetAnyCondition ? K::anyOf : K::allOf ), meetAnyCondition );
    }
}


// PropagationTerminationSettings

//! Create a `json` object from a shared pointer to a `PropagationTerminationSettings` object.
void to_json( nlohmann::json& jsonObject,
              const std::shared_ptr< PropagationTerminationSettings >& terminationSettings )
{
    if ( ! terminationSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Termination;

    // hybrid
    std::shared_ptr< PropagationHybridTerminationSettings > hybridTerminationSettings =
            std::dynamic_pointer_cast< PropagationHybridTerminationSettings >( terminationSettings );
    if ( hybridTerminationSettings )
    {
        jsonObject = hybridTerminationSettings;
        return;
    }

    // time
    std::shared_ptr< PropagationTimeTerminationSettings > timeTerminationSettings =
            std::dynamic_pointer_cast< PropagationTimeTerminationSettings >( terminationSettings );
    if ( timeTerminationSettings )
    {
        jsonObject[ K::variable ] = std::make_shared< VariableSettings >( independentVariable );
        jsonObject[ K::upperLimit ] = timeTerminationSettings->terminationTime_;
        return;
    }

    // cpu time
    std::shared_ptr< PropagationCPUTimeTerminationSettings > cpuTimeTerminationSettings =
            std::dynamic_pointer_cast< PropagationCPUTimeTerminationSettings >( terminationSettings );
    if ( cpuTimeTerminationSettings )
    {
        jsonObject[ K::variable ] = std::make_shared< VariableSettings >( cpuTimeVariable );
        jsonObject[ K::upperLimit ] = cpuTimeTerminationSettings->cpuTerminationTime_;
        return;
    }

    // dependent variable
    std::shared_ptr< PropagationDependentVariableTerminationSettings > dependentVariableTerminationSettings =
            std::dynamic_pointer_cast< PropagationDependentVariableTerminationSettings >(
                terminationSettings );
    assertNonnullptrPointer( dependentVariableTerminationSettings );
    jsonObject[ K::variable ] = dependentVariableTerminationSettings->dependentVariableSettings_;
    jsonObject[ dependentVariableTerminationSettings->useAsLowerLimit_ ? K::lowerLimit : K::upperLimit ] =
            dependentVariableTerminationSettings->limitValue_;
}

//! Create a shared pointer to a `PropagationTerminationSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject,
                std::shared_ptr< PropagationTerminationSettings >& terminationSettings )
{
    using namespace json_interface;
    using K = Keys::Termination;

    // Hybrid
    if ( isDefined( jsonObject, K::allOf ) || isDefined( jsonObject, K::anyOf ) )
    {
        terminationSettings = getAs< std::shared_ptr< PropagationHybridTerminationSettings > >( jsonObject );
        return;
    }

    const std::shared_ptr< VariableSettings > variable =
            getValue< std::shared_ptr< VariableSettings > >( jsonObject, K::variable );
    switch ( variable->variableType_ )
    {
    case independentVariable:
    {
        terminationSettings = std::make_shared< PropagationTimeTerminationSettings >(
                    getValue< double >( jsonObject, { K::upperLimit, SpecialKeys::root / Keys::finalEpoch } ) );
        return;
    }
    case cpuTimeVariable:
    {
        terminationSettings = std::make_shared< PropagationCPUTimeTerminationSettings >(
                    getValue< double >( jsonObject, K::upperLimit ) );
        return;
    }
    case dependentVariable:
    {
        // If both lower limit and upper limit, create hybrid satistying any of the two conditions
        if ( isDefined( jsonObject, K::lowerLimit ) && isDefined( jsonObject, K::upperLimit ) )
        {
            // Lower limit
            nlohmann::json lowerLimitObject = jsonObject;
            lowerLimitObject.erase( K::upperLimit );
            std::shared_ptr< PropagationTerminationSettings > lowerLimitCondition =
                    getAs< std::shared_ptr< PropagationTerminationSettings > >( lowerLimitObject );

            // Upper limit
            nlohmann::json upperLimitObject = jsonObject;
            upperLimitObject.erase( K::lowerLimit );
            std::shared_ptr< PropagationTerminationSettings > upperLimitCondition =
                    getAs< std::shared_ptr< PropagationTerminationSettings > >( upperLimitObject );

            // Combine
            terminationSettings = std::make_shared< PropagationHybridTerminationSettings >(
                        std::vector< std::shared_ptr< PropagationTerminationSettings > >(
            { lowerLimitObject, upperLimitObject } ), true );
            return;
        }

        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVar =
                std::dynamic_pointer_cast< SingleDependentVariableSaveSettings >( variable );
        assertNonnullptrPointer( dependentVar );

        // Limit value
        double limitValue;
        bool useLowerLimit;
        if ( isDefined( jsonObject, K::lowerLimit ) )
        {
            limitValue = getValue< double >( jsonObject, K::lowerLimit );
            useLowerLimit = true;
        }
        else
        {
            limitValue = getValue< double >( jsonObject, K::upperLimit );
            useLowerLimit = false;
        }

        terminationSettings = std::make_shared< PropagationDependentVariableTerminationSettings >(
                    dependentVar, limitValue, useLowerLimit );
        return;
    }
    default:
        throw std::runtime_error( "Termination settings cannot be defined as a function of variable: " +
                                  getVariableId( variable ) );
    }

}

} // namespace propagators

} // namespace tudat
