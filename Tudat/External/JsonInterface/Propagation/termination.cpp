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

#include "termination.h"

#include "variable.h"

namespace tudat
{

namespace propagators
{

// PropagationHybridTerminationSettings

//! Create a `json` object from a shared pointer to a `PropagationHybridTerminationSettings` object.
void to_json( json& jsonObject,
              const boost::shared_ptr< PropagationHybridTerminationSettings >& hybridTerminationSettings )
{
    if ( ! hybridTerminationSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Termination;

    jsonObject[ hybridTerminationSettings->fulFillSingleCondition_ ? K::anyOf : K::allOf ] =
            hybridTerminationSettings->terminationSettings_;
}

//! Create a shared pointer to a `PropagationHybridTerminationSettings` object from a `json` object.
void from_json( const json& jsonObject,
                boost::shared_ptr< PropagationHybridTerminationSettings >& hybridTerminationSettings )
{
    using namespace json_interface;
    using K = Keys::Termination;

    // Not-hybrid
    if ( ! defined( jsonObject, K::allOf ) || ! defined( jsonObject, K::anyOf ) )
    {
        hybridTerminationSettings = boost::make_shared< PropagationHybridTerminationSettings >(
                    std::vector< boost::shared_ptr< PropagationTerminationSettings > >(
        { getAs< boost::shared_ptr< PropagationTerminationSettings > >( jsonObject ) } ), true );
    }
    else
    {
        const bool satisfyAny = defined( jsonObject, K::anyOf );
        hybridTerminationSettings = boost::make_shared< PropagationHybridTerminationSettings >(
                    getValue< std::vector< boost::shared_ptr< PropagationTerminationSettings > > >(
                        jsonObject, satisfyAny ? K::anyOf : K::allOf ), satisfyAny );
    }
}


// PropagationTerminationSettings

//! Create a `json` object from a shared pointer to a `PropagationTerminationSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< PropagationTerminationSettings >& terminationSettings )
{
    if ( ! terminationSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Termination;

    boost::shared_ptr< PropagationTimeTerminationSettings > timeTerminationSettings =
            boost::dynamic_pointer_cast< PropagationTimeTerminationSettings >( terminationSettings );
    if ( timeTerminationSettings )  // time
    {
        jsonObject[ K::variable ] = boost::make_shared< VariableSettings >( independentVariable );
        jsonObject[ K::upperLimit ] = timeTerminationSettings->terminationTime_;
    }
    else  // dependent variable
    {
        boost::shared_ptr< PropagationDependentVariableTerminationSettings > dependentVariableTerminationSettings =
                boost::dynamic_pointer_cast< PropagationDependentVariableTerminationSettings >(
                    terminationSettings );
        enforceNonNullPointer( dependentVariableTerminationSettings );
        jsonObject[ K::variable ] = dependentVariableTerminationSettings->dependentVariableSettings_;
        jsonObject[ dependentVariableTerminationSettings->useAsLowerLimit_ ? K::lowerLimit : K::upperLimit ] =
                dependentVariableTerminationSettings->limitValue_;
    }
}

//! Create a shared pointer to a `PropagationTerminationSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< PropagationTerminationSettings >& terminationSettings )
{
    using namespace json_interface;
    using K = Keys::Termination;

    // Hybrid
    if ( defined( jsonObject, K::allOf ) || defined( jsonObject, K::anyOf ) )
    {
        terminationSettings = getAs< boost::shared_ptr< PropagationHybridTerminationSettings > >( jsonObject );
        return;
    }

    const boost::shared_ptr< VariableSettings > variable =
            getValue< boost::shared_ptr< VariableSettings > >( jsonObject, K::variable );
    switch ( variable->variableType_ )
    {
    case independentVariable:
    {
        terminationSettings = boost::make_shared< PropagationTimeTerminationSettings >(
                    getEpoch( jsonObject, K::upperLimit,
                              getEpoch< double >( jsonObject, SpecialKeys::root / Keys::endEpoch, TUDAT_NAN, true ) ) );
        return;
    }
    case dependentVariable:
    {
        // If both lower limit and upper limit, create hybrid satistying any of the two conditions
        if ( defined( jsonObject, K::lowerLimit ) || defined( jsonObject, K::upperLimit ) )
        {
            // Lower limit
            json lowerLimitObject = jsonObject;
            lowerLimitObject.erase( K::upperLimit );
            boost::shared_ptr< PropagationTerminationSettings > lowerLimitCondition =
                    getAs< boost::shared_ptr< PropagationTerminationSettings > >( lowerLimitObject );

            // Upper limit
            json upperLimitObject = jsonObject;
            upperLimitObject.erase( K::lowerLimit );
            boost::shared_ptr< PropagationTerminationSettings > upperLimitCondition =
                    getAs< boost::shared_ptr< PropagationTerminationSettings > >( upperLimitObject );

            // Combine
            terminationSettings = boost::make_shared< PropagationHybridTerminationSettings >(
                        std::vector< boost::shared_ptr< PropagationTerminationSettings > >(
            { lowerLimitObject, upperLimitObject } ), true );
            return;
        }

        const boost::shared_ptr< SingleDependentVariableSaveSettings > dependentVar =
                boost::dynamic_pointer_cast< SingleDependentVariableSaveSettings >( variable );
        enforceNonNullPointer( dependentVar );

        // Limit value
        double limitValue;
        bool useLowerLimit;
        if ( defined( jsonObject, K::lowerLimit ) )
        {
            limitValue = getNumeric< double >( jsonObject, K::lowerLimit );
            useLowerLimit = true;
        }
        else
        {
            limitValue = getNumeric< double >( jsonObject, K::upperLimit );
            useLowerLimit = false;
        }

        terminationSettings = boost::make_shared< PropagationDependentVariableTerminationSettings >(
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
