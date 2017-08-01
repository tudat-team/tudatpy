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

#include "dependentVariable.h"

namespace tudat
{

namespace propagators
{

//! Create a `json` object from a shared pointer to a `PropagationTerminationSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< PropagationTerminationSettings >& terminationSettings )
{
    if ( ! terminationSettings )
    {
        return;
    }
    using namespace json_interface;
    using K = Keys::Propagator::Termination;


}

//! Create a shared pointer to a `PropagationTerminationSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< PropagationTerminationSettings >& terminationSettings )
{
    using namespace json_interface;
    using K = Keys::Propagator::Termination;

    if ( defined( jsonObject, K::conditions ) )  // hybrid
    {
        PropagationHybridTerminationSettings defaults( { } );
        terminationSettings = boost::make_shared< PropagationHybridTerminationSettings >(
                    getValue< std::vector< boost::shared_ptr< PropagationTerminationSettings > > >(
                        jsonObject, K::conditions ),
                    getValue( jsonObject, K::stopIfSingleConditionMet, defaults.fulFillSingleCondition_ ) );
        return;
    }
    else
    {
        const json variable = getValue< json >( jsonObject, K::variable );
        if ( variable.is_string( ) )
        {
            if ( variable.get< std::string >( ) == "time" )  // time
            {
                terminationSettings = boost::make_shared< PropagationTimeTerminationSettings >(
                            getEpoch( jsonObject, K::limitValue, getEpoch< double >(
                                          jsonObject, SpecialKeys::root / KeyPaths::Simulation::endEpoch ) ) );
                return;
            }
        }
        // dependent variable
        terminationSettings = boost::make_shared< PropagationDependentVariableTerminationSettings >(
                    getAs< boost::shared_ptr< SingleDependentVariableSaveSettings > >( variable ),
                    getNumeric< double >( jsonObject, K::limitValue ),
                    getValue< bool >( jsonObject, K::useAsLowerLimit ) );
        return;
    }
}

} // namespace propagators

} // namespace tudat
