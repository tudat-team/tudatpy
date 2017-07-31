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


}

} // namespace propagators

} // namespace tudat
