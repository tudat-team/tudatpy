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

#include "valueAccess.h"

namespace tudat
{

namespace json_interface
{

//! -DOC
bool defined( const json& jsonObject, const KeyPath& keyPath )
{
    if ( getOptional< json >( jsonObject, keyPath ) )
    {
        return true;
    }
    else
    {
        return false;
    }
}

//! -DOC
boost::shared_ptr< json > getRootObject( const json& jsonObject )
{
    return getOptional< json >( jsonObject, SpecialKeys::rootObject );
}

//! -DOC
KeyPath getKeyPath( const json& jsonObject )
{
    return getValue< KeyPath >( jsonObject, SpecialKeys::keyPath, SpecialKeys::root );
}

//! -DOC
std::string getParentKey( const json& jsonObject, const std::string& errorMessage )
{
    try
    {
        return getKeyPath( jsonObject ).back( );
    }
    catch ( ... )
    {
        std::cerr << errorMessage << std::endl;
        throw;
    }
}

}  // namespace json_interface

}  // namespace tudat
