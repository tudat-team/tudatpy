/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_UNITTESTS_JSONINTERFACE_SUPPORT
#define TUDAT_UNITTESTS_JSONINTERFACE_SUPPORT

#include <fstream>

#include "json/src/json.hpp"
using json = nlohmann::json;

#include "Tudat/External/JsonInterface/Support/utilities.h"

namespace tudat
{

namespace unit_tests
{

template< typename T = json >
T readFile( const std::string& filename, const std::string& extension = "json" )
{
    const std::string filePath = ( boost::filesystem::path( __FILE__ ).parent_path( ) /
				   "inputs" / ( filename + "." + extension ) ).string( );
    return json::parse( std::ifstream( filePath ) ).get< T >( );
}

template< typename Enum >
bool isEnumConsistent( const std::vector< Enum >& values,
		       const std::map< Enum, std::string >& stringValues,
		       const std::vector< Enum >& usupportedValues )
{
    using namespace json_interface;

    // Create vector of supported values
    std::vector< Enum > supportedValues;
    for ( auto entry : stringValues )
    {
	Enum value = entry.first;
	if ( ! contains( usupportedValues, value ) )
	{
	    supportedValues.push_back( value );
	}
    }

    // Check that `values` contains all the elements in `supportedValues`, and viceversa
    return containsAllOf( values, supportedValues ) && containsAllOf( supportedValues, values );
}

template< typename Enum >
bool isEnumConsistent( const std::string& filename,
		       const std::map< Enum, std::string >& stringValues,
		       const std::vector< Enum >& usupportedValues )
{
    return isEnumConsistent( readFile< std::vector< Enum > >( filename ), stringValues, usupportedValues );
}

} // namespace unit_tests

} // namespace tudat

#endif // TUDAT_UNITTESTS_JSONINTERFACE_SUPPORT
