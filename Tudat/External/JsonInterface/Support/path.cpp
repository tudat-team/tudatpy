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

#include <boost/regex.hpp>

#include "path.h"

namespace tudat
{

namespace json_interface
{

//! Replace recognized paths with placeholders such as ${TUDAT_ROOT}.
std::string addPathPlaceholders( std::string path )
{
    using namespace boost;
    for ( std::vector< std::pair< std::string, std::string > >::reverse_iterator rit = pathPlaceholders.rbegin( );
          rit != pathPlaceholders.rend( ); ++rit )
    {
        const std::string placeholderId = rit->first;
        const std::string placeholderPath = filesystem::canonical( rit->second ).string( );
        path = regex_replace( path, regex( placeholderPath ), "${" + placeholderId + "}" );
    }
    filesystem::path currentPath = filesystem::current_path( );
    currentPath += filesystem::path::preferred_separator;
    return regex_replace( path, regex( currentPath.string( ) ), "" );
}

//! Replace recognized placeholders such as ${SRCROOT} with the actual paths.
std::string removePathPlaceholders( std::string path )
{
    using namespace boost;
    for ( std::vector< std::pair< std::string, std::string > >::reverse_iterator rit = pathPlaceholders.rbegin( );
          rit != pathPlaceholders.rend( ); ++rit )
    {
        const std::string placeholderId = rit->first;
        const std::string placeholderPath = filesystem::canonical( rit->second ).string( );
        path = regex_replace( path, regex( R"(\$\{)" + placeholderId + R"(\})" ), placeholderPath );
    }
    return path;
}

} // namespace json_interface

} // namespace tudat


namespace boost
{

namespace filesystem
{

//! Create a `json` object from a `path`.
void to_json( json& j, const path& p )
{
    const path absolutePath = p.is_absolute( ) ? p : current_path( ) / p;
    j = tudat::json_interface::addPathPlaceholders( weakly_canonical( absolutePath ).string( ) );
}

//! Create a path from a `json` object.
void from_json( const json& j, path& p )
{
    const path actualPath = tudat::json_interface::removePathPlaceholders( j.get< std::string >( ) );
    const path absolutePath = actualPath.is_absolute( ) ? actualPath : current_path( ) / actualPath;
    p = weakly_canonical( absolutePath );
}

}  // namespace filesystem

}  // namespace boost
