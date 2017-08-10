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

//! Return \p path with the recognized paths replaced by placeholders (such as ${TUDAT_ROOT_PATH}).
std::string pathAddingPlaceholders( std::string path )
{
    using namespace boost;
    for ( std::vector< std::pair< std::string, std::string > >::reverse_iterator rit = pathPlaceholders.rbegin( );
          rit != pathPlaceholders.rend( ); ++rit )
    {
        const std::string placeholderId = rit->first;
        const std::string placeholderPath = filesystem::canonical( rit->second ).string( );
        path = regex_replace( path, regex( placeholderPath ), "${" + placeholderId + "}" );
    }
    const std::string relativePath = filesystem::relative( path, filesystem::current_path( ) ).string( );
    if ( ! relativePath.empty( ) && relativePath.size( ) < path.size( ) )
    {
        path = relativePath;
    }
    return path;
}

//! Return \p path with the recognized path placeholders (such as ${TUDAT_ROOT_PATH}) replaced by the actual paths.
std::string pathRemovingPlaceholders( std::string path )
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
    j = tudat::json_interface::pathAddingPlaceholders( weakly_canonical( absolutePath ).string( ) );
}

//! Create a path from a `json` object.
void from_json( const json& j, path& p )
{
    const path actualPath = tudat::json_interface::pathRemovingPlaceholders( j.get< std::string >( ) );
    const path absolutePath = actualPath.is_absolute( ) ? actualPath : current_path( ) / actualPath;
    p = weakly_canonical( absolutePath );
}

}  // namespace filesystem

}  // namespace boost
