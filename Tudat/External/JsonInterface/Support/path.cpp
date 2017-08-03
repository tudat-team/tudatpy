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

//! Replace recognized paths with placeholders such as ${SRCROOT}.
std::string addPathPlaceholders( std::string path )
{
    for ( std::vector< std::pair< std::string, std::string > >::reverse_iterator rit = pathPlaceholders.rbegin( );
          rit != pathPlaceholders.rend( ); ++rit )
    {
        const std::string placeholderId = rit->first;
        const std::string placeholderPath = boost::filesystem::canonical( rit->second ).string( );
        path = boost::regex_replace( path, boost::regex( placeholderPath ), "${" + placeholderId + "}" );
    }
    return path;
}

//! Replace recognized placeholders such as ${SRCROOT} with the actual paths.
std::string removePathPlaceholders( std::string path )
{
    for ( std::vector< std::pair< std::string, std::string > >::reverse_iterator rit = pathPlaceholders.rbegin( );
          rit != pathPlaceholders.rend( ); ++rit )
    {
        const std::string placeholderId = rit->first;
        const std::string placeholderPath = boost::filesystem::canonical( rit->second ).string( );
        path = boost::regex_replace( path, boost::regex( R"(\$\{)" + placeholderId + R"(\})" ), placeholderPath );
    }
    return path;
}

} // namespace json_interface

} // namespace tudat
