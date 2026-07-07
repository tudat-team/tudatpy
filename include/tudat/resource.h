#ifndef TUDATRESOURCES_RESOURCE_H
#define TUDATRESOURCES_RESOURCE_H

#include <tudat/resource/config.hpp>

#define RESOURCE "/resource"
#define ATMOSPHERE_TABLES "/atmosphere_tables"
#define EARTH_ORIENTATION "/earth_orientation"
#define EARTH_DEFORMATION "/earth_deformation"
#define EPHEMERIS "/ephemeris"
#define GRAVITY_MODELS "/gravity_models"
#define QUADRATURE "/quadrature"
#define SPACE_WEATHER "/space_weather"
#define SPICE_KERNELS "/spice_kernels"
#define STAR_CATALOG_BIASES "/star_catalog_biases"
#define STATION_LOCATIONS "/station_locations"

#define MAX_PATH 255

#include <cstdlib>
#include <cstring>
#include <string>

// #define MAX_PREFIX_LEN 256
// #define MAX_RESOURCE_LEN strlen("/resource")
// #define MAX_RELATIVE_LEN 20

namespace tudat
{
namespace paths
{

// https://cboard.cprogramming.com/c-programming/164689-how-get-users-home-directory.html
static inline std::string get_homedir( void )
{
    char homedir[ MAX_PATH ];
#ifdef _WIN32
    snprintf( homedir, MAX_PATH, "%s%s", getenv( "HOMEDRIVE" ), getenv( "HOMEPATH" ) );
#else
    snprintf( homedir, MAX_PATH, "%s", getenv( "HOME" ) );
#endif
    return std::string( homedir );
}

static inline std::string get_prefix_path( )
{
    return std::string( TUDAT_RESOURCE_PREFIX ).c_str( );
}

static inline std::string get_hidden_path( )
{
    return std::string( get_homedir( ) + "/.tudat" ).c_str( );
}

static inline std::string get_resources_path( )
{
    if( const char* path_p = getenv( "TUDAT_RESOURCE_DIR" ) )
    {
        char resourcedir[ MAX_PATH ];
        snprintf( resourcedir, MAX_PATH, "%s", path_p );
        return std::string( resourcedir ).c_str( );
    }
    else  // TUDAT_RESOURCE_DIR is not set
    {
        return std::string( get_hidden_path( ) + RESOURCE ).c_str( );
    }
}

static inline std::string get_atmosphere_tables_path( )
{
    return std::string( get_hidden_path( ) + RESOURCE + ATMOSPHERE_TABLES ).c_str( );
}

static inline std::string get_earth_orientation_path( )
{
    return std::string( get_hidden_path( ) + RESOURCE + EARTH_ORIENTATION ).c_str( );
}

static inline std::string get_earth_deformation_path( )
{
    return std::string( get_hidden_path( ) + RESOURCE + EARTH_ORIENTATION ).c_str( );
}

static inline std::string get_ephemeris_path( )
{
    return std::string( get_hidden_path( ) + RESOURCE + EPHEMERIS ).c_str( );
}

static inline std::string get_gravity_models_path( )
{
    return std::string( get_hidden_path( ) + RESOURCE + GRAVITY_MODELS ).c_str( );
}

static inline std::string get_quadrature_path( )
{
    return std::string( get_hidden_path( ) + RESOURCE + QUADRATURE ).c_str( );
}

static inline std::string get_space_weather_path( )
{
    return std::string( get_hidden_path( ) + RESOURCE + SPACE_WEATHER ).c_str( );
}

static inline std::string get_spice_kernels_path( )
{
    return std::string( get_hidden_path( ) + RESOURCE + SPICE_KERNELS ).c_str( );
}

static inline std::string get_star_catalog_biases_path( )
{
    return std::string( get_hidden_path( ) + RESOURCE + STAR_CATALOG_BIASES ).c_str( );
}

static inline std::string get_station_locations_path( )
{
    return std::string( get_hidden_path( ) + RESOURCE + STATION_LOCATIONS ).c_str( );
}
}  // namespace paths
}  // namespace tudat

#endif  // TUDATRESOURCES_RESOURCE_H
