#ifndef TUDATRESOURCES_RESOURCE_H
#define TUDATRESOURCES_RESOURCE_H

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

static inline std::string get_resources_path( )
{
    if( const char* path_p = getenv( "TUDATPY_RESOURCE_DIR" ) )
    {
        char resourcedir[ MAX_PATH ];
        snprintf( resourcedir, MAX_PATH, "%s", path_p );
        return std::string( resourcedir ).c_str( );
    }
    else  // TUDATPY_RESOURCE_DIR is not set
    {
        return std::string( get_homedir( ) + "/.tudat" + RESOURCE ).c_str( );
    }
}

static inline std::string get_atmosphere_tables_path( )
{
    return std::string( get_resources_path( ) + ATMOSPHERE_TABLES ).c_str( );
}

static inline std::string get_earth_orientation_path( )
{
    return std::string( get_resources_path( ) + EARTH_ORIENTATION ).c_str( );
}

static inline std::string get_earth_deformation_path( )
{
    return std::string( get_resources_path( ) + EARTH_ORIENTATION ).c_str( );
}

static inline std::string get_ephemeris_path( )
{
    return std::string( get_resources_path( ) + EPHEMERIS ).c_str( );
}

static inline std::string get_gravity_models_path( )
{
    return std::string( get_resources_path( ) + GRAVITY_MODELS ).c_str( );
}

static inline std::string get_quadrature_path( )
{
    return std::string( get_resources_path( ) + QUADRATURE ).c_str( );
}

static inline std::string get_space_weather_path( )
{
    return std::string( get_resources_path( ) + SPACE_WEATHER ).c_str( );
}

static inline std::string get_spice_kernels_path( )
{
    return std::string( get_resources_path( ) + SPICE_KERNELS ).c_str( );
}

static inline std::string get_star_catalog_biases_path( )
{
    return std::string( get_resources_path( ) + STAR_CATALOG_BIASES ).c_str( );
}

static inline std::string get_station_locations_path( )
{
    return std::string( get_resources_path( ) + STATION_LOCATIONS ).c_str( );
}
}  // namespace paths
}  // namespace tudat

#endif  // TUDATRESOURCES_RESOURCE_H
