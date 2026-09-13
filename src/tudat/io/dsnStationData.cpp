/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/dsnStationData.h"

namespace tudat
{

namespace input_output
{

const std::map< std::string, Eigen::Vector3d >& getApproximateDsnGroundStationPositions( )
{
    // excluding DSS-23 and DSS-33-T19 since they are still in planning as of Rev. O (2024)
    static const std::map< std::string, Eigen::Vector3d > dsnStationPositionsItrf93 = {
        { "DSS-12",
          ( Eigen::Vector3d( ) << -2350443.812, -4651980.837, +3665630.988 )
                  .finished( ) },  // https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/stations/a_old_versions/dsnstns.cmt
        { "DSS-13", ( Eigen::Vector3d( ) << -2351112.659, -4655530.636, +3660912.728 ).finished( ) },
        { "DSS-14", ( Eigen::Vector3d( ) << -2353621.420, -4641341.472, +3677052.318 ).finished( ) },
        { "DSS-15", ( Eigen::Vector3d( ) << -2353538.958, -4641649.429, +3676669.984 ).finished( ) },
        { "DSS-24", ( Eigen::Vector3d( ) << -2354906.711, -4646840.095, +3669242.325 ).finished( ) },
        { "DSS-25", ( Eigen::Vector3d( ) << -2355022.014, -4646953.204, +3669040.567 ).finished( ) },
        { "DSS-26", ( Eigen::Vector3d( ) << -2354890.797, -4647166.328, +3668871.755 ).finished( ) },
        { "DSS-27", ( Eigen::Vector3d( ) << -2349915.428, -4656756.406, +3660096.469 ).finished( ) },
        { "DSS-34", ( Eigen::Vector3d( ) << -4461147.093, +2682439.239, -3674393.133 ).finished( ) },
        { "DSS-35", ( Eigen::Vector3d( ) << -4461273.090, +2682568.925, -3674152.093 ).finished( ) },
        { "DSS-36", ( Eigen::Vector3d( ) << -4461168.415, +2682814.657, -3674083.901 ).finished( ) },
        { "DSS-42",
          ( Eigen::Vector3d( ) << -4460981.016, +2682413.525, -3674582.072 )
                  .finished( ) },  // https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/stations/a_old_versions/dsnstns.cmt
        { "DSS-43", ( Eigen::Vector3d( ) << -4460894.917, +2682361.507, -3674748.152 ).finished( ) },
        { "DSS-45", ( Eigen::Vector3d( ) << -4460935.578, +2682765.661, -3674380.982 ).finished( ) },
        //{ "DSS-47",
        //( Eigen::Vector3d( ) << -4750915.837, +2792906.182, -3200483.747 )
        //          .finished( ) },  // W196:
        // https://www.narrabri.atnf.csiro.au/observing/users_guide/html_old_20090512/Cartesian_Coordinates.html
        // (may not be entirely accurate)
        { "DSS-53", ( Eigen::Vector3d( ) << +4849338.209, -360657.812, +4114746.173 ).finished( ) },
        { "DSS-54", ( Eigen::Vector3d( ) << +4849434.488, -360723.8999, +4114618.835 ).finished( ) },
        { "DSS-55", ( Eigen::Vector3d( ) << +4849525.256, -360606.0932, +4114495.084 ).finished( ) },
        { "DSS-56", ( Eigen::Vector3d( ) << +4849421.679, -360549.659, +4114646.987 ).finished( ) },
        { "DSS-61",
          ( Eigen::Vector3d( ) << +4849245.211, -360278.166, +4114884.445 )
                  .finished( ) },  // https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/stations/a_old_versions/dsnstns.cmt
        { "DSS-63", ( Eigen::Vector3d( ) << +4849092.518, -360180.3480, +4115109.251 ).finished( ) },
        { "DSS-65", ( Eigen::Vector3d( ) << +4849339.634, -360427.6637, +4114750.733 ).finished( ) },

    };

    return dsnStationPositionsItrf93;
}

int getDsnComplexId( const std::string& stationName )
{
    if( stationName.size( ) < 5 || stationName.substr( 0, 4 ) != "DSS-" )
    {
        throw std::runtime_error( "Error when retrieving DSN complex ID: station name " + stationName + " not recognized." );
    }

    switch( stationName[ 4 ] )
    {
        case '1':
        case '2':
            return 10;  // Goldstone
        case '3':
        case '4':
            return 40;  // Canberra
        case '5':
        case '6':
            return 60;  // Madrid
        default:
            throw std::runtime_error( "Error when retrieving DSN complex ID: station name " + stationName + " not recognized." );
    }
}

std::map< int, std::vector< std::string > > getDefaultDsnStationNamesPerComplex( )
{
    std::map< int, std::vector< std::string > > stations;
    for( const auto& station : getApproximateDsnGroundStationPositions( ) )
    {
        stations[ getDsnComplexId( station.first ) ].push_back( station.first );
    }
    return stations;
}

}  // namespace input_output

}  // namespace tudat
