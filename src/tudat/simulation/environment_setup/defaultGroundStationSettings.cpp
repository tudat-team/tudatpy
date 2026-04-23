/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/defaultGroundStationSettings.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/readSinexFile.h"

#include "tudat/simulation/environment_setup/createGroundStations.h"

namespace tudat
{
namespace simulation_setup
{

std::map< std::string, Eigen::Vector3d > getApproximateDsnGroundStationPositions( )
{
    // excluding DSS-23 and DSS-33-T19 since they are still in planning as of Rev. O (2024)
    std::map< std::string, Eigen::Vector3d > dsnStationPositionsItrf93 = {
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
          ( Eigen::Vector3d( ) << +4849245.211, -0360278.166, +4114884.445 )
                  .finished( ) },  // https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/stations/a_old_versions/dsnstns.cmt
        { "DSS-63", ( Eigen::Vector3d( ) << +4849092.518, -360180.3480, +4115109.251 ).finished( ) },
        { "DSS-65", ( Eigen::Vector3d( ) << +4849339.634, -360427.6637, +4114750.733 ).finished( ) },

    };

    return dsnStationPositionsItrf93;
}

Eigen::Vector3d getApproximateGroundStationPosition( std::string stationName )
{
    Eigen::Vector3d groundStationPosition;

    std::map< std::string, Eigen::Vector3d > dsnMap = getApproximateDsnGroundStationPositions( );
    if( dsnMap.count( stationName ) != 0 )
    {
        groundStationPosition = dsnMap.at( stationName );
    }
    else
    {
        throw std::runtime_error( "Error when retrieving approximate ground station position: station name " + stationName +
                                  "not recognized." );
    }

    return groundStationPosition;
}

const static std::string pysctrackGroundStationPosFile = tudat::paths::getStationLocationDataPath( ) + "/glo.sit";
const static std::string pysctrackGroundStationVelFile = tudat::paths::getStationLocationDataPath( ) + "/glo.vel";
const static std::string pysctrackGroundStationCodesFile = tudat::paths::getStationLocationDataPath( ) + "/ns_codes.dat";
const static std::string MPCGroundStationPosFile = tudat::paths::getStationLocationDataPath( ) + "/mpc.sit";
const static std::string MPCGroundStationVelFile = tudat::paths::getStationLocationDataPath( ) + "/mpc.vel";
const static std::string MPCGroundStationCodesFile = tudat::paths::getStationLocationDataPath( ) + "/mpc_codes.dat";

std::map< std::string, Eigen::Vector3d > getCombinedApproximateGroundStationPositions( )
{
    auto combinedMap = getApproximateDsnGroundStationPositions( );
    auto const& vlbiMap = getVlbiStationPositions( );
    combinedMap.insert( vlbiMap.begin( ), vlbiMap.end( ) );
    return combinedMap;
}

std::map< std::string, Eigen::Vector3d >& getVlbiStationPositions( )
{
    static std::map< std::string, Eigen::Vector3d > stationPositions =
            utilities::getMapFromFile< std::string, Eigen::Vector3d >( pysctrackGroundStationPosFile, '$', " \t" );
    return stationPositions;
}

std::map< std::string, Eigen::Vector3d >& getVlbiStationVelocities( )
{
    static std::map< std::string, Eigen::Vector3d > stationVelocities =
            utilities::getMapFromFile< std::string, Eigen::Vector3d >( pysctrackGroundStationVelFile, '$', " \t" );
    return stationVelocities;
}

std::map< std::string, Eigen::Vector3d >& getMPCStationPositions( )
{
    static std::map< std::string, Eigen::Vector3d > stationPositions =
            utilities::getMapFromFile< std::string, Eigen::Vector3d >( MPCGroundStationPosFile, '$', " \t" );
    return stationPositions;
}

std::map< std::string, Eigen::Vector3d >& getMPCStationVelocities( )
{
    static std::map< std::string, Eigen::Vector3d > stationVelocities =
            utilities::getMapFromFile< std::string, Eigen::Vector3d >( MPCGroundStationVelFile, '$', " \t" );
    return stationVelocities;
}

Eigen::Vector3d getDsnStationVelocity( std::string stationName )
{
    Eigen::Vector3d goldstoneStationVelocity( -0.0180, 0.0065, -0.0038 );
    goldstoneStationVelocity /= physical_constants::JULIAN_YEAR;
    Eigen::Vector3d canberraStationVelocity( -0.0335, -0.0041, 0.0392 );
    canberraStationVelocity /= physical_constants::JULIAN_YEAR;
    Eigen::Vector3d madridStationVelocity( -0.0100, 0.0242, 0.0156 );
    madridStationVelocity /= physical_constants::JULIAN_YEAR;
    Eigen::Vector3d stationVelocityItrf93 = Eigen::Vector3d::Constant( TUDAT_NAN );
    if( stationName[ 4 ] == '1' || stationName[ 4 ] == '2' )
    {
        stationVelocityItrf93 = goldstoneStationVelocity;
    }
    else if( stationName[ 4 ] == '3' || stationName[ 4 ] == '4' )
    {
        stationVelocityItrf93 = canberraStationVelocity;
    }
    else if( stationName[ 4 ] == '5' || stationName[ 4 ] == '6' )
    {
        stationVelocityItrf93 = madridStationVelocity;
    }
    else
    {
        throw std::runtime_error( "Error when retrieving approximate ground station velocity: station name " + stationName +
                                  "not recognized." );
    }

    return stationVelocityItrf93;
}

std::shared_ptr< GroundStationSettings > getDsnStationSetting( std::string stationName )
{
    // DSS positions: at 2003.0 with respect to ITRF93
    double stationPositionsReferenceEpoch = 3.0 * physical_constants::JULIAN_YEAR;
    std::map< std::string, Eigen::Vector3d > dsnStationPositionsItrf93 = getApproximateDsnGroundStationPositions( );

    if( dsnStationPositionsItrf93.find( stationName ) == dsnStationPositionsItrf93.end( ) )
    {
        throw std::runtime_error( "Error when retrieving approximate ground station position: station name " + stationName +
                                  "not recognized." );
    }

    // Get the station velocity in ITRF93
    Eigen::Vector3d stationVelocityItrf93 = getDsnStationVelocity( stationName );

    // Convert ground station state to ITRF2014
    Eigen::Vector6d stationStateItrf2014 = reference_frames::convertGroundStationStateArbitraryItrfToItrf2014(
            ( Eigen::Vector6d( ) << dsnStationPositionsItrf93[ stationName ], stationVelocityItrf93 ).finished( ),
            stationPositionsReferenceEpoch,
            "ITRF93" );

    std::shared_ptr< GroundStationMotionSettings > stationMotion =
            std::make_shared< LinearGroundStationMotionSettings >( stationStateItrf2014.segment( 3, 3 ), stationPositionsReferenceEpoch );

    std::shared_ptr< GroundStationSettings > stationSettings =
            std::make_shared< GroundStationSettings >( stationName, stationStateItrf2014.segment( 0, 3 ) );
    stationSettings->addStationMotionSettings( stationMotion );

    return stationSettings;
}

std::vector< std::shared_ptr< GroundStationSettings > > getDsnStationSettings( )
{
    std::map< std::string, Eigen::Vector3d > dsnStationPositionsItrf93 = getApproximateDsnGroundStationPositions( );

    std::vector< std::shared_ptr< GroundStationSettings > > stationSettingsList;

    for( auto it : dsnStationPositionsItrf93 )
    {
        std::shared_ptr< GroundStationSettings > stationSettings( getDsnStationSetting( it.first ) );

        stationSettingsList.push_back( stationSettings );
    }

    return stationSettingsList;
}

std::vector< std::shared_ptr< GroundStationSettings > > getEvnStationSettings( )
{
    std::vector< std::shared_ptr< GroundStationSettings > > stationSettingsList;

    std::map< std::string, Eigen::Vector3d > stationPositions = getVlbiStationPositions( );
    std::map< std::string, Eigen::Vector3d > stationVelocities = getVlbiStationVelocities( );

    std::vector< std::string > commonStationNames;
    for( auto it : stationPositions )
    {
        if( stationVelocities.count( it.first ) > 0 )
        {
            commonStationNames.push_back( it.first );
        }
    }

    for( unsigned int i = 0; i < commonStationNames.size( ); i++ )
    {
        std::shared_ptr< GroundStationMotionSettings > stationMotion = std::make_shared< LinearGroundStationMotionSettings >(
                stationVelocities.at( commonStationNames.at( i ) ) * 1.0E-3 / physical_constants::JULIAN_YEAR, 0.0 );

        std::shared_ptr< GroundStationSettings > stationSettings =
                std::make_shared< GroundStationSettings >( commonStationNames.at( i ), stationPositions.at( commonStationNames.at( i ) ) );
        stationSettings->addStationMotionSettings( stationMotion );
        stationSettingsList.push_back( stationSettings );
    }
    return stationSettingsList;
}

std::vector< std::shared_ptr< GroundStationSettings > > getMPCStationSettings( )
{
    std::vector< std::shared_ptr< GroundStationSettings > > stationSettingsList;

    std::map< std::string, Eigen::Vector3d > stationPositions = getMPCStationPositions( );
    std::map< std::string, Eigen::Vector3d > stationVelocities = getMPCStationVelocities( );

    std::vector< std::string > commonStationNames;
    for( auto it : stationPositions )
    {
        if( stationVelocities.count( it.first ) > 0 )
        {
            commonStationNames.push_back( it.first );
        }
    }

    for( unsigned int i = 0; i < commonStationNames.size( ); i++ )
    {
        std::shared_ptr< GroundStationMotionSettings > stationMotion = std::make_shared< LinearGroundStationMotionSettings >(
                stationVelocities.at( commonStationNames.at( i ) ) * 1.0E-3 / physical_constants::JULIAN_YEAR, 0.0 );

        std::shared_ptr< GroundStationSettings > stationSettings =
                std::make_shared< GroundStationSettings >( commonStationNames.at( i ), stationPositions.at( commonStationNames.at( i ) ) );
        stationSettings->addStationMotionSettings( stationMotion );
        stationSettingsList.push_back( stationSettings );
    }
    return stationSettingsList;
}

std::vector< std::shared_ptr< GroundStationSettings > > getRadioTelescopeStationSettings( )
{
    std::vector< std::shared_ptr< GroundStationSettings > > stations = getEvnStationSettings( );
    std::vector< std::shared_ptr< GroundStationSettings > > dsnStations = getDsnStationSettings( );
    stations.insert( stations.begin( ), dsnStations.begin( ), dsnStations.end( ) );
    return stations;
}

std::map< double, Eigen::Vector3d > getPiecewiseEccentricityDisplacementList(
        const std::vector< input_output::SinexStationEccentricity >& stationEccentricityHistory )
{
    std::map< double, Eigen::Vector3d > displacementList;
    if( stationEccentricityHistory.empty( ) )
    {
        return displacementList;
    }

    std::vector< input_output::SinexStationEccentricity > sortedEccentricityHistory = stationEccentricityHistory;
    std::sort( sortedEccentricityHistory.begin( ),
               sortedEccentricityHistory.end( ),
               []( const input_output::SinexStationEccentricity& firstEntry, const input_output::SinexStationEccentricity& secondEntry ) {
                   return firstEntry.startEpoch_ < secondEntry.startEpoch_;
               } );

    for( unsigned int i = 0; i + 1 < sortedEccentricityHistory.size( ); i++ )
    {
        if( sortedEccentricityHistory.at( i ).endEpoch_ > sortedEccentricityHistory.at( i + 1 ).startEpoch_ )
        {
            throw std::runtime_error( "Error when creating eccentricity displacement list from SINEX data: eccentricity arcs overlap." );
        }
    }

    const double minimumEpoch = -std::numeric_limits< double >::max( );
    if( sortedEccentricityHistory.at( 0 ).hasOpenEnd_ )
    {
        displacementList[ minimumEpoch ] = sortedEccentricityHistory.at( 0 ).eccentricity_;
    }

    for( unsigned int i = 0; i < sortedEccentricityHistory.size( ); i++ )
    {
        const input_output::SinexStationEccentricity& eccentricityEntry = sortedEccentricityHistory.at( i );
        displacementList[ eccentricityEntry.startEpoch_ ] = eccentricityEntry.eccentricity_;

        if( !eccentricityEntry.hasOpenEnd_ )
        {
            displacementList[ eccentricityEntry.endEpoch_ ] = Eigen::Vector3d::Zero( );
        }
    }

    return displacementList;
}

std::string getDefaultIlrsSinexStateFilePath( )
{
    return paths::getStationLocationDataPath( ) + "/SLRF2020_POS+VEL_2025.05.13.snx";
}

std::string getDefaultIlrsSinexEccentricityFilePath( )
{
    return paths::getStationLocationDataPath( ) + "/slrecc.250513.ILRS.xyz.snx";
}

std::vector< std::shared_ptr< GroundStationSettings > > getIlrsStationSettingsFromSinexDomes( const std::vector< std::string >& domesIds,
                                                                                              const std::string& sinexStateFile,
                                                                                              const std::string& sinexEccentricityFile,
                                                                                              const bool throwExceptionOnMissingData )
{
    const std::map< std::string, input_output::SinexStationState > sinexStateData = input_output::readSinexStationData( sinexStateFile );

    std::map< std::string, std::vector< input_output::SinexStationEccentricity > > sinexEccentricityData;
    if( !sinexEccentricityFile.empty( ) )
    {
        sinexEccentricityData = input_output::readSinexStationEccentricities( sinexEccentricityFile );
    }

    std::vector< std::shared_ptr< GroundStationSettings > > stationSettingsList;
    for( const std::string& domesId : domesIds )
    {
        if( sinexStateData.count( domesId ) == 0 )
        {
            if( throwExceptionOnMissingData )
            {
                throw std::runtime_error( "Error when creating ILRS station settings from SINEX: DOMES id " + domesId +
                                          " is not in SINEX state file." );
            }
            continue;
        }

        const input_output::SinexStationState currentState = sinexStateData.at( domesId );
        if( !( currentState.position_( 0 ) == currentState.position_( 0 ) ) )
        {
            if( throwExceptionOnMissingData )
            {
                throw std::runtime_error( "Error when creating ILRS station settings from SINEX: no Cartesian position for DOMES id " +
                                          domesId + "." );
            }
            continue;
        }

        std::shared_ptr< GroundStationSettings > stationSettings =
                std::make_shared< GroundStationSettings >( domesId, currentState.position_ );

        if( sinexEccentricityData.count( domesId ) > 0 )
        {
            const std::map< double, Eigen::Vector3d > displacementList =
                    getPiecewiseEccentricityDisplacementList( sinexEccentricityData.at( domesId ) );
            if( !displacementList.empty( ) )
            {
                stationSettings->addStationMotionSettings(
                        std::make_shared< PiecewiseConstantGroundStationMotionSettings >( displacementList ) );
            }
        }

        if( currentState.velocity_( 0 ) == currentState.velocity_( 0 ) )
        {
            const double referenceEpoch =
                    ( currentState.referenceEpoch_ == currentState.referenceEpoch_ ) ? currentState.referenceEpoch_ : 0.0;
            stationSettings->addStationMotionSettings(
                    std::make_shared< LinearGroundStationMotionSettings >( currentState.velocity_, referenceEpoch ) );
        }

        stationSettingsList.push_back( stationSettings );
    }

    return stationSettingsList;
}

}  // namespace simulation_setup

}  // namespace tudat
