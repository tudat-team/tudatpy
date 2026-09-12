/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <memory>
#include "tudat/simulation/environment_setup/defaultGroundStationSettings.h"
#include "tudat/io/dsnStationData.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/readSinexFile.h"

#include "tudat/simulation/environment_setup/createGroundStations.h"

namespace tudat
{
namespace simulation_setup
{

const std::map< std::string, Eigen::Vector3d >& getApproximateDsnGroundStationPositions( )
{
    return input_output::getApproximateDsnGroundStationPositions( );
}

int getDsnComplexId( const std::string& stationName )
{
    return input_output::getDsnComplexId( stationName );
}

std::map< int, std::vector< std::string > > getDefaultDsnStationNamesPerComplex( )
{
    return input_output::getDefaultDsnStationNamesPerComplex( );
}

Eigen::Vector3d getApproximateGroundStationPosition( const std::string& stationName )
{
    Eigen::Vector3d groundStationPosition;

    const std::map< std::string, Eigen::Vector3d >& dsnMap = getApproximateDsnGroundStationPositions( );
    if( dsnMap.count( stationName ) != 0 )
    {
        groundStationPosition = dsnMap.at( stationName );
    }
    else
    {
        throw std::runtime_error( "Error when retrieving approximate ground station position: station name " + stationName +
                                  " not recognized." );
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

const std::map< std::string, Eigen::Vector3d >& getVlbiStationPositions( )
{
    static const std::map< std::string, Eigen::Vector3d > stationPositions =
            utilities::getMapFromFile< std::string, Eigen::Vector3d >( pysctrackGroundStationPosFile, '$', " \t" );
    return stationPositions;
}

const std::map< std::string, Eigen::Vector3d >& getVlbiStationVelocities( )
{
    static const std::map< std::string, Eigen::Vector3d > stationVelocities =
            utilities::getMapFromFile< std::string, Eigen::Vector3d >( pysctrackGroundStationVelFile, '$', " \t" );
    return stationVelocities;
}

const std::map< std::string, Eigen::Vector3d >& getMPCStationPositions( )
{
    static const std::map< std::string, Eigen::Vector3d > stationPositions =
            utilities::getMapFromFile< std::string, Eigen::Vector3d >( MPCGroundStationPosFile, '$', " \t" );
    return stationPositions;
}

const std::map< std::string, Eigen::Vector3d >& getMPCStationVelocities( )
{
    static const std::map< std::string, Eigen::Vector3d > stationVelocities =
            utilities::getMapFromFile< std::string, Eigen::Vector3d >( MPCGroundStationVelFile, '$', " \t" );
    return stationVelocities;
}

Eigen::Vector3d getDsnStationVelocity( const std::string& stationName )
{
    Eigen::Vector3d goldstoneStationVelocity( -0.0180, 0.0065, -0.0038 );
    goldstoneStationVelocity /= physical_constants::JULIAN_YEAR;
    Eigen::Vector3d canberraStationVelocity( -0.0335, -0.0041, 0.0392 );
    canberraStationVelocity /= physical_constants::JULIAN_YEAR;
    Eigen::Vector3d madridStationVelocity( -0.0100, 0.0242, 0.0156 );
    madridStationVelocity /= physical_constants::JULIAN_YEAR;
    Eigen::Vector3d stationVelocityItrf93 = Eigen::Vector3d::Constant( TUDAT_NAN );

    int complexId = getDsnComplexId( stationName );
    switch( complexId )
    {
        case 10:
            stationVelocityItrf93 = goldstoneStationVelocity;
            break;
        case 40:
            stationVelocityItrf93 = canberraStationVelocity;
            break;
        case 60:
            stationVelocityItrf93 = madridStationVelocity;
            break;
        default:
            throw std::runtime_error( "Error when retrieving approximate ground station velocity: station complex of " + stationName +
                                      " not recognized." );
    }

    return stationVelocityItrf93;
}

std::shared_ptr< GroundStationSettings > getDsnStationSetting( const std::string& stationName )
{
    // DSS positions: at 2003.0 with respect to ITRF93
    double stationPositionsReferenceEpoch = 3.0 * physical_constants::JULIAN_YEAR;
    Eigen::Vector3d stationPositionItrf93 = getApproximateGroundStationPosition( stationName );

    // Get the station velocity in ITRF93
    Eigen::Vector3d stationVelocityItrf93 = getDsnStationVelocity( stationName );

    // Convert ground station state to ITRF2014
    Eigen::Vector6d stationStateItrf2014 = reference_frames::convertGroundStationStateArbitraryItrfToItrf2014(
            ( Eigen::Vector6d( ) << stationPositionItrf93, stationVelocityItrf93 ).finished( ), stationPositionsReferenceEpoch, "ITRF93" );

    auto stationSettings = std::make_shared< GroundStationSettings >( stationName, stationStateItrf2014.segment( 0, 3 ) );

    if( stationName == "DSS-65" )
    {
        // DSS-65 was moved in 2005 to its current new location
        // define piecewise-constant displacement according to
        // https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/stations/a_old_versions/earthstns_itrf93_050714.cmt
        const Eigen::Vector3d dss65Pre2005Position = ( Eigen::Vector3d( ) << 4849336.6176, -360488.6349, 4114748.9218 ).finished( );
        const std::map< double, Eigen::Vector3d > displacementList = {
            { basic_astrodynamics::DateTime( 1987, 1, 1, 0, 0, 0 ).epoch< double >( ),
              dss65Pre2005Position - stationStateItrf2014.segment( 0, 3 ) },
            { basic_astrodynamics::DateTime( 2005, 7, 3, 0, 0, 0 ).epoch< double >( ), Eigen::Vector3d::Zero( ) }
        };

        std::shared_ptr< GroundStationMotionSettings > dss65Displacement =
                std::make_shared< PiecewiseConstantGroundStationMotionSettings >( displacementList );
        stationSettings->addStationMotionSettings( dss65Displacement );
    }

    std::shared_ptr< GroundStationMotionSettings > stationMotion =
            std::make_shared< LinearGroundStationMotionSettings >( stationStateItrf2014.segment( 3, 3 ), stationPositionsReferenceEpoch );

    stationSettings->addStationMotionSettings( stationMotion );

    return stationSettings;
}

std::vector< std::shared_ptr< GroundStationSettings > > getDsnStationSettings( )
{
    auto const& dsnStationPositionsItrf93 = getApproximateDsnGroundStationPositions( );

    std::vector< std::shared_ptr< GroundStationSettings > > stationSettingsList;

    for( auto const& [ stationName, position ] : dsnStationPositionsItrf93 )
    {
        std::shared_ptr< GroundStationSettings > stationSettings( getDsnStationSetting( stationName ) );

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
    for( auto const& [ stationName, position ] : stationPositions )
    {
        if( stationVelocities.count( stationName ) > 0 )
        {
            commonStationNames.push_back( stationName );
        }
    }

    for( auto const& stationName : commonStationNames )
    {
        std::shared_ptr< GroundStationMotionSettings > stationMotion = std::make_shared< LinearGroundStationMotionSettings >(
                stationVelocities.at( stationName ) * 1.0E-3 / physical_constants::JULIAN_YEAR, 0.0 );

        std::shared_ptr< GroundStationSettings > stationSettings =
                std::make_shared< GroundStationSettings >( stationName, stationPositions.at( stationName ) );
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
    for( auto const& [ stationName, position ] : stationPositions )
    {
        if( stationVelocities.count( stationName ) > 0 )
        {
            commonStationNames.push_back( stationName );
        }
    }

    for( auto const& stationName : commonStationNames )
    {
        std::shared_ptr< GroundStationMotionSettings > stationMotion = std::make_shared< LinearGroundStationMotionSettings >(
                stationVelocities.at( stationName ) * 1.0E-3 / physical_constants::JULIAN_YEAR, 0.0 );

        std::shared_ptr< GroundStationSettings > stationSettings =
                std::make_shared< GroundStationSettings >( stationName, stationPositions.at( stationName ) );
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
