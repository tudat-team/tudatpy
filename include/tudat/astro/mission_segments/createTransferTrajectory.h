/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      Note that the exact implementation of Newton-Raphson as root finder should be updated if
 *      someone would want to use a different root-finding technique.
 *
 *      By default the eccentricity is used as the iteration procedure. This is because in
 *      optimizing a Cassini-like trajectory, the pericenter radius had about 2-4 NaN values in
 *      100000 times the gravity assist calculation. The eccentricity iteration had no NaN values
 *      for a similar run in which 100000 gravity assist calculations were done. Also the
 *      eccentricity seemed to require slightly less iterations (does not necessarily mean it is
 *      faster or more accurate).
 *
 */

#ifndef TUDAT_CREATE_TRANSFER_TRAJECTORY_H
#define TUDAT_CREATE_TRANSFER_TRAJECTORY_H

#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "tudat/astro/mission_segments/transferLeg.h"
#include "tudat/astro/mission_segments/transferNode.h"
#include "tudat/astro/mission_segments/transferTrajectory.h"
#include "tudat/simulation/environment/body.h"

namespace tudat
{
namespace mission_segments
{


class TransferLegSettings
{
public:
    TransferLegSettings(
            const std::string& departureBodyName,
            const std::string& arrivalBodyName,
            const std::string& centralBodyName,
            const TransferLegTypes legType ){ }

    const std::string departureBodyName_;
    const std::string arrivalBodyName_;
    const std::string centralBodyName_;

    const TransferLegTypes legType_;
};


class TransferNodeSettings
{
public:
    TransferNodeSettings(
            const std::string& nodeBodyName,
            const TransferNodeTypes nodeType ):
        nodeBodyName_( nodeBodyName ), nodeType_( nodeType ){ }

    const std::string nodeBodyName_;
    const TransferNodeTypes nodeType_;
};

class PoweredSwingbyNodeSettings: public TransferNodeSettings
{
public:
    PoweredSwingbyNodeSettings(
            const std::string& nodeBodyName,
            const TransferNodeTypes nodeType,
            const double minimumPeriapsisRadius ):
        TransferNodeSettings( nodeBodyName, nodeType ),
    minimumPeriapsisRadius_( minimumPeriapsisRadius ){ }

    const std::string nodeBodyName_;
    const TransferNodeTypes nodeType_;
    const double minimumPeriapsisRadius_;
};

class EscapeAndDepartureNodeSettings: public TransferNodeSettings
{
public:
    EscapeAndDepartureNodeSettings(
            const std::string& nodeBodyName,
            const TransferNodeTypes nodeType,
            const double departureSemiMajorAxis,
            const double departureEccentricity,
            const double maximumExcessVelocity ):
        TransferNodeSettings( nodeBodyName, nodeType ),
    departureSemiMajorAxis_( departureSemiMajorAxis ),
    departureEccentricity_( departureEccentricity ),
    maximumExcessVelocity_( maximumExcessVelocity ){ }

    const std::string nodeBodyName_;
    const TransferNodeTypes nodeType_;
    const double departureSemiMajorAxis_;
    const double departureEccentricity_;
    const double maximumExcessVelocity_;
};

class CaptureAndInsertionNodeSettings: public TransferNodeSettings
{
public:
    CaptureAndInsertionNodeSettings(
            const std::string& nodeBodyName,
            const TransferNodeTypes nodeType,
            const double captureSemiMajorAxis,
            const double captureEccentricity ):
        TransferNodeSettings( nodeBodyName, nodeType ),
    captureSemiMajorAxis_( captureSemiMajorAxis ),
    captureEccentricity_( captureEccentricity ){ }

    const std::string nodeBodyName_;
    const TransferNodeTypes nodeType_;
    const double captureSemiMajorAxis_;
    const double captureEccentricity_;
};

std::shared_ptr< TransferLeg > createTransferLeg(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< TransferLegSettings > legSettings,
        const double legStartTime,
        const double legEndTime,
        const Eigen::VectorXd& legFreeParameters )
{
    if( bodyMap.count( legSettings->centralBodyName_ ) == 0 )
    {
        throw std::runtime_error( "Error when making transfer leg, central body " + legSettings->centralBodyName_  + " not found." );
    }
    else if( bodyMap.at( legSettings->centralBodyName_ )->getGravityFieldModel( ) == nullptr )
    {
        throw std::runtime_error( "Error when making transfer leg, central body " + legSettings->centralBodyName_  + " has no gravity field." );
    }

    if( bodyMap.count( legSettings->departureBodyName_ ) == 0 )
    {
        throw std::runtime_error( "Error when making transfer leg, departure body " + legSettings->departureBodyName_  + " not found." );
    }
    else if( bodyMap.at( legSettings->departureBodyName_ )->getEphemeris( ) == nullptr )
    {
        throw std::runtime_error( "Error when making transfer leg, departure body " + legSettings->departureBodyName_  + " has no ephemeris." );
    }

    if( bodyMap.count( legSettings->arrivalBodyName_ ) == 0 )
    {
        throw std::runtime_error( "Error when making transfer leg, arrival body " + legSettings->arrivalBodyName_  + " not found." );
    }
    else if( bodyMap.at( legSettings->arrivalBodyName_ )->getEphemeris( ) == nullptr )
    {
        throw std::runtime_error( "Error when making transfer leg, arrival body " + legSettings->arrivalBodyName_  + " has no ephemeris." );
    }


    const double centralBodyGravitationalParameter =
            bodyMap.at( legSettings->centralBodyName_ )->getGravityFieldModel( )->getGravitationalParameter( );
    const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris =
            bodyMap.at( legSettings->departureBodyName_ )->getEphemeris( );
    const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris =
            bodyMap.at( legSettings->arrivalBodyName_ )->getEphemeris( );

    std::shared_ptr< TransferLeg >  transferLeg;
    switch( legSettings->legType_ )
    {
    case unpowered_unperturbed_leg:
    {
        if( legFreeParameters.rows( ) != 0 )
        {
            throw std::runtime_error(
                        "Error when making unpowered_unperturbed_leg, number of extra free parameters should be zero, was " +
                        std::to_string( legFreeParameters.rows( ) ) );
        }
        transferLeg = std::make_shared< UnpoweredUnperturbedTransferLeg >(
                    departureBodyEphemeris, arrivalBodyEphemeris,
                    ( Eigen::VectorXd( 2 )<< legStartTime, legEndTime ).finished( ),
                    centralBodyGravitationalParameter );
        break;
    }
    case dsm_position_based_leg:
    {
        if( legFreeParameters.rows( ) != 4 )
        {
            throw std::runtime_error(
                        "Error when making dsm_position_based_leg, number of extra free parameters should be four, was " +
                        std::to_string( legFreeParameters.rows( ) ) );
        }
        transferLeg = std::make_shared< DsmPositionBasedTransferLeg >(
                    departureBodyEphemeris, arrivalBodyEphemeris,
                    ( Eigen::VectorXd( 6 )<< legStartTime, legEndTime, legFreeParameters ).finished( ),
                    centralBodyGravitationalParameter );
        break;
    }
    case dsm_velocity_based_leg:
    {
        if( legFreeParameters.rows( ) != 3 )
        {
            throw std::runtime_error(
                        "Error when making dsm_velocity_based_leg, number of extra free parameters should be three, was " +
                        std::to_string( legFreeParameters.rows( ) ) );
        }
        if( previousTransferLeg == nullptr )
        {
            throw std::runtime_error( "Error when making dsm_velocity_based_leg, no previous transfer leg found" );
        }
        transferLeg = std::make_shared< DsmVelocityBasedTransferLeg >(
                    departureBodyEphemeris, arrivalBodyEphemeris,
                    ( Eigen::VectorXd( 5 )<< legStartTime, legEndTime, legFreeParameters ).finished( ),
                    centralBodyGravitationalParameter );
        break;
    }
    default:
        throw std::runtime_error( "Error when making transfer leg, leg type not recognized" );
    }
    return transferLeg;
}

std::shared_ptr< TransferNode > createTransferNode(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< TransferNodeSettings > nodeSettings,
        const double nodeTime,
        const Eigen::VectorXd& nodeFreeParameters,
        const std::shared_ptr< TransferLeg > incomingTransferLeg,
        const std::shared_ptr< TransferLeg > outgoingTransferLeg )
{

    if( bodyMap.count( nodeSettings->nodeBodyName_ ) == 0 )
    {
        throw std::runtime_error( "Error when making transfer node, body " + nodeSettings->nodeBodyName_  + " not found." );
    }
    else if( bodyMap.at( nodeSettings->nodeBodyName_ )->getGravityFieldModel( ) == nullptr )
    {
        throw std::runtime_error( "Error when making transfer node, body " + nodeSettings->nodeBodyName_  + " has no gravity field." );
    }
    else if( bodyMap.at( nodeSettings->nodeBodyName_ )->getEphemeris( ) == nullptr )
    {
        throw std::runtime_error( "Error when making transfer node, body " + nodeSettings->nodeBodyName_  + " has no ephemeris." );
    }

    const double centralBodyGravitationalParameter =
            bodyMap.at( nodeSettings->nodeBodyName_ )->getGravityFieldModel( )->getGravitationalParameter( );
    const std::shared_ptr< ephemerides::Ephemeris > centralBodyEphemeris =
            bodyMap.at( nodeSettings->nodeBodyName_ )->getEphemeris( );

    std::shared_ptr< TransferNode > transferNode;
    switch( nodeSettings->nodeType_ )
    {
    case powered_swingby:
    {
        if( nodeFreeParameters.rows( ) != 0 )
        {
            throw std::runtime_error(
                        "Error when making powered_swingby, number of extra free parameters should be zero, was " +
                        std::to_string( nodeFreeParameters.rows( ) ) );
        }

        std::shared_ptr< PoweredSwingbyNodeSettings > swingbySettings =
                std::dynamic_pointer_cast< PoweredSwingbyNodeSettings >( nodeSettings );
        if( swingbySettings == nullptr )
        {
            throw std::runtime_error( "Error when making powered_swingby node, type is inconsistent" );
        }

        std::function< Eigen::Vector3d( ) > incomingVelocityFunction =
                std::bind( &TransferLeg::getArrivalVelocity, incomingTransferLeg );
        std::function< Eigen::Vector3d( ) > outgoingVelocityFunction =
                std::bind( &TransferLeg::getDepartureVelocity, outgoingTransferLeg );
        transferNode = std::make_shared< PoweredSwingbyTransferNode >(
                    centralBodyEphemeris, ( Eigen::VectorXd( 1 )<< nodeTime ).finished( ),
                    centralBodyGravitationalParameter, swingbySettings->minimumPeriapsisRadius_,
                    incomingVelocityFunction, outgoingVelocityFunction );


        break;
    }
    case escape_and_departure:
    {
        if( nodeFreeParameters.rows( ) != 0 )
        {
            throw std::runtime_error(
                        "Error when making powered_swingby, number of extra free parameters should be zero, was " +
                        std::to_string( nodeFreeParameters.rows( ) ) );
        }

        std::shared_ptr< EscapeAndDepartureNodeSettings > escapeAndDepartureSettings =
                std::dynamic_pointer_cast< EscapeAndDepartureNodeSettings >( nodeSettings );
        if( escapeAndDepartureSettings == nullptr )
        {
            throw std::runtime_error( "Error when making escape_and_departure node, type is inconsistent" );
        }


        std::function< Eigen::Vector3d( ) > outgoingVelocityFunction =
                std::bind( &TransferLeg::getDepartureVelocity, outgoingTransferLeg );
        transferNode = std::make_shared< EscapeAndDepartureTransferNode >(
                    centralBodyEphemeris, ( Eigen::VectorXd( 1 )<< nodeTime ).finished( ),
                    centralBodyGravitationalParameter,
                    escapeAndDepartureSettings->departureSemiMajorAxis_,
                    escapeAndDepartureSettings->departureEccentricity_,
                    escapeAndDepartureSettings->maximumExcessVelocity_,
                    outgoingVelocityFunction );
        break;
    }
    case capture_and_insertion:
    {
        if( nodeFreeParameters.rows( ) != 0 )
        {
            throw std::runtime_error(
                        "Error when making powered_swingby, number of extra free parameters should be zero, was " +
                        std::to_string( nodeFreeParameters.rows( ) ) );
        }

        std::shared_ptr< CaptureAndInsertionNodeSettings > captureAndInsertionSettings =
                std::dynamic_pointer_cast< CaptureAndInsertionNodeSettings >( nodeSettings );
        if( captureAndInsertionSettings == nullptr )
        {
            throw std::runtime_error( "Error when making capture_and_insertion node, type is inconsistent" );
        }


        std::function< Eigen::Vector3d( ) > outgoingVelocityFunction =
                std::bind( &TransferLeg::getDepartureVelocity, outgoingTransferLeg );
        transferNode = std::make_shared< CaptureAndInsertionTransferNode >(
                    centralBodyEphemeris, ( Eigen::VectorXd( 1 )<< nodeTime ).finished( ),
                    centralBodyGravitationalParameter,
                    captureAndInsertionSettings->departureSemiMajorAxis_,
                    captureAndInsertionSettings->departureEccentricity_,
                    outgoingVelocityFunction );
        break;
    }
    default:
        throw std::runtime_error( "Error when making transfer node, node type not recognized" );
    }
}

void createTransferTrajectory(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        const std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
        const std::vector< double >& nodeTimes,
        const std::vector< Eigen::VectorXd >& legFreeParameters,
        const std::vector< Eigen::VectorXd >& nodeFreeParameters )
{
    std::vector< std::shared_ptr< TransferLeg > > legs;

    for( int i = 0; i < legSettings.size( ); i++ )
    {
        legs.push_back(
                    createTransferLeg(
                        bodyMap, legSettings.at( i ), nodeTimes.at( i ), nodeTimes.at( i + 1 ),
                        legFreeParameters.at( i ) ) );
    }

    std::vector< std::shared_ptr< TransferNode > > nodes;

    std::shared_ptr< TransferNode > incomingTransferLeg = nullptr;
    std::shared_ptr< TransferNode > outgoingTransferLeg = nullptr;

    for( int i = 0; i < nodeSettings.size( ); i++ )
    {
        if( i != 0 )
        {
            incomingTransferLeg = legs.at( i - 1 );
        }
        else
        {
            incomingTransferLeg = nullptr;
        }

        if( i != nodeSettings.size( ) - 1 )
        {
            outgoingTransferLeg = legs.at( i );
        }
        else
        {
            outgoingTransferLeg = nullptr;
        }


        nodes.push_back(
                    createTransferNode(
                        bodyMap, nodeSettings.at( i ), nodeTimes.at( i ),
                        nodeFreeParameters.at( i ), incomingTransferLeg, outgoingTransferLeg ) );
    }

    return std::make_shared< TransferTrajectory >( legs, nodes );

}


} // namespace mission_segments

} // namespace tudat

#endif // TUDAT_CREATE_TRANSFER_TRAJECTORY_H
