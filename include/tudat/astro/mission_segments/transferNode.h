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

#ifndef TUDAT_TRANSFER_NODE_H
#define TUDAT_TRANSFER_NODE_H

#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "tudat/astro/ephemerides/ephemeris.h"
#include "tudat/astro/mission_segments/escapeAndCapture.h"
#include "tudat/astro/mission_segments/gravityAssist.h"

namespace tudat
{
namespace mission_segments
{


enum TransferNodeTypes
{
    swingby_free_outgoing_velocity,
    swingby_fixed_outgoing_velocity,
    escape_and_departure,
    capture_and_insertion
};

class TransferNode
{
public:
    TransferNode(
            const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
            const Eigen::VectorXd nodeParameters ):
        nodeEphemeris_( nodeEphemeris ), nodeParameters_( nodeParameters_ ){ }

    void updateNodeParameters( const Eigen::VectorXd nodeParameters )
    {
        nodeParameters_ = nodeParameters;
        computeNode( );
    }

    virtual double getNodeDeltaV( ) = 0;

    virtual TransferNodeTypes getTransferNodeType( ) = 0;

    virtual Eigen::Vector3d getIncomingVelocity( ) = 0;

    virtual Eigen::Vector3d getOutgoingVelocity( ) = 0;

protected:

    virtual void computeNode( ) = 0;

    std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris_;

    Eigen::VectorXd nodeParameters_;

};

class EscapeAndDepartureTransferNode: public TransferNode
{
public:
    EscapeAndDepartureTransferNode(
            const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
            const Eigen::VectorXd nodeParameters,
            const double centralBodyGravitationalParameter,
            const double departureSemiMajorAxis,
            const double departureEccentricity,
            const double maximumExcessVelocity,
            const std::function< Eigen::Vector3d( ) > outgoingVelocityFunction ):
        TransferNode( nodeEphemeris, nodeParameters ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
        departureSemiMajorAxis_( departureSemiMajorAxis ),
        departureEccentricity_( departureEccentricity ),
        maximumExcessVelocity_( maximumExcessVelocity ),
        outgoingVelocityFunction_( outgoingVelocityFunction ){ }

    double getNodeDeltaV( )
    {
        return escapeDeltaV_;
    }

    TransferNodeTypes getTransferNodeType( )
    {
        return escape_and_departure;
    }

    Eigen::Vector3d getIncomingVelocity( )
    {
        throw std::runtime_error( "Error, no incoming velocity can be given for escape_and_departure leg" );
    }

    Eigen::Vector3d getOutgoingVelocity( )
    {
        return outgoingVelocity_;
    }
protected:

    void computeNode( )
    {
        if( nodeParameters_.rows( ) != 1 )
        {
            throw std::runtime_error( "Error when computing PoweredSwingbyTransferNode, incorrect input size" );
        }

        nodeTime_ = nodeParameters_( 0 );

        outgoingVelocity_ = outgoingVelocityFunction_( );
        centralBodyVelocity_ = nodeEphemeris_->getCartesianState( nodeTime_ ).segment( 3, 3 );

        double totalRequiredDeltaV = mission_segments::computeEscapeOrCaptureDeltaV(
                    centralBodyGravitationalParameter_,
                    departureSemiMajorAxis_, departureEccentricity_,
                    ( outgoingVelocity_ - centralBodyVelocity_ ).norm( ) );

        if( totalRequiredDeltaV < maximumExcessVelocity_ )
        {
            escapeDeltaV_ = 0.0;
        }
        else
        {
            escapeDeltaV_ = totalRequiredDeltaV - maximumExcessVelocity_;
        }
    }

    double nodeTime_;


    double centralBodyGravitationalParameter_;

    double departureSemiMajorAxis_;

    double departureEccentricity_;

    double maximumExcessVelocity_;

    std::function< Eigen::Vector3d( ) > outgoingVelocityFunction_;



    Eigen::Vector3d outgoingVelocity_;

    Eigen::Vector3d centralBodyVelocity_;

    double escapeDeltaV_;
};

class CaptureAndInsertionTransferNode: public TransferNode
{
    CaptureAndInsertionTransferNode(
            const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
            const Eigen::VectorXd nodeParameters,
            const double centralBodyGravitationalParameter,
            const double captureSemiMajorAxis,
            const double captureEccentricity,
            const std::function< Eigen::Vector3d( ) > incomingVelocityFunction ):
        TransferNode( nodeEphemeris, nodeParameters ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
        captureSemiMajorAxis_( captureSemiMajorAxis ),
        captureEccentricity_( captureEccentricity ),
        incomingVelocityFunction_( incomingVelocityFunction ){ }

    double getNodeDeltaV( )
    {
        return captureAndInsertionDeltaV_;
    }

    TransferNodeTypes getTransferNodeType( )
    {
        return capture_and_insertion;
    }

    Eigen::Vector3d getIncomingVelocity( )
    {
        return incomingVelocity_;
    }

    Eigen::Vector3d getOutgoingVelocity( )
    {
        throw std::runtime_error( "Error, no outgoing velocity can be given for capture_and_insertion leg" );
    }

protected:

    void computeNode( )
    {
        if( nodeParameters_.rows( ) != 1 )
        {
            throw std::runtime_error( "Error when computing PoweredSwingbyTransferNode, incorrect input size" );
        }

        nodeTime_ = nodeParameters_( 0 );

        incomingVelocity_ = incomingVelocityFunction_( );
        centralBodyVelocity_ = nodeEphemeris_->getCartesianState( nodeTime_ ).segment( 3, 3 );

        captureAndInsertionDeltaV_ = mission_segments::computeEscapeOrCaptureDeltaV(
                    centralBodyGravitationalParameter_,
                    captureSemiMajorAxis_, captureEccentricity_,
                    ( incomingVelocity_ - centralBodyVelocity_ ).norm( ) );
    }

    double nodeTime_;


    double centralBodyGravitationalParameter_;

    double captureSemiMajorAxis_;

    double captureEccentricity_;

    std::function< Eigen::Vector3d( ) > incomingVelocityFunction_;


    Eigen::Vector3d incomingVelocity_;

    Eigen::Vector3d centralBodyVelocity_;

    double captureAndInsertionDeltaV_;
};

class SwingbyWithFixedOutgoingVelocity: public TransferNode
{
public:
    SwingbyWithFixedOutgoingVelocity(
            const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
            const Eigen::VectorXd nodeParameters,
            const double centralBodyGravitationalParameter,
            const double minimumPeriapsisRadius,
            const std::function< Eigen::Vector3d( ) > incomingVelocityFunction,
            const std::function< Eigen::Vector3d( ) > outgoingVelocityFunction ):
        TransferNode( nodeEphemeris, nodeParameters ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
        minimumPeriapsisRadius_( minimumPeriapsisRadius ),
        incomingVelocityFunction_( incomingVelocityFunction ),
        outgoingVelocityFunction_( outgoingVelocityFunction )
    {
        computeNode( );
    }

    double getNodeDeltaV( )
    {
        return swingbyDeltaV_;
    }

    TransferNodeTypes getTransferNodeType( )
    {
        return swingby_fixed_outgoing_velocity;
    }

    Eigen::Vector3d getIncomingVelocity( )
    {
        return incomingVelocity_;
    }

    Eigen::Vector3d getOutgoingVelocity( )
    {
        return outgoingVelocity_;
    }

protected:

    void computeNode( )
    {
        if( nodeParameters_.rows( ) != 1 )
        {
            throw std::runtime_error( "Error when computing PoweredSwingbyTransferNode, incorrect input size" );
        }

        nodeTime_ = nodeParameters_( 0 );
        incomingVelocity_ = incomingVelocityFunction_( );
        outgoingVelocity_ = outgoingVelocityFunction_( );

        centralBodyVelocity_ = nodeEphemeris_->getCartesianState( nodeTime_ ).segment( 3, 3 );

        swingbyDeltaV_ = calculateGravityAssistDeltaV(
                    centralBodyGravitationalParameter_,
                    centralBodyVelocity_,
                    incomingVelocity_,
                    outgoingVelocity_,
                    minimumPeriapsisRadius_ );

    }

    double centralBodyGravitationalParameter_;

    double minimumPeriapsisRadius_;

    const std::function< Eigen::Vector3d( ) > incomingVelocityFunction_;

    const std::function< Eigen::Vector3d( ) > outgoingVelocityFunction_;


    double nodeTime_;

    Eigen::Vector3d incomingVelocity_;

    Eigen::Vector3d outgoingVelocity_;


    Eigen::Vector3d centralBodyVelocity_;

    double swingbyDeltaV_;

};

class SwingbyWithFreeOutgoingVelocity: public TransferNode
{
public:
    SwingbyWithFreeOutgoingVelocity(
            const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
            const Eigen::VectorXd nodeParameters,
            const double centralBodyGravitationalParameter,
            const std::function< Eigen::Vector3d( ) > incomingVelocityFunction ):
        PoweredSwingbyTransferNode( nodeEphemeris, nodeParameters ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
        incomingVelocityFunction_( incomingVelocityFunction )
    {
        computeNode( );
    }

    double getNodeDeltaV( )
    {
        return swingbyDeltaV_;
    }

    TransferNodeTypes getTransferNodeType( )
    {
        return swingby_free_outgoing_velocity;
    }

    Eigen::Vector3d getIncomingVelocity( )
    {
        return incomingVelocity_;
    }

    Eigen::Vector3d getOutgoingVelocity( )
    {
        return outgoingVelocity_;
    }

protected:

    void computeNode( )
    {
        if( nodeParameters_.rows( ) != 4 )
        {
            throw std::runtime_error( "Error when computing PoweredSwingbyTransferNode, incorrect input size" );
        }

        nodeTime_ = nodeParameters_( 0 );
        periapsisRadius_ = nodeParameters_( 1 );
        outgoingRotationAngle_ = nodeParameters_( 2 );
        swingbyDeltaV_ = nodeParameters_( 3 );

        incomingVelocity_ = incomingVelocityFunction_( );

        centralBodyVelocity_ = nodeEphemeris_->getCartesianState( nodeTime_ ).segment( 3, 3 );

        // Prepare the gravity assist propagator module.
        outgoingVelocity_ = mission_segments::calculatePoweredGravityAssistOutgoingVelocity(
                    centralBodyGravitationalParameter_,
                    centralBodyVelocity_, incomingVelocity_,
                    outgoingRotationAngle_, periapsisRadius_ );

    }

    double centralBodyGravitationalParameter_;

    const std::function< Eigen::Vector3d( ) > incomingVelocityFunction_;



    double nodeTime_;

    double periapsisRadius_;

    double swingbyDeltaV_;

    double outgoingRotationAngle_;


    Eigen::Vector3d incomingVelocity_;

    Eigen::Vector3d centralBodyVelocity_;


    Eigen::Vector3d outgoingVelocity_;

};

} // namespace mission_segments

} // namespace tudat

#endif // TUDAT_TRANSFER_NODE_H
