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
    swingby,
    escape_and_departure,
    capture_and_insertion
};


class TransferNode
{
public:
    TransferNode(
            const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
            const TransferNodeTypes nodeType );

    virtual ~TransferNode( ){ }

    void updateNodeParameters( const Eigen::VectorXd nodeParameters );

    double getNodeDeltaV( );

    TransferNodeTypes getTransferNodeType( );

    virtual bool nodeComputesOutgoingVelocity( ) = 0;

    virtual bool nodeComputesIncomingVelocity ( ) = 0;

    virtual Eigen::Vector3d getIncomingVelocity( );

    virtual Eigen::Vector3d getOutgoingVelocity( );

protected:

    void updateNodeState( const double nodeTime );

    virtual void computeNode( ) = 0;

    double nodeTime_;

    // Constant inputs
    std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris_;
    TransferNodeTypes nodeType_;

    // Input modified per iteration
    Eigen::VectorXd nodeParameters_;

    // Values computed per iteration
    Eigen::Vector3d incomingVelocity_;
    Eigen::Vector3d outgoingVelocity_;
    double totalNodeDeltaV_;

    Eigen::Vector6d nodeState_;

};


class DepartureWithFixedOutgoingVelocityNode: public TransferNode
{
public:
    DepartureWithFixedOutgoingVelocityNode(
            const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
            const double centralBodyGravitationalParameter,
            const double departureSemiMajorAxis,
            const double departureEccentricity,
            const std::function< Eigen::Vector3d( ) > outgoingVelocityFunction );

    virtual ~DepartureWithFixedOutgoingVelocityNode( ){ }

    Eigen::Vector3d getIncomingVelocity( );

    bool nodeComputesOutgoingVelocity( );

    bool nodeComputesIncomingVelocity( );

protected:

    void computeNode( );

    // Constant inputs
    double centralBodyGravitationalParameter_;
    double departureSemiMajorAxis_;
    double departureEccentricity_;
    std::function< Eigen::Vector3d( ) > outgoingVelocityFunction_;
};


class DepartureWithFreeOutgoingVelocityNode: public TransferNode
{
public:
    DepartureWithFreeOutgoingVelocityNode(
            const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
            const double centralBodyGravitationalParameter,
            const double departureSemiMajorAxis,
            const double departureEccentricity );

    virtual ~DepartureWithFreeOutgoingVelocityNode( ){ }

    Eigen::Vector3d getIncomingVelocity( );

    bool nodeComputesOutgoingVelocity( );

    bool nodeComputesIncomingVelocity( );

protected:

    void computeNode( );


    // Constant inputs
    double centralBodyGravitationalParameter_;
    double departureSemiMajorAxis_;
    double departureEccentricity_;

    // Input modified per iteration (extracted from nodeParameters_)
    double outgoingExcessVelocityMagnitude_;
    double outgoingExcessVelocityInPlaneAngle_;
    double outgoingExcessVelocityOutOfPlaneAngle_;

    // Values computed per iteration
    Eigen::Vector3d centralBodyPosition_;

};


class CaptureWithFixedIncomingVelocityNode: public TransferNode
{
public:
    CaptureWithFixedIncomingVelocityNode(
            const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
            const double centralBodyGravitationalParameter,
            const double captureSemiMajorAxis,
            const double captureEccentricity,
            const std::function< Eigen::Vector3d( ) > incomingVelocityFunction );

    virtual ~CaptureWithFixedIncomingVelocityNode( ){ }

    Eigen::Vector3d getOutgoingVelocity( );

    bool nodeComputesOutgoingVelocity( );

    bool nodeComputesIncomingVelocity( );

protected:

    void computeNode( );

    // Constant inputs
    double centralBodyGravitationalParameter_;
    double captureSemiMajorAxis_;
    double captureEccentricity_;
    std::function< Eigen::Vector3d( ) > incomingVelocityFunction_;
};


class CaptureWithFreeIncomingVelocityNode: public TransferNode
{
public:
    CaptureWithFreeIncomingVelocityNode(
            const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
            const double centralBodyGravitationalParameter,
            const double captureSemiMajorAxis,
            const double captureEccentricity);

    virtual ~CaptureWithFreeIncomingVelocityNode( ){ }

    Eigen::Vector3d getOutgoingVelocity( );

    bool nodeComputesOutgoingVelocity( );

    bool nodeComputesIncomingVelocity( );

protected:

    void computeNode( );

    // Constant inputs
    double centralBodyGravitationalParameter_;
    double captureSemiMajorAxis_;
    double captureEccentricity_;

    // Input modified per iteration (extracted from nodeParameters_)
    double incomingExcessVelocityMagnitude_;
    double incomingExcessVelocityInPlaneAngle_;
    double incomingExcessVelocityOutOfPlaneAngle_;
};


class SwingbyWithFixedIncomingFixedOutgoingVelocity: public TransferNode
{
public:
    SwingbyWithFixedIncomingFixedOutgoingVelocity(
            const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
            const double centralBodyGravitationalParameter,
            const double minimumPeriapsisRadius,
            const std::function< Eigen::Vector3d( ) > incomingVelocityFunction,
            const std::function< Eigen::Vector3d( ) > outgoingVelocityFunction );

    virtual ~SwingbyWithFixedIncomingFixedOutgoingVelocity( ){ }

    bool nodeComputesOutgoingVelocity( );

    bool nodeComputesIncomingVelocity( );

protected:

    void computeNode( );

    // Constant inputs
    double centralBodyGravitationalParameter_;
    double minimumPeriapsisRadius_;
    const std::function< Eigen::Vector3d( ) > incomingVelocityFunction_;
    const std::function< Eigen::Vector3d( ) > outgoingVelocityFunction_;
};


class SwingbyWithFixedIncomingFreeOutgoingVelocity: public TransferNode
{
public:
    SwingbyWithFixedIncomingFreeOutgoingVelocity(
            const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
            const double centralBodyGravitationalParameter,
            const std::function< Eigen::Vector3d( ) > incomingVelocityFunction );

    virtual ~SwingbyWithFixedIncomingFreeOutgoingVelocity( ){ }

    bool nodeComputesOutgoingVelocity( );

    bool nodeComputesIncomingVelocity( );

protected:

    void computeNode( );

    // Constant inputs
    double centralBodyGravitationalParameter_;
    const std::function< Eigen::Vector3d( ) > incomingVelocityFunction_;

    // Input modified per iteration (extracted from nodeParameters_)
    double periapsisRadius_;
    double swingbyDeltaV_;
    double outgoingRotationAngle_;
};


class SwingbyWithFreeIncomingFreeOutgoingVelocity: public TransferNode
{
public:
   SwingbyWithFreeIncomingFreeOutgoingVelocity(
           const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
           const double centralBodyGravitationalParameter);

   virtual ~SwingbyWithFreeIncomingFreeOutgoingVelocity( ){ }

   bool nodeComputesOutgoingVelocity( );

   bool nodeComputesIncomingVelocity( );

protected:

   void computeNode( );

   // Constant inputs
   double centralBodyGravitationalParameter_;
   const std::function< Eigen::Vector3d( ) > incomingVelocityFunction_;

   // Input modified per iteration (extracted from nodeParameters_)
   double incomingExcessVelocityMagnitude_;
   double incomingExcessVelocityInPlaneAngle_;
   double incomingExcessVelocityOutOfPlaneAngle_;
   double periapsisRadius_;
   double swingbyDeltaV_;
   double outgoingRotationAngle_;
};


class SwingbyWithFreeIncomingFixedOutgoingVelocity: public TransferNode
{
public:
   SwingbyWithFreeIncomingFixedOutgoingVelocity(
           const std::shared_ptr< ephemerides::Ephemeris > nodeEphemeris,
           const double centralBodyGravitationalParameter,
           const std::function< Eigen::Vector3d( ) > outgoingVelocityFunction );

   virtual ~SwingbyWithFreeIncomingFixedOutgoingVelocity( ){ }

   bool nodeComputesOutgoingVelocity( );

   bool nodeComputesIncomingVelocity( );

protected:

   void computeNode( );

   // Constant inputs
   double centralBodyGravitationalParameter_;
   const std::function< Eigen::Vector3d( ) > outgoingVelocityFunction_;

   // Input modified per iteration (extracted from nodeParameters_)
   double periapsisRadius_;
   double swingbyDeltaV_;
   double incomingRotationAngle_;
};


} // namespace mission_segments

} // namespace tudat

#endif // TUDAT_TRANSFER_NODE_H
