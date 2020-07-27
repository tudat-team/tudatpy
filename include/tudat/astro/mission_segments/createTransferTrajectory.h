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

static std::map< TransferLegTypes, bool > legRequiresInputFromNode =
{
    { unpowered_unperturbed_leg, false },
    { dsm_position_based_leg, false },
    { dsm_velocity_based_leg, true }
};


class TransferLegSettings
{
public:
    TransferLegSettings(
            const std::string& departureBodyName,
            const std::string& arrivalBodyName,
            const std::string& centralBodyName,
            const TransferLegTypes legType ):
        departureBodyName_( departureBodyName ),
        arrivalBodyName_( arrivalBodyName ),
        centralBodyName_( centralBodyName ),
        legType_( legType ){ }

    std::string departureBodyName_;
    std::string arrivalBodyName_;
    std::string centralBodyName_;

    TransferLegTypes legType_;
};


class TransferNodeSettings
{
public:
    TransferNodeSettings(
            const std::string& nodeBodyName,
            const TransferNodeTypes nodeType ):
        nodeBodyName_( nodeBodyName ), nodeType_( nodeType ){ }

    virtual ~TransferNodeSettings( ){ }

    std::string nodeBodyName_;
    TransferNodeTypes nodeType_;
};

class SwingbyNodeSettings: public TransferNodeSettings
{
public:
    SwingbyNodeSettings(
            const std::string& nodeBodyName,
            const TransferNodeTypes nodeType,
            const double minimumPeriapsisRadius ):
        TransferNodeSettings( nodeBodyName, nodeType ),
        minimumPeriapsisRadius_( minimumPeriapsisRadius ){ }

    std::string nodeBodyName_;
    TransferNodeTypes nodeType_;
    double minimumPeriapsisRadius_;
};

class EscapeAndDepartureNodeSettings: public TransferNodeSettings
{
public:
    EscapeAndDepartureNodeSettings(
            const std::string& nodeBodyName,
            const TransferNodeTypes nodeType,
            const double departureSemiMajorAxis,
            const double departureEccentricity ):
        TransferNodeSettings( nodeBodyName, nodeType ),
        departureSemiMajorAxis_( departureSemiMajorAxis ),
        departureEccentricity_( departureEccentricity ){ }

    std::string nodeBodyName_;
    TransferNodeTypes nodeType_;
    double departureSemiMajorAxis_;
    double departureEccentricity_;
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

    std::string nodeBodyName_;
    TransferNodeTypes nodeType_;
    double captureSemiMajorAxis_;
    double captureEccentricity_;
};

std::shared_ptr< TransferLeg > createTransferLeg(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< TransferLegSettings > legSettings,
        const double legStartTime,
        const double legEndTime,
        const Eigen::VectorXd& legFreeParameters,
        const std::shared_ptr< TransferNode > departureNode = nullptr );

std::shared_ptr< TransferNode > createTransferNode(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< TransferNodeSettings > nodeSettings,
        const double nodeTime,
        const Eigen::VectorXd& nodeFreeParameters,
        const std::shared_ptr< TransferLeg > incomingTransferLeg,
        const std::shared_ptr< TransferLeg > outgoingTransferLeg,
        const bool nodeComputesOutgoingVelocity );

std::shared_ptr< TransferTrajectory > createTransferTrajectory(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        const std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings,
        const std::vector< double >& nodeTimes,
        const std::vector< Eigen::VectorXd >& legFreeParameters,
        const std::vector< Eigen::VectorXd >& nodeFreeParameters );

void printFreeParameterDefinitions(
        const std::vector< std::shared_ptr< TransferLegSettings > >& legSettings,
        const std::vector< std::shared_ptr< TransferNodeSettings > >& nodeSettings );

} // namespace mission_segments

} // namespace tudat

#endif // TUDAT_CREATE_TRANSFER_TRAJECTORY_H
