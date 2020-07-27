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

#ifndef TUDAT_TRANSFER_LEG_H
#define TUDAT_TRANSFER_LEG_H

#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "tudat/astro/ephemerides/ephemeris.h"
#include "tudat/astro/mission_segments/lambertTargeterIzzo.h"


namespace tudat
{
namespace mission_segments
{


enum TransferLegTypes
{
    unpowered_unperturbed_leg,
    dsm_position_based_leg,
    dsm_velocity_based_leg
};


class TransferLeg
{
public:
    TransferLeg(
            const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
            const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
            const TransferLegTypes legType,
            const Eigen::VectorXd legParameters );

    void updateLegParameters( const Eigen::VectorXd legParameters );

    double getLegDeltaV( );

    TransferLegTypes getTransferLegType( );

//    virtual bool departureVelocityIsPredetermined( ) = 0;

    Eigen::Vector3d getDepartureVelocity( );

    Eigen::Vector3d getArrivalVelocity( );

protected:

    virtual void computeTransfer( ) = 0;

    void updateDepartureAndArrivalBodies(
            const double departureTime,
            const double arrivalTime );

    std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris_;
    std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris_;
    TransferLegTypes legType_;

    Eigen::VectorXd legParameters_;

    double departureTime_;
    double arrivalTime_;

    Eigen::Vector6d departureBodyState_;
    Eigen::Vector6d arrivalBodyState_;
    Eigen::Vector3d departureVelocity_;
    Eigen::Vector3d arrivalVelocity_;

    double legTotalDeltaV_;

};

class UnpoweredUnperturbedTransferLeg : public TransferLeg
{
public:
    UnpoweredUnperturbedTransferLeg(
            const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
            const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
            const Eigen::VectorXd legParameters,
            const double centralBodyGravitationalParameter );

protected:

    virtual void computeTransfer( );

    double centralBodyGravitationalParameter_;
};



class DsmPositionBasedTransferLeg : public TransferLeg
{
public:
    DsmPositionBasedTransferLeg(
            const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
            const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
            const Eigen::VectorXd legParameters,
            const double centralBodyGravitationalParameter );

protected:

    virtual void computeTransfer( );

    void calculateDsmLocation( );

    double centralBodyGravitationalParameter_;

    double dsmTimeOfFlightFraction_;
    double dimensionlessRadiusDsm_;
    double inPlaneAngle_;
    double outOfPlaneAngle_;


    Eigen::Vector3d dsmLocation_;
    Eigen::Vector3d velocityBeforeDsm_;
    Eigen::Vector3d velocityAfterDsm_;
    double dsmTime_;
    double dsmDeltaV_;
};


class DsmVelocityBasedTransferLeg : public TransferLeg
{
public:
    DsmVelocityBasedTransferLeg(
            const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
            const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
            const Eigen::VectorXd legParameters,
            const double centralBodyGravitationalParameter,
            const std::function< Eigen::Vector3d( ) > departureVelocityFunction );

protected:

    virtual void computeTransfer( );

    double centralBodyGravitationalParameter_;
    std::function< Eigen::Vector3d( ) > departureVelocityFunction_;

    double dsmTimeOfFlightFraction_;
    double excessVelocityMagnitude_;
    double excessVelocityInPlaneAngle_;
    double excessVelocityOutOfPlaneAngle_;

    Eigen::Vector3d dsmLocation_;
    Eigen::Vector3d velocityBeforeDsm_;
    Eigen::Vector3d velocityAfterDsm_;
    double dsmTime_;
    double dsmDeltaV_;
};


} // namespace mission_segments

} // namespace tudat

#endif // TUDAT_TRANSFER_LEG_H
