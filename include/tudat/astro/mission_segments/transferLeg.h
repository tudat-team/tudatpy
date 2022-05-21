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

#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/astro/ephemerides/ephemeris.h"
#include "tudat/astro/mission_segments/lambertTargeterIzzo.h"
#include "tudat/basics/utilities.h"

namespace tudat
{
namespace mission_segments
{


enum TransferLegTypes
{
    unpowered_unperturbed_leg,
    dsm_position_based_leg,
    dsm_velocity_based_leg,
    hodographic_low_thrust_leg,
    spherical_shaping_low_thrust_leg
};

struct TrajectoryManeuver
{
public:

    TrajectoryManeuver(
            const Eigen::Vector3d position,
            const double velocityChange,
            const double time ):
    position_( position ), velocityChange_( velocityChange ), time_( time ){ }

    TrajectoryManeuver( ):
        position_( Eigen::Vector3d::Constant( TUDAT_NAN ) ),
        velocityChange_( TUDAT_NAN ),
        time_( TUDAT_NAN ){ }

    Eigen::Vector3d getPosition( )
    {
        return position_;
    }

    double getVelocityChange( )
    {
        return velocityChange_;
    }

    double getManeuverTime( )
    {
        return time_;
    }
private:
    Eigen::Vector3d position_;

    double velocityChange_;

    double time_;
};

class TransferLeg
{
public:
    TransferLeg(
            const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
            const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
            const TransferLegTypes legType );

    virtual ~TransferLeg( ){ }

    void updateLegParameters( const Eigen::VectorXd legParameters );

    double getLegDeltaV( );

    TransferLegTypes getTransferLegType( );

    Eigen::Vector3d getDepartureVelocity( );

    Eigen::Vector3d getArrivalVelocity( );

    virtual int getNumberOfImpulsiveManeuvers( )
    {
        return 0;
    }

    virtual const TrajectoryManeuver& getTrajectoryManeuver( const int maneuverIndex )
    {
        throw std::runtime_error( "Error, no maneuvers present in current leg" );
    }

    double getLegTimeOfFlight( )
    {
        return arrivalTime_ - departureTime_;
    }

    virtual void getStateAlongTrajectory( Eigen::Vector6d& stateAlongTrajectory,
                                          const double time ) = 0;

    Eigen::Vector6d getStateAlongTrajectory( const double time )
    {
        Eigen::Vector6d stateAlongTrajectory;
        getStateAlongTrajectory( stateAlongTrajectory, time );
        return stateAlongTrajectory;
    }

    void getStatesAlongTrajectory( std::map< double, Eigen::Vector6d >& statesAlongTrajectory,
                                   const std::vector< double >& times )
    {
        statesAlongTrajectory.clear( );
        Eigen::Vector6d currentStateAlongTrajectory;
        for( unsigned int i = 0; i < times.size( ); i++ )
        {
            getStateAlongTrajectory( currentStateAlongTrajectory, times.at( i ) );
            statesAlongTrajectory[ times.at( i ) ] = currentStateAlongTrajectory;
        }
    }

    void getStatesAlongTrajectory( std::map< double, Eigen::Vector6d >& statesAlongTrajectory,
                                   const int numberOfDataPoints )
    {
        statesAlongTrajectory.clear( );
        std::vector< double > times = utilities::linspace( departureTime_, arrivalTime_, numberOfDataPoints );
        getStatesAlongTrajectory( statesAlongTrajectory, times );
    }

    //! Get single value of inertial cartesian thrust acceleration.
    virtual void getThrustAccelerationAlongTrajectory ( Eigen::Vector3d& thrustAccelerationAlongTrajectory,
                                                        const double time )
    {
        thrustAccelerationAlongTrajectory = ( Eigen::Vector3d() << 0.0, 0.0, 0.0 ).finished();
    }

    //! Get single value of cartesian thrust acceleration.
    Eigen::Vector3d getThrustAccelerationAlongTrajectory ( const double time )
    {
        Eigen::Vector3d thrustAccelerationAlongTrajectory;
        getThrustAccelerationAlongTrajectory( thrustAccelerationAlongTrajectory, time );
        return thrustAccelerationAlongTrajectory;
    }

    //! Get inertial cartesian thrust acceleration along transfer leg.
    void getThrustAccelerationsAlongTrajectory(std::map< double, Eigen::Vector3d >& thrustAccelerationsAlongTrajectory,
                                               const std::vector< double >& times )
    {
        thrustAccelerationsAlongTrajectory.clear( );
        Eigen::Vector3d currentThrustAccelerationAlongTrajectory;
        for( unsigned int i = 0; i < times.size( ); i++ )
        {
            getThrustAccelerationAlongTrajectory(currentThrustAccelerationAlongTrajectory, times.at(i));
            thrustAccelerationsAlongTrajectory[ times.at(i ) ] = currentThrustAccelerationAlongTrajectory;
        }
    }

    //! Get inertial cartesian thrust acceleration along transfer leg.
    void getThrustAccelerationsAlongTrajectory(std::map< double, Eigen::Vector3d >& thrustAccelerationsAlongTrajectory,
                                               const int numberOfDataPoints )
    {
        thrustAccelerationsAlongTrajectory.clear( );
        std::vector< double > times = utilities::linspace( departureTime_, arrivalTime_, numberOfDataPoints );
        getThrustAccelerationsAlongTrajectory(thrustAccelerationsAlongTrajectory, times);
    }

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
    double timeOfFlight_;

    Eigen::Vector6d departureBodyState_;
    Eigen::Vector6d arrivalBodyState_;
    Eigen::Vector3d departureVelocity_;
    Eigen::Vector3d arrivalVelocity_;

    double legTotalDeltaV_;

};

class UnpoweredUnperturbedTransferLeg : public TransferLeg
{
public:
    using TransferLeg::getStateAlongTrajectory;

    UnpoweredUnperturbedTransferLeg(
            const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
            const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
            const double centralBodyGravitationalParameter );

    virtual ~UnpoweredUnperturbedTransferLeg( ){ }

    void getStateAlongTrajectory( Eigen::Vector6d& stateAlongTrajectory,
                                  const double time );

protected:

    void computeTransfer( );

    double centralBodyGravitationalParameter_;

    Eigen::Vector6d constantKeplerianState_;
};



class DsmTransferLeg : public TransferLeg
{
public:

    DsmTransferLeg(
            const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
            const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
            const TransferLegTypes legType,
            const double centralBodyGravitationalParameter ):
        TransferLeg( departureBodyEphemeris, arrivalBodyEphemeris, legType ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter){ }

    virtual ~DsmTransferLeg( ){ }

    void getStateAlongTrajectory( Eigen::Vector6d& stateAlongTrajectory,
                                  const double time );

    double getDsmTime( )
    {
        return trajectoryManeuver_.getManeuverTime( );
    }

    Eigen::Vector3d getDsmLocation( )
    {
        return trajectoryManeuver_.getPosition( );
    }

    int getNumberOfImpulsiveManeuvers( )
    {
        return 1;
    }

    const TrajectoryManeuver& getTrajectoryManeuver( const int maneuverIndex )
    {
        if( maneuverIndex == 0 )
        {
            return trajectoryManeuver_;
        }
        else
        {
            throw std::runtime_error( "Error when retrieving maneuver from DSM leg, only one maneuver is available" );
        }
    }

protected:

    void calculateKeplerianElements( );

    double centralBodyGravitationalParameter_;
    Eigen::Vector3d velocityBeforeDsm_;
    Eigen::Vector3d velocityAfterDsm_;

    TrajectoryManeuver trajectoryManeuver_;

    Eigen::Vector6d constantKeplerianStateBeforeDsm_;
    Eigen::Vector6d constantKeplerianStateAfterDsm_;

};


class DsmPositionBasedTransferLeg : public DsmTransferLeg
{
public:
    using TransferLeg::getStateAlongTrajectory;

    DsmPositionBasedTransferLeg(
            const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
            const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
            const double centralBodyGravitationalParameter );

    virtual ~DsmPositionBasedTransferLeg( ){ }

protected:

    void computeTransfer( );

//    void calculateDsmLocation( );

    double dsmTimeOfFlightFraction_;
    double dimensionlessRadiusDsm_;
    double inPlaneAngle_;
    double outOfPlaneAngle_;
};

class DsmVelocityBasedTransferLeg : public DsmTransferLeg
{
public:
    using TransferLeg::getStateAlongTrajectory;

    DsmVelocityBasedTransferLeg(
            const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
            const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
            const double centralBodyGravitationalParameter,
            const std::function< Eigen::Vector3d( ) > departureVelocityFunction );

    virtual ~DsmVelocityBasedTransferLeg( ){ }

protected:

    void computeTransfer( );

    std::function< Eigen::Vector3d( ) > departureVelocityFunction_;

    double dsmTimeOfFlightFraction_;

};


} // namespace mission_segments

} // namespace tudat

#endif // TUDAT_TRANSFER_LEG_H
