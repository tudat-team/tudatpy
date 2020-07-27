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


enum TransferLegTypes
{
    unpowered_unperturbed_leg,
    dsm_position_based_leg,
    dsm_velocity_based_leg
};

namespace tudat
{
namespace mission_segments
{

class TransferLeg
{
public:
    TransferLeg(
            const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
            const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
            const Eigen::VectorXd legParameters ):
        departureBodyEphemeris_( departureBodyEphemeris ), arrivalBodyEphemeris_( arrivalBodyEphemeris ),
        legParameters_( legParameters )
    { }

    void updateLegParameters( const Eigen::VectorXd legParameters )
    {
        legParameters_ = legParameters;
        computeTransfer( );
    }

    virtual double getLegDeltaV( ) = 0;

    virtual TransferLegTypes getTransferLegType( ) = 0;

    virtual bool departureVelocityIsPredetermined( ) = 0;

    Eigen::Vector3d getDepartureVelocity( )
    {
        return departureVelocity_;
    }

    Eigen::Vector3d getArrivalVelocity( )
    {
        return arrivalVelocity_;
    }

protected:

    virtual void computeTransfer( ) = 0;

    void updateDepartureAndArrivalBodies(
            const double departureTime,
            const double arrivalTime )
    {
        departureTime_ = departureTime;
        arrivalTime_ = arrivalTime;
        departureBodyState_ = departureBodyEphemeris_->getCartesianState( departureTime_ );
        arrivalBodyState_ = arrivalBodyEphemeris_->getCartesianState( arrivalTime );
    }

    std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris_;

    std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris_;

    Eigen::VectorXd legParameters_;

    double departureTime_;

    double arrivalTime_;

    Eigen::Vector6d departureBodyState_;

    Eigen::Vector6d arrivalBodyState_;

    Eigen::Vector3d departureVelocity_;

    Eigen::Vector3d arrivalVelocity_;

};

class UnpoweredUnperturbedTransferLeg : public TransferLeg
{
public:
    UnpoweredUnperturbedTransferLeg(
            const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
            const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
            const Eigen::VectorXd legParameters,
            const double centralBodyGravitationalParameter ):
        TransferLeg( departureBodyEphemeris, arrivalBodyEphemeris, legParameters ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter )
    {
        computeTransfer( );
    }

    double getLegDeltaV( )
    {
        return 0.0;
    }

    TransferLegTypes getTransferLegType( )
    {
        return unpowered_unperturbed_leg;
    }

    bool departureVelocityIsPredetermined( )
    {
        return false;
    }

protected:

    virtual void computeTransfer( )
    {
        if( legParameters_.rows( ) != 2 )
        {
            throw std::runtime_error( "Error when computing UnpoweredUnperturbedTransferLeg, incorrect input size" );
        }
        updateDepartureAndArrivalBodies(
                    legParameters_( 0 ), legParameters_( 1 ) );

        // Calculate and set the spacecraft velocities after departure and before arrival.
        mission_segments::solveLambertProblemIzzo(
                    departureBodyState_.segment( 0, 3 ), arrivalBodyState_.segment( 0, 3 ),
                    legParameters_( 1 ) - legParameters_( 0 ), centralBodyGravitationalParameter_,
                    departureVelocity_, arrivalVelocity_ );
    }

    double centralBodyGravitationalParameter_;

};



class DsmPositionBasedTransferLeg : public TransferLeg
{
public:
    DsmPositionBasedTransferLeg(
            const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
            const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
            const Eigen::VectorXd legParameters,
            const double centralBodyGravitationalParameter ):
        TransferLeg( departureBodyEphemeris, arrivalBodyEphemeris, legParameters ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter )
    {
        computeTransfer( );
    }

    double getLegDeltaV( )
    {
        return dsmDeltaV_;
    }    

    TransferLegTypes getTransferLegType( )
    {
        return dsm_position_based_leg;
    }

    bool departureVelocityIsPredetermined( )
    {
        return false;
    }

protected:

    virtual void computeTransfer( )
    {
        if( legParameters_.rows( ) != 2 )
        {
            throw std::runtime_error( "Error when computing UnpoweredUnperturbedTransferLeg, incorrect input size" );
        }
        updateDepartureAndArrivalBodies(
                    legParameters_( 0 ), legParameters_( 1 ) );

        double timeOfFlight = arrivalTime_ - departureTime_;
        dsmTimeOfFlightFraction_ = legParameters_( 2 );
        dimensionlessRadiusDsm_ = legParameters_( 3 );
        inPlaneAngle_ = legParameters_( 4 );
        outOfPlaneAngle_ = legParameters_( 5 );

        calculateDsmLocation( );


        // Calculate the DSM time of application from the time of flight fraction.
        dsmTime_ = dsmTimeOfFlightFraction_ * timeOfFlight;

        // Calculate and set the spacecraft velocities after departure, before and after the DSM, and
        // before arrival using two lambert targeters and all the corresponding positions and flight
        // times.
        mission_segments::solveLambertProblemIzzo( departureBodyState_.segment( 0, 3 ), dsmLocation_, dsmTime_,
                                                   centralBodyGravitationalParameter_,
                                                   departureVelocity_, velocityBeforeDsm_ );
        mission_segments::solveLambertProblemIzzo( dsmLocation_, arrivalBodyState_.segment( 0, 3 ), timeOfFlight -
                                                   dsmTime_, centralBodyGravitationalParameter_,
                                                   velocityAfterDsm_, arrivalVelocity_ );


        //Calculate the deltaV originating from the DSM.
        dsmDeltaV_ = ( velocityAfterDsm_ - velocityBeforeDsm_ ).norm( );
    }

    //! Calculates the Dsm location
    void calculateDsmLocation( )
    {
        // Calculate the required unit vectors
        const Eigen::Vector3d unitVector1 = departureBodyState_.segment( 0, 3 ).normalized( );
        const Eigen::Vector3d unitVector3 = ( unitVector1.cross( departureBodyState_.segment( 3, 3 ) ) ).normalized( );
        const Eigen::Vector3d unitVector2 = unitVector3.cross( unitVector1 );

        // Calculate the absolute DSM radius.
        const double absoluteRadiusDsm = dimensionlessRadiusDsm_ * departureBodyState_.segment( 0, 3 ).norm( );

        // Calculate the radius in the central body reference frame.
        dsmLocation_ = std::cos( inPlaneAngle_ ) * std::cos( outOfPlaneAngle_ ) * absoluteRadiusDsm *
                unitVector1 +
                std::sin( inPlaneAngle_ ) * std::cos( outOfPlaneAngle_ ) * absoluteRadiusDsm *
                unitVector2 +
                std::sin( outOfPlaneAngle_ ) * absoluteRadiusDsm * unitVector3;
    }

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
            const double centralBodyGravitationalParameter ):
        TransferLeg( departureBodyEphemeris, arrivalBodyEphemeris, legParameters ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter )
    {
        computeTransfer( );
    }

    double getLegDeltaV( )
    {
        return dsmDeltaV_;
    }

    TransferLegTypes getTransferLegType( )
    {
        return dsm_velocity_based_leg;
    }

    bool departureVelocityIsPredetermined( )
    {
        return false;
    }

protected:

    virtual void computeTransfer( )
    {
        if( legParameters_.rows( ) != 2 )
        {
            throw std::runtime_error( "Error when computing UnpoweredUnperturbedTransferLeg, incorrect input size" );
        }
        updateDepartureAndArrivalBodies(
                    legParameters_( 0 ), legParameters_( 1 ) );

        double timeOfFlight = arrivalTime_ - departureTime_;
        dsmTimeOfFlightFraction_ = legParameters_( 2 );
        excessVelocityMagnitude_ = legParameters_( 3 );
        excessVelocityInPlaneAngle_ = legParameters_( 4 );
        excessVelocityOutOfPlaneAngle_ = legParameters_( 5 );

        // Calculate the DSM time of application from the time of flight fraction.
        dsmTime_ = dsmTimeOfFlightFraction_ * timeOfFlight;

        departureVelocity_ = departureBodyState_.segment( 3, 3, ) +
                excessVelocityMagnitude_ * std::cos( excessVelocityInPlaneAngle_ ) *
                std::cos( excessVelocityOutOfPlaneAngle_ ) * unitVector1 +
                excessVelocityMagnitude_ * std::sin( excessVelocityInPlaneAngle_ ) *
                std::cos( excessVelocityOutOfPlaneAngle_ ) * unitVector2 +
                excessVelocityMagnitude_ * std::sin( excessVelocityOutOfPlaneAngle_ ) * unitVector3;

        // Transfer the initial position and velocity into a vectorXd object with Cartesian
        // coordinates.
        Eigen::Vector6d departureCartesianElements, departureKeplerianElements;
        departureCartesianElements.segment( 0, 3 ) = departureBodyPosition_;
        departureCartesianElements.segment( 3, 3 ) = velocityAfterDeparture_;

        // Convert the cartesian elements into keplerian elements.
        departureKeplerianElements = orbital_element_conversions::convertCartesianToKeplerianElements(
                    departureCartesianElements, centralBodyGravitationalParameter_ );

        // Propagate the keplerian elements until the moment of application of the DSM.
        Eigen::Vector6d dsmArrivalCartesianElements, dsmArrivalKeplerianElements;
        dsmArrivalKeplerianElements = orbital_element_conversions::propagateKeplerOrbit( departureKeplerianElements,
                                                                                         dsmTime_, centralBodyGravitationalParameter_ );

        // Convert the keplerian elements back into Cartesian elements.
        dsmArrivalCartesianElements = orbital_element_conversions::convertKeplerianToCartesianElements(
                    dsmArrivalKeplerianElements, centralBodyGravitationalParameter_ );

        // Set the corresponding position and velocity vectors.
        dsmLocation_ = cartesianElements.segment( 0, 3 );
        velocityBeforeDsm_ = cartesianElements.segment( 3, 3 );

        // Calculate the velocities after the DSM and before the arrival body.
        mission_segments::solveLambertProblemIzzo( dsmLocation_, arrivalBodyState_.segment( 0, 3 ), timeOfFlight -
                                                   dsmTime_, centralBodyGravitationalParameter_,
                                                   velocityAfterDsm_, arrivalBodyState_.segment( 3, 3 ) );

        // Calculate the deltaV needed for the DSM.
        dsmDeltaV_ = ( velocityAfterDsm_ - velocityBeforeDsm_ ).norm( );

    }


    double centralBodyGravitationalParameter_;


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
