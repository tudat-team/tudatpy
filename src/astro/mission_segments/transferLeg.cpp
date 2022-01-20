#include "tudat/astro/mission_segments/transferLeg.h"
#include "tudat/astro/mission_segments/lambertRoutines.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/keplerPropagator.h"

namespace tudat
{
namespace mission_segments
{

TransferLeg::TransferLeg(
        const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
        const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
        const TransferLegTypes legType ):
    departureBodyEphemeris_( departureBodyEphemeris ), arrivalBodyEphemeris_( arrivalBodyEphemeris ),
    legType_( legType ), legParameters_( Eigen::VectorXd::Zero( 0 ) ){ }

void TransferLeg::updateLegParameters( const Eigen::VectorXd legParameters )
{
    legParameters_ = legParameters;
    computeTransfer( );
}

double TransferLeg::getLegDeltaV( )
{
    return legTotalDeltaV_;
}

TransferLegTypes TransferLeg::getTransferLegType( )
{
    return legType_;
}

Eigen::Vector3d TransferLeg::getDepartureVelocity( )
{
    return departureVelocity_;
}

Eigen::Vector3d TransferLeg::getArrivalVelocity( )
{
    return arrivalVelocity_;
}

void TransferLeg::updateDepartureAndArrivalBodies(
        const double departureTime,
        const double arrivalTime )
{
    departureTime_ = departureTime;
    arrivalTime_ = arrivalTime;
    departureBodyState_ = departureBodyEphemeris_->getCartesianState( departureTime_ );
    arrivalBodyState_ = arrivalBodyEphemeris_->getCartesianState( arrivalTime_ );
}



UnpoweredUnperturbedTransferLeg::UnpoweredUnperturbedTransferLeg(
        const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
        const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
        const double centralBodyGravitationalParameter ):
    TransferLeg( departureBodyEphemeris, arrivalBodyEphemeris, unpowered_unperturbed_leg ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter )
{
    legTotalDeltaV_ = 0.0;
}

void UnpoweredUnperturbedTransferLeg::getStateAlongTrajectory( Eigen::Vector6d& stateAlongTrajectory,
                                                               const double time )
{
    if( time < departureTime_ || time > arrivalTime_ )
    {
        throw std::runtime_error( "Error when requesting state along unpowered unperturbed leg, requested time is outside bounds" );
    }

    stateAlongTrajectory = orbital_element_conversions::convertKeplerianToCartesianElements(
                orbital_element_conversions::propagateKeplerOrbit(
                    constantKeplerianState_, time - departureTime_, centralBodyGravitationalParameter_ ),
                centralBodyGravitationalParameter_ );
}


void UnpoweredUnperturbedTransferLeg::computeTransfer( )
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

    Eigen::Vector6d initialState;
    initialState.segment( 0, 3 ) = departureBodyState_.segment( 0, 3 );
    initialState.segment( 3, 3 ) = departureVelocity_;

    constantKeplerianState_ = orbital_element_conversions::convertCartesianToKeplerianElements(
                initialState, centralBodyGravitationalParameter_ );
}

void DsmTransferLeg::calculateKeplerianElements( )
{
    Eigen::Vector6d initialState;
    initialState.segment( 0, 3 ) = departureBodyState_.segment( 0, 3 );
    initialState.segment( 3, 3 ) = departureVelocity_;

    constantKeplerianStateBeforeDsm_ = orbital_element_conversions::convertCartesianToKeplerianElements(
                initialState, centralBodyGravitationalParameter_ );

    Eigen::Vector6d stateAfterDsm;
    stateAfterDsm.segment( 0, 3 ) = trajectoryManeuver_.getPosition( );
    stateAfterDsm.segment( 3, 3 ) = velocityAfterDsm_;

    constantKeplerianStateAfterDsm_ = orbital_element_conversions::convertCartesianToKeplerianElements(
                stateAfterDsm, centralBodyGravitationalParameter_ );
}


void DsmTransferLeg::getStateAlongTrajectory( Eigen::Vector6d& stateAlongTrajectory,
                                              const double time )
{
    if( time < departureTime_ || time > arrivalTime_ )
    {
        throw std::runtime_error( "Error when requesting state along unpowered unperturbed leg, requested time is outside bounds" );
    }

    if( time < trajectoryManeuver_.getManeuverTime( ) )
    {
        stateAlongTrajectory = orbital_element_conversions::convertKeplerianToCartesianElements(
                    orbital_element_conversions::propagateKeplerOrbit(
                        constantKeplerianStateBeforeDsm_, time - departureTime_, centralBodyGravitationalParameter_ ),
                    centralBodyGravitationalParameter_ );
    }
    else
    {
        stateAlongTrajectory = orbital_element_conversions::convertKeplerianToCartesianElements(
                    orbital_element_conversions::propagateKeplerOrbit(
                        constantKeplerianStateAfterDsm_, time - trajectoryManeuver_.getManeuverTime( ), centralBodyGravitationalParameter_ ),
                    centralBodyGravitationalParameter_ );
    }
}


DsmPositionBasedTransferLeg::DsmPositionBasedTransferLeg(
        const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
        const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
        const double centralBodyGravitationalParameter ):
    DsmTransferLeg( departureBodyEphemeris, arrivalBodyEphemeris, dsm_position_based_leg, centralBodyGravitationalParameter )
{ }



//! Calculates the Dsm location
Eigen::Vector3d calculatePositionBasedDsmLocation(
        const Eigen::Vector6d& departureBodyState,
        const double dimensionlessRadiusDsm,
        const double inPlaneAngle,
        const double outOfPlaneAngle )
{
    // Calculate the required unit vectors
    const Eigen::Vector3d unitVector1 = departureBodyState.segment( 0, 3 ).normalized( );
    const Eigen::Vector3d unitVector3 =
            ( unitVector1.cross(
                  Eigen::Vector3d( departureBodyState.segment( 3, 3 ) ) ) ).normalized( );
    const Eigen::Vector3d unitVector2 = unitVector3.cross( unitVector1 );

    // Calculate the absolute DSM radius.
    const double absoluteRadiusDsm = dimensionlessRadiusDsm * departureBodyState.segment( 0, 3 ).norm( );

    // Calculate the radius in the central body reference frame.
    return std::cos( inPlaneAngle ) * std::cos( outOfPlaneAngle ) * absoluteRadiusDsm *
            unitVector1 +
            std::sin( inPlaneAngle ) * std::cos( outOfPlaneAngle ) * absoluteRadiusDsm *
            unitVector2 +
            std::sin( outOfPlaneAngle ) * absoluteRadiusDsm * unitVector3;
}

void DsmPositionBasedTransferLeg::computeTransfer( )
{
    if( legParameters_.rows( ) != 6 )
    {
        throw std::runtime_error( "Error when computing DsmPositionBasedTransferLeg, incorrect input size" );
    }
    updateDepartureAndArrivalBodies(
                legParameters_( 0 ), legParameters_( 1 ) );

    double timeOfFlight = arrivalTime_ - departureTime_;
    dsmTimeOfFlightFraction_ = legParameters_( 2 );
    dimensionlessRadiusDsm_ = legParameters_( 3 );
    inPlaneAngle_ = legParameters_( 4 );
    outOfPlaneAngle_ = legParameters_( 5 );

    Eigen::Vector3d dsmLocation = calculatePositionBasedDsmLocation(
                departureBodyState_, dimensionlessRadiusDsm_, inPlaneAngle_, outOfPlaneAngle_ );


    // Calculate the DSM time of application from the time of flight fraction.
    double dsmTime = dsmTimeOfFlightFraction_ * timeOfFlight;

    // Calculate and set the spacecraft velocities after departure, before and after the DSM, and
    // before arrival using two lambert targeters and all the corresponding positions and flight
    // times.
    mission_segments::solveLambertProblemIzzo( departureBodyState_.segment( 0, 3 ), dsmLocation, dsmTime,
                                               centralBodyGravitationalParameter_,
                                               departureVelocity_, velocityBeforeDsm_ );
    mission_segments::solveLambertProblemIzzo( dsmLocation, arrivalBodyState_.segment( 0, 3 ), timeOfFlight -
                                               dsmTime, centralBodyGravitationalParameter_,
                                               velocityAfterDsm_, arrivalVelocity_ );


    //Calculate the deltaV originating from the DSM.
    Eigen::Vector3d dsmManeuver = ( velocityAfterDsm_ - velocityBeforeDsm_ );
    legTotalDeltaV_ = dsmManeuver.norm( );
    trajectoryManeuver_ = TrajectoryManeuver( dsmLocation, dsmManeuver.norm( ), dsmTime + departureTime_ );

    calculateKeplerianElements( );
}


void computeVelocityBasedDsmState(
        const Eigen::Vector6d& departureCartesianElements,
        const double centralBodyGravitationalParameter,
        const double dsmTime,
        Eigen::Vector3d& dsmLocation,
        Eigen::Vector3d& velocityBeforeDsm )
{
    // Transfer the initial position and velocity into a vectorXd object with Cartesian
    // coordinates.
    Eigen::Vector6d departureKeplerianElements = orbital_element_conversions::convertCartesianToKeplerianElements(
                departureCartesianElements, centralBodyGravitationalParameter );

    // Propagate the keplerian elements until the moment of application of the DSM.
    Eigen::Vector6d dsmArrivalKeplerianElements = orbital_element_conversions::propagateKeplerOrbit(
                departureKeplerianElements, dsmTime, centralBodyGravitationalParameter );

    // Convert the keplerian elements back into Cartesian elements.
    Eigen::Vector6d dsmArrivalCartesianElements = orbital_element_conversions::convertKeplerianToCartesianElements(
                dsmArrivalKeplerianElements, centralBodyGravitationalParameter );

    // Set the corresponding position and velocity vectors.
    dsmLocation = dsmArrivalCartesianElements.segment( 0, 3 );
    velocityBeforeDsm = dsmArrivalCartesianElements.segment( 3, 3 );
}


DsmVelocityBasedTransferLeg::DsmVelocityBasedTransferLeg(
        const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
        const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
        const double centralBodyGravitationalParameter,
        const std::function< Eigen::Vector3d( ) > departureVelocityFunction ):
    DsmTransferLeg( departureBodyEphemeris, arrivalBodyEphemeris, dsm_velocity_based_leg, centralBodyGravitationalParameter ),
    departureVelocityFunction_( departureVelocityFunction )
{ }

void DsmVelocityBasedTransferLeg::computeTransfer( )
{
    if( legParameters_.rows( ) != 3 )
    {
        throw std::runtime_error( "Error when computing DsmVelocityBasedTransferLeg, incorrect input size" );
    }
    updateDepartureAndArrivalBodies(
                legParameters_( 0 ), legParameters_( 1 ) );

    double timeOfFlight = arrivalTime_ - departureTime_;
    dsmTimeOfFlightFraction_ = legParameters_( 2 );

    // Calculate the DSM time of application from the time of flight fraction.
    double dsmTime = dsmTimeOfFlightFraction_ * timeOfFlight;

    departureVelocity_ = departureVelocityFunction_( );

    // Transfer the initial position and velocity into a vectorXd object with Cartesian
    // coordinates.
    Eigen::Vector6d departureCartesianElements;
    departureCartesianElements.segment( 0, 3 ) = departureBodyState_.segment( 0, 3 );
    departureCartesianElements.segment( 3, 3 ) = departureVelocity_;

    Eigen::Vector3d dsmLocation;
    computeVelocityBasedDsmState(
                departureCartesianElements,
                centralBodyGravitationalParameter_, dsmTime,
                dsmLocation, velocityBeforeDsm_ );

    // Calculate the velocities after the DSM and before the arrival body.
    mission_segments::solveLambertProblemIzzo(
                dsmLocation, arrivalBodyState_.segment( 0, 3 ), timeOfFlight -  dsmTime,
                centralBodyGravitationalParameter_, velocityAfterDsm_, arrivalVelocity_ );

    // Calculate the deltaV needed for the DSM.
    Eigen::Vector3d dsmManeuver = ( velocityAfterDsm_ - velocityBeforeDsm_ );
    legTotalDeltaV_ = dsmManeuver.norm( );
    trajectoryManeuver_ = TrajectoryManeuver( dsmLocation, dsmManeuver.norm( ), dsmTime + departureTime_ );

    calculateKeplerianElements( );

}


} // namespace mission_segments

} // namespace tudat
