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
        const TransferLegTypes legType,
        const Eigen::VectorXd legParameters ):
    departureBodyEphemeris_( departureBodyEphemeris ), arrivalBodyEphemeris_( arrivalBodyEphemeris ),
    legType_( legType ), legParameters_( legParameters ){ }

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
        const Eigen::VectorXd legParameters,
        const double centralBodyGravitationalParameter ):
    TransferLeg( departureBodyEphemeris, arrivalBodyEphemeris, unpowered_unperturbed_leg, legParameters ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter )
{
//    computeTransfer( );
    legTotalDeltaV_ = 0.0;
}

void UnpoweredUnperturbedTransferLeg::getStateAlongTrajectory( std::map< double, Eigen::Vector6d >& statesAlongTrajectory,
                              const std::vector< double >& timePoints )
{
    Eigen::Vector6d initialState;
    initialState.segment( 0, 3 ) = departureBodyState_.segment( 0, 3 );
    initialState.segment( 3, 3 ) = departureVelocity_;

    Eigen::Vector6d initialKeplerianState = orbital_element_conversions::convertCartesianToKeplerianElements(
                initialState, centralBodyGravitationalParameter_ );

    statesAlongTrajectory = orbital_element_conversions::getKeplerOrbitCartesianStateHistory(
            initialKeplerianState, departureTime_, timePoints, centralBodyGravitationalParameter_ );
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
}

void DsmTransferLeg::getStateAlongTrajectory( std::map< double, Eigen::Vector6d >& statesAlongTrajectory,
                              const std::vector< double >& timePoints )
{
    Eigen::Vector6d initialState;
    initialState.segment( 0, 3 ) = departureBodyState_.segment( 0, 3 );
    initialState.segment( 3, 3 ) = departureVelocity_;

    Eigen::Vector6d initialKeplerianState = orbital_element_conversions::convertCartesianToKeplerianElements(
                initialState, centralBodyGravitationalParameter_ );

    auto lowerBound = std::lower_bound ( timePoints.begin( ), timePoints.end( ), dsmTime_ );
    std::vector< double > timesUpToDsm;
    std::vector< double > timesAfterDsm;

    timesUpToDsm.insert( timesUpToDsm.begin( ), timePoints.begin( ), lowerBound );
    timesAfterDsm.insert( timesAfterDsm.begin( ), lowerBound, timePoints.end( ) );

     std::map< double, Eigen::Vector6d > statesUpToDsm = orbital_element_conversions::getKeplerOrbitCartesianStateHistory(
                    initialKeplerianState, departureTime_, timesUpToDsm, centralBodyGravitationalParameter_ );

    Eigen::Vector6d stateAfterDsm;
    stateAfterDsm.segment( 0, 3 ) = dsmLocation_;
    stateAfterDsm.segment( 3, 3 ) = velocityAfterDsm_;;

    Eigen::Vector6d keplerianStateAfterDsm = orbital_element_conversions::convertCartesianToKeplerianElements(
                stateAfterDsm, centralBodyGravitationalParameter_ );

    std::map< double, Eigen::Vector6d > statesAfterDsm = orbital_element_conversions::getKeplerOrbitCartesianStateHistory(
                   keplerianStateAfterDsm, dsmTime_, timesAfterDsm, centralBodyGravitationalParameter_ );

    statesAlongTrajectory = std::move( statesUpToDsm );
    statesAlongTrajectory.insert( statesAfterDsm.begin( ), statesAfterDsm.end( ) );

}

DsmPositionBasedTransferLeg::DsmPositionBasedTransferLeg(
        const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
        const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
        const Eigen::VectorXd legParameters,
        const double centralBodyGravitationalParameter ):
    DsmTransferLeg( departureBodyEphemeris, arrivalBodyEphemeris, dsm_position_based_leg, legParameters, centralBodyGravitationalParameter )
{
//    computeTransfer( );
}



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

    dsmLocation_ = calculatePositionBasedDsmLocation(
                departureBodyState_, dimensionlessRadiusDsm_, inPlaneAngle_, outOfPlaneAngle_ );


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
    legTotalDeltaV_ = ( velocityAfterDsm_ - velocityBeforeDsm_ ).norm( );
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
        const Eigen::VectorXd legParameters,
        const double centralBodyGravitationalParameter,
        const std::function< Eigen::Vector3d( ) > departureVelocityFunction ):
    DsmTransferLeg( departureBodyEphemeris, arrivalBodyEphemeris, dsm_velocity_based_leg, legParameters, centralBodyGravitationalParameter ),
    departureVelocityFunction_( departureVelocityFunction )
{
//    computeTransfer( );
}

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
    dsmTime_ = dsmTimeOfFlightFraction_ * timeOfFlight;

    departureVelocity_ = departureVelocityFunction_( );

    // Transfer the initial position and velocity into a vectorXd object with Cartesian
    // coordinates.
    Eigen::Vector6d departureCartesianElements;
    departureCartesianElements.segment( 0, 3 ) = departureBodyState_.segment( 0, 3 );
    departureCartesianElements.segment( 3, 3 ) = departureVelocity_;

    computeVelocityBasedDsmState(
                departureCartesianElements,
                centralBodyGravitationalParameter_, dsmTime_,
                dsmLocation_, velocityBeforeDsm_ );

    // Calculate the velocities after the DSM and before the arrival body.
    mission_segments::solveLambertProblemIzzo(
                dsmLocation_, arrivalBodyState_.segment( 0, 3 ), timeOfFlight -  dsmTime_,
                centralBodyGravitationalParameter_, velocityAfterDsm_, arrivalVelocity_ );

    // Calculate the deltaV needed for the DSM.
    legTotalDeltaV_ = ( velocityAfterDsm_ - velocityBeforeDsm_ ).norm( );

}


} // namespace mission_segments

} // namespace tudat
