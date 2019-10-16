/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_SIMS_FLANAGAN_LEG_H
#define TUDAT_SIMS_FLANAGAN_LEG_H

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include <math.h>
#include <vector>
#include <Eigen/Dense>

namespace tudat
{
namespace low_thrust_trajectories
{

class SimsFlanaganLeg
{
public:

    //! Constructor.
    SimsFlanaganLeg( const Eigen::Vector6d& stateAtDeparture,
                     const Eigen::Vector6d& stateAtArrival,
                     const double maximumThrust,
                     const std::function< double ( const double ) > specificImpulseFunction,
                     const double timeOfFlight,
                     simulation_setup::NamedBodyMap& bodyMap,
                     std::vector< Eigen::Vector3d >& throttles,
                     const std::string bodyToPropagate,
                     const std::string centralBody ):
    stateAtDeparture_( stateAtDeparture ), stateAtArrival_( stateAtArrival ), maximumThrust_( maximumThrust ),
    specificImpulseFunction_( specificImpulseFunction ), timeOfFlight_( timeOfFlight ), bodyMap_( bodyMap ),
    throttles_( throttles ), bodyToPropagate_( bodyToPropagate ), centralBody_( centralBody )
    {

        // Retrieve number of segments from the size of the throttles vector.
        numberSegments_ = throttles_.size();

        // Calculate number of segments for both the forward propagation (from departure to match point)
        // and the backward propagation (from arrival to match point).
        numberSegmentsForwardPropagation_ = ( numberSegments_ + 1 ) / 2;
        numberSegmentsBackwardPropagation_ = numberSegments_ / 2;

        // Compute time at match point (half of the time of flight).
        timeAtMatchPoint_ = timeOfFlight_ / 2.0;

        // Compute the times at the different nodes of the leg.
        segmentDurationForwardPropagation_ = timeAtMatchPoint_ / numberSegmentsForwardPropagation_;
        segmentDurationBackwardPropagation_ = timeAtMatchPoint_ / numberSegmentsBackwardPropagation_;
        for ( int i = 0 ; i <= numberSegmentsForwardPropagation_ ; i++ )
        {
            timesAtNodes_.push_back( i * segmentDurationForwardPropagation_ );
        }
        for ( int i = 1 ; i <= numberSegmentsBackwardPropagation_ ; i++ )
        {
            timesAtNodes_.push_back( timeAtMatchPoint_ + i * segmentDurationBackwardPropagation_ );
        }

        // Initialise value of the total deltaV.
        totalDeltaV_ = 0.0;

        // Retrieve gravitational parameter of the central body.
        centralBodyGravitationalParameter_ = bodyMap_[ centralBody_ ]->getGravityFieldModel()->getGravitationalParameter();

        // Retrieve initial mass of the spacecraft.
        initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();
    }


    //! Default destructor.
    ~SimsFlanaganLeg( ) { }

    //! Propagate the spacecraft trajectory from departure to match point (forward propagation).
    void propagateForwardFromDepartureToMatchPoint( );

    //! Propagate the spacecraft trajectory from arrival to match point (backward propagation).
    void propagateBackwardFromArrivalToMatchPoint( );

    //! Returns initial state at leg departure.
    Eigen::VectorXd getStateAtLegDeparture( )
    {
        return stateAtDeparture_;
    }

    //! Returns final state at leg arrival.
    Eigen::VectorXd getStateAtLegArrival( )
    {
        return stateAtArrival_;
    }

    //! Returns maximum allowed thrust.
    double getMaximumThrustValue( )
    {
        return maximumThrust_;
    }

    //! Returns number of segments into which the leg is subdivided.
    int getNumberOfSegments( )
    {
        return numberSegments_;
    }

    //! Returns time of flight.
    double getTimeOfFlight( )
    {
        return timeOfFlight_;
    }

    //! Return vector of throttles.
    std::vector< Eigen::Vector3d > getThrottles( )
    {
        return throttles_;
    }

    //! Return state vector at match point from forward propagation.
    Eigen::Vector6d getStateAtMatchPointForwardPropagation( )
    {
        return stateAtMatchPointFromForwardPropagation_;
    }

    //! Return state vector at match point from backward propagation.
    Eigen::Vector6d getStateAtMatchPointBackwardPropagation( )
    {
        return stateAtMatchPointFromBackwardPropagation_;
    }

    //! Return spacecraft mass at match point from forward propagation.
    double getMassAtMatchPointForwardPropagation( )
    {
        return massAtMatchPointFromForwardPropagation_;
    }

    //! Return spacecraft mass at match point from backward propagation.
    double getMassAtMatchPointBackwardPropagation( )
    {
        return massAtMatchPointFromBackwardPropagation_;
    }

    //! Return total deltaV required by the trajectory.
    double getTotalDeltaV( )
    {
        totalDeltaV_ = 0.0;

        double currentMass = initialSpacecraftMass_;

        for ( int currentSegment = 0 ; currentSegment < numberSegments_ ; currentSegment++ )
        {
            double segmentDuration;
            if ( currentSegment < numberSegmentsForwardPropagation_ )
            {
                segmentDuration = segmentDurationForwardPropagation_;
            }
            else
            {
                segmentDuration = segmentDurationBackwardPropagation_;
            }

            Eigen::Vector3d deltaVvector = maximumThrust_ / currentMass
                                * segmentDuration * throttles_[ currentSegment ];

            currentMass = propagateMassToSegment( currentSegment );

            totalDeltaV_ += deltaVvector.norm( );
        }

        return totalDeltaV_;
    }

    basic_astrodynamics::AccelerationMap getLowThrustTrajectoryAccelerationMap( );

    //! Propagate the trajectory to given time.
    Eigen::Vector6d propagateTrajectoryForward( double initialTime, double finalTime, Eigen::Vector6d initialState, double segmentDuration );

    //! Propagate the trajectory to given time.
    Eigen::Vector6d propagateTrajectoryBackward( double initialTime, double finalTime, Eigen::Vector6d initialState, double segmentDuration );

    //! Propagate the trajectory to set of epochs (forward propagation).
    std::map< double, Eigen::Vector6d > propagateTrajectoryForward(
            std::vector< double > epochs, std::map< double, Eigen::Vector6d >& propagatedTrajectory, double segmentDuration );

    //! Propagate the trajectory to set of epochs (backward propagation).
    std::map< double, Eigen::Vector6d > propagateTrajectoryBackward(
            std::vector< double > epochs, std::map< double, Eigen::Vector6d >& propagatedTrajectory, double segmentDuration );

    //! Propagate the trajectory to set of epochs.
    void propagateTrajectory(
            std::vector< double > epochs, std::map< double, Eigen::Vector6d >& propagatedTrajectory );

    //! Propagate the trajectory inside one segment.
    Eigen::Vector6d propagateInsideForwardSegment( double initialTime, double finalTime, double segmentDuration, Eigen::Vector6d initialState );

    //! Propagate the trajectory inside one segment.
    Eigen::Vector6d propagateInsideBackwardSegment( double initialTime, double finalTime, double segmentDuration, Eigen::Vector6d initialState );

    //! Propagate the mass from departure to a given segment.
    double propagateMassToSegment( int indexSegment );

    int convertTimeToLegSegment( double currentTime );

    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > getThrustAccelerationSettingsFullLeg( );



protected:

    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > getConstantThrustAccelerationSettingsPerSegment(
            unsigned int indexSegment );

    basic_astrodynamics::AccelerationMap getAccelerationModelPerSegment( unsigned int indexSegment );

private:

    //! State vector of the vehicle at the leg departure.
    Eigen::Vector6d stateAtDeparture_;

    //! State vector of the vehicle at the leg arrival.
    Eigen::Vector6d stateAtArrival_;

    //! Maximum allowed thrust.
    double maximumThrust_;

    //! Specific impulse.
    std::function< double ( const double ) > specificImpulseFunction_;

    //! Number of segments into which the leg is subdivided.
    int numberSegments_;

    //! Time of flight for the leg.
    double timeOfFlight_;

    //! Pointer to the body map object.
    simulation_setup::NamedBodyMap bodyMap_;

    //! Vector of throttles for each segment of the leg.
    std::vector< Eigen::Vector3d > throttles_;

    //! Gravitational parameter of the central body of the 2-body problem.
    double centralBodyGravitationalParameter_;

    //! Name of the body to be propagated.
    std::string bodyToPropagate_;

    //! Name of the central body.
    std::string centralBody_;

    //! Number of segments for the forward propagation from departure to match point.
    int numberSegmentsForwardPropagation_;

    //! Number of segments for the backward propagation from arrival to match point.
    int numberSegmentsBackwardPropagation_;

    //! Segment duration for forward propagation.
    double segmentDurationForwardPropagation_;

    //! Segment duration for backward propagation.
    double segmentDurationBackwardPropagation_;

    //! Time defining the match point (half of the time of flight).
    double timeAtMatchPoint_;

    //! State at match point from forward propagation.
    Eigen::Vector6d stateAtMatchPointFromForwardPropagation_;

    //! State at match point from backward propagation.
    Eigen::Vector6d stateAtMatchPointFromBackwardPropagation_;

    //! Mass of the body to be propagated at match point from forward propagation.
    double massAtMatchPointFromForwardPropagation_;

    //! Mass of the body to be propagated at match point from backward propagation.
    double massAtMatchPointFromBackwardPropagation_;

    //! Total deltaV.
    double totalDeltaV_;

    //! Vector containing the time associated to each node of the leg.
    std::vector< double > timesAtNodes_;

    //! Initial mass of the spacecraft.
    double initialSpacecraftMass_;

};


} // namespace low_thrust_trajectories
} // namespace tudat

#endif // TUDAT_SIMS_FLANAGAN_LEG_H
