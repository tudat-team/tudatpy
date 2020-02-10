/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_SIMS_FLANAGAN_MODEL_H
#define TUDAT_SIMS_FLANAGAN_MODEL_H

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include <cmath>
#include <vector>
#include <Eigen/Dense>

namespace tudat
{
namespace low_thrust_trajectories
{

class SimsFlanaganModel
{
public:

    //! Constructor.
    SimsFlanaganModel( const Eigen::Vector6d& stateAtDeparture,
                       const Eigen::Vector6d& stateAtArrival,
                       const double centralBodyGravitationalParameter,
                       const double initialSpacecraftMass,
                       const double maximumThrust,
                       const std::function< double ( const double ) > specificImpulseFunction,
                       const double timeOfFlight,
                       std::vector< Eigen::Vector3d >& throttles ):
        stateAtDeparture_( stateAtDeparture ), stateAtArrival_( stateAtArrival ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
        initialSpacecraftMass_( initialSpacecraftMass ),
        maximumThrust_( maximumThrust ),
        specificImpulseFunction_( specificImpulseFunction ), timeOfFlight_( timeOfFlight ),
        throttles_( throttles )
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

        propagateMassToSegments( );

        // Initialise value of the total deltaV.
        totalDeltaV_ = 0.0;
    }


    //! Default destructor.
    ~SimsFlanaganModel( ) { }

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

    //! Return total deltaV required by the trajectory.
    double getTotalDeltaV( );

    basic_astrodynamics::AccelerationMap getLowThrustTrajectoryAccelerationMap(
            const simulation_setup::NamedBodyMap& bodyMapTest,
            const std::string& bodyToPropagate,
            const std::string& centralBody );

    //! Propagate the trajectory to given time.
    Eigen::Vector6d propagateTrajectoryForward(
            double initialTime, double finalTime, Eigen::Vector6d initialState, double segmentDuration );

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

    void propagateMassToSegments( );

    int convertTimeToLegSegment( double currentTime );

    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > getThrustAccelerationSettingsFullLeg(
            const simulation_setup::NamedBodyMap& bodyMapTest );

    double getMassAtSegment( const int segment )
    {
        return segmentMasses_.at( segment );
    }

protected:

    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > getConstantThrustAccelerationSettingsPerSegment(
            unsigned int indexSegment );

    basic_astrodynamics::AccelerationMap getAccelerationModelPerSegment(
            const unsigned int indexSegment,
            const simulation_setup::NamedBodyMap& bodyMapTest,
            const std::string& bodyToPropagate,
            const std::string& centralBody );

private:

    //! State vector of the vehicle at the leg departure.
    Eigen::Vector6d stateAtDeparture_;

    //! State vector of the vehicle at the leg arrival.
    Eigen::Vector6d stateAtArrival_;

    //! Gravitational parameter of the central body of the 2-body problem.
    double centralBodyGravitationalParameter_;

    double initialSpacecraftMass_;

    //! Maximum allowed thrust.
    double maximumThrust_;

    //! Specific impulse.
    std::function< double ( const double ) > specificImpulseFunction_;

    //! Number of segments into which the leg is subdivided.
    int numberSegments_;

    //! Time of flight for the leg.
    double timeOfFlight_;

    //! Vector of throttles for each segment of the leg.
    std::vector< Eigen::Vector3d > throttles_;

    std::vector< double > segmentMasses_;


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

    //! Total deltaV.
    double totalDeltaV_;

    //! Vector containing the time associated to each node of the leg.
    std::vector< double > timesAtNodes_;


};


} // namespace low_thrust_trajectories
} // namespace tudat

#endif // TUDAT_SIMS_FLANAGAN_MODEL_H
