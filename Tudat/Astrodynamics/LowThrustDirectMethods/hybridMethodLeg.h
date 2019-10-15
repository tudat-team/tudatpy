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

#ifndef HYBRIDMETHODLEG_H
#define HYBRIDMETHODLEG_H

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include <math.h>
#include <vector>
#include <Eigen/Dense>

namespace tudat
{
namespace low_thrust_direct_methods
{

class HybridMethodLeg
{
public:

    //! Constructor.
    HybridMethodLeg( const Eigen::Vector6d& stateAtDeparture,
                     const Eigen::Vector6d& stateAtArrival,
                     const Eigen::VectorXd& initialCoStates,
                     const Eigen::VectorXd& finalCoStates,
                     const double maximumThrust,
                     const double specificImpulse,
                     const double timeOfFlight,
                     simulation_setup::NamedBodyMap& bodyMap,
                     const std::string bodyToPropagate,
                     const std::string centralBody,
                     std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings ):
    stateAtDeparture_( stateAtDeparture ), stateAtArrival_( stateAtArrival ), initialCoStates_( initialCoStates ),
    finalCoStates_( finalCoStates ), maximumThrust_( maximumThrust ),
    specificImpulse_( specificImpulse ), timeOfFlight_( timeOfFlight ), bodyMap_( bodyMap ),
    bodyToPropagate_( bodyToPropagate ), centralBody_( centralBody ), integratorSettings_( integratorSettings )
    {
        // Initialise value of the total deltaV.
        totalDeltaV_ = 0.0;

        // Retrieve gravitational parameter of the central body.
        centralBodyGravitationalParameter_ = bodyMap_[ centralBody_ ]->getGravityFieldModel()->getGravitationalParameter();

        // Retrieve initial mass of the spacecraft.
        initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();

        // Define function returning the current MEE costates.
        costatesFunction_ = [ = ]( const double currentTime )
        {
            Eigen::VectorXd currentCostates;
            currentCostates.resize( 5 );

            for ( int i = 0 ; i < 5 ; i++ )
            {
                currentCostates[ i ] = initialCoStates_[ i ]
                        + ( currentTime / timeOfFlight_ ) * ( finalCoStates_[ i ] - initialCoStates_[ i ] );
            }
            return currentCostates;
        };

        // Initialise mass at time of flight (before propagation).
       Eigen::Vector6d propagatedStateAtTimeOfFlight = propagateTrajectory( );

    }


    //! Default destructor.
    ~HybridMethodLeg( ) { }

    //! Retrieve MEE costates-based thrust acceleration.
    std::shared_ptr< simulation_setup::AccelerationSettings > getMEEcostatesBasedThrustAccelerationSettings( );

    std::shared_ptr< simulation_setup::ThrustMagnitudeSettings > getMEEcostatesBasedThrustMagnitudeSettings( );
    std::shared_ptr< simulation_setup::ThrustDirectionGuidanceSettings > getMEEcostatesBasedThrustDirectionSettings( );

    //! Retrieve hybrid method acceleration model (including thrust and central gravity acceleration)
    basic_astrodynamics::AccelerationMap getLowThrustTrajectoryAccelerationMap( );

    //! Propagate the spacecraft trajectory to time of flight.
    Eigen::Vector6d propagateTrajectory( );

    //! Propagate the spacecraft trajectory to a given time.
    Eigen::Vector6d  propagateTrajectory( double initialTime, double finalTime, Eigen::Vector6d initialState, double initialMass );

    //! Propagate the trajectory to set of epochs.
    std::map< double, Eigen::Vector6d > propagateTrajectory(
            std::vector< double > epochs, std::map< double, Eigen::Vector6d >& propagatedTrajectory );

    //! Return the deltaV associated with the thrust profile of the trajectory.
    double computeDeltaV( );

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

    //! Returns propagated mass when the time of flight is reached.
    double getMassAtTimeOfFlight( )
    {
        return massAtTimeOfFlight_;
    }

    //! Returns maximum allowed thrust.
    double getMaximumThrustValue( )
    {
        return maximumThrust_;
    }

    //! Returns time of flight.
    double getTimeOfFlight( )
    {
        return timeOfFlight_;
    }

    //! Return total deltaV required by the trajectory.
    double getTotalDeltaV( )
    {
        return computeDeltaV( );//totalDeltaV_;
    }

    //! Return the current MEE co-states.
    std::function< Eigen::VectorXd( const double ) > getCostatesFunction_( )
    {
        return costatesFunction_;
    }

protected:

private:

    //! State vector of the vehicle at the leg departure.
    Eigen::Vector6d stateAtDeparture_;

    //! State vector of the vehicle at the leg arrival.
    Eigen::Vector6d stateAtArrival_;

    //! Modified equinoctial elements.

    //! Initial co-states vector.
    Eigen::Vector6d initialCoStates_;

    //! Final co-states vector.
    Eigen::Vector6d finalCoStates_;

    //! Function returning the current MEE co-states.
    std::function< Eigen::VectorXd( const double ) > costatesFunction_;

    //! Maximum allowed thrust.
    double maximumThrust_;

    //! Specific impulse.
    double specificImpulse_;

    //! Time of flight for the leg.
    double timeOfFlight_;

    //! Pointer to the body map object.
    simulation_setup::NamedBodyMap bodyMap_;

    //! Gravitational parameter of the central body of the 2-body problem.
    double centralBodyGravitationalParameter_;

    //! Name of the body to be propagated.
    std::string bodyToPropagate_;

    //! Name of the central body.
    std::string centralBody_;

    //! Integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings_;

    //! Total deltaV.
    double totalDeltaV_;

    //! Initial mass of the spacecraft.
    double initialSpacecraftMass_;

    //! Mass of the spacecraft at the end of the propagation.
    double massAtTimeOfFlight_;

};


} // namespace low_thrust_direct_methods
} // namespace tudat

#endif // HYBRIDMETHODLEG_H
