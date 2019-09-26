/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      More information on the trajectory design code and how the quantities
 *      are calculated can be found in [Musegaas, 2012], who is also the author
 *      of this code
 *
*/

#ifndef TUDAT_LOW_THRUST_LEG_H
#define TUDAT_LOW_THRUST_LEG_H

#include <vector>

#include <Eigen/Core>
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"

//#include "Tudat/Astrodynamics/TrajectoryDesign/missionLeg.h"

namespace tudat
{
namespace transfer_trajectories
{

//! Low-thrust leg base class.
/*!
 * A base class for calculating a low-thrust trajectory leg.
 */
class LowThrustLeg /*: public MissionLeg*/
{
public:
    //! Constructor with immediate definition of parameters.
    /*!
     *  Constructor, sets objects and functions from which relevant environment and state variables
     *  are retrieved.
     *  \param departureBodyPosition location of the departure body.
     *  \param arrivalBodyPosition position of the target body.
     *  \param timeOfFlight Length of the leg.
     *  \param departureBodyVelocity velocity of the departure body.
     *  \param centralBodyGravitationalParameter gravitational parameter of the central body (most cases the Sun).
     */
    LowThrustLeg( const Eigen::Vector6d& stateAtDeparture,
                  const Eigen::Vector6d& stateAtArrival,
                  const double timeOfFlight,
                  simulation_setup::NamedBodyMap& bodyMap,
                  const std::string bodyToPropagate,
                  const std::string centralBody
/*              const Eigen::Vector3d& arrivalBodyPosition,
              const double timeOfFlight,
              const Eigen::Vector3d& departureBodyVelocity,
              const double centralBodyGravitationalParameter*/ ):
        stateAtDeparture_( stateAtDeparture ), stateAtArrival_( stateAtArrival ), timeOfFlight_( timeOfFlight ),
        bodyMap_( bodyMap ), bodyToPropagate_( bodyToPropagate ), centralBody_( centralBody )
    {
        // Retrieve initial mass of the spacecraft.
        initialMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass( );
    }
//        MissionLeg( departureBodyPosition, timeOfFlight, departureBodyVelocity, centralBodyGravitationalParameter),
//        arrivalBodyPosition_( arrivalBodyPosition ){ }

    virtual ~LowThrustLeg( ){ }

    //! Convert time to independent variable.
    virtual double convertTimeToIndependentVariable( const double time ) = 0;
//    {
//        return time;
//    }

    //! Convert independent variable to time.
    virtual double convertIndependentVariableToTime( const double independentVariable ) = 0;
//    {
//        return independentVariable;
//    }

    //! Compute direction thrust acceleration in cartesian coordinates.
    virtual Eigen::Vector3d computeCurrentThrustAccelerationDirection( double currentTime ) = 0;

    //! Compute magnitude thrust acceleration.
    virtual double computeCurrentThrustAccelerationMagnitude( double currentTime ) = 0;

//    //! Propagate state to a given time.
//    Eigen::Vector6d propagateTrajectory(
//            const Eigen::Vector6d initialState,
//            const double initialMass,
//            const double initialTime,
//            const double finalTime );

    //! Compute state history.
    virtual void getTrajectory(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::Vector6d >& propagatedTrajectory/*,
            const Eigen::Vector6d initialState,
            const double initialMass,
            const double initialTime*/ ) = 0;

//    //! Propagate mass to a given time.
//    double propagateMass(
//            const double initialMass,
//            const double initialTime,
//            const double finalTime );

    //! Compute current mass of the spacecraft.
    double computeCurrentMass( const double independentVariable,
                               std::function< double ( const double ) > specificImpulseFunction,
                               std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Return mass profile.
    void getMassProfile(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::VectorXd >& massProfile,
            std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Retrieve acceleration model (thrust).
    std::shared_ptr< propulsion::ThrustAcceleration > getLowThrustAccelerationModel(
            std::function< double( const double ) > specificImpulseFunction/*,
            std::function< Eigen::Vector3d( const double ) > thrustAccelerationDirectionFunction,
            std::function< double ( const double ) > thrustAccelerationMagnitudeFunction*/ );

    //! Retrieve acceleration map (thrust and central gravity accelerations).
    virtual basic_astrodynamics::AccelerationMap retrieveLowThrustAccelerationMap(
            std::function< double ( const double ) > specificImpulseFunction ) = 0;

    //! Compute current thrust vector.
    Eigen::Vector3d computeCurrentThrust(
            double time,
            std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Return thrust profile.
    void getThrustProfile( std::vector< double >& epochsVector,
                           std::map< double, Eigen::VectorXd >& thrustProfile,
                           std::function< double ( const double ) > specificImpulseFunction,
                           std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Compute current thrust vector.
    Eigen::Vector3d computeCurrentThrustAcceleration( double time )
    {
        double independentVariable = convertTimeToIndependentVariable( time );
        return computeCurrentThrustAccelerationMagnitude( independentVariable ) * computeCurrentThrustAccelerationDirection( independentVariable );
    }

    //! Return thrust acceleration profile.
    void getThrustAccelerationProfile(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::VectorXd >& thrustAccelerationProfile );


    //! Compute total deltaV.
    virtual double computeDeltaV( ) = 0;

//    //! Full propagation.
//    void computeSemiAnalyticalAndFullPropagation(
////            simulation_setup::NamedBodyMap& bodyMap,
//            std::function< double ( const double ) > specificImpulseFunction,
//            const std::shared_ptr<numerical_integrators::IntegratorSettings<double> > integratorSettings,
//            std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
//                    std::shared_ptr< propagators::PropagatorSettings< double > > >& propagatorSettings,
//            std::map< double, Eigen::VectorXd >& fullPropagationResults,
//            std::map< double, Eigen::VectorXd >& semiAnalyticalResults,
//            std::map< double, Eigen::VectorXd>& dependentVariablesHistory/*,
//            const bool isMassPropagated*/ );


////    //! Update the ephemeris.
////    /*!
////     * Sets the positions and the velocities to the newly specified values. Required for re-using
////     * the class, without re-initializing it.
////     *  \param departureBodyPosition sets the new departure body position.
////     *  \param arrivalBodyPosition sets the new arrival body position.
////     *  \param departureBodyVelocity sets the new departure body velocity.
////     */
////    void updateEphemeris( const Eigen::Vector3d& departureBodyPosition,
////                          const Eigen::Vector3d& arrivalBodyPosition,
////                          const Eigen::Vector3d& departureBodyVelocity )
////    {
////        departureBodyPosition_ = departureBodyPosition;
////        departureBodyVelocity_ = departureBodyVelocity;
////        arrivalBodyPosition_ = arrivalBodyPosition;
////    }

protected:

    //! Cartesian state at departure.
    Eigen::Vector6d stateAtDeparture_;

    //! Cartesian state at arrival.
    Eigen::Vector6d stateAtArrival_;

    //! Time of flight of the trajectory leg.
    double timeOfFlight_;

    //! Body map object.
    simulation_setup::NamedBodyMap bodyMap_;

    //! Name of the body to be propagated.
    std::string bodyToPropagate_;

    //! Name of the central body.
    std::string centralBody_;

    //! Initial mass of the spacecraft.
    double initialMass_;

//    //! The arrival body position.
//    /*!
//     * The position of the arrival body at the arrival time.
//     */
//    Eigen::Vector3d arrivalBodyPosition_;

    //! Compute current mass of the spacecraft between two epochs.
    double computeCurrentMass(
            const double timeInitialEpoch,
            const double timeFinalEpoch,
            const double massInitialEpoch,
            std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

private:


};
} // namespace transfer_trajectories
} // namespace tudat

#endif // TUDAT_LOW_THRUST_LEG_H
