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
 *      More information on the trajectory design code and how the quantities
 *      are calculated can be found in [Musegaas, 2012], who is also the author
 *      of this code
 *
*/

#ifndef TUDAT_LOW_THRUST_LEG_H
#define TUDAT_LOW_THRUST_LEG_H

#include <iostream>
#include <vector>
#include <map>

#include <Eigen/Core>
#include "tudat/basics/basicTypedefs.h"

namespace tudat
{
namespace low_thrust_trajectories
{

//! Low-thrust leg base class.
class LowThrustLeg
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
                  const bool legModelIsForceBased ):
        stateAtDeparture_( stateAtDeparture ), stateAtArrival_( stateAtArrival ), timeOfFlight_( timeOfFlight ),
        legModelIsForceBased_( legModelIsForceBased ){ }

    virtual ~LowThrustLeg( ){ }


    bool getLegModelIsForceBased( )
    {
        return legModelIsForceBased_;
    }

    double getTimeOfFlight( )
    {
        return timeOfFlight_;
    }


    //! Convert time to independent variable.
    virtual double convertTimeToIndependentVariable( const double time ) = 0;

    //! Convert independent variable to time.
    virtual double convertIndependentVariableToTime( const double independentVariable ) = 0;

    //! Compute current cartesian state.
    virtual Eigen::Vector6d computeCurrentStateVector( const double currentTime ) = 0;




    //! Compute state history.
    virtual void getTrajectory( std::vector< double >& epochsVector, std::map< double, Eigen::Vector6d >& propagatedTrajectory ) = 0;

    std::map< double, Eigen::Vector6d > getTrajectory( std::vector< double >& epochsVector )
    {
       std::map< double, Eigen::Vector6d > propagatedTrajectory;
       getTrajectory( epochsVector, propagatedTrajectory );
       return propagatedTrajectory;
    }

    //! TODO-REFACTOR: REMOVE
    Eigen::Vector6d getStateAtEpoch( const double epoch )
    {
        std::vector< double > epochs = { epoch };
        return getTrajectory( epochs ).begin( )->second;
    }




//    //! Compute direction thrust acceleration in cartesian coordinates.
//    virtual Eigen::Vector3d computeCurrentThrustAccelerationDirection( double currentIndependentVariable ) = 0;

//    //! Compute magnitude thrust acceleration.
//    virtual double computeCurrentThrustAccelerationMagnitude( double currentIndependentVariable ) = 0;



    double getForceBasedThrustMagnitude(
            const double currentTime,
            const double timeOffset )
    {
        return computeCurrentThrustForce( currentTime, timeOffset ).norm( );
    }

    double getAccelerationBasedThrustMagnitude(
            const double currentTime,
            const double timeOffset,
            const double currentVehicleMass )
    {
        return computeCurrentThrustAcceleration( currentTime, timeOffset ).norm( ) * currentVehicleMass;
    }


    Eigen::Vector3d getThrustDirection(
            const double currentTime,
            const double timeOffset )
    {
        if( legModelIsForceBased_ )
        {
            return computeCurrentThrustForce( currentTime, timeOffset ).normalized( );
        }
        else
        {
            return computeCurrentThrustAcceleration( currentTime, timeOffset ).normalized( );
        }
    }

    //! Compute current thrust vector.
   virtual Eigen::Vector3d computeCurrentThrustForce(
           const double currentTime,
           const double timeOffset )
   {
        if( legModelIsForceBased_ )
        {
            throw std::runtime_error( "Error, thrust force function not yet overridden" );
        }
        else
        {
            throw std::runtime_error( "Error, requested thrust force from acceleration based model" );
        }
   }

    //! Compute current thrust vector.
    virtual Eigen::Vector3d computeCurrentThrustAcceleration(
            const double currentTime,
            const double timeOffset  )
    {
         if( !legModelIsForceBased_ )
         {
             throw std::runtime_error( "Error, thrust acceleration function not yet overridden" );
         }
         else
         {
             throw std::runtime_error( "Error, requested thrust acceleration from force based model" );
         }
    }



//    //! Compute current mass of the spacecraft.
//    double computeCurrentMass( const double independentVariable,
//                               std::function< double ( const double ) > specificImpulseFunction,
//                               std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

//    //! Return mass profile.
//    // LOW THRUST TODO: MOVE TO FREE FUNCTION
//    void getMassProfile(
//            const std::vector< double >& epochsVector,
//            std::map< double, Eigen::VectorXd >& massProfile,
//            const std::function< double ( const double ) > specificImpulseFunction,
//            const std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

//    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > getLowThrustAccelerationSettings(
//            const simulation_setup::SystemOfBodies& bodies,
//            const std::string& bodyToPropagate,
//            const std::function< double( const double ) > specificImpulseFunction,
//            const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
//            const double timeOffset = 0.0 );



//    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > getLowThrustAccelerationSettings(
//            const simulation_setup::SystemOfBodies& bodies,
//            const std::string& bodyToPropagate,
//            const double specificImpulse,
//            const double timeOffset );

//    //! Retrieve acceleration model (thrust).
//    std::shared_ptr< propulsion::ThrustAcceleration > getLowThrustAccelerationModel(
//            const simulation_setup::SystemOfBodies& bodies,
//            const std::string& bodyToPropagate,
//            const std::function< double( const double ) > specificImpulseFunction,
//            const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
//            const double timeOffset = 0.0 );

//    //! Retrieve acceleration map (thrust and central gravity accelerations).
//    virtual basic_astrodynamics::AccelerationMap retrieveLowThrustAccelerationMap(
//            const simulation_setup::SystemOfBodies& bodies,
//            const std::string& bodyToPropagate,
//            const std::string& centralBody,
//            const std::function< double ( const double ) > specificImpulseFunction,
//            const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings ) = 0;


//    //! Return thrust profile.
//    // LOW THRUST TODO: MOVE TO FREE FUNCTION
//    virtual void getThrustForceProfile( std::vector< double >& epochsVector,
//                           std::map< double, Eigen::VectorXd >& thrustProfile,
//                           std::function< double ( const double ) > specificImpulseFunction,
//                           std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );


//    //! Return thrust acceleration profile.
//    virtual void getThrustAccelerationProfile(
//            std::vector< double >& epochsVector,
//            std::map< double, Eigen::VectorXd >& thrustAccelerationProfile,
//            std::function< double ( const double ) > specificImpulseFunction,
//            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

//    std::map< double, Eigen::VectorXd > getThrustAccelerationProfile(
//            std::vector< double >& epochsVector,
//            std::function< double ( const double ) > specificImpulseFunction,
//            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
//    {
//        std::map< double, Eigen::VectorXd > thrustProfile;
//        getThrustAccelerationProfile( epochsVector, thrustProfile, specificImpulseFunction, integratorSettings );
//        return thrustProfile;
//    }


    //! Compute total deltaV.
    virtual double computeDeltaV( ) = 0;


//    //! Full propagation.
//    void computeSemiAnalyticalAndFullPropagation(
//            const simulation_setup::SystemOfBodies& bodies,
//            const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
//            const std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
//                    std::shared_ptr< propagators::PropagatorSettings< double > > >& propagatorSettings,
//            std::map< double, Eigen::VectorXd >& fullPropagationResults,
//            std::map< double, Eigen::Vector6d >& semiAnalyticalResults,
//            std::map< double, Eigen::VectorXd >& dependentVariablesHistory );

//    //! Define appropriate propagator settings for the full propagation.
//    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
//    std::shared_ptr< propagators::PropagatorSettings< double > > > createLowThrustPropagatorSettings(
//            const simulation_setup::SystemOfBodies& bodies,
//            const std::string& bodyToPropagate,
//            const std::string& centralBody,
//            const std::function< double( const double ) > specificImpulseFunction,
//            const basic_astrodynamics::AccelerationMap perturbingAccelerationsMap,
//            const std::shared_ptr< numerical_integrators::IntegratorSettings< double > >& integratorSettings,
//            const std::shared_ptr< propagators::DependentVariableSaveSettings >& dependentVariablesToSave );

//    //! Define appropriate propagator settings for the full propagation.
//    virtual std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
//    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > createLowThrustTranslationalStatePropagatorSettings(
//            const std::string& bodyToPropagate,
//            const std::string& centralBody,
//            const basic_astrodynamics::AccelerationMap& accelerationModelMap,
//            const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave ) = 0;


protected:

    //! Cartesian state at departure.
    Eigen::Vector6d stateAtDeparture_;

    //! Cartesian state at arrival.
    Eigen::Vector6d stateAtArrival_;

    //! Time of flight of the trajectory leg.
    double timeOfFlight_;

//    double initialMass_;

    bool legModelIsForceBased_;

//    //! Compute current mass of the spacecraft between two epochs.
//    double computeCurrentMass(
//            const double timeInitialEpoch,
//            const double timeFinalEpoch,
//            const double massInitialEpoch,
//            const std::function< double ( const double ) > specificImpulseFunction,
//            const std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

private:


};
} // namespace low_thrust_trajectories
} // namespace tudat

#endif // TUDAT_LOW_THRUST_LEG_H
