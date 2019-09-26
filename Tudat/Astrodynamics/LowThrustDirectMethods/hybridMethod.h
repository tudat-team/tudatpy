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

#ifndef HYBRIDMETHOD_H
#define HYBRIDMETHOD_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <map>
#include "pagmo/algorithm.hpp"

#include "Tudat/Astrodynamics/LowThrustDirectMethods/hybridMethodLeg.h"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/lowThrustLeg.h"

namespace tudat
{
namespace low_thrust_direct_methods
{

//! Transform thrust model as a function of time into hybrid method thrust model.
Eigen::Matrix< double, 10 , 1 > convertToHybridMethodThrustModel( std::function< Eigen::Vector3d( const double ) > thrustModelWrtTime );

class HybridMethod : public transfer_trajectories::LowThrustLeg
{
public:

    //! Constructor.
    HybridMethod(
            const Eigen::Vector6d& stateAtDeparture,
            const Eigen::Vector6d& stateAtArrival,
            const double maximumThrust,
            const double specificImpulse,
            const double timeOfFlight,
            simulation_setup::NamedBodyMap& bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
            pagmo::algorithm optimisationAlgorithm,
            const int numberOfGenerations,
            const int numberOfIndividualsPerPopulation,
            const double relativeToleranceConstraints = 1.0e-6,
            std::pair< std::function< Eigen::Vector3d( const double ) >, double > initialGuessThrustModel = std::make_pair( nullptr, 0.0 ) ) :
        LowThrustLeg( stateAtDeparture, stateAtArrival, timeOfFlight, bodyMap, bodyToPropagate, centralBody ),
//        stateAtDeparture_( stateAtDeparture ),
//        stateAtArrival_( stateAtArrival ),
        maximumThrust_( maximumThrust ),
        specificImpulse_( specificImpulse ),
//        timeOfFlight_( timeOfFlight ),
//        bodyMap_( bodyMap ),
//        bodyToPropagate_( bodyToPropagate ),
//        centralBody_( centralBody ),
        integratorSettings_( integratorSettings ),
        optimisationAlgorithm_( optimisationAlgorithm ),
        numberOfGenerations_( numberOfGenerations ),
        numberOfIndividualsPerPopulation_( numberOfIndividualsPerPopulation ),
        relativeToleranceConstraints_( relativeToleranceConstraints ) //,
//        initialGuessThrustModel_( initialGuessThrustModel )
    {

        // Store initial spacecraft mass.
        initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();

        // Convert the thrust model proposed as initial guess into simplified thrust model adapted to the hybrid method.
        if ( initialGuessThrustModel.first != nullptr )
        {
            initialGuessThrustModel_.first = convertToHybridMethodThrustModel( initialGuessThrustModel.first );
        }
        else
        {
            Eigen::VectorXd emptyVector;
            initialGuessThrustModel_.first = emptyVector; //Eigen::VectorXd::Zero( 10 ); // emptyVector;
        }
        initialGuessThrustModel_.second = initialGuessThrustModel.second;

        // Perform optimisation
        std::pair< std::vector< double >, std::vector< double > > bestIndividual = performOptimisation( );
        championFitness_ = bestIndividual.first;
        championDesignVariables_ = bestIndividual.second;


        // Transform vector of design variables into 3D vector of throttles.
        Eigen::VectorXd initialCostates; initialCostates.resize( 5 );
        Eigen::VectorXd finalCostates; finalCostates.resize( 5 );
        for ( unsigned int i = 0 ; i < 5 ; i++ )
        {
            initialCostates[ i ] = championDesignVariables_[ i ];
            finalCostates[ i ] = championDesignVariables_[ i + 5 ];
        }

        bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );

        // Create Sims-Flanagan leg from the best optimisation individual.
        hybridMethodLeg_ = std::make_shared< HybridMethodLeg >( stateAtDeparture_, stateAtArrival_, initialCostates, finalCostates, maximumThrust_,
                                                                specificImpulse_, timeOfFlight_, bodyMap_, bodyToPropagate_, centralBody_, integratorSettings );

    }

    //! Default destructor.
    ~HybridMethod( ) { }

    //! Convert time to independent variable.
    double convertTimeToIndependentVariable( const double time )
    {
        return time;
    }

    //! Convert independent variable to time.
    double convertIndependentVariableToTime( const double independentVariable )
    {
        return independentVariable;
    }

    //! Perform optimisation.
    std::pair< std::vector< double >, std::vector< double > > performOptimisation( );

    //! Compute DeltaV.
    double computeDeltaV( )
    {
        return championFitness_[ 0 ];
    }

    //! Compute current cartesian state.
    Eigen::Vector6d computeCurrentStateVector( const double currentTime );

    //! Compute state history.
    void getTrajectory(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::Vector6d >& propagatedTrajectory )
    {
        propagatedTrajectory = hybridMethodLeg_->propagateTrajectory( epochsVector, propagatedTrajectory );
    }

//    //! Return mass profile.
//    void getMassProfile(
//            std::vector< double >& epochsVector,
//            std::map< double, Eigen::VectorXd >& massProfile,
//            std::function< double ( const double ) > specificImpulseFunction,
//            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Compute direction thrust acceleration in cartesian coordinates.
    Eigen::Vector3d computeCurrentThrustAccelerationDirection( double currentTime );

    //! Compute magnitude thrust acceleration.
    double computeCurrentThrustAccelerationMagnitude( double currentTime );

//    //! Compute current thrust vector.
//    Eigen::Vector3d computeCurrentThrustAcceleration( double time );

//    //! Return thrust acceleration profile.
//    void getThrustAccelerationProfile(
//            std::vector< double >& epochsVector,
//            std::map< double, Eigen::VectorXd >& thrustAccelerationProfile );

//    //! Compute current thrust vector.
//    Eigen::Vector3d computeCurrentThrust(
//            double time,
//            std::function< double ( const double ) > specificImpulseFunction,
//            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

//    //! Return thrust profile.
//    void getThrustProfile( std::vector< double >& epochsVector,
//           std::map< double, Eigen::VectorXd >& thrustProfile,
//           std::function< double ( const double ) > specificImpulseFunction,
//           std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );


    //! Return best individual.
    std::vector< double > getBestIndividual( )
    {
        return championDesignVariables_;
    }

    //! Return fitness of best individual.
    std::vector< double > getBestIndividualFitness( )
    {
        return championFitness_;
    }


    //! Retrieve acceleration map (thrust and central gravity accelerations).
    basic_astrodynamics::AccelerationMap retrieveLowThrustAccelerationMap( std::function< double ( const double ) > specificImpulseFunction );

//    void computeHybridMethodTrajectoryAndFullPropagation(
//         std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
//            std::shared_ptr< propagators::PropagatorSettings< double > > >& propagatorSettings,
//         std::map< double, Eigen::VectorXd >& fullPropagationResults,
//         std::map< double, Eigen::Vector6d >& hybridMethodResults,
//         std::map< double, Eigen::VectorXd>& dependentVariablesHistory );

    //! Define appropriate translational state propagator settings for the full propagation.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > createLowThrustTranslationalStatePropagatorSettings(
            basic_astrodynamics::AccelerationMap accelerationModelMap,
            std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave );


protected:

//    //! Compute current mass of the spacecraft between two epochs.
//    double computeCurrentMass(
//            const double timeInitialEpoch,
//            const double timeFinalEpoch,
//            const double massInitialEpoch,
//            std::function< double ( const double ) > specificImpulseFunction,
//            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

private:

//    //! State vector of the vehicle at the leg departure.
//    Eigen::Vector6d stateAtDeparture_;

//    //! State vector of the vehicle at the leg arrival.
//    Eigen::Vector6d stateAtArrival_;

    //! Maximum allowed thrust.
    double maximumThrust_;

    //! Specific impulse.
    double specificImpulse_;

//    //! Time of flight for the leg.
//    double timeOfFlight_;

//    //! Body map object.
//    simulation_setup::NamedBodyMap bodyMap_;

//    //! Name of the body to be propagated.
//    std::string bodyToPropagate_;

//    //! Name of the central body.
//    std::string centralBody_;

    //! Integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings_;

    //! Optimisation algorithm to be used to solve the Sims-Flanagan problem.
    pagmo::algorithm optimisationAlgorithm_;

    //! Number of generations for the optimisation algorithm.
    int numberOfGenerations_;

    //! Number of individuals per population for the optimisation algorithm.
    int numberOfIndividualsPerPopulation_;

    //! Relative tolerance for optimisation constraints.
    double relativeToleranceConstraints_;

    //! Initial guess for the optimisation.
    //! The first element contains the thrust throttles corresponding to the initial guess for the thrust model.
    //! The second element defines the bounds around the initial time (in percentage).
    std::pair< Eigen::VectorXd, double > initialGuessThrustModel_;

    //! Fitness vector of the optimisation best individual.
    std::vector< double > championFitness_;

    //! Design variables vector corresponding to the optimisation best individual.
    std::vector< double > championDesignVariables_;

    //! Initial mass of the spacecraft.
    double initialSpacecraftMass_;

    //! Hybrid method leg corresponding to the best optimisation output.
    std::shared_ptr< HybridMethodLeg > hybridMethodLeg_;

};


} // namespace low_thrust_direct_methods
} // namespace tudat

#endif // HYBRIDMETHOD_H
