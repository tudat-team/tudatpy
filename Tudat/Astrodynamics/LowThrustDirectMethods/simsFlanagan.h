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

#ifndef SIMSFLANAGAN_H
#define SIMSFLANAGAN_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <map>
//#include "pagmo/algorithm.hpp"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/simsFlanaganLeg.h"
//#include "Tudat/Astrodynamics/LowThrustDirectMethods/lowThrustLeg.h"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/optimisationSettings.h"

//#include "tudatExampleApplications/libraryExamples/PaGMOEx/Problems/applicationOutput.h"
//#include "tudatExampleApplications/libraryExamples/PaGMOEx/Problems/getAlgorithm.h"
//#include "tudatExampleApplications/libraryExamples/PaGMOEx/Problems/saveOptimizationResults.h"

namespace tudat
{
namespace low_thrust_direct_methods
{

//! Transform thrust model as a function of time into Sims Flanagan thrust model.
std::vector< Eigen::Vector3d > convertToSimsFlanaganThrustModel( std::function< Eigen::Vector3d( const double ) > thrustModelWrtTime,
                                                                 const double maximumThrust,
                                                                 const double timeOfFlight, const int numberSegmentsForwardPropagation,
                                                                 const int numberSegmentsBackwardPropagation );


class SimsFlanagan : public transfer_trajectories::LowThrustLeg
{
public:

    //! Constructor.
    SimsFlanagan(
            const Eigen::Vector6d& stateAtDeparture,
            const Eigen::Vector6d& stateAtArrival,
            const double maximumThrust,
            const std::function< double ( const double ) > specificImpulseFunction,
            const int numberSegments,
            const double timeOfFlight,
            simulation_setup::NamedBodyMap bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
//            pagmo::algorithm optimisationAlgorithm,
//            const int numberOfGenerations,
//            const int numberOfIndividualsPerPopulation,
            std::shared_ptr< transfer_trajectories::OptimisationSettings > optimisationSettings ) :
//            const double relativeToleranceConstraints = 1.0e-6,
//            std::pair< std::function< Eigen::Vector3d( const double ) >, double > initialGuessThrustModel = std::make_pair( nullptr, 0.0 ) ) :
        LowThrustLeg( stateAtDeparture, stateAtArrival, timeOfFlight, bodyMap, bodyToPropagate, centralBody ),
//        stateAtDeparture_( stateAtDeparture ),
//        stateAtArrival_( stateAtArrival ),
        maximumThrust_( maximumThrust ),
        specificImpulseFunction_( specificImpulseFunction ),
        numberSegments_( numberSegments ),
//        timeOfFlight_( timeOfFlight ),
//        bodyMap_( bodyMap ),
//        bodyToPropagate_( bodyToPropagate ),
//        centralBody_( centralBody ),
//        optimisationAlgorithm_( optimisationAlgorithm ),
//        numberOfGenerations_( numberOfGenerations ),
//        numberOfIndividualsPerPopulation_( numberOfIndividualsPerPopulation ),
//        relativeToleranceConstraints_( relativeToleranceConstraints ),
        optimisationSettings_( optimisationSettings )//,
//        initialGuessThrustModel_( initialGuessThrustModel )
    {

        // Store initial spacecraft mass.
        initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();

        // Calculate number of segments for both the forward propagation (from departure to match point)
        // and the backward propagation (from arrival to match point).
        numberSegmentsForwardPropagation_ = ( numberSegments_ + 1 ) / 2;
        numberSegmentsBackwardPropagation_ = numberSegments_ / 2;

        // Convert the thrust model proposed as initial guess into simplified thrust model adapted to the Sims-Flanagan method.
        if ( optimisationSettings_->initialGuessThrustModel_.first /*initialGuessThrustModel.first*/ != nullptr )
        {
            initialGuessThrustModel_.first = convertToSimsFlanaganThrustModel( optimisationSettings_->initialGuessThrustModel_.first /* initialGuessThrustModel.first*/,
                                                                               maximumThrust_, timeOfFlight_,
                                                                               numberSegmentsForwardPropagation_, numberSegmentsBackwardPropagation_ );
        }
        else
        {
            initialGuessThrustModel_.first = std::vector< Eigen::Vector3d >( );
        }
        initialGuessThrustModel_.second = optimisationSettings_->initialGuessThrustModel_.second; // initialGuessThrustModel.second;

        // Perform optimisation
        std::pair< std::vector< double >, std::vector< double > > bestIndividual = performOptimisation( );
        championFitness_ = bestIndividual.first;
        championDesignVariables_ = bestIndividual.second;

        // Transform best design variables into vector of throttles.
        std::vector< Eigen::Vector3d > throttles;
        for ( int i = 0 ; i < numberSegments ; i++ )
        {
            throttles.push_back( ( Eigen::Vector3d( ) << championDesignVariables_[ i * 3 ], championDesignVariables_[ i * 3 + 1 ],
                    championDesignVariables_[ i * 3 + 2 ] ).finished( ) );
        }

        bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );

        // Create Sims-Flanagan leg from the best optimisation individual.
        simsFlanaganLeg_ = std::make_shared< SimsFlanaganLeg >( stateAtDeparture_, stateAtArrival_, maximumThrust_, specificImpulseFunction_,
                                            timeOfFlight_, bodyMap_, throttles, bodyToPropagate_, centralBody_ );

    }

    //! Default destructor.
    ~SimsFlanagan( ) { }

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
        return simsFlanaganLeg_->getTotalDeltaV( ); // deltaV_; //championFitness_[ 0 ];
    }

    //! Compute current cartesian state.
    Eigen::Vector6d computeCurrentStateVector( const double currentTime );

    //! Compute state history.
    void getTrajectory(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::Vector6d >& propagatedTrajectory )
    {
        simsFlanaganLeg_->propagateTrajectory( epochsVector, propagatedTrajectory );
    }

    Eigen::Vector3d computeCurrentThrust( double time,
                                          std::function< double ( const double ) > specificImpulseFunction,
                                          std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Compute direction thrust acceleration in cartesian coordinates.
    Eigen::Vector3d computeCurrentThrustAccelerationDirection(
            double currentTime, std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Compute magnitude thrust acceleration.
    double computeCurrentThrustAccelerationMagnitude(
            double currentTime, std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Return thrust acceleration profile.
    void getThrustAccelerationProfile(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::VectorXd >& thrustAccelerationProfile,
            std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

//    //! Return mass profile.
//    void getMassProfile(
//            std::vector< double >& epochsVector,
//            std::map< double, Eigen::VectorXd >& massProfile,
//            std::function< double ( const double ) > specificImpulseFunction,
//            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

//    //! Return thrust profile.
//    void getThrustProfile(
//            std::vector< double >& epochsVector,
//            std::map< double, Eigen::VectorXd >& thrustProfile,
//            std::function< double ( const double ) > specificImpulseFunction,
//            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
//    {

//    }

//    //! Return thrust acceleration profile.
//    void getThrustAccelerationProfile(
//            std::vector< double >& epochsVector,
//            std::map< double, Eigen::VectorXd >& thrustAccelerationProfile );


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


    //! Return best Sims-Flanagan leg after optimisation.
    std::shared_ptr< SimsFlanaganLeg > getOptimalSimsFlanaganLeg( )
    {
        return simsFlanaganLeg_;
    }


    //! Retrieve acceleration map (thrust and central gravity accelerations).
    basic_astrodynamics::AccelerationMap retrieveLowThrustAccelerationMap(
            std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings );


    //! Define appropriate translational state propagator settings for the full propagation.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > createLowThrustTranslationalStatePropagatorSettings(
            basic_astrodynamics::AccelerationMap accelerationModelMap,
            std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave );


protected:

    int convertTimeToLegSegment( double currentTime )
    {
        return simsFlanaganLeg_->convertTimeToLegSegment( currentTime );
    }

private:

//    //! State vector of the vehicle at the leg departure.
//    Eigen::Vector6d stateAtDeparture_;

//    //! State vector of the vehicle at the leg arrival.
//    Eigen::Vector6d stateAtArrival_;

    //! Maximum allowed thrust.
    double maximumThrust_;

    //! Specific impulse function.
    std::function< double ( const double ) > specificImpulseFunction_;

    //! Number of segments into which the leg is subdivided.
    int numberSegments_;

//    //! Time of flight for the leg.
//    double timeOfFlight_;

//    //! Body map object.
//    simulation_setup::NamedBodyMap bodyMap_;

//    //! Name of the body to be propagated.
//    std::string bodyToPropagate_;

//    //! Name of the central body.
//    std::string centralBody_;

//    //! Optimisation algorithm to be used to solve the Sims-Flanagan problem.
//    pagmo::algorithm optimisationAlgorithm_;

//    //! Number of generations for the optimisation algorithm.
//    int numberOfGenerations_;

//    //! Number of individuals per population for the optimisation algorithm.
//    int numberOfIndividualsPerPopulation_;

//    //! Relative tolerance for optimisation constraints.
//    double relativeToleranceConstraints_;

    //! Initial guess for the optimisation.
    //! The first element contains the thrust throttles corresponding to the initial guess for the thrust model.
    //! The second element defines the bounds around the initial time (in percentage).
    std::pair< std::vector< Eigen::Vector3d >, double > initialGuessThrustModel_;

    //! Optimisation settings.
    std::shared_ptr< transfer_trajectories::OptimisationSettings > optimisationSettings_;

    //! Fitness vector of the optimisation best individual.
    std::vector< double > championFitness_;

    //! Design variables vector corresponding to the optimisation best individual.
    std::vector< double > championDesignVariables_;

    //! Initial mass of the spacecraft.
    double initialSpacecraftMass_;

    //! Number of segments for the forward propagation from departure to match point.
    int numberSegmentsForwardPropagation_;

    //! Number of segments for the backward propagation from arrival to match point.
    int numberSegmentsBackwardPropagation_;

    //! DeltaV corresponding to best trajectory.
    double deltaV_;

    //! Sims-Flanagan leg corresponding to the best optimisation output.
    std::shared_ptr< SimsFlanaganLeg > simsFlanaganLeg_;

};


} // namespace low_thrust_direct_methods
} // namespace tudat

#endif // SIMSFLANAGAN_H
