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
#include "pagmo/algorithm.hpp"

namespace tudat
{
namespace low_thrust_direct_methods
{

//! Transform thrust model as a function of time into Sims Flanagan thrust model.
std::vector< Eigen::Vector3d > convertToSimsFlanaganThrustModel( std::function< Eigen::Vector3d( const double ) > thrustModelWrtTime,
                                                                 const double maximumThrust,
                                                                 const double timeOfFlight, const int numberSegmentsForwardPropagation,
                                                                 const int numberSegmentsBackwardPropagation );


class SimsFlanagan
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
            pagmo::algorithm optimisationAlgorithm,
            const int numberOfGenerations,
            const int numberOfIndividualsPerPopulation,
            const double relativeToleranceConstraints = 1.0e-6,
            std::pair< std::function< Eigen::Vector3d( const double ) >, double > initialGuessThrustModel = std::make_pair( nullptr, 0.0 ) ) :
        stateAtDeparture_( stateAtDeparture ),
        stateAtArrival_( stateAtArrival ),
        maximumThrust_( maximumThrust ),
        specificImpulseFunction_( specificImpulseFunction ),
        numberSegments_( numberSegments ),
        timeOfFlight_( timeOfFlight ),
        bodyMap_( bodyMap ),
        bodyToPropagate_( bodyToPropagate ),
        centralBody_( centralBody ),
        optimisationAlgorithm_( optimisationAlgorithm ),
        numberOfGenerations_( numberOfGenerations ),
        numberOfIndividualsPerPopulation_( numberOfIndividualsPerPopulation ),
        relativeToleranceConstraints_( relativeToleranceConstraints ) //,
//        initialGuessThrustModel_( initialGuessThrustModel )
    {

        // Store initial spacecraft mass.
        initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();

        // Calculate number of segments for both the forward propagation (from departure to match point)
        // and the backward propagation (from arrival to match point).
        numberSegmentsForwardPropagation_ = ( numberSegments_ + 1 ) / 2;
        numberSegmentsBackwardPropagation_ = numberSegments_ / 2;

        // Convert the thrust model proposed as initial guess into simplified thrust model adapted to the Sims-Flanagan method.
        if ( initialGuessThrustModel.first != nullptr )
        {
            initialGuessThrustModel_.first = convertToSimsFlanaganThrustModel( initialGuessThrustModel.first, maximumThrust_, timeOfFlight_,
                                                                               numberSegmentsForwardPropagation_, numberSegmentsBackwardPropagation_ );
        }
        else
        {
            initialGuessThrustModel_.first = std::vector< Eigen::Vector3d >( );
        }
        initialGuessThrustModel_.second = initialGuessThrustModel.second;

        // Perform optimisation
        std::pair< std::vector< double >, std::vector< double > > bestIndividual = performOptimisation( );
        championFitness_ = bestIndividual.first;
        championDesignVariables_ = bestIndividual.second;

    }

    //! Default destructor.
    ~SimsFlanagan( ) { }

    //! Perform optimisation.
    std::pair< std::vector< double >, std::vector< double > > performOptimisation( );

    //! Compute DeltaV.
    double computeDeltaV( )
    {
        return deltaV_; //championFitness_[ 0 ];
    }

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



    //! Function to compute the Sims Flanagan trajectory and the propagation fo the full problem.
    void computeSimsFlanaganTrajectoryAndFullPropagation(
         std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
         std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > >& propagatorSettings,
         std::map< double, Eigen::VectorXd >& fullPropagationResults,
         std::map< double, Eigen::Vector6d >& SimsFlanaganResults,
         std::map< double, Eigen::VectorXd>& dependentVariablesHistory );


protected:

private:

    //! State vector of the vehicle at the leg departure.
    Eigen::Vector6d stateAtDeparture_;

    //! State vector of the vehicle at the leg arrival.
    Eigen::Vector6d stateAtArrival_;

    //! Maximum allowed thrust.
    double maximumThrust_;

    //! Specific impulse function.
    std::function< double ( const double ) > specificImpulseFunction_;

    //! Number of segments into which the leg is subdivided.
    int numberSegments_;

    //! Time of flight for the leg.
    double timeOfFlight_;

    //! Body map object.
    simulation_setup::NamedBodyMap bodyMap_;

    //! Name of the body to be propagated.
    std::string bodyToPropagate_;

    //! Name of the central body.
    std::string centralBody_;

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
    std::pair< std::vector< Eigen::Vector3d >, double > initialGuessThrustModel_;

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

};


} // namespace low_thrust_direct_methods
} // namespace tudat

#endif // SIMSFLANAGAN_H
