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

namespace tudat
{
namespace low_thrust_direct_methods
{


class HybridMethod
{
public:

    //! Constructor.
    HybridMethod(
            const Eigen::Vector6d& stateAtDeparture,
            const Eigen::Vector6d& stateAtArrival,
            const double maximumThrust,
            const std::function< double ( const double ) > specificImpulseFunction,
            const double timeOfFlight,
            simulation_setup::NamedBodyMap bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
            pagmo::algorithm optimisationAlgorithm,
            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
            const double relativeToleranceConstraints = 1.0e-6 ) :
        stateAtDeparture_( stateAtDeparture ),
        stateAtArrival_( stateAtArrival ),
        maximumThrust_( maximumThrust ),
        specificImpulseFunction_( specificImpulseFunction ),
        timeOfFlight_( timeOfFlight ),
        bodyMap_( bodyMap ),
        bodyToPropagate_( bodyToPropagate ),
        centralBody_( centralBody ),
        optimisationAlgorithm_( optimisationAlgorithm ),
        integratorSettings_( integratorSettings ),
        relativeToleranceConstraints_( relativeToleranceConstraints )
    {

        // Store initial spacecraft mass.
        initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();

//        // Perform optimisation
//        std::pair< std::vector< double >, std::vector< double > > bestIndividual = performOptimisation( );
//        championFitness_ = bestIndividual.first;
//        championDesignVariables_ = bestIndividual.second;

    }

    //! Default destructor.
    ~HybridMethod( ) { }

    //! Perform optimisation.
    std::pair< std::vector< double >, std::vector< double > > performOptimisation( );

    //! Compute DeltaV.
    double computeDeltaV( )
    {
        return championFitness_[ 0 ];
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


    void computeHybridMethodTrajectoryAndFullPropagation(
         std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > >& propagatorSettings,
         std::map< double, Eigen::VectorXd >& fullPropagationResults,
         std::map< double, Eigen::Vector6d >& hybridMethodResults,
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

    //! Integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings_;

    //! Relative tolerance for optimisation constraints.
    double relativeToleranceConstraints_;

    //! Fitness vector of the optimisation best individual.
    std::vector< double > championFitness_;

    //! Design variables vector corresponding to the optimisation best individual.
    std::vector< double > championDesignVariables_;

    //! Initial mass of the spacecraft.
    double initialSpacecraftMass_;

};


} // namespace low_thrust_direct_methods
} // namespace tudat

#endif // HYBRIDMETHOD_H
