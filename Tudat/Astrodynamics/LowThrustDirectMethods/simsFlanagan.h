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

namespace tudat
{
namespace low_thrust_direct_methods
{


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
            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
            const propagators::TranslationalPropagatorType propagatorType = propagators::cowell,
            const bool useHighOrderSolution = false ) :
        stateAtDeparture_( stateAtDeparture ),
        stateAtArrival_( stateAtArrival ),
        maximumThrust_( maximumThrust ),
        specificImpulseFunction_( specificImpulseFunction ),
        numberSegments_( numberSegments ),
        timeOfFlight_( timeOfFlight ),
        bodyMap_( bodyMap ),
        bodyToPropagate_( bodyToPropagate ),
        centralBody_( centralBody ),
        integratorSettings_( integratorSettings ),
        propagatorType_( propagatorType ),
        useHighOrderSolution_( useHighOrderSolution )
    {

        // Store initial spacecraft mass.
        initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();

//        // Perform optimisation
//        std::pair< std::vector< double >, std::vector< double > > bestIndividual = performOptimisation( );
//        championFitness_ = bestIndividual.first;
//        championDesignVariables_ = bestIndividual.second;

    }

    //! Default destructor.
    ~SimsFlanagan( ) { }

    //! Perform optimisation.
    std::pair< std::vector< double >, std::vector< double > > performOptimisation( );

    //! Compute DeltaV.
    double computeDeltaV( )
    {
        return championFitness_[ 0 ];
    }

    //! Function to compute the Sims Flanagan trajectory and the propagation fo the full problem.
    void computeSimsFlanaganTrajectoryAndFullPropagation(
         std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > >& propagatorSettings,
         std::map< double, Eigen::VectorXd >& fullPropagationResults,
         std::map< double, Eigen::VectorXd >& SimsFlanaganResults,
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

    //! Integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings_;

    //! Propagator type.
    propagators::TranslationalPropagatorType propagatorType_;

    //! Boolean defining which of the low or high order solutions is used.
    bool useHighOrderSolution_;

    //! Fitness vector of the optimisation best individual.
    std::vector< double > championFitness_;

    //! Design variables vector corresponding to the optimisation best individual.
    std::vector< double > championDesignVariables_;

    //! Initial mass of the spacecraft.
    double initialSpacecraftMass_;

};


} // namespace low_thrust_direct_methods
} // namespace tudat

#endif // SIMSFLANAGAN_H
