/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SIMS_FLANAGAN_OPTIMISATION_SETUP_H
#define TUDAT_SIMS_FLANAGAN_OPTIMISATION_SETUP_H

#include <vector>
#include <utility>
#include <limits>

#include <Eigen/Core>

#include "pagmo/island.hpp"
#include "pagmo/io.hpp"
#include "pagmo/problem.hpp"

#include "tudat/simulation/simulation.h"


namespace tudat
{
namespace low_thrust_trajectories
{

using namespace pagmo;

struct SimsFlanaganProblem
{

    typedef Eigen::Matrix< double, 6, 1 > StateType;

    //! Default constructor, required for Pagmo compatibility
    SimsFlanaganProblem( ){ }

    //! Constructor.
    SimsFlanaganProblem( const Eigen::Vector6d& stateAtDeparture,
                         const Eigen::Vector6d& stateAtArrival,
                         const double centralBodyGravitationalParameter,
                         const double initialSpacecraftMass,
                         const double maximumThrust,
                         const std::function< double ( const double ) > specificImpulseFunction,
                         const int numberSegments,
                         const double timeOfFlight,
                         const std::pair< std::vector< double >, double > initialGuessThrustModel,
                         const double relativeToleranceConstraints = 1.0e-6 );

    //! Calculate the fitness as a function of the parameter vector x
    std::vector< double > fitness( const std::vector< double > &x ) const;

    //! Retrieve the allowable limits of the parameter vector x: pair containing minima and maxima of parameter values
    std::pair< std::vector< double >, std::vector< double > > get_bounds() const;

    //! Retrieve the name of the problem
    std::string get_name( ) const;

    //! Retrieve the number of objectives in problem, e.g. the size of the vector returned by the fitness function
    vector_double::size_type get_nobj() const
    {
        return 1u;
    }

    vector_double::size_type get_nic() const
    {
        return 0u;
    }

    vector_double::size_type get_nec() const
    {
        return 0u;
    }


protected:

private:

    //! State vector of the vehicle at the leg departure.
    Eigen::Vector6d stateAtDeparture_;

    //! State vector of the vehicle at the leg arrival.
    Eigen::Vector6d stateAtArrival_;

    double centralBodyGravitationalParameter_;

    double initialSpacecraftMass_;

    //! Maximum allowed thrust.
    double maximumThrust_;

    //! Specific impulse function.
    std::function< double ( const double ) > specificImpulseFunction_;

    //! Number of segments into which the leg is subdivided.
    int numberSegments_;

    //! Time of flight for the leg.
    double timeOfFlight_;

    //! Initial guess for the optimisation.
    //! The first element contains the thrust throttles corresponding to the initial guess for the thrust model.
    //! The second element defines the bounds around the initial time (in percentage).
    std::pair< std::vector< double >, double > initialGuessThrustModel_;

    //! Thrust throttles for the thrust model initial guess.
    std::vector< double > initialGuessThrottles_;

    //! Relative margin w.r.t. initial guess.
    double relativeMarginWrtInitialGuess_;

    //! Relative tolerance for optimisation constraints.
    double relativeToleranceConstraints_;


};

} // namespace low_thrust_trajectories

} // namespace tudat

#endif // TUDAT_SIMS_FLANAGAN_OPTIMISATION_SETUP_H
