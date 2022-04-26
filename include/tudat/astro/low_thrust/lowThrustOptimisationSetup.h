/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *
 *    Kluever (2010), Low-Thrust Trajectory Optimization Using Orbital Averaging and Control Parameterization, In: Conway,
 *    (editor) Spacecraft trajectory optimization. Cambridge University Press, 2010.
 *    Boudestijn (2014), DEVELOPMENT OF A LOW -THRUST EARTH-CENTERED TRANSFER OPTIMIZER FOR THE PRELIMINARY MISSION DESIGN PHASE,
 *    M.Sc. Thesis, Delft University of Technology
 */

#ifndef TUDAT_LOW_THRUST_OPTIMISATION_SETUP_H
#define TUDAT_LOW_THRUST_OPTIMISATION_SETUP_H

#include <Eigen/Geometry>

#include <boost/bind/bind.hpp>
#include <functional>

#include "tudat/astro/low_thrust/simsFlanagan.h"
#include "tudat/astro/low_thrust/lowThrustLegSettings.h"
#include "pagmo/island.hpp"
#include "pagmo/io.hpp"
#include "pagmo/problem.hpp"

using namespace boost::placeholders;

namespace tudat
{

namespace low_thrust_trajectories
{

/*!
 *  The class defined in this file is to be used in a Pagmo optimization.
 */

using namespace pagmo;

//! Definition optimisation problem to minimise the total deltaV required by the trajectory.
struct TrajectoryOptimisationProblem
{

    typedef Eigen::Matrix< double, 6, 1 > StateType;

    //! Default constructor, required for Pagmo compatibility
    TrajectoryOptimisationProblem( ){ }

    //! Constructor.
    TrajectoryOptimisationProblem(
            std::function< Eigen::Vector6d( const double ) > departureStateFunction,
            std::function< Eigen::Vector6d( const double ) > arrivalStateFunction,
            std::pair< double, double > departureTimeBounds,
            std::pair< double, double > timeOfFlightBounds,
            const std::shared_ptr< low_thrust_trajectories::LowThrustLegSettings >& lowThrustLegSettings );

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

    //! Function returning the state at departure as a function of the departure time.
    std::function< Eigen::Vector6d( const double ) > departureStateFunction_;

    //! Function returning the state at arrival as a function of the arrival time.
    std::function< Eigen::Vector6d( const double ) > arrivalStateFunction_;

    //! Initial spacecraft mass.
    double initialSpacecraftMass_;

    //! Bounds for departure time.
    std::pair< double, double > departureTimeBounds_;

    //! Bounds for time of flight.
    std::pair< double, double > timeOfFlightBounds_;

    //! Pointer to settings for low-thrust trajectory leg.
    std::shared_ptr< low_thrust_trajectories::LowThrustLegSettings > lowThrustLegSettings_;


};

} // namespace low_thrust_trajectories

} // namespace tudat

#endif // TUDAT_LOW_THRUST_OPTIMISATION_SETUP_H
