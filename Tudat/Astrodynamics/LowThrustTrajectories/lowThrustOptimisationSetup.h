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

#include <boost/bind.hpp>
#include <functional>

#include "Tudat/Astrodynamics/LowThrustTrajectories/hybridMethod.h"
//#include "Tudat/Astrodynamics/LowThrustTrajectories/hybridMethodLeg.h"
//#include "Tudat/Astrodynamics/LowThrustTrajectories/hybridOptimisationSetup.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/simsFlanagan.h"
//#include "Tudat/Astrodynamics/LowThrustTrajectories/simsFlanaganLeg.h"
//#include "Tudat/Astrodynamics/LowThrustTrajectories/simsFlanaganOptimisationSetup.h"
//#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsHodographicShaping.h"
//#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/hodographicShaping.h"
//#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsSphericalShaping.h"
//#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionSphericalShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/sphericalShaping.h"

//#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustLeg.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustLegSettings.h"

//#include "pagmo/algorithm.hpp"
#include "pagmo/island.hpp"
#include "pagmo/io.hpp"
#include "pagmo/serialization.hpp"
#include "pagmo/problem.hpp"


namespace tudat
{

namespace low_thrust_trajectories
{

/*!
 *  The class defined in this file is to be used in a Pagmo optimization.
 */

using namespace pagmo;
//using namespace tudat;

//! Definition optimisation problem to minimise the total deltaV required by the trajectory.
struct TrajectoryOptimisationProblem
{

    typedef Eigen::Matrix< double, 6, 1 > StateType;

    //! Default constructor, required for Pagmo compatibility
    TrajectoryOptimisationProblem( ){ }

    //! Constructor.
    TrajectoryOptimisationProblem(
            simulation_setup::NamedBodyMap bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
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
        return 0u; //numberSegments_;
    }

    vector_double::size_type get_nec() const
    {
        return 0u; //6u;
    }

//    //! Serialization function for Pagmo compatibility
//    template <typename Archive>
//    void serialize(Archive &ar)
//    {
//        ar(problemBounds_);
//    }

protected:

private:

//    //! State vector of the vehicle at the leg departure.
//    Eigen::Vector6d stateAtDeparture_;

//    //! State vector of the vehicle at the leg arrival.
//    Eigen::Vector6d stateAtArrival_;

//    //! Maximum allowed thrust.
//    double maximumThrust_;

//    //! Specific impulse function.
//    std::function< double ( const double ) > specificImpulseFunction_;

//    //! Number of segments into which the leg is subdivided.
//    int numberSegments_;

//    //! Time of flight for the leg.
//    double timeOfFlight_;

    //! Body map.
    mutable simulation_setup::NamedBodyMap bodyMap_;

    //! Name of the body to be propagated.
    const std::string bodyToPropagate_;

    //! Name of the central body.
    const std::string centralBody_;

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

//    //! Initial guess for the optimisation.
//    //! The first element contains the thrust throttles corresponding to the initial guess for the thrust model.
//    //! The second element defines the bounds around the initial time (in percentage).
//    std::pair< std::vector< Eigen::Vector3d >, double > initialGuessThrustModel_;

//    //! Thrust throttles for the thrust model initial guess.
//    std::vector< Eigen::Vector3d > initialGuessThrottles_;

//    //! Relative margin w.r.t. initial guess.
//    double relativeMarginWrtInitialGuess_;

//    //! Relative tolerance for optimisation constraints.
//    double relativeToleranceConstraints_;


};

} // namespace low_thrust_trajectories

} // namespace tudat

#endif // TUDAT_LOW_THRUST_OPTIMISATION_SETUP_H
