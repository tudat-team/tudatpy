/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/LowThrustTrajectories/lowThrustOptimisationSetup.h"

namespace tudat
{
namespace low_thrust_trajectories
{

TrajectoryOptimisationProblem::TrajectoryOptimisationProblem(
        std::function< Eigen::Vector6d( const double ) > departureStateFunction,
        std::function< Eigen::Vector6d( const double ) > arrivalStateFunction,
        std::pair< double, double > departureTimeBounds,
        std::pair< double, double > timeOfFlightBounds,
        const std::shared_ptr< low_thrust_trajectories::LowThrustLegSettings >& lowThrustLegSettings ) :
    departureStateFunction_( departureStateFunction ),
    arrivalStateFunction_( arrivalStateFunction ),
    departureTimeBounds_( departureTimeBounds ),
    timeOfFlightBounds_( timeOfFlightBounds ),
    lowThrustLegSettings_( lowThrustLegSettings )
{
    initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();

}


//! Descriptive name of the problem
std::string TrajectoryOptimisationProblem::get_name() const {
    return "Low-thrust trajectory leg optimisation to minimise the required deltaV.";
}

//! Get bounds
std::pair< std::vector< double >, std::vector< double > > TrajectoryOptimisationProblem::get_bounds() const {

    // Define lower bounds.
    std::vector< double > lowerBounds;
    lowerBounds.push_back( departureTimeBounds_.first );
    lowerBounds.push_back( timeOfFlightBounds_.first );

    // Define upper bounds.
    std::vector< double > upperBounds;
    upperBounds.push_back( departureTimeBounds_.second );
    upperBounds.push_back( timeOfFlightBounds_.second );

    return { lowerBounds, upperBounds };
}


//! Fitness function.
std::vector< double > TrajectoryOptimisationProblem::fitness( const std::vector< double > &designVariables ) const{

    std::vector< double > fitness;

    double departureTime = designVariables[ 0 ];
    double timeOfFlight = designVariables[ 1 ];

    Eigen::Vector6d stateAtDeparture = departureStateFunction_( departureTime );
    Eigen::Vector6d stateAtArrival = arrivalStateFunction_( departureTime + timeOfFlight );

    std::shared_ptr< LowThrustLeg > lowThrustLeg = createLowThrustLeg(
                lowThrustLegSettings_, stateAtDeparture, stateAtArrival, timeOfFlight );

    fitness.push_back( lowThrustLeg->computeDeltaV( ) );

    return fitness;
}

} // namespace low_thrust_trajectories

} // namespace tudat
