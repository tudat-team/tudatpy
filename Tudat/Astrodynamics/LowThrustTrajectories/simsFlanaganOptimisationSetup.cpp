/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/LowThrustTrajectories/simsFlanaganOptimisationSetup.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/simsFlanaganModel.h"

namespace tudat
{
namespace low_thrust_trajectories
{

SimsFlanaganProblem::SimsFlanaganProblem(
        const Eigen::Vector6d &stateAtDeparture,
        const Eigen::Vector6d &stateAtArrival,
        const double centralBodyGravitationalParameter,
        const double initialSpacecraftMass,
        const double maximumThrust,
        const std::function< double ( const double ) > specificImpulseFunction,
        const int numberSegments,
        const double timeOfFlight,
        const std::pair< std::vector< double >, double > initialGuessThrustModel,
        const double relativeToleranceConstraints ) :
    stateAtDeparture_( stateAtDeparture ),
    stateAtArrival_( stateAtArrival ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    initialSpacecraftMass_( initialSpacecraftMass ),
    maximumThrust_( maximumThrust ),
    specificImpulseFunction_( specificImpulseFunction ),
    numberSegments_( numberSegments ),
    timeOfFlight_( timeOfFlight ),
    initialGuessThrustModel_( initialGuessThrustModel ),
    relativeToleranceConstraints_( relativeToleranceConstraints )
{
    // Retrieve initial guess.
    initialGuessThrottles_ = initialGuessThrustModel_.first;
    relativeMarginWrtInitialGuess_ = initialGuessThrustModel_.second;

    if ( ( relativeMarginWrtInitialGuess_ != TUDAT_NAN ) && ( ( relativeMarginWrtInitialGuess_ < 0.0) || ( relativeMarginWrtInitialGuess_ > 1.0 ) ) )
    {
        std::string errorMessage = "Error when using initial guess to create Sims-Flanagan optimisation problem, the relative margin w.r.t. the "
                                   "initial guess should be constrained between [0,1], the current value is "
                + std::to_string( relativeMarginWrtInitialGuess_ ) + ".";

        throw std::runtime_error( errorMessage );
    }
}


//! Descriptive name of the problem
std::string SimsFlanaganProblem::get_name() const {
    return "Sims-Flanagan direct method to compute a low-thrust trajectory";
}

//! Get bounds
std::pair< std::vector< double >, std::vector< double > > SimsFlanaganProblem::get_bounds() const {

    // Define lower bounds.
    std::vector< double > lowerBounds;
    std::vector< double > upperBounds;

    if ( initialGuessThrottles_.size( ) != 0 )
    {
        if ( initialGuessThrottles_.size( ) != static_cast< unsigned int >( 3 * numberSegments_ ) )
        {
            throw std::runtime_error( "Error when providing an initial guess for Sims-Flanagan, vector size unconsistent with"
                                      " number of segments into which the trajectory is subdivided." );
        }
        else
        {
            for ( int i = 0 ; i < numberSegments_ ; i++ )
            {
                for ( int k = 0 ; k < 3 ; k++ )
                {
                    double lowerBoundFromInitialGuess;
                    double upperBoundFromInitialGuess;
                    if ( initialGuessThrottles_[ i * 3 + k ] >= 0 )
                    {
                        lowerBoundFromInitialGuess = ( 1.0 - relativeMarginWrtInitialGuess_ ) * initialGuessThrottles_[ i * 3 + k ];
                        upperBoundFromInitialGuess = ( 1.0 + relativeMarginWrtInitialGuess_ ) * initialGuessThrottles_[ i * 3 + k ];
                    }
                    else
                    {
                        lowerBoundFromInitialGuess = ( 1.0 + relativeMarginWrtInitialGuess_ ) * initialGuessThrottles_[ i * 3 + k ];
                        upperBoundFromInitialGuess = ( 1.0 - relativeMarginWrtInitialGuess_ ) * initialGuessThrottles_[ i * 3 + k ];
                    }

                    if ( std::fabs( lowerBoundFromInitialGuess <= 1.0 ) )
                    {
                        lowerBounds.push_back( lowerBoundFromInitialGuess );
                    }
                    else
                    {
                        lowerBounds.push_back( - 1.0 );
                    }

                    if ( std::fabs( upperBoundFromInitialGuess <= 1.0 ) )
                    {
                        upperBounds.push_back( upperBoundFromInitialGuess );
                    }
                    else
                    {
                        upperBounds.push_back( 1.0 );
                    }
                }
            }
        }
    }
    else
    {
        for ( int i = 0 ; i < 3 * numberSegments_ ; i++ )
        {
            lowerBounds.push_back( - 1.0 );
            upperBounds.push_back( 1.0 );
        }
    }

    return { lowerBounds, upperBounds };
}

//! Fitness function.
std::vector< double > SimsFlanaganProblem::fitness( const std::vector< double > &designVariables ) const
{
    // Transform vector of design variables into 3D vector of throttles.
    std::vector< Eigen::Vector3d > throttles;

    // Check consistency of the size of the design variables vector.
    if ( designVariables.size( ) != static_cast< unsigned int >( 3 * numberSegments_ ) )
    {
        throw std::runtime_error( "Error, size of the design variables vector unconsistent with number of segments." );
    }

    for ( int i = 0 ; i < numberSegments_ ; i++ )
    {
        throttles.push_back( ( Eigen::Vector3d( ) << designVariables[ i * 3 ],
                             designVariables[ i * 3 + 1 ], designVariables[ i * 3 + 2 ] ).finished( ) );
    }

    std::vector< double > fitness;

    // Create Sims Flanagan trajectory leg.
    low_thrust_trajectories::SimsFlanaganModel currentLeg = low_thrust_trajectories::SimsFlanaganModel(
                stateAtDeparture_, stateAtArrival_, centralBodyGravitationalParameter_, initialSpacecraftMass_,
                maximumThrust_, specificImpulseFunction_, timeOfFlight_, throttles );

    // Forward propagation from departure to match point.
    currentLeg.propagateForwardFromDepartureToMatchPoint();

    // Backward propagation from arrival to match point.
    currentLeg.propagateBackwardFromArrivalToMatchPoint();

    // Fitness
    double deltaV = currentLeg.getTotalDeltaV( );

    // Equality constraints (must be ... = 0 )
    std::vector< double > equalityConstraints;

    // Spacecraft position and velocity must be continuous at the match point.
    for ( int i = 0 ; i < 3 ; i++ )
    {
        equalityConstraints.push_back( std::fabs( currentLeg.getStateAtMatchPointForwardPropagation( )[ i ]
                                       - currentLeg.getStateAtMatchPointBackwardPropagation( )[ i ] ) /
                                       physical_constants::ASTRONOMICAL_UNIT );
        equalityConstraints.push_back( std::fabs( currentLeg.getStateAtMatchPointForwardPropagation( )[ i + 3 ]
                                       - currentLeg.getStateAtMatchPointBackwardPropagation( )[ i + 3 ] ) /
                                       stateAtDeparture_.segment(3 ,3 ).norm() );
    }


    // Inequality constraints.
    std::vector< double > inequalityConstraints;

    // Magnitude of the normalised thrust vector must be inferior or equal to 1 )
    for ( unsigned int currentThrottle = 0 ; currentThrottle < throttles.size() ; currentThrottle++ )
    {
        if ( throttles[ currentThrottle ].norm( ) > 1.0 )
        {
            inequalityConstraints.push_back( throttles[ currentThrottle ][ 0 ] * throttles[ currentThrottle ][ 0 ]
                    + throttles[ currentThrottle ][ 1 ] * throttles[ currentThrottle ][ 1 ]
                    + throttles[ currentThrottle ][ 2 ] * throttles[ currentThrottle ][ 2 ] - 1.0 );
        }
        else
        {
            inequalityConstraints.push_back( 0.0 );
        }
    }

    // Compute auxiliary variables.
    Eigen::VectorXd c; c.resize( equalityConstraints.size( ) + inequalityConstraints.size( ) );
    Eigen::VectorXd r; r.resize( equalityConstraints.size( ) + inequalityConstraints.size( ) );
    Eigen::VectorXd epsilon; epsilon.resize( equalityConstraints.size( ) + inequalityConstraints.size( ) );
    for ( unsigned int i = 0 ; i < equalityConstraints.size( ) ; i++ )
    {
        c[ i ] = 1.0 / relativeToleranceConstraints_;
        r[ i ] = 1.0 - relativeToleranceConstraints_ * c[ i ];
        epsilon[ i ] = equalityConstraints[ i ] * c[ i ] + r[ i ];
    }
    for ( unsigned int i = 0 ; i < inequalityConstraints.size( ) ; i++ )
    {
        c[ equalityConstraints.size( ) + i ] = 1.0 / relativeToleranceConstraints_;
        r[ equalityConstraints.size( ) + i ] = 1.0 - relativeToleranceConstraints_ * c[ equalityConstraints.size( ) + i ];
        epsilon[ equalityConstraints.size( ) + i ] = inequalityConstraints[ i ] * c[ equalityConstraints.size( ) + i ] + r[ equalityConstraints.size( ) + i ];
    }

    double weightDeltaV = 1.0;
    double weightConstraints = 10.0;
    double optimisationObjective = weightDeltaV * deltaV + weightConstraints * ( epsilon.norm( ) * epsilon.norm( ) );


    // Optimisation objectives
    fitness.push_back( optimisationObjective );


    return fitness;
}

} // namespace low_thrust_trajectories

} // namespace tudat


