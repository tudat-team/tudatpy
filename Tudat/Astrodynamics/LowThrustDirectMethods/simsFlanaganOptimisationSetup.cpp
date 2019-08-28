/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/LowThrustDirectMethods/simsFlanaganOptimisationSetup.h"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/simsFlanaganLeg.h"
//#include "Tudat/SimulationSetup/tudatSimulationHeader.h"

using namespace tudat;

SimsFlanaganProblem::SimsFlanaganProblem(
        const Eigen::Vector6d &stateAtDeparture,
        const Eigen::Vector6d &stateAtArrival,
        const double maximumThrust,
        const std::function< double ( const double ) > specificImpulseFunction,
        const int numberSegments,
        const double timeOfFlight,
        simulation_setup::NamedBodyMap bodyMap,
        const std::string bodyToPropagate,
        const std::string centralBody,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const propagators::TranslationalPropagatorType propagatorType,
//        const bool useHighOrderSolution,
        const bool optimiseTimeOfFlight,
        const std::pair< double, double > timeOfFlightBounds ) :
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
    optimiseTimeOfFlight_( optimiseTimeOfFlight ),
    timeOfFlightBounds_( timeOfFlightBounds )
{
    initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();
}


//! Descriptive name of the problem
std::string SimsFlanaganProblem::get_name() const {
    return "Sims-Flanagan direct method to compute a low-thrust trajectory";
}

//! Get bounds
std::pair< std::vector< double >, std::vector< double > > SimsFlanaganProblem::get_bounds() const {

    // Define lower bounds.
    std::vector< double > lowerBounds;
    for ( int i = 0 ; i < numberSegments_ ; i++ )
    {
        // Lower bound for the 3 components of the thrust vector.
        lowerBounds.push_back( - 1.0 );
        lowerBounds.push_back( - 1.0 );
        lowerBounds.push_back( - 1.0 );
    }

    // Define upper bounds.
    std::vector< double > upperBounds;
    for ( int i = 0 ; i < numberSegments_ ; i++ )
    {
        // Upper bound for the 3 components of the thrust vector.
        upperBounds.push_back( 1.0 );
        upperBounds.push_back( 1.0 );
        upperBounds.push_back( 1.0 );
    }

    // Add lower and upper bounds for time of flight optimisation if necessary.
    if ( optimiseTimeOfFlight_ )
    {
        if ( timeOfFlightBounds_.first == TUDAT_NAN || timeOfFlightBounds_.second == TUDAT_NAN )
        {
            throw std::runtime_error( "Error when trying to optimise time of flight in Sims-Flanagan method, boundaries for time of flight"
                                      "are not defined." );
        }

        lowerBounds.push_back( timeOfFlightBounds_.first );
        upperBounds.push_back( timeOfFlightBounds_.second );
    }

    return { lowerBounds, upperBounds };
}

//! Fitness function.
std::vector< double > SimsFlanaganProblem::fitness( const std::vector< double > &designVariables ) const{

    // Re-initialise mass of the spacecraft.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );

    // Transform vector of design variables into 3D vector of throttles.
    std::vector< Eigen::Vector3d > throttles;

    // Check consistency of the size of the design variables vector.
    if ( ( ( designVariables.size( ) != 3 * numberSegments_ ) && ( !optimiseTimeOfFlight_ ) )
         || ( ( designVariables.size() != 3 * numberSegments_ + 1 ) && ( optimiseTimeOfFlight_ ) ) ) /*( designVariables.size() / 3 != numberSegments_ ) || ( designVariables.size() % 3 != 0 ) )*/
    {
        throw std::runtime_error( "Error, size of the design variables vector unconsistent with number of segments." );
    }

    for ( unsigned int i = 0 ; i < numberSegments_ ; i++ )
    {
        throttles.push_back( ( Eigen::Vector3d( ) << designVariables[ i * 3 ],
                             designVariables[ i * 3 + 1 ], designVariables[ i * 3 + 2 ] ).finished( ) );
    }

    // Update timeOfFlight_ if minimising the time of flight is one of the optimisation objectives.
    double timeOfFlight = timeOfFlight_;
    if ( optimiseTimeOfFlight_ )
    {
        timeOfFlight = designVariables[ 3 * numberSegments_ ];
    }

    std::vector< double > fitness;

    // Create Sims Flanagan trajectory leg.
    low_thrust_direct_methods::SimsFlanaganLeg currentLeg = low_thrust_direct_methods::SimsFlanaganLeg( stateAtDeparture_,
                                                                                                        stateAtArrival_,
                                                                                                        maximumThrust_,
                                                                                                        specificImpulseFunction_,
                                                                                                        timeOfFlight,
                                                                                                        bodyMap_,
                                                                                                        throttles,
                                                                                                        bodyToPropagate_,
                                                                                                        centralBody_ );

    // Forward propagation from departure to match point.
        currentLeg.propagateForwardFromDepartureToMatchPoint();

//    // Backward propagation from arrival to match point.
        currentLeg.propagateBackwardFromArrivalToMatchPoint();


    // Fitness -> here total deltaV (can be updated -> choice left to the user (deltaV, mass, TOF,... ?) )
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

//    for ( int i = 0 ; i < 6 ; i++ )
//    {
//        std::cout << "difference in extended state at match point: " << equalityConstraints[ i ] << "\n\n";
//    }


    // Inequality constraints
    std::vector< double > inequalityConstraints;
    // Magnitude of the normalised thrust vector must be inferior or equal to 1 )
    for ( unsigned int currentThrottle = 0 ; currentThrottle < throttles.size() ; currentThrottle++ )
    {
        inequalityConstraints.push_back( throttles[ currentThrottle ][ 0 ] * throttles[ currentThrottle ][ 0 ]
                + throttles[ currentThrottle ][ 1 ] * throttles[ currentThrottle ][ 1 ]
                + throttles[ currentThrottle ][ 2 ] * throttles[ currentThrottle ][ 2 ] - 1.0 );
    }

//    std::cout << "number inequality constraints: " << inequalityConstraints.size() << "\n\n";
//    for ( int i = 0 ; i < inequalityConstraints.size() ; i++ )
//    {
//        std::cout << "inequality constraints: " << inequalityConstraints[ i ] << "\n\n";
//    }


    // Output of the fitness function.

    // Optimisation objectives
    fitness.push_back( deltaV );
    if ( optimiseTimeOfFlight_ )
    {
        fitness.push_back( timeOfFlight );
    }

    // Return equality constraints
    for ( unsigned int i = 0 ; i < equalityConstraints.size() ; i++ )
    {
        fitness.push_back( equalityConstraints[ i ] );
    }

    // Return inequality constraints.
    for ( unsigned int i = 0 ; i < inequalityConstraints.size() ; i++ )
    {
        fitness.push_back( inequalityConstraints[ i ] );
    }

    return fitness;
}


