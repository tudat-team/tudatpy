/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/LowThrustTrajectories/hybridOptimisationSetup.h"
#include "tudat/astro/LowThrustTrajectories/hybridMethodModel.h"

namespace tudat
{
namespace low_thrust_trajectories
{

HybridMethodProblem::HybridMethodProblem(
        const Eigen::Vector6d &stateAtDeparture,
        const Eigen::Vector6d &stateAtArrival,
        const double maximumThrust,
        const double specificImpulse,
        const double timeOfFlight,
        simulation_setup::NamedBodyMap bodyMap,
        const std::string bodyToPropagate,
        const std::string centralBody,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::pair< std::vector< double >, double > initialGuessThrustModel,
        const std::pair< double, double > initialAndFinalMEEcostatesBounds,
        const double relativeToleranceConstraints ) :
    stateAtDeparture_( stateAtDeparture ),
    stateAtArrival_( stateAtArrival ),
    maximumThrust_( maximumThrust ),
    specificImpulse_( specificImpulse ),
    timeOfFlight_( timeOfFlight ),
    bodyMap_( bodyMap ),
    bodyToPropagate_( bodyToPropagate ),
    centralBody_( centralBody ),
    integratorSettings_( integratorSettings ),
    initialGuessThrustModel_( initialGuessThrustModel ),
    initialAndFinalMEEcostatesBounds_( initialAndFinalMEEcostatesBounds ),
    relativeToleranceConstraints_( relativeToleranceConstraints )
{
    initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();

    // Retrieve initial guess.
    guessInitialAndFinalCostates_ = initialGuessThrustModel_.first;
    relativeMarginWrtInitialGuess_ = initialGuessThrustModel_.second;

    if ( ( relativeMarginWrtInitialGuess_ != TUDAT_NAN ) && ( ( relativeMarginWrtInitialGuess_ < 0.0) || ( relativeMarginWrtInitialGuess_ > 1.0 ) ) )
    {
        std::string errorMessage = "Error when using initial guess to create hybrid method optimisation problem, the relative margin w.r.t. the "
                                   "initial guess should be constrained between [0,1], the current value is "
                + std::to_string( relativeMarginWrtInitialGuess_ ) + ".";

        throw std::runtime_error( errorMessage );
    }
}


//! Descriptive name of the problem
std::string HybridMethodProblem::get_name() const {
    return "Hybrid method to compute a low-thrust trajectory";
}

//! Get bounds
std::pair< std::vector< double >, std::vector< double > > HybridMethodProblem::get_bounds() const {

    // Define lower bounds.
    std::vector< double > lowerBounds;

    // Define upper bounds.
    std::vector< double > upperBounds;

    if ( guessInitialAndFinalCostates_.size( ) != 0 )
    {
        if ( guessInitialAndFinalCostates_.size() != 10 )
        {
            throw std::runtime_error( "Error when providing an initial guess for hybrid method, size of the vector unconsistent"
                                      "with the expected 5 initial and 5 final MEE costate values." );
        }
        else
        {
            for ( int i = 0 ; i < 10 ; i++ )
            {
                double lowerBoundsFromInitialGuess = ( 1.0 - relativeMarginWrtInitialGuess_ ) * guessInitialAndFinalCostates_[ i ];
                lowerBounds.push_back( lowerBoundsFromInitialGuess );
            }

            for ( int i = 0 ; i < 10 ; i++ )
            {
                double upperBoundsFromInitialGuess = ( 1.0 + relativeMarginWrtInitialGuess_ ) * guessInitialAndFinalCostates_[ i ];
                upperBounds.push_back( upperBoundsFromInitialGuess );
            }
        }
    }
    else
    {
        for ( int i = 0 ; i < 10 ; i ++ )
        {
            lowerBounds.push_back( initialAndFinalMEEcostatesBounds_.first );
            upperBounds.push_back( initialAndFinalMEEcostatesBounds_.second );
        }
    }

    return { lowerBounds, upperBounds };
}

//! Fitness function.
std::vector< double > HybridMethodProblem::fitness( const std::vector< double > &designVariables ) const{

    // Re-initialise mass of the spacecraft.
    bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );

    // Transform vector of design variables into 3D vector of throttles.
    Eigen::VectorXd initialCostates; initialCostates.resize( 5 );
    Eigen::VectorXd finalCostates; finalCostates.resize( 5 );

    // Check consistency of the size of the design variables vector.
    if ( designVariables.size( ) != 10 )
    {
        throw std::runtime_error( "Error, size of the design variables vector unconsistent with initial and final "
                                  "MEE costates sizes." );
    }

    for ( unsigned int i = 0 ; i < 5 ; i++ )
    {
        initialCostates[ i ] = designVariables[ i ];
        finalCostates[ i ] = designVariables[ i + 5 ];
    }

    std::vector< double > fitness;

    // Create hybrid method leg.
    low_thrust_trajectories::HybridMethodModel currentLeg = low_thrust_trajectories::HybridMethodModel(
                stateAtDeparture_, stateAtArrival_, initialCostates, finalCostates, maximumThrust_,
                specificImpulse_, timeOfFlight_, bodyMap_, bodyToPropagate_, centralBody_, integratorSettings_ );

    // Propagate until time of flight is reached.
    Eigen::Vector6d finalPropagatedState = currentLeg.propagateTrajectory( );

    // Convert final propagated state to MEE.
    Eigen::Vector6d finalPropagatedMEEstate = orbital_element_conversions::convertCartesianToModifiedEquinoctialElements(
                finalPropagatedState, bodyMap_[ centralBody_ ]->getGravityFieldModel()->getGravitationalParameter() );

    // Convert targeted final state to MEE.
    Eigen::Vector6d finalTargetedMEEstate = orbital_element_conversions::convertCartesianToModifiedEquinoctialElements(
                stateAtArrival_, bodyMap_[ centralBody_ ]->getGravityFieldModel()->getGravitationalParameter() );


    // Fitness
    double deltaV = currentLeg.getTotalDeltaV( );

    // Equality constraints (must be ... = 0 )
    std::vector< double > equalityConstraints;

    // Differences in MEE at arrival.
    equalityConstraints.push_back( std::fabs( finalPropagatedMEEstate[ orbital_element_conversions::semiParameterIndex ]
                                   - finalTargetedMEEstate[ orbital_element_conversions::semiParameterIndex ] ) );
    equalityConstraints.push_back( std::fabs( finalPropagatedMEEstate[ orbital_element_conversions::fElementIndex ]
                                   - finalTargetedMEEstate[ orbital_element_conversions::fElementIndex ] ) );
    equalityConstraints.push_back( std::fabs( finalPropagatedMEEstate[ orbital_element_conversions::gElementIndex ]
                                   - finalTargetedMEEstate[ orbital_element_conversions::gElementIndex ] ) );
    equalityConstraints.push_back( std::fabs( finalPropagatedMEEstate[ orbital_element_conversions::hElementIndex ]
                                   - finalTargetedMEEstate[ orbital_element_conversions::hElementIndex ] ) );
    equalityConstraints.push_back( std::fabs( finalPropagatedMEEstate[ orbital_element_conversions::kElementIndex ]
                                   - finalTargetedMEEstate[ orbital_element_conversions::kElementIndex ] ) );

    // Compute auxiliary variables.
    Eigen::Vector6d c;
    Eigen::Vector6d r;
    Eigen::Vector6d epsilon;
    for ( int i = 0 ; i < 6 ; i++ )
    {
        c[ i ] = 1.0 / relativeToleranceConstraints_;
        r[ i ] = 1.0 - relativeToleranceConstraints_ * c[ i ];
        epsilon[ i ] = equalityConstraints[ i ] * c[ i ] + r[ i ];
    }

    double weightDeltaV = 1.0;
    double weightConstraints = 10.0;
    double optimisationObjective =  weightDeltaV * deltaV + weightConstraints * ( epsilon.norm( ) * epsilon.norm( ) );


    // Output of the fitness function..
    fitness.push_back( optimisationObjective );


    return fitness;
}

} // namespace low_thrust_trajectories
} // namespace tudat


