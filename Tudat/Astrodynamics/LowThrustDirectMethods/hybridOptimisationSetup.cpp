/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/LowThrustDirectMethods/hybridOptimisationSetup.h"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/hybridMethodLeg.h"
//#include "Tudat/SimulationSetup/tudatSimulationHeader.h"

using namespace tudat;

HybridMethodProblem::HybridMethodProblem(
        const Eigen::Vector6d &stateAtDeparture,
        const Eigen::Vector6d &stateAtArrival,
        const double maximumThrust,
        const std::function< double ( const double ) > specificImpulseFunction,
        const double timeOfFlight,
        simulation_setup::NamedBodyMap bodyMap,
        const std::string bodyToPropagate,
        const std::string centralBody,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const double relativeToleranceConstraints ) :
    stateAtDeparture_( stateAtDeparture ),
    stateAtArrival_( stateAtArrival ),
    maximumThrust_( maximumThrust ),
    specificImpulseFunction_( specificImpulseFunction ),
    timeOfFlight_( timeOfFlight ),
    bodyMap_( bodyMap ),
    bodyToPropagate_( bodyToPropagate ),
    centralBody_( centralBody ),
    integratorSettings_( integratorSettings ),
    relativeToleranceConstraints_( relativeToleranceConstraints )
{
    initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();
}


//! Descriptive name of the problem
std::string HybridMethodProblem::get_name() const {
    return "Hybrid method to compute a low-thrust trajectory";
}

//! Get bounds
std::pair< std::vector< double >, std::vector< double > > HybridMethodProblem::get_bounds() const {

    // Define lower bounds.
    std::vector< double > lowerBounds;
    for ( int i = 0 ; i < 10 ; i++ )
    {
        // Lower bound for the 5 MEE costates.
        lowerBounds.push_back( - 10.0 );
    }

    // Define upper bounds.
    std::vector< double > upperBounds;
    for ( int i = 0 ; i < 10 ; i++ )
    {
        // Upper bound for the 5 MEE costates.
        upperBounds.push_back( 10.0 );
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
    low_thrust_direct_methods::HybridMethodLeg currentLeg = low_thrust_direct_methods::HybridMethodLeg( stateAtDeparture_,
                                                                                                        stateAtArrival_,
                                                                                                        initialCostates,
                                                                                                        finalCostates,
                                                                                                        maximumThrust_,
                                                                                                        specificImpulseFunction_,
                                                                                                        timeOfFlight_,
                                                                                                        bodyMap_,
                                                                                                        bodyToPropagate_,
                                                                                                        centralBody_ );

    // Propagate until time of flight is reached.
    Eigen::Vector6d finalPropagatedState = currentLeg.propagateTrajectory( integratorSettings_ );

    // Convert final propagated state to keplerian elements.
    Eigen::Vector6d finalPropagatedKeplerianElements = orbital_element_conversions::convertCartesianToKeplerianElements(
                finalPropagatedState, bodyMap_[ centralBody_ ]->getGravityFieldModel()->getGravitationalParameter() );

    // Convert targeted final state to keplerian elements.
    Eigen::Vector6d finalTargetedKeplerianElements = orbital_element_conversions::convertCartesianToKeplerianElements(
                stateAtArrival_, bodyMap_[ centralBody_ ]->getGravityFieldModel()->getGravitationalParameter() );


    // Fitness -> here total deltaV (can be updated -> choice left to the user (deltaV, mass, TOF,... ?) )
    double deltaV = currentLeg.getTotalDeltaV( );

    // Equality constraints (must be ... = 0 )
    std::vector< double > equalityConstraints;

//    // Spacecraft position and velocity must be continuous at the match point.
//    for ( int i = 0 ; i < 3 ; i++ )
//    {
//        // Differences in cartesian elements.
//        equalityConstraints.push_back( std::fabs( finalPropagatedState[ i ] - stateAtArrival_[ i ] ) ); /*/
//                                       physical_constants::ASTRONOMICAL_UNIT );*/
//        equalityConstraints.push_back( std::fabs( finalPropagatedState[ i + 3 ] - stateAtArrival_[ i + 3 ] ) ); /*/
//                                       stateAtDeparture_.segment(3 ,3 ).norm() );*/
//    }

    // Differences in keplerian elements.
    equalityConstraints.push_back( std::fabs( finalPropagatedKeplerianElements[ orbital_element_conversions::semiMajorAxisIndex ]
                                   - finalTargetedKeplerianElements[ orbital_element_conversions::semiMajorAxisIndex ] ) /
                                   finalTargetedKeplerianElements[ orbital_element_conversions::semiMajorAxisIndex ]
              * mathematical_constants::PI );
    equalityConstraints.push_back( std::fabs( finalPropagatedKeplerianElements[ orbital_element_conversions::eccentricityIndex ]
                                   - finalTargetedKeplerianElements[ orbital_element_conversions::eccentricityIndex ] )
            / finalTargetedKeplerianElements[ orbital_element_conversions::eccentricityIndex ] * mathematical_constants::PI );
    equalityConstraints.push_back( std::fabs( finalPropagatedKeplerianElements[ orbital_element_conversions::inclinationIndex ]
                                   - finalTargetedKeplerianElements[ orbital_element_conversions::inclinationIndex ] ) );
    equalityConstraints.push_back( std::fabs( finalPropagatedKeplerianElements[ orbital_element_conversions::longitudeOfAscendingNodeIndex ]
                                   - finalTargetedKeplerianElements[ orbital_element_conversions::longitudeOfAscendingNodeIndex ] ) );
    equalityConstraints.push_back( std::fabs( finalPropagatedKeplerianElements[ orbital_element_conversions::argumentOfPeriapsisIndex ]
                                   - finalTargetedKeplerianElements[ orbital_element_conversions::argumentOfPeriapsisIndex ] ) );

    std::cout << "targeted keplerian elements: " << finalTargetedKeplerianElements.transpose() << "\n\n"; // [ orbital_element_conversions::semiMajorAxisIndex ] << "\n\n";
    std::cout << "propagated keplerian elements: " << finalPropagatedKeplerianElements.transpose() << "\n\n"; //[ orbital_element_conversions::semiMajorAxisIndex ] << "\n\n";

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

    double weightTimeOfFlight = 1.0;
    double weightMass = 1.0;
    double weightConstraints = 10.0;
    double optimisationObjective = equalityConstraints[ 0 ] + equalityConstraints[ 1 ] + equalityConstraints[ 2 ]
            + equalityConstraints[ 3 ] + equalityConstraints[ 4 ]; /*weightTimeOfFlight * timeOfFlight + weightMass *
            ( 1.0 - currentLeg.getMassAtTimeOfFlight( ) / initialSpacecraftMass_ )*/
//            + weightConstraints * ( epsilon.norm( ) * epsilon.norm( ) );

//    for ( int i = 0 ; i < 6 ; i++ )
//    {
//        std::cout << "difference in extended state at match point: " << equalityConstraints[ i ] << "\n\n";
//    }


    // Output of the fitness function.

//    // Optimisation objectives
//    fitness.push_back( deltaV );

//    // Return equality constraints
//    for ( unsigned int i = 0 ; i < equalityConstraints.size() ; i++ )
//    {
//        fitness.push_back( equalityConstraints[ i ] );
//    }

    // Output of the fitness function.
    fitness.push_back( optimisationObjective );


    return fitness;
}


