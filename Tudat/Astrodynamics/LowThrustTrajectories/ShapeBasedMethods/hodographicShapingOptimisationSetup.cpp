/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "hodographicShapingOptimisationSetup.h"

using namespace tudat::ephemerides;
using namespace tudat::basic_astrodynamics;
using namespace tudat::orbital_element_conversions;
using namespace tudat::basic_mathematics;
using namespace tudat::input_output;
using namespace tudat::low_thrust_trajectories;
using namespace pagmo;

namespace tudat
{
namespace shape_based_methods
{

// Calculates the fitness
std::vector< double > FixedTimeHodographicShapingOptimisationProblem::fitness( const std::vector< double > &x ) const
{

    int numberFreeCoefficientsRadialFunction = radialVelocityFunctionComponents_.size( ) - 3;
    int numberFreeCoefficientsNormalFunction = normalVelocityFunctionComponents_.size( ) - 3;
    int numberFreeCoefficientsAxialFunction = axialVelocityFunctionComponents_.size( ) - 3;

    if( numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction !=
            static_cast< int >( x.size( ) ) )
    {
        throw std::runtime_error( "Error, size of design variables vector unconsistent with number of base function components"
                                  "when making a hodographic shaping optimisation problem." );
    }

    Eigen::VectorXd freeCoefficientsRadialVelocityFunction( numberFreeCoefficientsRadialFunction );
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction( numberFreeCoefficientsNormalFunction );
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction( numberFreeCoefficientsAxialFunction );

    for ( int i = 0 ; i < numberFreeCoefficientsRadialFunction ; i++ )
    {
        freeCoefficientsRadialVelocityFunction[ i ] = x[ i ];
    }
    for( int i = 0 ; i < numberFreeCoefficientsNormalFunction ; i++ )
    {
        freeCoefficientsNormalVelocityFunction[ i ] = x[ i + numberFreeCoefficientsRadialFunction ];
    }
    for ( int i = 0 ; i < numberFreeCoefficientsAxialFunction ; i++ )
    {
        freeCoefficientsAxialVelocityFunction[ i ] = x[ i + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction ];
    }

    HodographicShaping hodographicShaping = HodographicShaping(
                initialState_, finalState_, timeOfFlight_, centralBodyGravitationalParameter_,
                 numberOfRevolutions_, radialVelocityFunctionComponents_,
                normalVelocityFunctionComponents_, axialVelocityFunctionComponents_,
                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction );

    std::vector< double > fitnessVector;
    fitnessVector.push_back( hodographicShaping.computeDeltaV( ) );

    return fitnessVector;
}

// Calculates the fitness
std::vector< double > HodographicShapingOptimisationProblem::fitness( const std::vector< double > &x ) const
{
    double departureTime = x.at( 0 );
    double timeOfFlight = x.at( 1 );
    double arrivalTime = departureTime + timeOfFlight;

    std::vector< BaseFunctionVector > basisFunctions = basisFunctionsFunction_( timeOfFlight );

    BaseFunctionVector radialVelocityFunctionComponents = basisFunctions.at( 0 );
    BaseFunctionVector normalVelocityFunctionComponents = basisFunctions.at( 1 );
    BaseFunctionVector axialVelocityFunctionComponents = basisFunctions.at( 2 );

    int numberFreeCoefficientsRadialFunction = radialVelocityFunctionComponents.size( ) - 3;
    int numberFreeCoefficientsNormalFunction = normalVelocityFunctionComponents.size( ) - 3;
    int numberFreeCoefficientsAxialFunction = axialVelocityFunctionComponents.size( ) - 3;

    if( numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction + 5 !=
            static_cast< int >( x.size( ) ) )
    {
        throw std::runtime_error( "Error, size of design variables vector unconsistent with number of base function components"
                                  "when making a hodographic shaping optimisation problem." );
    }

    Eigen::VectorXd freeCoefficientsRadialVelocityFunction( numberFreeCoefficientsRadialFunction );
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction( numberFreeCoefficientsNormalFunction );
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction( numberFreeCoefficientsAxialFunction );

    for ( int i = 0 ; i < numberFreeCoefficientsRadialFunction ; i++ )
    {
        freeCoefficientsRadialVelocityFunction[ i ] = x[ i + 2 ];
    }
    for( int i = 0 ; i < numberFreeCoefficientsNormalFunction ; i++ )
    {
        freeCoefficientsNormalVelocityFunction[ i ] = x[ i + 2 + numberFreeCoefficientsRadialFunction ];
    }
    for ( int i = 0 ; i < numberFreeCoefficientsAxialFunction ; i++ )
    {
        freeCoefficientsAxialVelocityFunction[ i ] = x[ i + 2 + numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction ];
    }

    Eigen::Vector6d departureState = initialStateFunction_( departureTime );
    departureState( 3 ) += x.at( 7 );
    departureState( 4 ) += x.at( 8 );
    departureState( 5 ) += x.at( 9 );

    HodographicShaping hodographicShaping = HodographicShaping(
                departureState, finalStateFunction_( arrivalTime ),
                timeOfFlight, centralBodyGravitationalParameter_,
                numberOfRevolutions_, radialVelocityFunctionComponents,
                normalVelocityFunctionComponents, axialVelocityFunctionComponents,
                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction,
                freeCoefficientsAxialVelocityFunction, initialMass_ );

    std::vector< double > fitnessVector;
    fitnessVector.push_back( hodographicShaping.computeDeltaV( ) );

    if( minimizeMaximumThrust_ )
    {
        Eigen::VectorXd epochsVector = Eigen::VectorXd::LinSpaced( 100, 0, timeOfFlight );

        double maximumAcceleration = 0.0;
        for( int i = 0; i < epochsVector.rows( ); i++ )
        {
            double currentAcceleration =
                    hodographicShaping.computeCurrentThrustAccelerationMagnitude( epochsVector( i ) );
            if( maximumAcceleration < currentAcceleration )
            {
                maximumAcceleration = currentAcceleration;
            }
        }
        fitnessVector.push_back( maximumAcceleration );

    }

    return fitnessVector;
}

} // namespace shape_based_methods
} // namespace tudat




