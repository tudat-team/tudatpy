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
using namespace tudat::low_thrust_trajectories; //NEED TO CHANGE THIS TO: transfer_trajectories
//using namespace tudat::shape_based_methods;
//using namespace tudat;
using namespace pagmo;

namespace tudat
{
namespace shape_based_methods
{

// Calculates the fitness
std::vector< double > HodographicShapingOptimisationProblem::fitness( const std::vector< double > &x ) const
{

    int numberFreeCoefficientsRadialFunction = radialVelocityFunctionComponents_.size( ) - 3;
    int numberFreeCoefficientsNormalFunction = normalVelocityFunctionComponents_.size( ) - 3;
    int numberFreeCoefficientsAxialFunction = axialVelocityFunctionComponents_.size( ) - 3;

    if ( numberFreeCoefficientsRadialFunction + numberFreeCoefficientsNormalFunction + numberFreeCoefficientsAxialFunction != x.size( ) )
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

    HodographicShaping hodographicShaping = HodographicShaping( initialState_, finalState_, timeOfFlight_, numberOfRevolutions_,
                                                                bodyMap_, bodyToPropagate_, centralBody_, radialVelocityFunctionComponents_,
                                                                normalVelocityFunctionComponents_, axialVelocityFunctionComponents_,
                                                                freeCoefficientsRadialVelocityFunction, freeCoefficientsNormalVelocityFunction,
                                                                freeCoefficientsAxialVelocityFunction );

    std::vector< double > fitnessVector;
    fitnessVector.push_back( hodographicShaping.computeDeltaV( ) );

    return fitnessVector;
}


} // namespace shape_based_methods
} // namespace tudat




