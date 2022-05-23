/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_RECOMMENDED_BASE_FUNCTIONS_HODOGRAPHIC_SHAPING_H
#define TUDAT_RECOMMENDED_BASE_FUNCTIONS_HODOGRAPHIC_SHAPING_H

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <cmath>
#include <boost/filesystem.hpp>

#include "tudat/astro/low_thrust/shape_based/baseFunctionsHodographicShaping.h"
#include "tudat/astro/low_thrust/shape_based/createBaseFunctionHodographicShaping.h"

namespace tudat
{
namespace shape_based_methods
{

void getRecommendedBaseFunctions(
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& radialVelocityFunctionComponents,
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& normalVelocityFunctionComponents,
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& axialVelocityFunctionComponents,
        Eigen::VectorXd& freeCoefficientsRadialVelocityFunction,
        Eigen::VectorXd& freeCoefficientsNormalVelocityFunction,
        Eigen::VectorXd& freeCoefficientsAxialVelocityFunction,
        const double timeOfFlight,
        const int numberOfRevolutions );


void getRecommendedRadialVelocityBaseFunctions(
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& radialVelocityFunctionComponents,
        Eigen::VectorXd& freeCoefficientsRadialVelocityFunction,
        const double timeOfFlight );

void getRecommendedNormalBaseFunctions(
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& normalVelocityFunctionComponents,
        Eigen::VectorXd& freeCoefficientsNormalVelocityFunction,
        const double timeOfFlight );

void getRecommendedAxialVelocityBaseFunctions(
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& axialVelocityFunctionComponents,
        Eigen::VectorXd& freeCoefficientsAxialVelocityFunction,
        const double timeOfFlight,
        const int numberOfRevolutions );


std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > getRecommendedRadialVelocityBaseFunctions(
        const double timeOfFlight );

std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > getRecommendedNormalBaseFunctions(
        const double timeOfFlight );

std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > getRecommendedAxialVelocityBaseFunctions(
        const double timeOfFlight,
        const int numberOfRevolutions );

} // namespace shape_based_methods
} // namespace tudat

#endif // TUDAT_RECOMMENDED_BASE_FUNCTIONS_HODOGRAPHIC_SHAPING_H
