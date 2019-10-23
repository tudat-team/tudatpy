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
#include <math.h>
#include <boost/filesystem.hpp>

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/hodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/baseFunctionsHodographicShaping.h"
#include "Tudat/Astrodynamics/LowThrustTrajectories/ShapeBasedMethods/createBaseFunctionHodographicShaping.h"

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
        const int numberOfRevolutions )
{
    double frequency = 2.0 * mathematical_constants::PI / timeOfFlight;
    double scaleFactor = 1.0 / timeOfFlight;

    // Create base function settings for the components of the radial velocity composite function.
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the radial velocity composite function.
    radialVelocityFunctionComponents.clear( );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, thirdRadialVelocityBaseFunctionSettings ) );

    // Create base function settings for the components of the normal velocity composite function.
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the normal velocity composite function.
    normalVelocityFunctionComponents.clear( );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, thirdNormalVelocityBaseFunctionSettings ) );


    // Create base function settings for the components of the axial velocity composite function.
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions + 0.5 ) * frequency );
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );

    // Set components for the axial velocity function.
    axialVelocityFunctionComponents.clear( );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );

    // Initialize free coefficients vector for radial velocity function.
    freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for normal velocity function.
    freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );

    // Initialize free coefficients vector for axial velocity function.
    freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );

}


void getRecommendedRadialVelocityBaseFunctions(
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& radialVelocityFunctionComponents,
        Eigen::VectorXd& freeCoefficientsRadialVelocityFunction,
        const double timeOfFlight )
{
    double scaleFactor = 1.0 / timeOfFlight;

    // Create base function settings for the components of the radial velocity composite function.
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstRadialVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondRadialVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdRadialVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the radial velocity composite function.
    radialVelocityFunctionComponents.clear( );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::constant, firstRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, secondRadialVelocityBaseFunctionSettings ) );
    radialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, thirdRadialVelocityBaseFunctionSettings ) );

    // Initialize free coefficients vector for radial velocity function.
    freeCoefficientsRadialVelocityFunction = Eigen::VectorXd::Zero( 0 );

}


void getRecommendedNormalAxialBaseFunctions(
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& normalVelocityFunctionComponents,
        Eigen::VectorXd& freeCoefficientsNormalVelocityFunction,
        const double timeOfFlight )
{
    double scaleFactor = 1.0 / timeOfFlight;

    // Create base function settings for the components of the normal velocity composite function.
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstNormalVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::BaseFunctionHodographicShapingSettings >( );
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondNormalVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 1.0, scaleFactor );
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdNormalVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerFunctionHodographicShapingSettings >( 2.0, scaleFactor );

    // Create components of the normal velocity composite function.
    normalVelocityFunctionComponents.clear( );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::constant, firstNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, secondNormalVelocityBaseFunctionSettings ) );
    normalVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPower, thirdNormalVelocityBaseFunctionSettings ) );


    // Initialize free coefficients vector for normal velocity function.
    freeCoefficientsNormalVelocityFunction = Eigen::VectorXd::Zero( 0 );

}

void getRecommendedAxialVelocityBaseFunctions(
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& axialVelocityFunctionComponents,
        Eigen::VectorXd& freeCoefficientsAxialVelocityFunction,
        const double timeOfFlight,
        const int numberOfRevolutions )
{
    double frequency = 2.0 * mathematical_constants::PI / timeOfFlight;
    double scaleFactor = 1.0 / timeOfFlight;

    // Create base function settings for the components of the axial velocity composite function.
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > firstAxialVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::TrigonometricFunctionHodographicShapingSettings >( ( numberOfRevolutions + 0.5 ) * frequency );
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > secondAxialVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >
            ( 3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );
    std::shared_ptr< shape_based_methods::BaseFunctionHodographicShapingSettings > thirdAxialVelocityBaseFunctionSettings =
            std::make_shared< shape_based_methods::PowerTimesTrigonometricFunctionHodographicShapingSettings >(
                3.0, ( numberOfRevolutions + 0.5 ) * frequency, scaleFactor );

    // Set components for the axial velocity function.
    axialVelocityFunctionComponents.clear( );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::cosine, firstAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerCosine, secondAxialVelocityBaseFunctionSettings ) );
    axialVelocityFunctionComponents.push_back(
                createBaseFunctionHodographicShaping( shape_based_methods::scaledPowerSine, thirdAxialVelocityBaseFunctionSettings ) );

    // Initialize free coefficients vector for axial velocity function.
    freeCoefficientsAxialVelocityFunction = Eigen::VectorXd::Zero( 0 );

}
} // namespace shape_based_methods
} // namespace tudat

#endif // TUDAT_RECOMMENDED_BASE_FUNCTIONS_HODOGRAPHIC_SHAPING_H
