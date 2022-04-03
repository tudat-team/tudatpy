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


#include <cmath>
#include <iostream>
#include "tudat/basics/timeType.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/low_thrust/shape_based/hodographicShaping.h"


namespace tudat
{
namespace shape_based_methods
{


//! Constructor which sets radial, normal and axial velocity functions and boundary conditions.
HodographicShaping::HodographicShaping(
        const Eigen::Vector6d& initialState,
        const Eigen::Vector6d& finalState,
        const double timeOfFlight,
        const double centralBodyGravitationalParameter,
        const int numberOfRevolutions,
        const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& radialVelocityFunctionComponents,
        const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& normalVelocityFunctionComponents,
        const std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& axialVelocityFunctionComponents,
        const Eigen::VectorXd& freeCoefficientsRadialVelocityFunction,
        const Eigen::VectorXd& freeCoefficientsNormalVelocityFunction,
        const Eigen::VectorXd& freeCoefficientsAxialVelocityFunction ) :
    ShapeBasedMethod( initialState, finalState, timeOfFlight ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    numberOfRevolutions_( numberOfRevolutions ),
    freeCoefficientsRadialVelocityFunction_( freeCoefficientsRadialVelocityFunction ),
    freeCoefficientsNormalVelocityFunction_( freeCoefficientsNormalVelocityFunction ),
    freeCoefficientsAxialVelocityFunction_( freeCoefficientsAxialVelocityFunction ),
    radialPositionHasBeenNegative_( true )
{
    // Define composite function in radial direction.
    Eigen::VectorXd radialVelocityCoefficients;
    radialVelocityCoefficients.resize( 3 + freeCoefficientsRadialVelocityFunction_.size() );
    radialVelocityCoefficients.segment( 0, 3 ) = Eigen::Vector3d::Zero();
    radialVelocityCoefficients.segment( 3, freeCoefficientsRadialVelocityFunction_.size() ) = freeCoefficientsRadialVelocityFunction_;

    radialVelocityFunction_ = std::make_shared< shape_based_methods::CompositeFunctionHodographicShaping >(
                radialVelocityFunctionComponents, radialVelocityCoefficients );

    // Define composite function in normal direction.
    Eigen::VectorXd normalVelocityCoefficients;
    normalVelocityCoefficients.resize( 3 + freeCoefficientsNormalVelocityFunction_.size() );
    normalVelocityCoefficients.segment( 0, 3 ) = Eigen::Vector3d::Zero();
    normalVelocityCoefficients.segment( 3, freeCoefficientsNormalVelocityFunction_.size() ) = freeCoefficientsNormalVelocityFunction_;

    normalVelocityFunction_ = std::make_shared< shape_based_methods::CompositeFunctionHodographicShaping >(
                normalVelocityFunctionComponents, normalVelocityCoefficients );

    // Define composite function in axial direction.
    Eigen::VectorXd axialVelocityCoefficients;
    axialVelocityCoefficients.resize( 3 + freeCoefficientsAxialVelocityFunction_.size() );
    axialVelocityCoefficients.segment( 0, 3 ) = Eigen::Vector3d::Zero();
    axialVelocityCoefficients.segment( 3, freeCoefficientsAxialVelocityFunction_.size() ) = freeCoefficientsAxialVelocityFunction_;

    axialVelocityFunction_ = std::make_shared< shape_based_methods::CompositeFunctionHodographicShaping >(
                axialVelocityFunctionComponents, axialVelocityCoefficients );


    // Compute initial cylindrical state.
    Eigen::Vector6d initialCylindricalState = coordinate_conversions::convertCartesianToCylindricalState( stateAtDeparture_ );

    // Compute final cylindrical state.
    Eigen::Vector6d finalCylindricalState = coordinate_conversions::convertCartesianToCylindricalState( stateAtArrival_ );

    // Compute number of free coefficients for each velocity function.
    numberOfFreeCoefficientsRadialVelocityFunction_ = radialVelocityFunction_->getNumberOfCompositeFunctionComponents() - 3;
    numberOfFreeCoefficientsNormalVelocityFunction_ = normalVelocityFunction_->getNumberOfCompositeFunctionComponents() - 3;
    numberOfFreeCoefficientsAxialVelocityFunction = axialVelocityFunction_->getNumberOfCompositeFunctionComponents() - 3;
    std::cout<<"Radial: "<<numberOfFreeCoefficientsRadialVelocityFunction_<<" "<<radialVelocityFunctionComponents.size( ) - 3<<std::endl;
    if ( numberOfFreeCoefficientsRadialVelocityFunction_ < 0 || numberOfFreeCoefficientsNormalVelocityFunction_ < 0 ||
         numberOfFreeCoefficientsAxialVelocityFunction < 0 )
    {
        std::cerr << "The number of base functions in one of the velocity functions is smaller than "
                  << "the number of corresponding boundary conditions. The boundary conditions cannot be set!\n";
    }
    else
    {
        // Define numerical quadrature settings, required to compute the current polar angle and final deltaV.
        quadratureSettings_ = std::make_shared< numerical_quadrature::GaussianQuadratureSettings< double > >( 0.0, 64 );

        // Set boundary conditions in the radial direction.
        radialBoundaryConditions_.push_back( initialCylindricalState[ 0 ] );
        radialBoundaryConditions_.push_back( finalCylindricalState[ 0 ] );
        radialBoundaryConditions_.push_back( initialCylindricalState[ 3 ] );
        radialBoundaryConditions_.push_back( finalCylindricalState[ 3 ] );

        // Set boundary conditions in the normal direction.
        normalBoundaryConditions_.push_back( initialCylindricalState[ 4 ] );
        normalBoundaryConditions_.push_back( finalCylindricalState[ 4 ] );

        // Final value polar angle
        normalBoundaryConditions_.push_back(
                    numberOfRevolutions_ * 2.0 * mathematical_constants::PI
                    + basic_mathematics::computeModulo( ( finalCylindricalState[ 1 ] - initialCylindricalState[ 1 ] ),
                2.0 * mathematical_constants::PI ) );

        // Set boundary conditions in the axial direction.
        axialBoundaryConditions_.push_back( initialCylindricalState[ 2 ] );
        axialBoundaryConditions_.push_back( finalCylindricalState[ 2 ] );
        axialBoundaryConditions_.push_back( initialCylindricalState[ 5 ] );
        axialBoundaryConditions_.push_back( finalCylindricalState[ 5 ] );

        // Compute inverse of matrices containing boundary values.
        inverseMatrixRadialBoundaryValues_ = computeInverseMatrixRadialOrAxialBoundaries( radialVelocityFunction_ );
        inverseMatrixNormalBoundaryValues_ = computeInverseMatrixNormalBoundaries( normalVelocityFunction_ );
        inverseAxialMatrixBoundaryValues_ = computeInverseMatrixRadialOrAxialBoundaries( axialVelocityFunction_ );

        // Satisfy boundary conditions.
        satisfyRadialBoundaryConditions( freeCoefficientsRadialVelocityFunction_ );
        satisfyNormalBoundaryConditions( freeCoefficientsNormalVelocityFunction_ );
        satisfyAxialBoundaryConditions( freeCoefficientsAxialVelocityFunction_ );
    }


}


//! Compute inverse of matrix filled with boundary conditions in radial or axial direction.
Eigen::Matrix3d HodographicShaping::computeInverseMatrixRadialOrAxialBoundaries(
        std::shared_ptr< CompositeFunctionHodographicShaping > velocityFunction )
{
    Eigen::Matrix3d matrixBoundaryValues, inverseMatrixBoundaryValues;
    matrixBoundaryValues <<
                         velocityFunction->getComponentFunctionIntegralCurrentValue(0, timeOfFlight_)
                            - velocityFunction->getComponentFunctionIntegralCurrentValue(0, 0.0),
            velocityFunction->getComponentFunctionIntegralCurrentValue(1, timeOfFlight_)
            - velocityFunction->getComponentFunctionIntegralCurrentValue(1, 0.0),
            velocityFunction->getComponentFunctionIntegralCurrentValue(2, timeOfFlight_)
            - velocityFunction->getComponentFunctionIntegralCurrentValue(2, 0.0),
            velocityFunction->getComponentFunctionCurrentValue( 0, 0.0 ),
            velocityFunction->getComponentFunctionCurrentValue( 1, 0.0 ),
            velocityFunction->getComponentFunctionCurrentValue( 2, 0.0 ),
            velocityFunction->getComponentFunctionCurrentValue( 0, timeOfFlight_ ),
            velocityFunction->getComponentFunctionCurrentValue( 1, timeOfFlight_ ),
            velocityFunction->getComponentFunctionCurrentValue( 2, timeOfFlight_ );

    // Compute inverse of boundary-value matrix.
    inverseMatrixBoundaryValues = matrixBoundaryValues.inverse();

    // Return inverse of boundary-value matrix.
    return inverseMatrixBoundaryValues;
}

//! Compute inverse of matrix filled with boundary conditions in normal direction.
Eigen::Matrix2d HodographicShaping::computeInverseMatrixNormalBoundaries(
        std::shared_ptr< CompositeFunctionHodographicShaping > velocityFunction )
{
    Eigen::Matrix2d matrixBoundaryValues, inverseMatrixBoundaryValues;

    matrixBoundaryValues <<
                            velocityFunction->getComponentFunctionCurrentValue( 0, 0.0 ),
            velocityFunction->getComponentFunctionCurrentValue( 1, 0.0 ),
            velocityFunction->getComponentFunctionCurrentValue( 0, timeOfFlight_ ),
            velocityFunction->getComponentFunctionCurrentValue( 1, timeOfFlight_ );

    // Compute inverse of boundary-value matrix.
    inverseMatrixBoundaryValues = matrixBoundaryValues.inverse();

    // Return inverse of boundary-value matrix.
    return inverseMatrixBoundaryValues;
}


void HodographicShaping::satisfyRadialBoundaryConditions( Eigen::VectorXd freeCoefficients ){

    // Vector containing boundary conditions on radial distance and initial and final radial velocity.
    Eigen::Vector3d vectorBoundaryConditionsRadial;
    vectorBoundaryConditionsRadial[ 0 ] = radialBoundaryConditions_[ 1 ] - radialBoundaryConditions_[ 0 ];
    vectorBoundaryConditionsRadial[ 1 ] = radialBoundaryConditions_[ 2 ];
    vectorBoundaryConditionsRadial[ 2 ] = radialBoundaryConditions_[ 3 ];

    // Subtract boundary values of free components of velocity function from corresponding boundary conditions.
    for ( int i = 0 ; i < numberOfFreeCoefficientsRadialVelocityFunction_ ; i++ )
    {
        vectorBoundaryConditionsRadial[ 0 ] -=
                freeCoefficients( i ) * ( radialVelocityFunction_->getComponentFunctionIntegralCurrentValue(i + 3,
                                                                                                            timeOfFlight_)
                                          -
                        radialVelocityFunction_->getComponentFunctionIntegralCurrentValue(i + 3, 0.0) );
        vectorBoundaryConditionsRadial[ 1 ] -= freeCoefficients( i ) *
                radialVelocityFunction_->getComponentFunctionCurrentValue( i + 3, 0.0 );
        vectorBoundaryConditionsRadial[ 2 ] -= freeCoefficients( i ) *
                radialVelocityFunction_->getComponentFunctionCurrentValue( i + 3, timeOfFlight_ );
    }

    // Compute fixed coefficients.
    Eigen::Vector3d fixedCoefficientsRadial = inverseMatrixRadialBoundaryValues_ * vectorBoundaryConditionsRadial;

    // Create vector containing all radial velocity function coefficients.
    Eigen::VectorXd radialVelocityFunctionCoefficients =
            Eigen::VectorXd::Zero( fixedCoefficientsRadial.size() + numberOfFreeCoefficientsRadialVelocityFunction_ );

    // Check whether the radial velocity function has free coefficients.
    if ( numberOfFreeCoefficientsRadialVelocityFunction_ == 0 )
    {
        radialVelocityFunctionCoefficients << fixedCoefficientsRadial;
    }
    else
    {
        radialVelocityFunctionCoefficients << fixedCoefficientsRadial,
                freeCoefficients.segment( 0, numberOfFreeCoefficientsRadialVelocityFunction_);
    }

    // Set the coefficients of the radial velocity function.
    radialVelocityFunction_->resetCompositeFunctionCoefficients( radialVelocityFunctionCoefficients );

}

void HodographicShaping::satisfyAxialBoundaryConditions( Eigen::VectorXd freeCoefficients ){

    // Vector containing boundary conditions on axial distance and initial and final axial velocity.
    Eigen::Vector3d vectorBoundaryConditionsAxial;
    vectorBoundaryConditionsAxial[ 0 ] = axialBoundaryConditions_[ 1 ] - axialBoundaryConditions_[ 0 ];
    vectorBoundaryConditionsAxial[ 1 ] = axialBoundaryConditions_[ 2 ];
    vectorBoundaryConditionsAxial[ 2 ] = axialBoundaryConditions_[ 3 ];

    // Subtract boundary values of free components of velocity function from corresponding boundary conditions.
    for ( int i = 0 ; i < numberOfFreeCoefficientsAxialVelocityFunction ; i++ )
    {
        vectorBoundaryConditionsAxial[ 0 ] -= freeCoefficients( i )
                * ( axialVelocityFunction_->getComponentFunctionIntegralCurrentValue(i + 3, timeOfFlight_)
                    - axialVelocityFunction_->getComponentFunctionIntegralCurrentValue(i + 3, 0.0) );
        vectorBoundaryConditionsAxial[ 1 ] -= freeCoefficients( i )
                * axialVelocityFunction_->getComponentFunctionCurrentValue( i + 3, 0.0 );
        vectorBoundaryConditionsAxial[ 2 ] -= freeCoefficients( i )
                * axialVelocityFunction_->getComponentFunctionCurrentValue( i + 3, timeOfFlight_ );
    }

    // Compute fixed coefficients.
    Eigen::Vector3d fixedCoefficientsAxial = inverseAxialMatrixBoundaryValues_ * vectorBoundaryConditionsAxial;

    // Create vector containing all axial velocity function coefficients.
    Eigen::VectorXd axialVelocityFunctionCoefficients =
            Eigen::VectorXd::Zero( fixedCoefficientsAxial.size() + numberOfFreeCoefficientsAxialVelocityFunction );

    // Check whether the axial velocity function has free coefficients.
    if ( numberOfFreeCoefficientsAxialVelocityFunction == 0 )
    {
        axialVelocityFunctionCoefficients << fixedCoefficientsAxial;
    }
    else
    {
        axialVelocityFunctionCoefficients << fixedCoefficientsAxial,
                freeCoefficients.segment( 0, numberOfFreeCoefficientsAxialVelocityFunction);
    }
    // Set the coefficients of the axial velocity function.
    axialVelocityFunction_->resetCompositeFunctionCoefficients( axialVelocityFunctionCoefficients );

}

void HodographicShaping::satisfyNormalBoundaryConditions( Eigen::VectorXd freeCoefficients ){

    // Compute coefficient of the third component of the composite function, from the required value of the final polar angle.
    double C3 = computeThirdFixedCoefficientAxialVelocityFromFinalPolarAngle( freeCoefficients );

    // Vector containing boundary conditions on initial and final normal velocity.
    Eigen::Vector2d vectorBoundaryConditionsNormal;
    vectorBoundaryConditionsNormal( 0 ) = normalBoundaryConditions_[ 0 ];
    vectorBoundaryConditionsNormal( 1 ) = normalBoundaryConditions_[ 1 ];

    // Subtract boundary values of velocity function component used to solve for final polar angle from corresponding boundary conditions.
    vectorBoundaryConditionsNormal[ 0 ] -= C3 * normalVelocityFunction_->getComponentFunctionCurrentValue( 2, 0.0 );
    vectorBoundaryConditionsNormal[ 1 ] -= C3 * normalVelocityFunction_->getComponentFunctionCurrentValue( 2, timeOfFlight_ );

    // Subtract boundary values of free velocity function components from corresponding boundary conditions.
    for ( int i = 0 ; i < numberOfFreeCoefficientsNormalVelocityFunction_ ; i++ )
    {
        vectorBoundaryConditionsNormal[ 0 ] -= freeCoefficients[ i ]
                * normalVelocityFunction_->getComponentFunctionCurrentValue( i + 3, 0.0 );
        vectorBoundaryConditionsNormal[ 1 ] -= freeCoefficients[ i ]
                * normalVelocityFunction_->getComponentFunctionCurrentValue( i + 3, timeOfFlight_ );
    }

    // Compute fixed coefficients by satisfying the boundary conditions.
    Eigen::Vector2d fixedCoefficientsNormal = inverseMatrixNormalBoundaryValues_ * vectorBoundaryConditionsNormal;

    // Create vector containing all normal velocity function coefficients.
    Eigen::VectorXd normalVelocityFunctionCoefficients =
            Eigen::VectorXd::Zero( fixedCoefficientsNormal.size() + 1 + numberOfFreeCoefficientsNormalVelocityFunction_ );

    // Check whether the normal velocity function has free coefficients.
    if ( numberOfFreeCoefficientsNormalVelocityFunction_ == 0 )
    {
        normalVelocityFunctionCoefficients << fixedCoefficientsNormal, C3;
    }
    else
    {
        normalVelocityFunctionCoefficients << fixedCoefficientsNormal, C3,
                freeCoefficients.segment( 0, numberOfFreeCoefficientsNormalVelocityFunction_);
    }

    // Set the coefficients of the normal velocity function.
    normalVelocityFunction_->resetCompositeFunctionCoefficients( normalVelocityFunctionCoefficients );

}

// Define angular velocity due to the third component of the composite function only.
double HodographicShaping::computeDerivativePolarAngleDueToThirdComponent(
        const double timeSinceDeparture, const Eigen::Vector2d& matrixK )
{

    double angularVelocityDueToThirdComponent = (
                matrixK( 0 ) * normalVelocityFunction_->getComponentFunctionCurrentValue( 0, timeSinceDeparture )
                + matrixK( 1 ) * normalVelocityFunction_->getComponentFunctionCurrentValue( 1, timeSinceDeparture )
                + normalVelocityFunction_->getComponentFunctionCurrentValue( 2, timeSinceDeparture ) ) /
            computeCurrentRadialDistance( timeSinceDeparture );

    return angularVelocityDueToThirdComponent;

}


// Define the angular velocity due to all the other components of the composite function, once combined.
double HodographicShaping::computeDerivativePolarAngleDueToOtherComponents(
        const double timeSinceDeparture, const Eigen::Vector2d& matrixL, const Eigen::VectorXd& freeCoefficients )
{

    double angularVelocityDueToFreeCoefficients = 0.0;
    for( int j = 0 ; j < numberOfFreeCoefficientsNormalVelocityFunction_ ; j++ )
    {
        angularVelocityDueToFreeCoefficients += freeCoefficients( j )
                * normalVelocityFunction_->getComponentFunctionCurrentValue( j + 3 , timeSinceDeparture );
    }

    double angularVelocityDueToOtherComponents = ( matrixL( 0 ) * normalVelocityFunction_->getComponentFunctionCurrentValue( 0 , timeSinceDeparture )
                                                   + matrixL( 1 ) * normalVelocityFunction_->getComponentFunctionCurrentValue( 1 , timeSinceDeparture )
                                                   + angularVelocityDueToFreeCoefficients )
            / computeCurrentRadialDistance( timeSinceDeparture );

    return angularVelocityDueToOtherComponents;

}

double HodographicShaping::computeThirdFixedCoefficientAxialVelocityFromFinalPolarAngle( Eigen::VectorXd freeCoefficients ){

    // Compute the third fixed coefficient of the normal velocity composite function, so that the condition on the final
    // polar angle is fulfilled.
    // The calculation is based on Equation (16) in Gondelach D., and Noomen R.
    // "Hodographic-shaping method for low-thrust interplanetary trajectory design." Journal of Spacecraft and Rockets 52.3 (2015): 728-738.

    // Vector containing boundary conditions on initial and final normal velocity.
    Eigen::Vector2d vectorBoundaryConditionsNormal;
    vectorBoundaryConditionsNormal[ 0 ] = normalBoundaryConditions_[ 0 ];
    vectorBoundaryConditionsNormal[ 1 ] = normalBoundaryConditions_[ 1 ];

    // Subtract boundary values of free velocity function components from corresponding boundary conditions.
    for ( int i = 0 ; i < numberOfFreeCoefficientsNormalVelocityFunction_; i++ )
    {
        vectorBoundaryConditionsNormal[ 0 ] -= freeCoefficients[ i ]
                * normalVelocityFunction_->getComponentFunctionCurrentValue( i + 3, 0.0 );
        vectorBoundaryConditionsNormal[ 1 ] -= freeCoefficients[ i ]
                * normalVelocityFunction_->getComponentFunctionCurrentValue( i + 3, timeOfFlight_ );
    }

    // Define matrix L, as proposed in ... (ADD PROPER REFERENCE AND EQUATION NUMBER)
    Eigen::Vector2d matrixL;
    matrixL = inverseMatrixNormalBoundaryValues_ * vectorBoundaryConditionsNormal;

    Eigen::Vector2d initialAndFinalValuesThirdComponentFunction(
                - normalVelocityFunction_->getComponentFunctionCurrentValue( 2, 0.0 ),
                - normalVelocityFunction_->getComponentFunctionCurrentValue( 2, timeOfFlight_ ) );

    // Define matrix K, as proposed in ... (ADD PROPER REFERENCE AND EQUATION NUMBER)
    Eigen::Vector2d matrixK;
    matrixK = inverseMatrixNormalBoundaryValues_ * initialAndFinalValuesThirdComponentFunction;


    // Define angular velocity due to the third component of the composite function only.
    std::function< double( const double ) > derivativePolarAngleDueToThirdComponent =
            std::bind( &HodographicShaping::computeDerivativePolarAngleDueToThirdComponent, this, std::placeholders::_1, matrixK );

    // Define the angular velocity due to all the other components of the composite function, once combined.
    std::function< double( const double ) > derivativePolarAngleDueToOtherComponents  =
            std::bind( &HodographicShaping::computeDerivativePolarAngleDueToOtherComponents, this, std::placeholders::_1, matrixL, freeCoefficients );

    // Define numerical quadratures.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > integratorPolarAngleDueToThirdComponentTest =
            numerical_quadrature::createQuadrature( derivativePolarAngleDueToThirdComponent, quadratureSettings_, timeOfFlight_ );

    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > integratorPolarAngleDueToOtherComponentsTest =
            numerical_quadrature::createQuadrature( derivativePolarAngleDueToOtherComponents, quadratureSettings_, timeOfFlight_ );

    return ( normalBoundaryConditions_[ 2 ] - integratorPolarAngleDueToOtherComponentsTest->getQuadrature() )
            / integratorPolarAngleDueToThirdComponentTest->getQuadrature();

}


//! Compute velocity components.
Eigen::Vector3d HodographicShaping::computeVelocityVectorInCylindricalCoordinates( double timeSinceDeparture )
{
    Eigen::Vector3d velocityComponents;

    // Compute radial velocity.
    velocityComponents[ 0 ] = radialVelocityFunction_->evaluateCompositeFunctionCurrentValue(timeSinceDeparture);

    // Compute normal velocity.
    velocityComponents[ 1 ] = normalVelocityFunction_->evaluateCompositeFunctionCurrentValue(timeSinceDeparture);

    // Compute axial velocity.
    velocityComponents[ 2 ] = axialVelocityFunction_->evaluateCompositeFunctionCurrentValue(timeSinceDeparture);

    return velocityComponents;
}


//! Compute angular velocity.
double HodographicShaping::evaluateDerivativePolarAngleWrtTime( const double timeSinceDeparture )
{
    // Return derivative of the polar angle w.r.t. time, i.e. angular velocity.
    return normalVelocityFunction_->evaluateCompositeFunctionCurrentValue(timeSinceDeparture)
            / ( radialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentValue(timeSinceDeparture)
                - radialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentValue(0.0) + radialBoundaryConditions_[0] );
}


//! Compute radial distance.
double HodographicShaping::computeCurrentRadialDistance( const double timeSinceDeparture )
{
    // Compute radial distance from the central body.
    return radialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentValue(timeSinceDeparture)
           - radialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentValue(0.0) + radialBoundaryConditions_[0];
}

//! Compute current polar angle.
double HodographicShaping::computeCurrentPolarAngle( double timeSinceDeparture )
{
    // Define the derivative of the polar angle, ie angular velocity function, as a function of time.
    std::function< double( const double ) > derivativeFunctionPolarAngle =
            std::bind( &HodographicShaping::evaluateDerivativePolarAngleWrtTime, this, std::placeholders::_1 );

    // Define numerical quadrature.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeFunctionPolarAngle, quadratureSettings_, timeSinceDeparture );

    double currentPolarAngle = quadrature->getQuadrature( ) + coordinate_conversions::convertCartesianToCylindricalState( stateAtDeparture_ )[ 1 ];

    return currentPolarAngle;

}

double HodographicShaping::computeCurrentAxialDistance( const double timeSinceDeparture )
{
    // Compute axial distance from the central body.
    return axialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentValue(timeSinceDeparture)
           - axialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentValue(0.0) + axialBoundaryConditions_[0];
}

//! Compute current cylindrical state.
Eigen::Vector6d HodographicShaping::computeStateVectorInCylindricalCoordinates( const double timeSinceDeparture )
{
    Eigen::Vector3d velocityCylindricalCoordinates = computeVelocityVectorInCylindricalCoordinates( timeSinceDeparture );

    Eigen::Vector6d cylindricalState = ( Eigen::Vector6d() <<
                                         computeCurrentRadialDistance( timeSinceDeparture ),
                                         computeCurrentPolarAngle( timeSinceDeparture ),
                                         computeCurrentAxialDistance( timeSinceDeparture ),
                                         velocityCylindricalCoordinates[ 0 ],
            velocityCylindricalCoordinates[ 1 ],
            velocityCylindricalCoordinates[ 2 ] ).finished();
    return cylindricalState;
}


//! Compute current cartesian state.
Eigen::Vector6d HodographicShaping::computeCurrentStateVector( const double timeSinceDeparture ){

    return coordinate_conversions::convertCylindricalToCartesianState( computeStateVectorInCylindricalCoordinates( timeSinceDeparture ) );

}


//! Compute thrust acceleration components.
Eigen::Vector3d HodographicShaping::computeThrustAccelerationInCylindricalCoordinates( double timeSinceDeparture )
{

    Eigen::Vector3d thrustAccelerationComponents;

    // Compute radial distance from the central body.
    double radialDistance = radialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentValue(timeSinceDeparture)
                            - radialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentValue(0.0) + radialBoundaryConditions_[0];
    // Compute axial distance from the central body.
    double axialDistance = axialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentValue(timeSinceDeparture)
                           - axialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentValue(0.0) + axialBoundaryConditions_[0];

    // Compute distance from the central body.
    double distanceFromCentralBody = sqrt( pow( radialDistance, 2.0 ) + pow( axialDistance, 2.0 ) );

    // Computation of normal velocity.
    double normalVelocity = normalVelocityFunction_->evaluateCompositeFunctionCurrentValue(timeSinceDeparture);

    // Compute angular velocity.
    double angularVelocity = normalVelocity / radialDistance;

    // Compute radial thrust acceleration.
    thrustAccelerationComponents[ 0 ] =
            radialVelocityFunction_->evaluateCompositeFunctionDerivativeCurrentValue(timeSinceDeparture)
            - angularVelocity * normalVelocity + centralBodyGravitationalParameter_ / std::pow( distanceFromCentralBody, 3.0 ) * radialDistance;

    // Compute normal thrust acceleration.
    thrustAccelerationComponents[ 1 ] =
            normalVelocityFunction_->evaluateCompositeFunctionDerivativeCurrentValue(timeSinceDeparture)
            + angularVelocity * radialVelocityFunction_->evaluateCompositeFunctionCurrentValue(timeSinceDeparture);

    // Compute axial thrust acceleration.
    thrustAccelerationComponents[ 2 ] =
            axialVelocityFunction_->evaluateCompositeFunctionDerivativeCurrentValue(timeSinceDeparture)
            + centralBodyGravitationalParameter_ / std::pow( distanceFromCentralBody, 3.0 ) * axialDistance;

    // Return total thrust acceleration.
    return thrustAccelerationComponents;
}


double HodographicShaping::computeCurrentThrustAccelerationMagnitude(
        const double timeSinceDeparture )
{
    // Return total thrust acceleration.
    return computeCurrentThrustAcceleration( timeSinceDeparture ).norm();
}

Eigen::Vector3d HodographicShaping::computeCurrentThrustAcceleration( double timeSinceDeparture )
{
    // Prevent out-of-range errors in quadrature
    if( timeSinceDeparture < 0.0 )
    {
        timeSinceDeparture = 0.0;
    }
    else if( timeSinceDeparture > timeOfFlight_ )
    {
        timeSinceDeparture = timeOfFlight_;
    }

    if( thrustAccelerationVectorCache_.count( timeSinceDeparture ) == 0 )
    {
        Eigen::Vector3d cylindricalAcceleration = computeThrustAccelerationInCylindricalCoordinates( timeSinceDeparture );
        Eigen::Vector3d cylindricalState = computeStateVectorInCylindricalCoordinates( timeSinceDeparture ).segment(0,3);
        Eigen::Vector3d cartesianState = computeCurrentStateVector( timeSinceDeparture ).segment(0,3);
        if( cylindricalState( 0 ) < 0.0 )
        {
            radialPositionHasBeenNegative_ = true;
        }
        Eigen::Vector3d cartesianAcceleration;

        cartesianAcceleration[ 0 ] = ( 1.0 / ( cartesianState[ 0 ] + ( std::pow( cartesianState[ 1 ], 2 ) / cartesianState[ 0 ] ) ) )
                * ( cylindricalState[ 0 ] * cylindricalAcceleration[ 0 ]
                - ( cartesianState[ 1 ] / cartesianState[ 0 ] ) * std::pow( cylindricalState[ 0 ], 1 ) * cylindricalAcceleration[ 1 ] );

        cartesianAcceleration[ 1 ] = ( std::pow( cylindricalState[ 0 ], 1 ) * cylindricalAcceleration[ 1 ]
                + cartesianState[ 1 ] * cartesianAcceleration[ 0 ] ) / cartesianState[ 0 ];

        cartesianAcceleration[ 2 ] = cylindricalAcceleration[ 2 ];
        thrustAccelerationVectorCache_[ timeSinceDeparture ] = cartesianAcceleration;
    }
    return thrustAccelerationVectorCache_[ timeSinceDeparture ];

}

Eigen::Vector3d HodographicShaping::computeCurrentThrustAcceleration( const double currentTime,
                                                                      const double timeOffset )
{
    return computeCurrentThrustAcceleration( currentTime - timeOffset );

}

//! Compute direction cartesian acceleration.
Eigen::Vector3d HodographicShaping::computeCurrentThrustAccelerationDirection(
        double timeSinceDeparture )
{
    return computeCurrentThrustAcceleration( timeSinceDeparture ).normalized();
}


//! Compute DeltaV.
double HodographicShaping::computeDeltaV( )
{
    radialPositionHasBeenNegative_ = false;

    // Define the derivative of the deltaV, ie thrust acceleration function, as a function of time.
    std::function< double( const double ) > derivativeFunctionDeltaVtest =
            std::bind( &HodographicShaping::computeCurrentThrustAccelerationMagnitude, this,
                       std::placeholders::_1 );

    // Define numerical quadrature.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeFunctionDeltaVtest, quadratureSettings_, timeOfFlight_ );

    double deltaV = quadrature->getQuadrature( );
    if( radialPositionHasBeenNegative_ )
    {
        deltaV *= 1000.0;
    }
    return deltaV;
}


} // namespace shape_based_methods
} // namespace tudat
