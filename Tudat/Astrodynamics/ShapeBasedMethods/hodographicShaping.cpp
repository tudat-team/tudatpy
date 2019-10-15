/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */


#include "hodographicShaping.h"
#include <math.h>
#include <iostream>
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include "Tudat/Basics/timeType.h"
#include "Tudat/Astrodynamics/Propulsion/thrustAccelerationModel.h"


namespace tudat
{
namespace shape_based_methods
{


//! Constructor which sets radial, normal and axial velocity functions and boundary conditions.
HodographicShaping::HodographicShaping(
        const Eigen::Vector6d initialState,
        const Eigen::Vector6d finalState,
        const double timeOfFlight,
        const int numberOfRevolutions,
        simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToPropagate,
        const std::string centralBody,
//        const double centralBodyGravitationalParameter,
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& radialVelocityFunctionComponents,
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& normalVelocityFunctionComponents,
        std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& axialVelocityFunctionComponents,
        const Eigen::VectorXd freeCoefficientsRadialVelocityFunction,
        const Eigen::VectorXd freeCoefficientsNormalVelocityFunction,
        const Eigen::VectorXd freeCoefficientsAxialVelocityFunction,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings ) :
    ShapeBasedMethodLeg( initialState, finalState, timeOfFlight, bodyMap, bodyToPropagate, centralBody, integratorSettings ),
    initialState_( initialState ),
    finalState_( finalState ),
    timeOfFlight_( timeOfFlight ),
    numberOfRevolutions_( numberOfRevolutions ),
//    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    freeCoefficientsRadialVelocityFunction_( freeCoefficientsRadialVelocityFunction ),
    freeCoefficientsNormalVelocityFunction_( freeCoefficientsNormalVelocityFunction ),
    freeCoefficientsAxialVelocityFunction_( freeCoefficientsAxialVelocityFunction )
{
    // Retrieve gravitational parameter of the central body.
    centralBodyGravitationalParameter_ = bodyMap_[ centralBody_ ]->getGravityFieldModel()->getGravitationalParameter( );

    // Define composite function in radial direction.
    Eigen::VectorXd radialVelocityCoefficients;
    radialVelocityCoefficients.resize( 3 + freeCoefficientsRadialVelocityFunction_.size() );
    radialVelocityCoefficients.segment( 0, 3 ) = Eigen::Vector3d::Zero();
    radialVelocityCoefficients.segment( 3, freeCoefficientsRadialVelocityFunction_.size() ) = freeCoefficientsRadialVelocityFunction_;

    radialVelocityFunction_ = std::make_shared< shape_based_methods::CompositeFunctionHodographicShaping >( radialVelocityFunctionComponents, radialVelocityCoefficients );

    // Define composite function in normal direction.
    Eigen::VectorXd normalVelocityCoefficients;
    normalVelocityCoefficients.resize( 3 + freeCoefficientsNormalVelocityFunction_.size() );
    normalVelocityCoefficients.segment( 0, 3 ) = Eigen::Vector3d::Zero();
    normalVelocityCoefficients.segment( 3, freeCoefficientsNormalVelocityFunction_.size() ) = freeCoefficientsNormalVelocityFunction_;

    normalVelocityFunction_ = std::make_shared< shape_based_methods::CompositeFunctionHodographicShaping >( normalVelocityFunctionComponents, normalVelocityCoefficients );

    // Define composite function in axial direction.
    Eigen::VectorXd axialVelocityCoefficients;
    axialVelocityCoefficients.resize( 3 + freeCoefficientsAxialVelocityFunction_.size() );
    axialVelocityCoefficients.segment( 0, 3 ) = Eigen::Vector3d::Zero();
    axialVelocityCoefficients.segment( 3, freeCoefficientsAxialVelocityFunction_.size() ) = freeCoefficientsAxialVelocityFunction_;

    axialVelocityFunction_ = std::make_shared< shape_based_methods::CompositeFunctionHodographicShaping >( axialVelocityFunctionComponents, axialVelocityCoefficients );


    // Compute initial cylindrical state.
    Eigen::Vector6d initialCylindricalState = coordinate_conversions::convertCartesianToCylindricalState( initialState_ );

    // Compute final cylindrical state.
    Eigen::Vector6d finalCylindricalState = coordinate_conversions::convertCartesianToCylindricalState( finalState_ );

    // Compute number of free coefficients for each velocity function.
    numberOfFreeCoefficientsRadialVelocityFunction_ = radialVelocityFunction_->getNumberOfCompositeFunctionComponents() - 3;
    numberOfFreeCoefficientsNormalVelocityFunction_ = normalVelocityFunction_->getNumberOfCompositeFunctionComponents() - 3;
    numberOfFreeCoefficientsAxialVelocityFunction = axialVelocityFunction_->getNumberOfCompositeFunctionComponents() - 3;
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
        radialBoundaryConditions_.push_back( initialCylindricalState[ 0 ] );   // initial radial distance
        radialBoundaryConditions_.push_back( finalCylindricalState[ 0 ] );     // final radial distance
        radialBoundaryConditions_.push_back( initialCylindricalState[ 3 ] );   // initial radial velocity
        radialBoundaryConditions_.push_back( finalCylindricalState[ 3 ] );     // final radial velocity

        // Set boundary conditions in the normal direction.
        normalBoundaryConditions_.push_back( initialCylindricalState[ 4 ] );   // initial normal velocity
        normalBoundaryConditions_.push_back( finalCylindricalState[ 4 ] );     // final normal velocity
        normalBoundaryConditions_.push_back( numberOfRevolutions_ * 2.0 * mathematical_constants::PI
                                             + basic_mathematics::computeModulo( ( finalCylindricalState[ 1 ] - initialCylindricalState[ 1 ] ),
                2.0 * mathematical_constants::PI ) ); // final polar angle

        // Set boundary conditions in the axial direction.
        axialBoundaryConditions_.push_back( initialCylindricalState[ 2 ] );    // initial axial distance
        axialBoundaryConditions_.push_back( finalCylindricalState[ 2 ] );      // final axial distance
        axialBoundaryConditions_.push_back( initialCylindricalState[ 5 ] );    // initial axial velocity
        axialBoundaryConditions_.push_back( finalCylindricalState[ 5 ] );      // final axial velocity

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
Eigen::Matrix3d HodographicShaping::computeInverseMatrixRadialOrAxialBoundaries( std::shared_ptr< CompositeFunctionHodographicShaping > velocityFunction )
{
    Eigen::Matrix3d matrixBoundaryValues, inverseMatrixBoundaryValues;
    matrixBoundaryValues <<
            velocityFunction->getComponentFunctionIntegralCurrentTime( 0, timeOfFlight_ )
                            - velocityFunction->getComponentFunctionIntegralCurrentTime( 0, 0.0 ),
            velocityFunction->getComponentFunctionIntegralCurrentTime( 1, timeOfFlight_ )
            - velocityFunction->getComponentFunctionIntegralCurrentTime( 1, 0.0 ),
            velocityFunction->getComponentFunctionIntegralCurrentTime( 2, timeOfFlight_ )
            - velocityFunction->getComponentFunctionIntegralCurrentTime( 2, 0.0 ),
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
Eigen::Matrix2d HodographicShaping::computeInverseMatrixNormalBoundaries( std::shared_ptr< CompositeFunctionHodographicShaping > velocityFunction )
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
                freeCoefficients( i ) * ( radialVelocityFunction_->getComponentFunctionIntegralCurrentTime( i + 3, timeOfFlight_ )
                                          - radialVelocityFunction_->getComponentFunctionIntegralCurrentTime( i + 3, 0.0 ) );
        vectorBoundaryConditionsRadial[ 1 ] -= freeCoefficients( i ) * radialVelocityFunction_->getComponentFunctionCurrentValue( i + 3, 0.0 );
        vectorBoundaryConditionsRadial[ 2 ] -= freeCoefficients( i ) * radialVelocityFunction_->getComponentFunctionCurrentValue( i + 3, timeOfFlight_ );
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
        radialVelocityFunctionCoefficients << fixedCoefficientsRadial, freeCoefficients.segment( 0, numberOfFreeCoefficientsRadialVelocityFunction_);
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
                * ( axialVelocityFunction_->getComponentFunctionIntegralCurrentTime( i + 3 , timeOfFlight_ )
                    - axialVelocityFunction_->getComponentFunctionIntegralCurrentTime( i + 3, 0.0 ) );
        vectorBoundaryConditionsAxial[ 1 ] -= freeCoefficients( i )
                * axialVelocityFunction_->getComponentFunctionCurrentValue( i + 3, 0.0 );
        vectorBoundaryConditionsAxial[ 2 ] -= freeCoefficients( i )
                * axialVelocityFunction_->getComponentFunctionCurrentValue( i + 3, timeOfFlight_ );
    }

    // Compute fixed coefficients.
    Eigen::Vector3d fixedCoefficientsAxial = inverseAxialMatrixBoundaryValues_ * vectorBoundaryConditionsAxial;

    // Create vector containing all axial velocity function coefficients.
    Eigen::VectorXd axialVelocityFunctionCoefficients = Eigen::VectorXd::Zero( fixedCoefficientsAxial.size() + numberOfFreeCoefficientsAxialVelocityFunction );

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

    Eigen::Vector2d initialAndFinalValuesThirdComponentFunction( - normalVelocityFunction_->getComponentFunctionCurrentValue( 2, 0.0 ),
                                  - normalVelocityFunction_->getComponentFunctionCurrentValue( 2, timeOfFlight_ ) );

    // Define matrix K, as proposed in ... (ADD PROPER REFERENCE AND EQUATION NUMBER)
    Eigen::Vector2d matrixK;
    matrixK = inverseMatrixNormalBoundaryValues_ * initialAndFinalValuesThirdComponentFunction;


    // Define angular velocity due to the third component of the composite function only.
    std::function< double( const double ) > derivativePolarAngleDueToThirdComponent = [ = ] ( const double currentTime ){

        double angularVelocityDueToThirdComponent = ( matrixK( 0 ) * normalVelocityFunction_->getComponentFunctionCurrentValue( 0, currentTime )
                + matrixK( 1 ) * normalVelocityFunction_->getComponentFunctionCurrentValue( 1, currentTime )
                + normalVelocityFunction_->getComponentFunctionCurrentValue( 2, currentTime ) ) / computeCurrentRadialDistance( currentTime );

        return angularVelocityDueToThirdComponent;

    };


    // Define the angular velocity due to all the other components of the composite function, once combined.
    std::function< double( const double ) > derivativePolarAngleDueToOtherComponents = [ = ] ( const double currentTime ){

        double angularVelocityDueToFreeCoefficients = 0.0;
        for( int j = 0 ; j < numberOfFreeCoefficientsNormalVelocityFunction_ ; j++ )
        {
            angularVelocityDueToFreeCoefficients += freeCoefficients( j )
                             * normalVelocityFunction_->getComponentFunctionCurrentValue( j + 3 , currentTime );
        }

        double angularVelocityDueToOtherComponents = ( matrixL( 0 ) * normalVelocityFunction_->getComponentFunctionCurrentValue( 0 , currentTime )
                                                     + matrixL( 1 ) * normalVelocityFunction_->getComponentFunctionCurrentValue( 1 , currentTime )
                                                     + angularVelocityDueToFreeCoefficients )
                / computeCurrentRadialDistance( currentTime );

        return angularVelocityDueToOtherComponents;

    };

    // Define numerical quadratures.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > integratorPolarAngleDueToThirdComponentTest =
            numerical_quadrature::createQuadrature( derivativePolarAngleDueToThirdComponent, quadratureSettings_, timeOfFlight_ );

    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > integratorPolarAngleDueToOtherComponentsTest =
            numerical_quadrature::createQuadrature( derivativePolarAngleDueToOtherComponents, quadratureSettings_, timeOfFlight_ );

    return ( normalBoundaryConditions_[ 2 ] - integratorPolarAngleDueToOtherComponentsTest->getQuadrature() )
            / integratorPolarAngleDueToThirdComponentTest->getQuadrature();

}


//! Compute velocity components.
Eigen::Vector3d HodographicShaping::computeVelocityVectorInCylindricalCoordinates( double currentTime )
{
    Eigen::Vector3d velocityComponents;

    // Compute radial velocity.
    velocityComponents[ 0 ] =  radialVelocityFunction_->evaluateCompositeFunctionCurrentTime( currentTime );

    // Compute normal velocity.
    velocityComponents[ 1 ] = normalVelocityFunction_->evaluateCompositeFunctionCurrentTime( currentTime );

    // Compute axial velocity.
    velocityComponents[ 2 ] = axialVelocityFunction_->evaluateCompositeFunctionCurrentTime( currentTime );

    return velocityComponents;
}


//! Compute angular velocity.
double HodographicShaping::evaluateDerivativePolarAngleWrtTime( const double currentTime )
{
    // Return derivative of the polar angle w.r.t. time, i.e. angular velocity.
    return normalVelocityFunction_->evaluateCompositeFunctionCurrentTime( currentTime )
           / ( radialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentTime( currentTime )
               - radialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentTime( 0.0 ) + radialBoundaryConditions_[0] );
}


//! Compute radial distance.
double HodographicShaping::computeCurrentRadialDistance( const double currentTime )
{
    // Compute radial distance from the central body.
    return radialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentTime( currentTime )
            - radialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentTime( 0.0 ) + radialBoundaryConditions_[0];
}

//! Compute current polar angle.
double HodographicShaping::computeCurrentPolarAngle( double currentTime )
{

    // Define the derivative of the polar angle, ie angular velocity function, as a function of time.
    std::function< double( const double ) > derivativeFunctionPolarAngle = [ = ] ( const double currentTime ){

        return evaluateDerivativePolarAngleWrtTime( currentTime );

    } ;

    // Define numerical quadrature.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeFunctionPolarAngle, quadratureSettings_, currentTime );

    return quadrature->getQuadrature( ) + coordinate_conversions::convertCartesianToCylindricalState( initialState_ )[ 1 ];

}

double HodographicShaping::computeCurrentAxialDistance( const double currentTime )
{
    // Compute axial distance from the central body.
    return axialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentTime( currentTime )
            - axialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentTime( 0.0 ) + axialBoundaryConditions_[0];
}

//! Compute current cylindrical state.
Eigen::Vector6d HodographicShaping::computeStateVectorInCylindricalCoordinates( const double currentTime ){

    Eigen::Vector3d velocityCylindricalCoordinates = computeVelocityVectorInCylindricalCoordinates( currentTime );

    Eigen::Vector6d cylindricalState = ( Eigen::Vector6d() <<
                                         computeCurrentRadialDistance( currentTime ),
                                         computeCurrentPolarAngle( currentTime ),
                                         computeCurrentAxialDistance( currentTime ),
                                         velocityCylindricalCoordinates[ 0 ],
                                         velocityCylindricalCoordinates[ 1 ],
                                         velocityCylindricalCoordinates[ 2 ] ).finished();

    return cylindricalState;
}


//! Compute current cartesian state.
Eigen::Vector6d HodographicShaping::computeCurrentStateVector( const double currentTime ){

    return coordinate_conversions::convertCylindricalToCartesianState( computeStateVectorInCylindricalCoordinates( currentTime ) );

}


//! Compute thrust acceleration components.
Eigen::Vector3d HodographicShaping::computeThrustAccelerationInCylindricalCoordinates( double currentTime )
{

    Eigen::Vector3d thrustAccelerationComponents;

    // Compute radial distance from the central body.
    double radialDistance = radialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentTime( currentTime )
            - radialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentTime( 0.0 ) + radialBoundaryConditions_[0];

    // Compute axial distance from the central body.
    double axialDistance = axialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentTime( currentTime )
            - axialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentTime( 0.0 ) + axialBoundaryConditions_[0];

    // Compute distance from the central body.
    double distanceFromCentralBody = sqrt( pow( radialDistance, 2.0 ) + pow( axialDistance, 2.0 ) );

    // Computation of normal velocity.
    double normalVelocity = normalVelocityFunction_->evaluateCompositeFunctionCurrentTime( currentTime );

    // Compute angular velocity.
    double angularVelocity = normalVelocity / radialDistance;

    // Compute radial thrust acceleration.
    thrustAccelerationComponents[ 0 ] = radialVelocityFunction_->evaluateCompositeFunctionDerivativeCurrentTime( currentTime )
            - angularVelocity * normalVelocity + centralBodyGravitationalParameter_ / std::pow( distanceFromCentralBody, 3.0 ) * radialDistance;

    // Compute normal thrust acceleration.
    thrustAccelerationComponents[ 1 ] = normalVelocityFunction_->evaluateCompositeFunctionDerivativeCurrentTime( currentTime )
            + angularVelocity * radialVelocityFunction_->evaluateCompositeFunctionCurrentTime( currentTime );

    // Compute axial thrust acceleration.
    thrustAccelerationComponents[ 2 ] = axialVelocityFunction_->evaluateCompositeFunctionDerivativeCurrentTime( currentTime )
            + centralBodyGravitationalParameter_ / std::pow( distanceFromCentralBody, 3.0 ) * axialDistance;

    // Return total thrust acceleration.
    return thrustAccelerationComponents;
}


double HodographicShaping::computeCurrentThrustAccelerationMagnitude(
        const double currentTime, std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings ){

//    // Compute radial distance from the central body.
//    double radialDistance = radialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentTime( currentTime )
//            - radialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentTime( 0.0 ) + radialBoundaryConditions_[0];

//    // Compute axial distance from the central body.
//    double axialDistance = axialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentTime( currentTime )
//            - axialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentTime( 0.0 ) + axialBoundaryConditions_[0];

//    // Compute distance from central body.
//    double distanceFromCentralBody = sqrt( pow( radialDistance, 2.0 ) + pow( axialDistance, 2.0 ) );

//    // Computation of normal velocity.
//    double normalVelocity = normalVelocityFunction_->evaluateCompositeFunctionCurrentTime( currentTime );

//    // Compute angular velocity.
//    double angularVelocity = normalVelocity / radialDistance;

//    // Compute radial thrust acceleration.
//    double radialThrustAcceleration =
//            radialVelocityFunction_->evaluateCompositeFunctionDerivativeCurrentTime( currentTime ) - angularVelocity * normalVelocity
//            + centralBodyGravitationalParameter_ / pow( distanceFromCentralBody, 3.0 ) * radialDistance;

//    // Compute normal thrust acceleration.
//    double normalThrustAcceleration = normalVelocityFunction_->evaluateCompositeFunctionDerivativeCurrentTime( currentTime )
//            + angularVelocity * radialVelocityFunction_->evaluateCompositeFunctionCurrentTime( currentTime );

//    // Compute thrust acceleration.
//    double axialThrustAcceleration = axialVelocityFunction_->evaluateCompositeFunctionDerivativeCurrentTime( currentTime )
//            + centralBodyGravitationalParameter_ / pow( distanceFromCentralBody, 3.0 ) * axialDistance;

//    // Return total thrust acceleration.
    return computeThrustAccelerationInCylindricalCoordinates( currentTime ).norm();
}


//! Compute cartesian acceleration.
Eigen::Vector3d HodographicShaping::computeCurrentThrustAccelerationVector( double currentTime ){

    Eigen::Vector3d cylindricalAcceleration = computeThrustAccelerationInCylindricalCoordinates( currentTime );
    Eigen::Vector3d cylindricalState = computeStateVectorInCylindricalCoordinates( currentTime ).segment(0,3);
    Eigen::Vector3d cartesianState = computeCurrentStateVector( currentTime ).segment(0,3);

    Eigen::Vector3d cartesianAcceleration;

    cartesianAcceleration[ 0 ] = ( 1.0 / ( cartesianState[ 0 ] + ( std::pow( cartesianState[ 1 ], 2 ) / cartesianState[ 0 ] ) ) )
            * ( cylindricalState[ 0 ] * cylindricalAcceleration[ 0 ]
            - ( cartesianState[ 1 ] / cartesianState[ 0 ] ) * std::pow( cylindricalState[ 0 ], 1 ) * cylindricalAcceleration[ 1 ] );

    cartesianAcceleration[ 1 ] = ( std::pow( cylindricalState[ 0 ], 1 ) * cylindricalAcceleration[ 1 ]
            + cartesianState[ 1 ] * cartesianAcceleration[ 0 ] ) / cartesianState[ 0 ];

    cartesianAcceleration[ 2 ] = cylindricalAcceleration[ 2 ];

    return cartesianAcceleration;

}

//! Compute direction cartesian acceleration.
Eigen::Vector3d HodographicShaping::computeCurrentThrustAccelerationDirection(
        double currentTime, std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    return computeCurrentThrustAccelerationVector( currentTime ).normalized();
}


//! Compute DeltaV.
double HodographicShaping::computeDeltaV( )
{

    // Define the derivative of the deltaV, ie thrust acceleration function, as a function of time.
    std::function< Eigen::Vector1d( const double, const Eigen::Vector1d& ) > derivativeFunctionDeltaV = [ = ]
            ( const double currentTime, const Eigen::Vector1d& currentState ){

        std::function< double( const double ) > specificImpulseFunction = [ = ]( double time )
        {
            return 0.0;
        };

        Eigen::Vector1d thrustAcceleration;
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings;
        thrustAcceleration[ 0 ] = computeCurrentThrustAccelerationMagnitude( currentTime, specificImpulseFunction, integratorSettings );
        return thrustAcceleration;

    } ;

    // Define the derivative of the deltaV, ie thrust acceleration function, as a function of time.
    std::function< double( const double ) > derivativeFunctionDeltaVtest = [ = ] ( const double currentTime ){

        return computeThrustAccelerationInCylindricalCoordinates( currentTime ).norm();

    } ;

    // Define numerical quadrature.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeFunctionDeltaVtest, quadratureSettings_, timeOfFlight_ );

    return quadrature->getQuadrature( );

}


////! Get low-thrust acceleration model from shaping method.
//std::shared_ptr< propulsion::ThrustAcceleration > HodographicShaping::getLowThrustAccelerationModel(
//        simulation_setup::NamedBodyMap& bodyMap,
//        const std::string& bodyToPropagate,
//        std::function< double( const double ) > specificImpulseFunction )
//{
//    // Define thrust direction settings from the direction of thrust acceleration retrieved from the shaping method.
//    std::shared_ptr< simulation_setup::CustomThrustDirectionSettings > thrustDirectionSettings =
//            std::make_shared< simulation_setup::CustomThrustDirectionSettings >(
//                std::bind( &HodographicShaping::computeCurrentThrustAccelerationDirection, this, std::placeholders::_1 ) );

//    std::shared_ptr< simulation_setup::Body > vehicle = bodyMap[ bodyToPropagate ];

//    // Define thrust magnitude function from the shaped trajectory.
//    std::function< double( const double ) > thrustMagnitudeFunction = [ = ]( const double currentTime )
//    {
//        // Compute current acceleration.
//        double currentAcceleration = computeCurrentThrustAccelerationMagnitude( currentTime );

//        // Compute current mass of the vehicle.
//        double currentMass = vehicle->getBodyMass();

//        // Compute and return magnitude of the low-thrust force.
//        return currentAcceleration * currentMass;
//    };

//    // Define thrust magnitude settings from thrust magnitude function.
//    std::shared_ptr< simulation_setup::FromFunctionThrustMagnitudeSettings > thrustMagnitudeSettings
//            = std::make_shared< simulation_setup::FromFunctionThrustMagnitudeSettings >(
//                thrustMagnitudeFunction, specificImpulseFunction );

//    // Define thrust acceleration settings.
//    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > thrustAccelerationSettings =
//            std::make_shared< simulation_setup::ThrustAccelerationSettings >(
//                thrustDirectionSettings, thrustMagnitudeSettings );

//    // Create low thrust acceleration model.
//    std::shared_ptr< propulsion::ThrustAcceleration > lowThrustAccelerationModel = createThrustAcceleratioModel( thrustAccelerationSettings, bodyMap, bodyToPropagate );

//    return lowThrustAccelerationModel;
//}


//void HodographicShaping::computeShapedTrajectoryAndFullPropagation(
//        simulation_setup::NamedBodyMap& bodyMap,
//        std::function< double( const double ) > specificImpulseFunction,
//        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
//        std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
//                std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > >& propagatorSettings,
//        std::map< double, Eigen::VectorXd >& fullPropagationResults,
//        std::map< double, Eigen::VectorXd >& shapingMethodResults,
//        std::map< double, Eigen::VectorXd >& dependentVariablesHistory,
//        const bool isMassPropagated ){

//    fullPropagationResults.clear();
//    shapingMethodResults.clear();
//    dependentVariablesHistory.clear();

//    std::string bodyToPropagate = propagatorSettings.first->bodiesToIntegrate_[ 0 ];

//    // Compute halved time of flight.
//    double halvedTimeOfFlight = timeOfFlight_ / 2.0;


//    // Compute state at half of the time of flight.
//    Eigen::Vector6d initialStateAtHalvedTimeOfFlight = computeCurrentStateVector( halvedTimeOfFlight );


//    // Create low thrust acceleration model.
//    std::shared_ptr< propulsion::ThrustAcceleration > lowThrustAccelerationModel =
//            getLowThrustAccelerationModel( /*bodyMap, bodyToPropagate,*/ specificImpulseFunction );

//    basic_astrodynamics::AccelerationMap accelerationMap = propagators::getAccelerationMapFromPropagatorSettings(
//                std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< double > >( propagatorSettings.first ) );

//    accelerationMap[ bodyToPropagate ][ bodyToPropagate ].push_back( lowThrustAccelerationModel );


//    // Create complete propagation settings (backward and forward propagations).
//    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
//            std::shared_ptr< propagators::PropagatorSettings< double > > > completePropagatorSettings;



//    // Define translational state propagation settings
//    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
//            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > translationalStatePropagatorSettings;

//    // Define backward translational state propagation settings.
//    translationalStatePropagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
//                        ( propagatorSettings.first->centralBodies_, accelerationMap, propagatorSettings.first->bodiesToIntegrate_,
//                          initialStateAtHalvedTimeOfFlight, propagatorSettings.first->getTerminationSettings(),
//                          propagatorSettings.first->propagator_, propagatorSettings.first->getDependentVariablesToSave() );

//    // Define forward translational state propagation settings.
//    translationalStatePropagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
//                        ( propagatorSettings.second->centralBodies_, accelerationMap, propagatorSettings.second->bodiesToIntegrate_,
//                          initialStateAtHalvedTimeOfFlight, propagatorSettings.second->getTerminationSettings(),
//                          propagatorSettings.second->propagator_, propagatorSettings.second->getDependentVariablesToSave() );



//    // If translational state and mass are propagated concurrently.
//    if ( isMassPropagated )
//    {
//        // Create mass rate models
//        std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
//        massRateModels[ bodyToPropagate ] = createMassRateModel( bodyToPropagate, std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
//                                                           bodyMap, accelerationMap );


////        // Propagate mass until half of the time of flight.
////        std::shared_ptr< propagators::PropagatorSettings< double > > massPropagatorSettingsToHalvedTimeOfFlight =
////                std::make_shared< propagators::MassPropagatorSettings< double > >( std::vector< std::string >{ bodyToPropagate }, massRateModels,
////                    ( Eigen::Vector1d() << bodyMap[ bodyToPropagate ]->getBodyMass() ).finished(),
////                    std::make_shared< propagators::PropagationTimeTerminationSettings >( halvedTimeOfFlight, true ) );

////        integratorSettings->initialTime_ = 0.0;

////        // Create dynamics simulation object.
////        propagators::SingleArcDynamicsSimulator< double, double > dynamicsSimulator(
////                    bodyMap, integratorSettings, massPropagatorSettingsToHalvedTimeOfFlight, true, false, false );

//        // Propagate spacecraft mass until half of the time of flight.
////        std::map< double, Eigen::VectorXd > propagatedMass = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );
//        double massAtHalvedTimeOfFlight = computeCurrentMass( 0.0, halvedTimeOfFlight, initialMass_, specificImpulseFunction, integratorSettings ); // propagatedMass.rbegin()->second[ 0 ];

//        // Create settings for propagating the mass of the vehicle.
//        std::pair< std::shared_ptr< propagators::MassPropagatorSettings< double > >,
//                std::shared_ptr< propagators::MassPropagatorSettings< double > > > massPropagatorSettings;

//        // Define backward mass propagation settings.
//        massPropagatorSettings.first = std::make_shared< propagators::MassPropagatorSettings< double > >(
//                    std::vector< std::string >{ bodyToPropagate }, massRateModels, ( Eigen::Matrix< double, 1, 1 >( ) << massAtHalvedTimeOfFlight ).finished( ),
//                                                                      propagatorSettings.first->getTerminationSettings() );

//        // Define forward mass propagation settings.
//        massPropagatorSettings.second = std::make_shared< propagators::MassPropagatorSettings< double > >(
//                    std::vector< std::string >{ bodyToPropagate }, massRateModels, ( Eigen::Matrix< double, 1, 1 >( ) << massAtHalvedTimeOfFlight ).finished( ),
//                                                                      propagatorSettings.second->getTerminationSettings() );


//        // Create list of propagation settings.
//        std::pair< std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > >,
//                std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > > > propagatorSettingsVector;

//        // Backward propagator settings vector.
//        propagatorSettingsVector.first.push_back( translationalStatePropagatorSettings.first );
//        propagatorSettingsVector.first.push_back( massPropagatorSettings.first );

//        // Forward propagator settings vector.
//        propagatorSettingsVector.second.push_back( translationalStatePropagatorSettings.second );
//        propagatorSettingsVector.second.push_back( massPropagatorSettings.second );


//        // Backward hybrid propagation settings.
//        completePropagatorSettings.first = std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector.first,
//                    propagatorSettings.first->getTerminationSettings(), propagatorSettings.first->getDependentVariablesToSave() );

//        // Forward hybrid propagation settings.
//        completePropagatorSettings.second = std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector.second,
//                    propagatorSettings.second->getTerminationSettings(), propagatorSettings.second->getDependentVariablesToSave() );


//    }

//    // If only translational state is propagated.
//    else
//    {
//        // Backward hybrid propagation settings.
//        completePropagatorSettings.first = translationalStatePropagatorSettings.first;

//        // Forward hybrid propagation settings.
//        completePropagatorSettings.second =  translationalStatePropagatorSettings.second;
//    }


//    // Initialise integrator at half of the time of flight.
//    integratorSettings->initialTime_ = halvedTimeOfFlight;

//    // Perform forward propagation.
//    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationForwards( bodyMap, integratorSettings, completePropagatorSettings.second );
//    std::map< double, Eigen::VectorXd > stateHistoryFullProblemForwardPropagation = dynamicsSimulatorIntegrationForwards.getEquationsOfMotionNumericalSolution( );
//    std::map< double, Eigen::VectorXd > dependentVariableHistoryForwardPropagation = dynamicsSimulatorIntegrationForwards.getDependentVariableHistory( );

//    // Compute and save full propagation and shaping method results along the forward propagation direction.
//    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemForwardPropagation.begin( );
//         itr != stateHistoryFullProblemForwardPropagation.end( ); itr++ )
//    {
//        shapingMethodResults[ itr->first ] = computeCurrentStateVector( itr->first );
//        fullPropagationResults[ itr->first ] = itr->second;
//        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryForwardPropagation[ itr->first ];
//    }


//    // Define backward propagator settings variables.
//    integratorSettings->initialTimeStep_ = - integratorSettings->initialTimeStep_;
//    integratorSettings->initialTime_ = halvedTimeOfFlight;

//    // Perform the backward propagation.
//    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationBackwards( bodyMap, integratorSettings, completePropagatorSettings.first );
//    std::map< double, Eigen::VectorXd > stateHistoryFullProblemBackwardPropagation = dynamicsSimulatorIntegrationBackwards.getEquationsOfMotionNumericalSolution( );
//    std::map< double, Eigen::VectorXd > dependentVariableHistoryBackwardPropagation = dynamicsSimulatorIntegrationBackwards.getDependentVariableHistory( );

//    // Compute and save full propagation and shaping method results along the backward propagation direction
//    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemBackwardPropagation.begin( );
//         itr != stateHistoryFullProblemBackwardPropagation.end( ); itr++ )
//    {
//        shapingMethodResults[ itr->first ] = computeCurrentStateVector( itr->first );
//        fullPropagationResults[ itr->first ] = itr->second;
//        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryBackwardPropagation[ itr->first ];
//    }

//    // Reset initial integrator settings.
//    integratorSettings->initialTimeStep_ = - integratorSettings->initialTimeStep_;

//}

} // namespace shape_based_methods
} // namespace tudat
