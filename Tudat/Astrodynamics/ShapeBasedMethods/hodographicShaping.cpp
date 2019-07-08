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
HodographicShaping::HodographicShaping( CompositeFunction& radialVelocityFunction,
                                        CompositeFunction& normalVelocityFunction,
                                        CompositeFunction& axialVelocityFunction,
                                        Eigen::Vector6d initialState,
                                        std::vector< double > boundaryConditionsRadial,
                                        std::vector< double > boundaryConditionsNormal,
                                        std::vector< double > boundaryConditionsAxial,
                                        std::vector< double > boundaryValuesTime,
                                        double initialTime,
                                        double finalTime ):
    radialVelocityFunction_( radialVelocityFunction ),
    normalVelocityFunction_( normalVelocityFunction ),
    axialVelocityFunction_( axialVelocityFunction ),
    initialState_( initialState ),
    initialTime_( initialTime ),
    finalTime_( finalTime )
{

    // Compute number of free coefficients for each velocity function.
    numberOfFreeCoefficientsRadial = radialVelocityFunction.getNumberOfCompositeFunctionComponents()
                                    - std::min( 3, static_cast< int >( boundaryConditionsRadial.size() ) );
    numberOfFreeCoefficientsNormal = normalVelocityFunction.getNumberOfCompositeFunctionComponents()
                                    - std::min( 3, static_cast< int >( boundaryConditionsNormal.size() ) );
    numberOfFreeCoefficientsAxial = axialVelocityFunction.getNumberOfCompositeFunctionComponents()
                                    - std::min( 3, static_cast< int >( boundaryConditionsAxial.size() ) );
    if ( numberOfFreeCoefficientsRadial < 0 || numberOfFreeCoefficientsNormal < 0 || numberOfFreeCoefficientsAxial < 0 )
    {
        std::cerr << "The number of base functions in one of the velocity functions is smaller than "
             << "the number of corresponding boundary conditions. The boundary conditions cannot be set!\n";
    }
    else
    {
        // Set boundary conditions.
        boundaryConditionsRadial_ = boundaryConditionsRadial;
        boundaryConditionsNormal_ = boundaryConditionsNormal;
        boundaryConditionsAxial_ = boundaryConditionsAxial;
        boundaryValuesTime_ = boundaryValuesTime;

        // Compute inverse of matrices containing boundary values.
        inverseMatrixBoundaryValuesRadial = computeInverseMatrixRadialOrAxialBoundaries( radialVelocityFunction );
        inverseMatrixBoundaryValuesNormal = computeInverseMatrixNormalBoundaries( normalVelocityFunction );
        inverseMatrixBoundaryValuesAxial = computeInverseMatrixRadialOrAxialBoundaries( axialVelocityFunction );
    }

    // Define numerical quadrature settings, required to compute the current polar angle and final deltaV.
    quadratureSettings_ = std::make_shared< numerical_quadrature::GaussianQuadratureSettings< double > >( initialTime_, 64 );
}

//! Set boundary conditions.
void HodographicShaping::setBoundaryConditions( std::vector< double > boundaryConditionsRadial,
                                                std::vector< double > boundaryConditionsNormal,
                                                std::vector< double > boundaryConditionsAxial,
                                                std::vector< double > boundaryValuesTime )
{
    // Clear any previously set boundary conditions.
    boundaryConditionsRadial_.clear();
    boundaryConditionsNormal_.clear();
    boundaryConditionsAxial_.clear();
    boundaryValuesTime_.clear();

    // Compute number of free coefficients for each velocity function.
    numberOfFreeCoefficientsRadial = radialVelocityFunction_.getNumberOfCompositeFunctionComponents()
                                    - std::min( 3, static_cast< int >( boundaryConditionsRadial.size() ) );
    numberOfFreeCoefficientsNormal = normalVelocityFunction_.getNumberOfCompositeFunctionComponents()
                                    - std::min( 3, static_cast< int >( boundaryConditionsNormal.size() ) );
    numberOfFreeCoefficientsAxial = axialVelocityFunction_.getNumberOfCompositeFunctionComponents()
                                    - std::min( 3, static_cast< int >( boundaryConditionsAxial.size() ) );

    // Check consistency between number of composite function components and number of free parameters.
    if ( numberOfFreeCoefficientsRadial < 0 || numberOfFreeCoefficientsNormal < 0 || numberOfFreeCoefficientsAxial < 0 )
    {
        std::cerr << "The number of base functions in one of the velocity functions is smaller than "
             << "the number of corresponding boundary conditions. The boundary conditions cannot be reset.";
    }
    else
    {
        // Set boundary conditions.
        boundaryConditionsRadial_ = boundaryConditionsRadial;
        boundaryConditionsNormal_ = boundaryConditionsNormal;
        boundaryConditionsAxial_ = boundaryConditionsAxial;
        boundaryValuesTime_ = boundaryValuesTime;

        // Compute inverse of matrices containing boundary values.
        inverseMatrixBoundaryValuesRadial = computeInverseMatrixRadialOrAxialBoundaries( radialVelocityFunction_ );
        inverseMatrixBoundaryValuesNormal = computeInverseMatrixNormalBoundaries( normalVelocityFunction_ );
        inverseMatrixBoundaryValuesAxial = computeInverseMatrixRadialOrAxialBoundaries( axialVelocityFunction_ );
    }
}


//! Compute inverse of matrix filled with boundary conditions in radial or axial direction.
Eigen::Matrix3d HodographicShaping::computeInverseMatrixRadialOrAxialBoundaries( CompositeFunction& velocityFunction )
{
    Eigen::Matrix3d matrixBoundaryValues, inverseMatrixBoundaryValues;
    matrixBoundaryValues <<
            velocityFunction.getComponentFunctionIntegralCurrentTime( 0, finalTime_ )
                            - velocityFunction.getComponentFunctionIntegralCurrentTime( 0, 0.0 ),
            velocityFunction.getComponentFunctionIntegralCurrentTime( 1, finalTime_ )
            - velocityFunction.getComponentFunctionIntegralCurrentTime( 1, 0.0 ),
            velocityFunction.getComponentFunctionIntegralCurrentTime( 2, finalTime_ )
            - velocityFunction.getComponentFunctionIntegralCurrentTime( 2, 0.0 ),
            velocityFunction.getComponentFunctionCurrentValue( 0, 0.0 ),
            velocityFunction.getComponentFunctionCurrentValue( 1, 0.0 ),
            velocityFunction.getComponentFunctionCurrentValue( 2, 0.0 ),
            velocityFunction.getComponentFunctionCurrentValue( 0, finalTime_ ),
            velocityFunction.getComponentFunctionCurrentValue( 1, finalTime_ ),
            velocityFunction.getComponentFunctionCurrentValue( 2, finalTime_ );

    // Compute inverse of boundary-value matrix.
    inverseMatrixBoundaryValues = matrixBoundaryValues.inverse();

    // Return inverse of boundary-value matrix.
    return inverseMatrixBoundaryValues;
}

//! Compute inverse of matrix filled with boundary conditions in normal direction.
Eigen::Matrix2d HodographicShaping::computeInverseMatrixNormalBoundaries( CompositeFunction& velocityFunction )
{
    Eigen::Matrix2d matrixBoundaryValues, inverseMatrixBoundaryValues;

    matrixBoundaryValues <<
            velocityFunction.getComponentFunctionCurrentValue( 0, 0.0 ),
            velocityFunction.getComponentFunctionCurrentValue( 1, 0.0 ),
            velocityFunction.getComponentFunctionCurrentValue( 0, finalTime_ ),
            velocityFunction.getComponentFunctionCurrentValue( 1, finalTime_ );

    // Compute inverse of boundary-value matrix.
    inverseMatrixBoundaryValues = matrixBoundaryValues.inverse();

    // Return inverse of boundary-value matrix.
    return inverseMatrixBoundaryValues;
}


void HodographicShaping::satisfyRadialBoundaryConditions( Eigen::VectorXd freeCoefficients ){

    // Vector containing boundary conditions on radial distance and initial and final radial velocity.
    Eigen::Vector3d vectorBoundaryConditionsRadial;
    vectorBoundaryConditionsRadial[ 0 ] = boundaryConditionsRadial_[ 1 ] - boundaryConditionsRadial_[ 0 ];
    vectorBoundaryConditionsRadial[ 1 ] = boundaryConditionsRadial_[ 2 ];
    vectorBoundaryConditionsRadial[ 2 ] = boundaryConditionsRadial_[ 3 ];

    // Subtract boundary values of free components of velocity function from corresponding boundary conditions.
    for ( int i = 0 ; i < numberOfFreeCoefficientsRadial ; i++ )
    {
        vectorBoundaryConditionsRadial[ 0 ] -=
                freeCoefficients( i ) * ( radialVelocityFunction_.getComponentFunctionIntegralCurrentTime( i + 3, finalTime_ )
                                          - radialVelocityFunction_.getComponentFunctionIntegralCurrentTime( i + 3, 0.0 ) );
        vectorBoundaryConditionsRadial[ 1 ] -= freeCoefficients( i ) * radialVelocityFunction_.getComponentFunctionCurrentValue( i + 3, 0.0 );
        vectorBoundaryConditionsRadial[ 2 ] -= freeCoefficients( i ) * radialVelocityFunction_.getComponentFunctionCurrentValue( i + 3, finalTime_ );
    }

    // Compute fixed coefficients.
    Eigen::Vector3d fixedCoefficientsRadial = inverseMatrixBoundaryValuesRadial * vectorBoundaryConditionsRadial;

    // Create vector containing all radial velocity function coefficients.
    Eigen::VectorXd radialVelocityFunctionCoefficients =
            Eigen::VectorXd::Zero( fixedCoefficientsRadial.size() + numberOfFreeCoefficientsRadial );

    // Check whether the radial velocity function has free coefficients.
    if ( numberOfFreeCoefficientsRadial == 0 )
    {
        radialVelocityFunctionCoefficients << fixedCoefficientsRadial;
    }
    else
    {
        radialVelocityFunctionCoefficients << fixedCoefficientsRadial, freeCoefficients.segment(0, numberOfFreeCoefficientsRadial);
    }

    // Set the coefficients of the radial velocity function.
    radialVelocityFunction_.resetCompositeFunctionCoefficients( radialVelocityFunctionCoefficients );

}


void HodographicShaping::satisfyNormalBoundaryConditions( Eigen::VectorXd freeCoefficients ){

    // Vector containing boundary conditions on initial and final normal velocity.
    Eigen::Vector2d vectorBoundaryConditionsNormal;
    vectorBoundaryConditionsNormal( 0 ) = boundaryConditionsNormal_[ 0 ];
    vectorBoundaryConditionsNormal( 1 ) = boundaryConditionsNormal_[ 1 ];

    // Subtract boundary values of free velocity function components from corresponding boundary conditions.
    for ( int i = 0 ; i < numberOfFreeCoefficientsNormal ; i++ )
    {
        vectorBoundaryConditionsNormal[ 0 ] -= freeCoefficients( numberOfFreeCoefficientsRadial + i )
                * normalVelocityFunction_.getComponentFunctionCurrentValue( i + 2, 0.0 );
        vectorBoundaryConditionsNormal[ 1 ] -= freeCoefficients( numberOfFreeCoefficientsRadial + i )
                * normalVelocityFunction_.getComponentFunctionCurrentValue( i + 2, finalTime_ );
    }

    // Compute fixed coefficients by satisfying the boundary conditions.
    Eigen::Vector2d fixedCoefficientsNormal = inverseMatrixBoundaryValuesNormal * vectorBoundaryConditionsNormal;

    // Create vector containing all normal velocity function coefficients.
    Eigen::VectorXd normalVelocityFunctionCoefficients =
            Eigen::VectorXd::Zero( fixedCoefficientsNormal.size() + numberOfFreeCoefficientsNormal );

    // Check whether the normal velocity function has free coefficients.
    if ( numberOfFreeCoefficientsNormal == 0 )
    {
        normalVelocityFunctionCoefficients << fixedCoefficientsNormal;
    }
    else
    {
        normalVelocityFunctionCoefficients << fixedCoefficientsNormal,
                freeCoefficients.segment( numberOfFreeCoefficientsRadial, numberOfFreeCoefficientsNormal );
    }

    // Set the coefficients of the normal velocity function.
    normalVelocityFunction_.resetCompositeFunctionCoefficients( normalVelocityFunctionCoefficients );

}


void HodographicShaping::satisfyAxialBoundaryConditions( Eigen::VectorXd freeCoefficients ){

    // Vector containing boundary conditions on axial distance and initial and final axial velocity.
    Eigen::Vector3d vectorBoundaryConditionsAxial;
    vectorBoundaryConditionsAxial( 0 ) = boundaryConditionsAxial_[ 1 ] - boundaryConditionsAxial_[ 0 ];
    vectorBoundaryConditionsAxial( 1 ) = boundaryConditionsAxial_[ 2 ];
    vectorBoundaryConditionsAxial( 2 ) = boundaryConditionsAxial_[ 3 ];

    // Subtract boundary values of free components of velocity function from corresponding boundary conditions.
    for ( int i = 0 ; i < numberOfFreeCoefficientsAxial ; i++ )
    {
        vectorBoundaryConditionsAxial[ 0 ] -= freeCoefficients( numberOfFreeCoefficientsRadial + numberOfFreeCoefficientsNormal + i )
                * ( axialVelocityFunction_.getComponentFunctionIntegralCurrentTime( i + 3 , finalTime_ )
                    - axialVelocityFunction_.getComponentFunctionIntegralCurrentTime( i + 3, 0.0 ) );
        vectorBoundaryConditionsAxial[ 1 ] -= freeCoefficients( numberOfFreeCoefficientsRadial + numberOfFreeCoefficientsNormal + i )
                * axialVelocityFunction_.getComponentFunctionCurrentValue( i + 3, 0.0 );
        vectorBoundaryConditionsAxial[ 2 ] -= freeCoefficients( numberOfFreeCoefficientsRadial + numberOfFreeCoefficientsNormal + i )
                * axialVelocityFunction_.getComponentFunctionCurrentValue( i + 3, finalTime_ );
    }

    // Compute fixed coefficients.
    Eigen::Vector3d fixedCoefficientsAxial = inverseMatrixBoundaryValuesAxial * vectorBoundaryConditionsAxial;

    // Create vector containing all axial velocity function coefficients.
    Eigen::VectorXd axialVelocityFunctionCoefficients = Eigen::VectorXd::Zero( fixedCoefficientsAxial.size() + numberOfFreeCoefficientsAxial );

    // Check whether the axial velocity function has free coefficients.
    if ( numberOfFreeCoefficientsAxial == 0 )
    {
        axialVelocityFunctionCoefficients << fixedCoefficientsAxial;
    }
    else
    {
        axialVelocityFunctionCoefficients << fixedCoefficientsAxial,
                freeCoefficients.segment( numberOfFreeCoefficientsRadial + numberOfFreeCoefficientsNormal, numberOfFreeCoefficientsAxial);
    }
    // Set the coefficients of the axial velocity function.
    axialVelocityFunction_.resetCompositeFunctionCoefficients( axialVelocityFunctionCoefficients );

}

void HodographicShaping::satisfyNormalBoundaryConditionsWithFinalPolarAngle( Eigen::VectorXd freeCoefficients ){

    // Compute coefficient of the third component of the composite function, from the required value of the final polar angle.
    double C3 = computeThirdFixedCoefficientAxialVelocityFromFinalPolarAngle( freeCoefficients );

    // Vector containing boundary conditions on initial and final normal velocity.
    Eigen::Vector2d vectorBoundaryConditionsNormal;
    vectorBoundaryConditionsNormal( 0 ) = boundaryConditionsNormal_[ 0 ];
    vectorBoundaryConditionsNormal( 1 ) = boundaryConditionsNormal_[ 1 ];

    // Subtract boundary values of velocity function component used to solve for final polar angle from corresponding boundary conditions.
    vectorBoundaryConditionsNormal[ 0 ] -= C3 * normalVelocityFunction_.getComponentFunctionCurrentValue( 2, 0.0 );
    vectorBoundaryConditionsNormal[ 1 ] -= C3 * normalVelocityFunction_.getComponentFunctionCurrentValue( 2, finalTime_ );

    // Subtract boundary values of free velocity function components from corresponding boundary conditions.
    for ( int i = 0 ; i < numberOfFreeCoefficientsNormal ; i++ )
    {
        vectorBoundaryConditionsNormal[ 0 ] -= freeCoefficients[ numberOfFreeCoefficientsRadial + i ]
                * normalVelocityFunction_.getComponentFunctionCurrentValue( i + 3, 0.0 );
        vectorBoundaryConditionsNormal[ 1 ] -= freeCoefficients[ numberOfFreeCoefficientsRadial + i ]
                * normalVelocityFunction_.getComponentFunctionCurrentValue( i + 3, finalTime_ );
    }

    // Compute fixed coefficients by satisfying the boundary conditions.
    Eigen::Vector2d fixedCoefficientsNormal = inverseMatrixBoundaryValuesNormal * vectorBoundaryConditionsNormal;

    // Create vector containing all normal velocity function coefficients.
    Eigen::VectorXd normalVelocityFunctionCoefficients =
            Eigen::VectorXd::Zero( fixedCoefficientsNormal.size() + 1 + numberOfFreeCoefficientsNormal );

    // Check whether the normal velocity function has free coefficients.
    if ( numberOfFreeCoefficientsNormal == 0 )
    {
        normalVelocityFunctionCoefficients << fixedCoefficientsNormal, C3;
    }
    else
    {
        normalVelocityFunctionCoefficients << fixedCoefficientsNormal, C3,
                                        freeCoefficients.segment( numberOfFreeCoefficientsRadial, numberOfFreeCoefficientsNormal);
    }

    // Set the coefficients of the normal velocity function.
    normalVelocityFunction_.resetCompositeFunctionCoefficients( normalVelocityFunctionCoefficients );

}

double HodographicShaping::computeThrustAccelerationCurrentTime( const double currentTime ){

    // Computation of radial distance from the Sun by integrating radial velocity.
    double radialDistance = radialVelocityFunction_.evaluateCompositeFunctionIntegralCurrentTime( currentTime )
            - radialVelocityFunction_.evaluateCompositeFunctionIntegralCurrentTime( 0.0 ) + boundaryConditionsRadial_[0];

    // Computation of axial distance from the Sun by integrating axial velocity.
    double axialDistance = axialVelocityFunction_.evaluateCompositeFunctionIntegralCurrentTime( currentTime )
            - axialVelocityFunction_.evaluateCompositeFunctionIntegralCurrentTime( 0.0 ) + boundaryConditionsAxial_[0];

    double distanceFromSun = sqrt( pow( radialDistance, 2.0 ) + pow( axialDistance, 2.0 ) );

    // Computation of normal velocity.
    double normalVelocity = normalVelocityFunction_.evaluateCompositeFunctionCurrentTime( currentTime );

    // Computation of angular velocity: dTheta/dt = V_theta/r
    double angularVelocity = normalVelocity / radialDistance;

    //Computation of radial thrust acceleration using equations of motion:
    double radialThrustAcceleration =
            radialVelocityFunction_.evaluateCompositeFunctionDerivativeCurrentTime( currentTime ) - angularVelocity * normalVelocity
            + celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER / pow( distanceFromSun, 3.0 ) * radialDistance;

    // Computation of normal thrust acceleration using equations of motion:
    double normalThrustAcceleration = normalVelocityFunction_.evaluateCompositeFunctionDerivativeCurrentTime( currentTime )
            + angularVelocity * radialVelocityFunction_.evaluateCompositeFunctionCurrentTime( currentTime );

    //Computation of axial thrust acceleration using equations of motion:
    double axialThrustAcceleration = axialVelocityFunction_.evaluateCompositeFunctionDerivativeCurrentTime( currentTime )
            + celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER / pow( distanceFromSun, 3.0 ) * axialDistance;

    // Return total thrust acceleration: f_tot = sqrt( f_r^2 + f_theta^2 + f_z^2 )
    return ( Eigen::Vector3d() << radialThrustAcceleration, normalThrustAcceleration, axialThrustAcceleration ).finished().norm();
}


//! Compute cartesian acceleration.
Eigen::Vector3d HodographicShaping::computeCartesianAcceleration( double currentTime ){

    Eigen::Vector3d cylindricalAcceleration = computeThrustAccelerationComponents( currentTime );
    Eigen::Vector3d cylindricalState = computeCurrentCylindricalState( currentTime ).segment(0,3);
    Eigen::Vector3d cartesianState = computeCurrentCartesianState( currentTime ).segment(0,3);

    Eigen::Vector3d cartesianAcceleration;

    cartesianAcceleration[ 0 ] = ( 1.0 / ( cartesianState[ 0 ] + ( std::pow( cartesianState[ 1 ], 2 ) / cartesianState[ 0 ] ) ) )
            * ( cylindricalState[ 0 ] * cylindricalAcceleration[ 0 ]
            - ( cartesianState[ 1 ] / cartesianState[ 0 ] ) * std::pow( cylindricalState[ 0 ], 1 ) * cylindricalAcceleration[ 1 ] );

    cartesianAcceleration[ 1 ] = ( std::pow( cylindricalState[ 0 ], 1 ) * cylindricalAcceleration[ 1 ]
            + cartesianState[ 1 ] * cartesianAcceleration[ 0 ] ) / cartesianState[ 0 ];

    cartesianAcceleration[ 2 ] = cylindricalAcceleration[ 2 ];

    return cartesianAcceleration;

}

//! Convert cartesian acceleration into cylindrical one.
Eigen::Vector3d HodographicShaping::convertCartesianToCylindricalAcceleration( Eigen::Vector3d cartesianAcceleration, Eigen::Vector3d cartesianState )
{
    Eigen::Vector3d cylindricalAcceleration;
    cylindricalAcceleration[ 0 ] = ( cartesianState[ 0 ] * cartesianAcceleration[ 0 ] + cartesianState[ 1 ] * cartesianAcceleration[ 1 ] )
            / std::sqrt( std::pow( cartesianState[ 0 ], 2 ) + std::pow( cartesianState[ 1 ], 2 ) );
    cylindricalAcceleration[ 1 ] = ( cartesianState[ 0 ] * cartesianAcceleration[ 1 ] - cartesianState[ 1 ] * cartesianAcceleration[ 0 ] )
            / std::pow( std::sqrt( std::pow( cartesianState[ 0 ], 2 ) + std::pow( cartesianState[ 1 ], 2 ) ), 1 );
    cylindricalAcceleration[ 2 ] = cartesianAcceleration[ 2 ];

    return cylindricalAcceleration;


}

//! Compute magnitude cartesian acceleration.
double  HodographicShaping::computeMagnitudeCartesianAcceleration( double currentTime )
{
    return computeCartesianAcceleration( currentTime ).norm();
}

//! Compute direction cartesian acceleration.
Eigen::Vector3d HodographicShaping::computeDirectionCartesianAcceleration( double currentTime )
{
    return computeCartesianAcceleration( currentTime ).normalized();
}


//! Compute DeltaV.
double HodographicShaping::computeDeltaV( )
{

    // Define the derivative of the deltaV, ie thrust acceleration function, as a function of time.
    std::function< Eigen::Vector1d( const double, const Eigen::Vector1d& ) > derivativeFunctionDeltaV = [ = ]
            ( const double currentTime, const Eigen::Vector1d& currentState ){

        Eigen::Vector1d thrustAcceleration;
        thrustAcceleration[ 0 ] = computeThrustAccelerationCurrentTime( currentTime );
        return thrustAcceleration;

    } ;

    // Define the derivative of the deltaV, ie thrust acceleration function, as a function of time.
    std::function< double( const double ) > derivativeFunctionDeltaVtest = [ = ] ( const double currentTime ){

        return computeThrustAccelerationCurrentTime( currentTime );

    } ;

    // Define numerical quadrature.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeFunctionDeltaVtest, quadratureSettings_, finalTime_ );

    return quadrature->getQuadrature( );
}


//! Compute angular velocity.
double HodographicShaping::computeAngularVelocityCurrentTime( const double currentTime )
{
    // Return angular velocity.
    return normalVelocityFunction_.evaluateCompositeFunctionCurrentTime( currentTime )
           / ( radialVelocityFunction_.evaluateCompositeFunctionIntegralCurrentTime( currentTime )
               - radialVelocityFunction_.evaluateCompositeFunctionIntegralCurrentTime( 0.0 ) + boundaryConditionsRadial_[0] );
}

//! Compute final polar angle.
double HodographicShaping::computeFinalPolarAngle( )
{

    // Define the derivative of the polar angle, ie angular velocity function, as a function of time.
    std::function< Eigen::Vector1d( const double, const Eigen::Vector1d& ) > derivativeFunctionPolarAngle = [ = ]
            ( const double currentTime, const Eigen::Vector1d& currentState ){

        Eigen::Vector1d angularVelocity;
        angularVelocity[ 0 ] = computeAngularVelocityCurrentTime( currentTime );
        return angularVelocity;

    } ;


    // Define the derivative of the polar angle, ie angular velocity function, as a function of time.
    std::function< double( const double ) > derivativeFunctionPolarAngleTest = [ = ]
            ( const double currentTime ){

        return computeAngularVelocityCurrentTime( currentTime );

    } ;

    // Define numerical quadrature.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeFunctionPolarAngleTest, quadratureSettings_, finalTime_ );

    return quadrature->getQuadrature() + coordinate_conversions::convertCartesianToCylindricalState( initialState_ )[ 1 ];
}

//! Compute polar angle.
double HodographicShaping::computePolarAngle( double currentTime )
{

    // Define the derivative of the polar angle, ie angular velocity function, as a function of time.
    std::function< Eigen::Vector1d( const double, const Eigen::Vector1d& ) > derivativeFunctionPolarAngle = [ = ]
            ( const double currentTime, const Eigen::Vector1d& currentState ){

        Eigen::Vector1d angularVelocity;
        angularVelocity[ 0 ] = computeAngularVelocityCurrentTime( currentTime );
        return angularVelocity;

    } ;

    // Define the derivative of the polar angle, ie angular velocity function, as a function of time.
    std::function< double( const double ) > derivativeFunctionPolarAngleTest = [ = ] ( const double currentTime ){

        return computeAngularVelocityCurrentTime( currentTime );

    } ;

    // Define numerical quadrature.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeFunctionPolarAngleTest, quadratureSettings_, currentTime );

    return quadrature->getQuadrature( ) + coordinate_conversions::convertCartesianToCylindricalState( initialState_ )[ 1 ];

}


//! Compute thrust acceleration components.
Eigen::Vector3d HodographicShaping::computeThrustAccelerationComponents( double currentTime )
{
    double radialDistance, axialDistance, distanceFromSun, normalVelocity, angularVelocity;
    Eigen::Vector3d thrustAccelerationComponents;

    // Computation of radial distance from the central body by integrating radial velocity.
    radialDistance = radialVelocityFunction_.evaluateCompositeFunctionIntegralCurrentTime( currentTime )
            - radialVelocityFunction_.evaluateCompositeFunctionIntegralCurrentTime( 0.0 ) + boundaryConditionsRadial_[0];

    // Computation of axial distance from the central body by integrating axial velocity.
    axialDistance = axialVelocityFunction_.evaluateCompositeFunctionIntegralCurrentTime( currentTime )
            - axialVelocityFunction_.evaluateCompositeFunctionIntegralCurrentTime( 0.0 ) + boundaryConditionsAxial_[0];

    // Computation of distance from the central body.
    distanceFromSun = sqrt( pow( radialDistance, 2.0 ) + pow( axialDistance, 2.0 ) );

    // Computation of normal velocity.
    normalVelocity = normalVelocityFunction_.evaluateCompositeFunctionCurrentTime( currentTime );

    // Computation of angular velocity: dTheta/dt = V_theta/r
    angularVelocity = normalVelocity / radialDistance;

    //Computation of radial thrust acceleration using equations of motion:
    thrustAccelerationComponents[ 0 ] = radialVelocityFunction_.evaluateCompositeFunctionDerivativeCurrentTime( currentTime )
            - angularVelocity * normalVelocity + celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER
            / std::pow( distanceFromSun, 3.0 ) * radialDistance;

    // Computation of normal thrust acceleration using equations of motion:
    thrustAccelerationComponents[ 1 ] = normalVelocityFunction_.evaluateCompositeFunctionDerivativeCurrentTime( currentTime )
            + angularVelocity * radialVelocityFunction_.evaluateCompositeFunctionCurrentTime( currentTime );

    //Computation of axial thrust acceleration using equations of motion:
    thrustAccelerationComponents[ 2 ] = axialVelocityFunction_.evaluateCompositeFunctionDerivativeCurrentTime( currentTime )
            + celestial_body_constants::SUN_GRAVITATIONAL_PARAMETER / std::pow( distanceFromSun, 3.0 ) * axialDistance;

    // Return total thrust acceleration.
    return thrustAccelerationComponents;
}

//! Compute thrust acceleration direction.
Eigen::Vector3d HodographicShaping::computeThrustAccelerationDirection( double currentTime )
{

    Eigen::Vector3d thrustAccelerationDirection;
    Eigen::Vector3d thrustAcceleration = computeThrustAccelerationComponents( currentTime );

    thrustAccelerationDirection = ( Eigen::Vector3d() << thrustAcceleration[ 0 ], thrustAcceleration[ 1 ],
            thrustAcceleration[ 2 ] ).finished();

    // Return total thrust acceleration.
    return coordinate_conversions::convertCylindricalToCartesian( thrustAccelerationDirection ).normalized();
}


//! Compute velocity components.
std::vector< double > HodographicShaping::computeVelocityComponents( double currentTime )
{
    std::vector< double > velocityComponents;

    // Compute radial velocity.
    velocityComponents.push_back( radialVelocityFunction_.evaluateCompositeFunctionCurrentTime( currentTime ) );

    // Compute normal velocity.
    velocityComponents.push_back( normalVelocityFunction_.evaluateCompositeFunctionCurrentTime( currentTime ) );

    // Compute axial velocity.
    velocityComponents.push_back( axialVelocityFunction_.evaluateCompositeFunctionCurrentTime( currentTime ) );

    return velocityComponents;
}



//! Compute radial distance.
double HodographicShaping::computeRadialDistanceCurrentTime( const double currentTime )
{
    // Computation of radial distance from the Sun by integrating radial velocity.
    return radialVelocityFunction_.evaluateCompositeFunctionIntegralCurrentTime( currentTime )
            - radialVelocityFunction_.evaluateCompositeFunctionIntegralCurrentTime( 0.0 ) + boundaryConditionsRadial_[0];
}

double HodographicShaping::computeAxialDistanceCurrentTime( const double currentTime )
{
    return axialVelocityFunction_.evaluateCompositeFunctionIntegralCurrentTime( currentTime )
            - axialVelocityFunction_.evaluateCompositeFunctionIntegralCurrentTime( 0.0 ) + boundaryConditionsAxial_[0];
}

//! Compute current cylindrical state.
Eigen::Vector6d HodographicShaping::computeCurrentCylindricalState( const double currentTime ){

    std::vector< double > velocityCylindricalCoordinates= computeVelocityComponents( currentTime );

    Eigen::Vector6d cylindricalState = ( Eigen::Vector6d() <<
                                         computeRadialDistanceCurrentTime( currentTime ),
                                         computePolarAngle( currentTime ),
                                         computeAxialDistanceCurrentTime( currentTime ),
                                         velocityCylindricalCoordinates[0],
                                         velocityCylindricalCoordinates[1],
                                         velocityCylindricalCoordinates[2] ).finished();

    return cylindricalState;
}


//! Compute current cartesian state.
Eigen::Vector6d HodographicShaping::computeCurrentCartesianState( const double currentTime ){

    return coordinate_conversions::convertCylindricalToCartesianState( computeCurrentCylindricalState( currentTime ) );

}



double HodographicShaping::computeThirdFixedCoefficientAxialVelocityFromFinalPolarAngle( Eigen::VectorXd freeCoefficients ){

    // Compute the third fixed coefficient of the normal velocity composite function, so that the condition on the final
    // polar angle is fulfilled.
    // The calculation is based on Equation (16) in Gondelach D., and Noomen R.
    // "Hodographic-shaping method for low-thrust interplanetary trajectory design." Journal of Spacecraft and Rockets 52.3 (2015): 728-738.

    // Vector containing boundary conditions on initial and final normal velocity.
    Eigen::Vector2d vectorBoundaryConditionsNormal;
    vectorBoundaryConditionsNormal( 0 ) = boundaryConditionsNormal_[ 0 ];
    vectorBoundaryConditionsNormal( 1 ) = boundaryConditionsNormal_[ 1 ];

    // Subtract boundary values of free velocity function components from corresponding boundary conditions.
    for ( int i = 0 ; i < numberOfFreeCoefficientsNormal; i++ )
    {
        vectorBoundaryConditionsNormal[ 0 ] -= freeCoefficients[ numberOfFreeCoefficientsRadial + i ]
                * normalVelocityFunction_.getComponentFunctionCurrentValue( i + 3, 0.0 );
        vectorBoundaryConditionsNormal[ 1 ] -= freeCoefficients[ numberOfFreeCoefficientsRadial + i ]
                * normalVelocityFunction_.getComponentFunctionCurrentValue( i + 3, finalTime_ );
    }

    Eigen::Vector2d matrixL;
    matrixL = inverseMatrixBoundaryValuesNormal * vectorBoundaryConditionsNormal;

    Eigen::Vector2d initialAndFinalValuesThirdComponentFunction( - normalVelocityFunction_.getComponentFunctionCurrentValue( 2, 0.0 ),
                                  - normalVelocityFunction_.getComponentFunctionCurrentValue( 2, finalTime_ ) );
    Eigen::Vector2d matrixK;
    matrixK = inverseMatrixBoundaryValuesNormal * initialAndFinalValuesThirdComponentFunction;


    // Define angular velocity due to the third component of the composite function only.
    std::function< double( const double ) > derivativePolarAngleDueToThirdComponent = [ = ] ( const double currentTime ){

        double angularVelocityDueToThirdComponent = ( matrixK( 0 ) * normalVelocityFunction_.getComponentFunctionCurrentValue( 0, currentTime )
                + matrixK( 1 ) * normalVelocityFunction_.getComponentFunctionCurrentValue( 1, currentTime )
                + normalVelocityFunction_.getComponentFunctionCurrentValue( 2, currentTime ) ) / computeRadialDistanceCurrentTime( currentTime );

        return angularVelocityDueToThirdComponent;

    };


    // Define the angular velocity due to all the other components of the composite function, once combined.
    std::function< double( const double ) > derivativePolarAngleDueToOtherComponents = [ = ] ( const double currentTime ){

        double angularVelocityDueToFreeCoefficients = 0.0;
        for( int j = 0; j < numberOfFreeCoefficientsNormal ; j++ )
        {
            angularVelocityDueToFreeCoefficients += freeCoefficients( numberOfFreeCoefficientsRadial + j )
                             * normalVelocityFunction_.getComponentFunctionCurrentValue( j + 3 , currentTime );
        }

        double angularVelocityDueToOtherComponents = ( matrixL( 0 ) * normalVelocityFunction_.getComponentFunctionCurrentValue( 0 , currentTime )
                                                     + matrixL( 1 ) * normalVelocityFunction_.getComponentFunctionCurrentValue( 1 , currentTime )
                                                     + angularVelocityDueToFreeCoefficients )
                / computeRadialDistanceCurrentTime( currentTime );

        return angularVelocityDueToOtherComponents;

    };

    // Define numerical quadratures.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > integratorPolarAngleDueToThirdComponentTest =
            numerical_quadrature::createQuadrature( derivativePolarAngleDueToThirdComponent, quadratureSettings_, finalTime_ );

    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > integratorPolarAngleDueToOtherComponentsTest =
            numerical_quadrature::createQuadrature( derivativePolarAngleDueToOtherComponents, quadratureSettings_, finalTime_ );

    return ( boundaryConditionsNormal_[ 2 ] - integratorPolarAngleDueToOtherComponentsTest->getQuadrature() )
            / integratorPolarAngleDueToThirdComponentTest->getQuadrature();;

}


//! Get low-thrust acceleration model from shaping method.
std::shared_ptr< propulsion::ThrustAcceleration > HodographicShaping::getLowThrustAccelerationModel(
        simulation_setup::NamedBodyMap& bodyMap,
        const std::string& bodyToPropagate,
        std::function< double( const double ) > specificImpulseFunction )
{
    // Define thrust direction settings from the direction of thrust acceleration retrieved from the shaping method.
    std::shared_ptr< simulation_setup::CustomThrustDirectionSettings > thrustDirectionSettings =
            std::make_shared< simulation_setup::CustomThrustDirectionSettings >(
                std::bind( &HodographicShaping::computeDirectionCartesianAcceleration, this, std::placeholders::_1 ) );

    std::shared_ptr< simulation_setup::Body > vehicle = bodyMap[ bodyToPropagate ];

    // Define thrust magnitude function from the shaped trajectory.
    std::function< double( const double ) > thrustMagnitudeFunction = [ = ]( const double currentTime )
    {
        // Compute current acceleration.
        double currentAcceleration = computeMagnitudeCartesianAcceleration( currentTime );

        // Compute current mass of the vehicle.
        double currentMass = vehicle->getBodyMass();

        // Compute and return magnitude of the low-thrust force.
        return currentAcceleration * currentMass;
    };

    // Define thrust magnitude settings from thrust magnitude function.
    std::shared_ptr< simulation_setup::FromFunctionThrustMagnitudeSettings > thrustMagnitudeSettings
            = std::make_shared< simulation_setup::FromFunctionThrustMagnitudeSettings >(
                thrustMagnitudeFunction, specificImpulseFunction );

    // Define thrust acceleration settings.
    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > thrustAccelerationSettings =
            std::make_shared< simulation_setup::ThrustAccelerationSettings >(
                thrustDirectionSettings, thrustMagnitudeSettings );

    // Create low thrust acceleration model.
    std::shared_ptr< propulsion::ThrustAcceleration > lowThrustAccelerationModel = createThrustAcceleratioModel( thrustAccelerationSettings, bodyMap, bodyToPropagate );

    return lowThrustAccelerationModel;
}


void HodographicShaping::computeShapingTrajectoryAndFullPropagation(
        simulation_setup::NamedBodyMap& bodyMap,
        basic_astrodynamics::AccelerationMap& accelerationMap,
        const std::string& centralBody,
        const std::string& bodyToPropagate,
        std::function< double( const double ) > specificImpulseFunction,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
        std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettings,
        std::map< double, Eigen::VectorXd >& fullPropagationResults,
        std::map< double, Eigen::VectorXd >& shapingMethodResults,
        std::map< double, Eigen::VectorXd >& dependentVariablesHistory,
        propagators::TranslationalPropagatorType propagatorType,
        const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave ){

    fullPropagationResults.clear();
    shapingMethodResults.clear();
    dependentVariablesHistory.clear();

    // Compute halved time of flight.
    double halvedTimeOfFlight = finalTime_ / 2.0;


    // Compute state at half of the time of flight.
    Eigen::Vector6d initialStateAtHalvedTimeOfFlight = computeCurrentCartesianState( halvedTimeOfFlight );


    // Create low thrust acceleration model.
    std::shared_ptr< propulsion::ThrustAcceleration > lowThrustAccelerationModel = getLowThrustAccelerationModel( bodyMap, bodyToPropagate, specificImpulseFunction );

    accelerationMap[ bodyToPropagate ][ bodyToPropagate ].push_back( lowThrustAccelerationModel );


    std::vector< std::string > centralBodies;
    centralBodies.push_back( centralBody );

    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( bodyToPropagate );

    // Define forward propagator settings variables.
    integratorSettings->initialTime_ = initialTime_ + halvedTimeOfFlight;

    // Define propagation settings
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsForwardPropagation;
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > propagatorSettingsBackwardPropagation;

    // Define forward propagation settings.
    propagatorSettingsForwardPropagation = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( centralBodies, accelerationMap, bodiesToPropagate, initialStateAtHalvedTimeOfFlight, terminationSettings.second,
                          propagatorType, dependentVariablesToSave );

    // Define backward propagation settings.
    propagatorSettingsBackwardPropagation = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( centralBodies, accelerationMap, bodiesToPropagate, initialStateAtHalvedTimeOfFlight, terminationSettings.first,
                          propagatorType, dependentVariablesToSave );

    // Perform forward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationForwards( bodyMap, integratorSettings, propagatorSettingsForwardPropagation );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemForwardPropagation = dynamicsSimulatorIntegrationForwards.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistoryForwardPropagation = dynamicsSimulatorIntegrationForwards.getDependentVariableHistory( );

    // Compute and save full propagation and shaping method results along the forward propagation direction.
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemForwardPropagation.begin( );
         itr != stateHistoryFullProblemForwardPropagation.end( ); itr++ )
    {
        shapingMethodResults[ itr->first ] = computeCurrentCartesianState( itr->first );
        fullPropagationResults[ itr->first ] = itr->second;
        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryForwardPropagation[ itr->first ];
    }


    // Define backward propagator settings variables.
    integratorSettings->initialTimeStep_ = - integratorSettings->initialTimeStep_;
    integratorSettings->initialTime_ = initialTime_ + halvedTimeOfFlight;

    // Perform the backward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationBackwards( bodyMap, integratorSettings, propagatorSettingsBackwardPropagation );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemBackwardPropagation = dynamicsSimulatorIntegrationBackwards.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistoryBackwardPropagation = dynamicsSimulatorIntegrationBackwards.getDependentVariableHistory( );

    // Compute and save full propagation and shaping method results along the backward propagation direction
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemBackwardPropagation.begin( );
         itr != stateHistoryFullProblemBackwardPropagation.end( ); itr++ )
    {
        shapingMethodResults[ itr->first ] = computeCurrentCartesianState( itr->first );
        fullPropagationResults[ itr->first ] = itr->second;
        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryBackwardPropagation[ itr->first ];
    }

    // Reset initial integrator settings.
    integratorSettings->initialTimeStep_ = -integratorSettings->initialTimeStep_;

}

} // namespace shape_based_methods
} // namespace tudat
