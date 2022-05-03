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


#include "tudat/astro/low_thrust/shape_based/hodographicShapingLeg.h"
#include "tudat/math/basic/coordinateConversions.h"

namespace tudat
{
namespace shape_based_methods
{


//! Constructor which sets radial, normal and axial velocity functions and boundary conditions.
HodographicShapingLeg::HodographicShapingLeg(
        const std::shared_ptr< ephemerides::Ephemeris > departureBodyEphemeris,
        const std::shared_ptr< ephemerides::Ephemeris > arrivalBodyEphemeris,
        const double centralBodyGravitationalParameter,
        const std::function< Eigen::Vector3d( ) > departureVelocityFunction,
        const std::function< Eigen::Vector3d( ) > arrivalVelocityFunction,
        const HodographicBasisFunctionList& radialVelocityFunctionComponents,
        const HodographicBasisFunctionList& normalVelocityFunctionComponents,
        const HodographicBasisFunctionList& axialVelocityFunctionComponents ) :
    mission_segments::TransferLeg( departureBodyEphemeris, arrivalBodyEphemeris, mission_segments::hodographic_low_thrust_leg ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    departureVelocityFunction_( departureVelocityFunction ),
    arrivalVelocityFunction_( arrivalVelocityFunction ),
    numberOfFreeRadialCoefficients_( radialVelocityFunctionComponents.size( ) - 3 ),
    numberOfFreeNormalCoefficients_( normalVelocityFunctionComponents.size( ) - 3 ),
    numberOfFreeAxialCoefficients_( axialVelocityFunctionComponents.size( ) - 3 )
{
    if( numberOfFreeRadialCoefficients_ < 0 || numberOfFreeNormalCoefficients_ < 0 || numberOfFreeAxialCoefficients_ < 0 )
    {
        throw std::runtime_error(
                    "Error, input of number of free coefficients for hodographic shaping is smaller than 0. Each direction must have at least three components. Total number of coefficients is: " +
                    std::to_string( radialVelocityFunctionComponents.size( ) ) + " for radial direction, " +
                    std::to_string( normalVelocityFunctionComponents.size( ) ) + " for normal direction, " +
                    std::to_string( axialVelocityFunctionComponents.size( ) ) + " for axial direction, " );

    }

    fullCoefficientsRadialVelocityFunction_ = Eigen::VectorXd::Zero(numberOfFreeRadialCoefficients_ + 3 );
    fullCoefficientsNormalVelocityFunction_ = Eigen::VectorXd::Zero(numberOfFreeNormalCoefficients_ + 3 );
    fullCoefficientsAxialVelocityFunction_ = Eigen::VectorXd::Zero(numberOfFreeAxialCoefficients_ + 3 );

    radialVelocityFunction_ = std::make_shared< shape_based_methods::CompositeFunctionHodographicShaping >(
                radialVelocityFunctionComponents, fullCoefficientsRadialVelocityFunction_ );
    normalVelocityFunction_ = std::make_shared< shape_based_methods::CompositeFunctionHodographicShaping >(
                normalVelocityFunctionComponents, fullCoefficientsNormalVelocityFunction_ );
    axialVelocityFunction_ = std::make_shared< shape_based_methods::CompositeFunctionHodographicShaping >(
                axialVelocityFunctionComponents, fullCoefficientsAxialVelocityFunction_ );

    // Define numerical quadrature settings, required to compute the current polar angle and final deltaV.
    quadratureSettings_ = std::make_shared< numerical_quadrature::GaussianQuadratureSettings< double > >( 0.0, 64 );

}

void HodographicShapingLeg::computeTransfer( )
{
    updateDepartureAndArrivalBodies( legParameters_( 0 ), legParameters_( 1 ) );

    // Update number of revolutions, after testing if value is valid
    if ( legParameters_(2) < 0 )
    {
        throw std::runtime_error( "Error when updating hodographic shaping object, number of revolutions should be equal to or larger than 0" );
    }
    else if ( std::floor( legParameters_(2) ) != legParameters_(2) )
    {
        throw std::runtime_error( "Error when updating hodographic shaping object, number of revolutions should be an integer" );
    }
    else
    {
        numberOfRevolutions_ = int(legParameters_(2));
    }

    // Check whether number of remaining parameters is valid
    if( legParameters_.rows( ) - 3 != numberOfFreeRadialCoefficients_ + numberOfFreeNormalCoefficients_ + numberOfFreeAxialCoefficients_ )
    {
        throw std::runtime_error( "Error when updating hodographic shaping object, number of inputs is inconsistent" );
    }

    updateFreeCoefficients( );
    satisfyBoundaryConditions( );
    legTotalDeltaV_ = computeDeltaV( );
}

void HodographicShapingLeg::updateFreeCoefficients( )
{
    fullCoefficientsRadialVelocityFunction_.segment( 0, 3 ).setZero( );
    fullCoefficientsRadialVelocityFunction_.segment(3, numberOfFreeRadialCoefficients_ ) = legParameters_.segment(
            3, numberOfFreeRadialCoefficients_ );
    radialVelocityFunction_->resetCompositeFunctionCoefficients( fullCoefficientsRadialVelocityFunction_ );

    fullCoefficientsNormalVelocityFunction_.segment( 0, 3 ).setZero( );
    fullCoefficientsNormalVelocityFunction_.segment(3, numberOfFreeNormalCoefficients_ ) = legParameters_.segment(
            3 + numberOfFreeRadialCoefficients_, numberOfFreeNormalCoefficients_ );
    normalVelocityFunction_->resetCompositeFunctionCoefficients( fullCoefficientsNormalVelocityFunction_ );

    fullCoefficientsAxialVelocityFunction_.segment( 0, 3 ).setZero( );
    fullCoefficientsAxialVelocityFunction_.segment(3, numberOfFreeAxialCoefficients_ ) = legParameters_.segment(
            3 + numberOfFreeRadialCoefficients_ + numberOfFreeNormalCoefficients_, numberOfFreeAxialCoefficients_ );
    axialVelocityFunction_->resetCompositeFunctionCoefficients( fullCoefficientsAxialVelocityFunction_ );

}

void HodographicShapingLeg::satisfyBoundaryConditions( )
{
    departureVelocity_ = departureVelocityFunction_( );
    arrivalVelocity_ = arrivalVelocityFunction_( );

    // Select initial and final cartesian state
    Eigen::Vector6d initialCartesianState;
    initialCartesianState.segment(0, 3) = departureBodyState_.segment(0, 3);
    initialCartesianState.segment(3, 3 ) = departureVelocity_;

    // Normalize the final state.
    Eigen::Vector6d finalCartesianState;
    finalCartesianState.segment(0, 3 ) = arrivalBodyState_.segment( 0, 3 );
    finalCartesianState.segment(3, 3 ) = arrivalVelocity_;

    // Compute initial and final cylindrical state which is to be met by the hodographic shaping transfer
    Eigen::Vector6d initialCylindricalState = coordinate_conversions::convertCartesianToCylindricalState( initialCartesianState );
    Eigen::Vector6d finalCylindricalState = coordinate_conversions::convertCartesianToCylindricalState( finalCartesianState );

    // Set boundary conditions in the radial direction.
    radialBoundaryConditions_.clear( );
    radialBoundaryConditions_.push_back( initialCylindricalState[ 0 ] );
    radialBoundaryConditions_.push_back( finalCylindricalState[ 0 ] );
    radialBoundaryConditions_.push_back( initialCylindricalState[ 3 ] );
    radialBoundaryConditions_.push_back( finalCylindricalState[ 3 ] );

    // Set boundary conditions in the normal direction.
    normalBoundaryConditions_.clear( );
    normalBoundaryConditions_.push_back( initialCylindricalState[ 4 ] );
    normalBoundaryConditions_.push_back( finalCylindricalState[ 4 ] );
    // Final value polar angle
    normalBoundaryConditions_.push_back(
                numberOfRevolutions_ * 2.0 * mathematical_constants::PI
                + basic_mathematics::computeModulo( ( finalCylindricalState[ 1 ] - initialCylindricalState[ 1 ] ),
            2.0 * mathematical_constants::PI ) );

    // Set boundary conditions in the axial direction.
    axialBoundaryConditions_.clear( );
    axialBoundaryConditions_.push_back( initialCylindricalState[ 2 ] );
    axialBoundaryConditions_.push_back( finalCylindricalState[ 2 ] );
    axialBoundaryConditions_.push_back( initialCylindricalState[ 5 ] );
    axialBoundaryConditions_.push_back( finalCylindricalState[ 5 ] );

    // Compute inverse of matrices containing boundary values.
    inverseMatrixRadialBoundaryValues_ = computeInverseMatrixRadialOrAxialBoundaries( radialVelocityFunction_ );
    inverseMatrixNormalBoundaryValues_ = computeInverseMatrixNormalBoundaries( normalVelocityFunction_ );
    inverseMatrixAxialBoundaryValues_ = computeInverseMatrixRadialOrAxialBoundaries(axialVelocityFunction_ );

    // Satisfy boundary conditions.
    satisfyRadialBoundaryConditions( fullCoefficientsRadialVelocityFunction_.segment(3, numberOfFreeRadialCoefficients_ ) );
    satisfyNormalBoundaryConditions( fullCoefficientsNormalVelocityFunction_.segment(3, numberOfFreeNormalCoefficients_ ) );
    satisfyAxialBoundaryConditions( fullCoefficientsAxialVelocityFunction_.segment(3, numberOfFreeAxialCoefficients_ ) );
}


//! Compute inverse of matrix filled with boundary conditions in radial or axial direction.
Eigen::Matrix3d HodographicShapingLeg::computeInverseMatrixRadialOrAxialBoundaries(
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
Eigen::Matrix2d HodographicShapingLeg::computeInverseMatrixNormalBoundaries(
        std::shared_ptr< CompositeFunctionHodographicShaping > velocityFunction )
{
    Eigen::Matrix2d matrixBoundaryValues, inverseMatrixBoundaryValues;

    matrixBoundaryValues << velocityFunction->getComponentFunctionCurrentValue( 0, 0.0 ),
            velocityFunction->getComponentFunctionCurrentValue( 1, 0.0 ),
            velocityFunction->getComponentFunctionCurrentValue( 0, timeOfFlight_ ),
            velocityFunction->getComponentFunctionCurrentValue( 1, timeOfFlight_ );

    // Compute inverse of boundary-value matrix.
    inverseMatrixBoundaryValues = matrixBoundaryValues.inverse();

    // Return inverse of boundary-value matrix.
    return inverseMatrixBoundaryValues;
}


void HodographicShapingLeg::satisfyRadialBoundaryConditions( const Eigen::VectorXd& freeCoefficients){

    // Vector containing boundary conditions on radial distance and initial and final radial velocity.
    Eigen::Vector3d vectorBoundaryConditionsRadial;
    vectorBoundaryConditionsRadial[ 0 ] = radialBoundaryConditions_[ 1 ] - radialBoundaryConditions_[ 0 ];
    vectorBoundaryConditionsRadial[ 1 ] = radialBoundaryConditions_[ 2 ];
    vectorBoundaryConditionsRadial[ 2 ] = radialBoundaryConditions_[ 3 ];

    // Subtract boundary values of free components of velocity function from corresponding boundary conditions.
    for (int i = 0 ; i < numberOfFreeRadialCoefficients_ ; i++ )
    {
        vectorBoundaryConditionsRadial[ 0 ] -=
                freeCoefficients( i ) * ( radialVelocityFunction_->getComponentFunctionIntegralCurrentValue(i + 3, timeOfFlight_)
                - radialVelocityFunction_->getComponentFunctionIntegralCurrentValue(i + 3, 0.0) );
        vectorBoundaryConditionsRadial[ 1 ] -= freeCoefficients( i ) *
                radialVelocityFunction_->getComponentFunctionCurrentValue( i + 3, 0.0 );
        vectorBoundaryConditionsRadial[ 2 ] -= freeCoefficients( i ) *
                radialVelocityFunction_->getComponentFunctionCurrentValue( i + 3, timeOfFlight_ );
    }

    // Compute fixed coefficients.
    Eigen::Vector3d fixedCoefficientsRadial = inverseMatrixRadialBoundaryValues_ * vectorBoundaryConditionsRadial;

    // Create vector containing all radial velocity function coefficients.
    Eigen::VectorXd radialVelocityFunctionCoefficients =
            Eigen::VectorXd::Zero(fixedCoefficientsRadial.size() + numberOfFreeRadialCoefficients_ );

    // Check whether the radial velocity function has free coefficients.
    if (numberOfFreeRadialCoefficients_ == 0 )
    {
        radialVelocityFunctionCoefficients << fixedCoefficientsRadial;
    }
    else
    {
        radialVelocityFunctionCoefficients << fixedCoefficientsRadial, freeCoefficients;
    }

    // Set the coefficients of the radial velocity function.
    radialVelocityFunction_->resetCompositeFunctionCoefficients( radialVelocityFunctionCoefficients );

}

void HodographicShapingLeg::satisfyAxialBoundaryConditions( const Eigen::VectorXd& freeCoefficients ){

    // Vector containing boundary conditions on axial distance and initial and final axial velocity.
    Eigen::Vector3d vectorBoundaryConditionsAxial;
    vectorBoundaryConditionsAxial[ 0 ] = axialBoundaryConditions_[ 1 ] - axialBoundaryConditions_[ 0 ];
    vectorBoundaryConditionsAxial[ 1 ] = axialBoundaryConditions_[ 2 ];
    vectorBoundaryConditionsAxial[ 2 ] = axialBoundaryConditions_[ 3 ];

    // Subtract boundary values of free components of velocity function from corresponding boundary conditions.
    for (int i = 0 ; i < numberOfFreeAxialCoefficients_ ; i++ )
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
    Eigen::Vector3d fixedCoefficientsAxial = inverseMatrixAxialBoundaryValues_ * vectorBoundaryConditionsAxial;

    // Create vector containing all axial velocity function coefficients.
    Eigen::VectorXd axialVelocityFunctionCoefficients =
            Eigen::VectorXd::Zero(fixedCoefficientsAxial.size() + numberOfFreeAxialCoefficients_ );

    // Check whether the axial velocity function has free coefficients.
    if (numberOfFreeAxialCoefficients_ == 0 )
    {
        axialVelocityFunctionCoefficients << fixedCoefficientsAxial;
    }
    else
    {
        axialVelocityFunctionCoefficients << fixedCoefficientsAxial, freeCoefficients;
    }
    // Set the coefficients of the axial velocity function.
    axialVelocityFunction_->resetCompositeFunctionCoefficients( axialVelocityFunctionCoefficients );

}

void HodographicShapingLeg::satisfyNormalBoundaryConditions( const Eigen::VectorXd& freeCoefficients ){

    // Compute coefficient of the third component of the composite function, from the required value of the final polar angle.
    double C3 = computeThirdFixedCoefficientAxialVelocity(freeCoefficients);

    // Vector containing boundary conditions on initial and final normal velocity.
    Eigen::Vector2d vectorBoundaryConditionsNormal;
    vectorBoundaryConditionsNormal( 0 ) = normalBoundaryConditions_[ 0 ];
    vectorBoundaryConditionsNormal( 1 ) = normalBoundaryConditions_[ 1 ];

    // Subtract boundary values of velocity function component used to solve for final polar angle from corresponding boundary conditions.
    vectorBoundaryConditionsNormal[ 0 ] -= C3 * normalVelocityFunction_->getComponentFunctionCurrentValue( 2, 0.0 );
    vectorBoundaryConditionsNormal[ 1 ] -= C3 * normalVelocityFunction_->getComponentFunctionCurrentValue( 2, timeOfFlight_ );

    // Subtract boundary values of free velocity function components from corresponding boundary conditions.
    for (int i = 0 ; i < numberOfFreeNormalCoefficients_ ; i++ )
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
            Eigen::VectorXd::Zero(fixedCoefficientsNormal.size() + 1 + numberOfFreeNormalCoefficients_ );

    // Check whether the normal velocity function has free coefficients.
    if (numberOfFreeNormalCoefficients_ == 0 )
    {
        normalVelocityFunctionCoefficients << fixedCoefficientsNormal, C3;
    }
    else
    {
        normalVelocityFunctionCoefficients << fixedCoefficientsNormal, C3, freeCoefficients;
    }

    // Set the coefficients of the normal velocity function.
    normalVelocityFunction_->resetCompositeFunctionCoefficients( normalVelocityFunctionCoefficients );

}

// Define angular velocity due to the third component of the composite function only.
double HodographicShapingLeg::computeDerivativePolarAngleDueToThirdComponent(
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
double HodographicShapingLeg::computeDerivativePolarAngleDueToOtherComponents(
        const double timeSinceDeparture, const Eigen::Vector2d& matrixL, const Eigen::VectorXd& freeCoefficients )
{

    // Compute angular velocity associated with free coefficients
    double angularVelocityDueToFreeCoefficients = 0.0;
    for(int j = 0 ; j < numberOfFreeNormalCoefficients_ ; j++ )
    {
        angularVelocityDueToFreeCoefficients += freeCoefficients( j )
                * normalVelocityFunction_->getComponentFunctionCurrentValue( j + 3 , timeSinceDeparture );
    }
    angularVelocityDueToFreeCoefficients = angularVelocityDueToFreeCoefficients / computeCurrentRadialDistance( timeSinceDeparture );

    // Compute angular velocity associated with fixed coefficients (first and second coefficient)
    double angularVelocityDueToFirstAndSecondCoefficients =
            ( matrixL( 0 ) * normalVelocityFunction_->getComponentFunctionCurrentValue( 0 , timeSinceDeparture )
            + matrixL( 1 ) * normalVelocityFunction_->getComponentFunctionCurrentValue( 1 , timeSinceDeparture ) ) /
            computeCurrentRadialDistance( timeSinceDeparture );

    return angularVelocityDueToFreeCoefficients + angularVelocityDueToFirstAndSecondCoefficients;

}

double HodographicShapingLeg::computeThirdFixedCoefficientAxialVelocity ( const Eigen::VectorXd& freeCoefficients ){

    // Compute the third fixed coefficient of the normal velocity composite function, so that the condition on the final
    // polar angle is fulfilled.
    // The calculation is based on Equation (16) in Gondelach D., and Noomen R.
    // "Hodographic-shaping method for low-thrust interplanetary trajectory design." Journal of Spacecraft and Rockets 52.3 (2015): 728-738.

    // Vector containing boundary conditions on initial and final normal velocity.
    Eigen::Vector2d vectorBoundaryConditionsNormal;
    vectorBoundaryConditionsNormal[ 0 ] = normalBoundaryConditions_[ 0 ];
    vectorBoundaryConditionsNormal[ 1 ] = normalBoundaryConditions_[ 1 ];

    // Subtract boundary values of free velocity function components from corresponding boundary conditions.
    for (int i = 0 ; i < numberOfFreeNormalCoefficients_; i++ )
    {
        vectorBoundaryConditionsNormal[ 0 ] -= freeCoefficients[ i ]
                * normalVelocityFunction_->getComponentFunctionCurrentValue( i + 3, 0.0 );
        vectorBoundaryConditionsNormal[ 1 ] -= freeCoefficients[ i ]
                * normalVelocityFunction_->getComponentFunctionCurrentValue( i + 3, timeOfFlight_ );
    }

    // Define matrix L, according to equation (15) of Gondelach and Noomen
    Eigen::Vector2d matrixL;
    matrixL = inverseMatrixNormalBoundaryValues_ * vectorBoundaryConditionsNormal;

    Eigen::Vector2d initialAndFinalValuesThirdComponentFunction(
                - normalVelocityFunction_->getComponentFunctionCurrentValue( 2, 0.0 ),
                - normalVelocityFunction_->getComponentFunctionCurrentValue( 2, timeOfFlight_ ) );

    // Define matrix K, according to equation (15) of Gondelach and Noomen
    Eigen::Vector2d matrixK;
    matrixK = inverseMatrixNormalBoundaryValues_ * initialAndFinalValuesThirdComponentFunction;


    // Define angular velocity due to the third component of the composite function only.
    std::function< double( const double ) > derivativePolarAngleDueToThirdComponent =
            std::bind( &HodographicShapingLeg::computeDerivativePolarAngleDueToThirdComponent, this, std::placeholders::_1, matrixK );

    // Define the angular velocity due to all the other components of the composite function, once combined.
    std::function< double( const double ) > derivativePolarAngleDueToOtherComponents  =
            std::bind( &HodographicShapingLeg::computeDerivativePolarAngleDueToOtherComponents, this, std::placeholders::_1, matrixL, freeCoefficients );

    // Define numerical quadratures.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > integratorPolarAngleDueToThirdComponent =
            numerical_quadrature::createQuadrature( derivativePolarAngleDueToThirdComponent, quadratureSettings_, timeOfFlight_ );

    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > integratorPolarAngleDueToOtherComponents =
            numerical_quadrature::createQuadrature( derivativePolarAngleDueToOtherComponents, quadratureSettings_, timeOfFlight_ );

    return ( normalBoundaryConditions_[ 2 ] - integratorPolarAngleDueToOtherComponents->getQuadrature() )
           / integratorPolarAngleDueToThirdComponent->getQuadrature();

}


//! Compute velocity components.
Eigen::Vector3d HodographicShapingLeg::computeVelocityVectorInCylindricalCoordinates( const double timeSinceDeparture )
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
double HodographicShapingLeg::evaluateDerivativePolarAngleWrtTime( const double timeSinceDeparture )
{
    // Return derivative of the polar angle w.r.t. time, i.e. angular velocity.
    return normalVelocityFunction_->evaluateCompositeFunctionCurrentValue(timeSinceDeparture)
            / computeCurrentRadialDistance( timeSinceDeparture );
}


//! Compute radial distance.
double HodographicShapingLeg::computeCurrentRadialDistance( const double timeSinceDeparture )
{
    // Compute radial distance from the central body.
    return radialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentValue(timeSinceDeparture)
           - radialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentValue(0.0) + radialBoundaryConditions_[0];
}

//! Compute current polar angle.
double HodographicShapingLeg::computeCurrentPolarAngle( const double timeSinceDeparture )
{
    // Define the derivative of the polar angle, ie angular velocity function, as a function of time.
    std::function< double( const double ) > derivativeFunctionPolarAngle =
            std::bind( &HodographicShapingLeg::evaluateDerivativePolarAngleWrtTime, this, std::placeholders::_1 );

    // Define numerical quadrature.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeFunctionPolarAngle, quadratureSettings_, timeSinceDeparture );

    double currentPolarAngle = quadrature->getQuadrature( ) + radialBoundaryConditions_[ 1 ];

    return currentPolarAngle;

}

double HodographicShapingLeg::computeCurrentAxialDistance( const double timeSinceDeparture )
{
    // Compute axial distance from the central body.
    return axialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentValue(timeSinceDeparture)
           - axialVelocityFunction_->evaluateCompositeFunctionIntegralCurrentValue(0.0) + axialBoundaryConditions_[0];
}

//! Compute current cylindrical state.
Eigen::Vector6d HodographicShapingLeg::computeStateVectorInCylindricalCoordinates( const double timeSinceDeparture )
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
Eigen::Vector6d HodographicShapingLeg::computeCurrentCartesianState(const double timeSinceDeparture )
{
    if ( timeSinceDeparture < 0.0 || timeSinceDeparture > timeOfFlight_ )
    {
        throw std::runtime_error( "Error when computing state vector, requested time is outside bounds" );
    }

    return coordinate_conversions::convertCylindricalToCartesianState( computeStateVectorInCylindricalCoordinates( timeSinceDeparture ) );
}


//! Compute thrust acceleration components.
Eigen::Vector3d HodographicShapingLeg::computeThrustAccelerationInCylindricalCoordinates( const double timeSinceDeparture )
{

    Eigen::Vector3d thrustAccelerationComponents;

    // Compute radial distance from the central body.
    double radialDistance = computeCurrentRadialDistance( timeSinceDeparture );

    // Compute axial distance from the central body.
    double axialDistance = computeCurrentAxialDistance( timeSinceDeparture );

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


double HodographicShapingLeg::computeThrustAccelerationMagnitude(
        const double timeSinceDeparture )
{
    // Return total thrust acceleration.
    return computeThrustAcceleration(timeSinceDeparture).norm();
}

Eigen::Vector3d HodographicShapingLeg::computeThrustAcceleration(const double timeSinceDeparture )
{
    if ( timeSinceDeparture < 0.0 || timeSinceDeparture > timeOfFlight_ )
    {
        throw std::runtime_error( "Error when computing acceleration vector, requested time is outside bounds" );
    }

    if( thrustAccelerationVectorCache_.count( timeSinceDeparture ) == 0 )
    {
        Eigen::Vector3d cylindricalAcceleration = computeThrustAccelerationInCylindricalCoordinates( timeSinceDeparture );
        Eigen::Vector3d cylindricalPosition = computeStateVectorInCylindricalCoordinates(timeSinceDeparture ).segment(0, 3);
        Eigen::Vector3d cartesianPosition = computeCurrentCartesianState(timeSinceDeparture).segment(0, 3);

        Eigen::Vector3d cartesianAcceleration;

        cartesianAcceleration[ 0 ] = ( 1.0 / ( cartesianPosition[ 0 ] + ( std::pow(cartesianPosition[ 1 ], 2 ) / cartesianPosition[ 0 ] ) ) )
                * ( cylindricalPosition[ 0 ] * cylindricalAcceleration[ 0 ]
                    - ( cartesianPosition[ 1 ] / cartesianPosition[ 0 ] ) * cylindricalPosition[ 0 ] * cylindricalAcceleration[ 1 ] );

        cartesianAcceleration[ 1 ] = ( cylindricalPosition[ 0 ] * cylindricalAcceleration[ 1 ]
                                       + cartesianPosition[ 1 ] * cartesianAcceleration[ 0 ] ) / cartesianPosition[ 0 ];

        cartesianAcceleration[ 2 ] = cylindricalAcceleration[ 2 ];
        thrustAccelerationVectorCache_[ timeSinceDeparture ] = cartesianAcceleration;
    }

    return thrustAccelerationVectorCache_[ timeSinceDeparture ];

}

//! Compute direction cartesian acceleration.
Eigen::Vector3d HodographicShapingLeg::computeThrustAccelerationDirection(
        const double timeSinceDeparture )
{
    return computeThrustAcceleration(timeSinceDeparture).normalized();
}


//! Compute DeltaV.
double HodographicShapingLeg::computeDeltaV( )
{
    // Define the derivative of the deltaV, ie thrust acceleration function, as a function of time.
    std::function< double( const double ) > derivativeFunctionDeltaV =
            std::bind(&HodographicShapingLeg::computeThrustAccelerationMagnitude, this,
                      std::placeholders::_1 );

    // Define numerical quadrature.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature(derivativeFunctionDeltaV, quadratureSettings_, timeOfFlight_ );

    return quadrature->getQuadrature( );
}


} // namespace shape_based_methods
} // namespace tudat
