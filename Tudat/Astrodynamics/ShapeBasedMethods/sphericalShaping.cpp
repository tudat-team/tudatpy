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


#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Basics/timeType.h"
#include "Tudat/Astrodynamics/Propulsion/thrustAccelerationModel.h"
#include "sphericalShaping.h"

namespace tudat
{
namespace shape_based_methods
{


//! Constructur for spherical shaping.
SphericalShaping::SphericalShaping( Eigen::Vector6d initialState,
                                    Eigen::Vector6d finalState,
                                    double requiredTimeOfFlight,
                                    double initialValueCoefficientRadialInversePolynomial,
                                    Eigen::VectorXd freeCoefficientsRadialFunction,
                                    Eigen::VectorXd freeCoefficientsElevationFunction,
                                    double centralBodyGravitationalParameter,
                                    root_finders::RootFinderType rootFinderType,
                                    const double lowerBoundFreeCoefficient,
                                    const double upperBoundFreeCoefficient,
                                    const double initialGuessForFreeCoefficient,
                                    const int maxNumberOfIterations,
                                    const double requiredToleranceForTimeOfFlight ):
    initialState_( initialState ), finalState_( finalState ),
    freeCoefficientsRadialFunction_( freeCoefficientsRadialFunction ),
    freeCoefficientsElevationFunction_( freeCoefficientsElevationFunction ),
    requiredTimeOfFlight_( requiredTimeOfFlight ),
    initialValueCoefficientRadialInversePolynomial_( initialValueCoefficientRadialInversePolynomial ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
    rootFinderType_( rootFinderType ),
    lowerBoundFreeCoefficient_( lowerBoundFreeCoefficient ),
    upperBoundFreeCoefficient_( upperBoundFreeCoefficient ),
    initialGuessForFreeCoefficient_( initialGuessForFreeCoefficient ),
    maxNumberOfIterations_( maxNumberOfIterations ),
    requiredToleranceForTimeOfFlight_( requiredToleranceForTimeOfFlight )
{

    initialStateSphericalCoordinates_ = coordinate_conversions::convertCartesianToSphericalState( initialState_ );

    finalStateSphericalCoordinates_ = coordinate_conversions::convertCartesianToSphericalState( finalState_ );
    std::cout << "initial state spherical coordinates: " << initialStateSphericalCoordinates_ << "\n\n";
    std::cout << "final state spherical coordinates: " << finalStateSphericalCoordinates_ << "\n\n";

    if ( initialStateSphericalCoordinates_( 1 ) < 0.0 )
    {
        initialStateSphericalCoordinates_( 1 ) += 2.0 * mathematical_constants::PI;
    }
    if ( finalStateSphericalCoordinates_( 1 ) < 0.0 )
    {
        finalStateSphericalCoordinates_( 1 ) += 2.0 * mathematical_constants::PI;
    }

    double numberOfRevolutions = 1.0;
    if ( ( finalStateSphericalCoordinates_( 1 ) - initialStateSphericalCoordinates_( 1 ) ) < 0.0 )
    {
        finalAzimuthalAngle_ = finalStateSphericalCoordinates_( 1 ) + 2.0 * mathematical_constants::PI * ( numberOfRevolutions + 1.0 );
    }
    else
    {
        finalAzimuthalAngle_ = finalStateSphericalCoordinates_( 1 ) + 2.0 * mathematical_constants::PI * ( numberOfRevolutions );
    }


    double initialDerivativeAzimuthalAngle = initialStateSphericalCoordinates_[ 4 ]
            / ( initialStateSphericalCoordinates_[ 0 ] * std::cos( initialStateSphericalCoordinates_[ 2 ] ) );
    double finalDerivativeAzimuthalAngle = finalStateSphericalCoordinates_[ 4 ]
            / ( finalStateSphericalCoordinates_[ 0 ] * std::cos( finalStateSphericalCoordinates_[ 2 ] ) );

    initialStateThetaParametrized_ = ( Eigen::Vector6d() << initialStateSphericalCoordinates_[ 0 ],
            initialStateSphericalCoordinates_[ 1 ],
            initialStateSphericalCoordinates_[ 2 ],
            initialStateSphericalCoordinates_[ 3 ] / initialDerivativeAzimuthalAngle,
            initialStateSphericalCoordinates_[ 4 ] / initialDerivativeAzimuthalAngle,
            initialStateSphericalCoordinates_[ 5 ] / initialDerivativeAzimuthalAngle ).finished();

    finalStateThetaParametrized_ = ( Eigen::Vector6d() << finalStateSphericalCoordinates_[ 0 ],
            finalStateSphericalCoordinates_[ 1 ],
            finalStateSphericalCoordinates_[ 2 ],
            finalStateSphericalCoordinates_[ 3 ] / finalDerivativeAzimuthalAngle,
            finalStateSphericalCoordinates_[ 4 ] / finalDerivativeAzimuthalAngle,
            finalStateSphericalCoordinates_[ 5 ] / finalDerivativeAzimuthalAngle ).finished();

    // Retrieve the initial value of the azimuth angle.
    initialAzimuthalAngle_ = initialStateSphericalCoordinates_[ 1 ];

    // Define coefficients for radial distance and elevation angle composite functions.
    compositeRadialFunction_ = std::make_shared< CompositeRadialFunctionSphericalShaping >( freeCoefficientsRadialFunction_ );
    compositeElevationFunction_ = std::make_shared< CompositeElevationFunctionSphericalShaping >( freeCoefficientsElevationFunction_ );

    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
    std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
            std::make_shared< root_finders::RootFinderSettings >( rootFinderType_, requiredToleranceForTimeOfFlight_, maxNumberOfIterations_ );


    // Define settings for numerical quadrature, to be used to compute time of flight and final deltaV.
    quadratureSettings_ = std::make_shared< numerical_quadrature::GaussianQuadratureSettings < double > >( initialAzimuthalAngle_, 64 );

    // Iterate on the free coefficient value until the time of flight matches its required value.
    iterateToMatchRequiredTimeOfFlight( rootFinderSettings, lowerBoundFreeCoefficient_, upperBoundFreeCoefficient_, initialGuessForFreeCoefficient_ );

}

Eigen::MatrixXd SphericalShaping::computeInverseMatrixBoundaryConditions( )
{
    Eigen::MatrixXd matrixBoundaryConditions = Eigen::MatrixXd::Zero( 10, 10 );

    for ( int i = 0 ; i < 6 ; i++ )
    {
        int index = i;
        if ( i >= 2 )
        {
            index = i + 1;
        }
        matrixBoundaryConditions( 0, i ) = compositeRadialFunction_->getComponentFunctionCurrentValue( index, initialAzimuthalAngle_ );
        matrixBoundaryConditions( 1, i ) = compositeRadialFunction_->getComponentFunctionCurrentValue( index, finalAzimuthalAngle_ );
        matrixBoundaryConditions( 2, i ) = compositeRadialFunction_->getComponentFunctionFirstDerivative( index, initialAzimuthalAngle_ );
        matrixBoundaryConditions( 3, i ) = compositeRadialFunction_->getComponentFunctionFirstDerivative( index, finalAzimuthalAngle_ );
        matrixBoundaryConditions( 4, i ) = - std::pow( initialStateSphericalCoordinates_[ 0 ], 2.0 )
                * compositeRadialFunction_->getComponentFunctionSecondDerivative( index, initialAzimuthalAngle_ );
        matrixBoundaryConditions( 5, i ) = - std::pow( finalStateSphericalCoordinates_[ 0 ], 2.0 )
                * compositeRadialFunction_->getComponentFunctionSecondDerivative( index, finalAzimuthalAngle_ );
    }

    // Compute value of variable alpha at initial time.
    double initialValueAlpha = computeInitialAlphaValue();

    // Compute value of variable alpha at final time.
    double finalValueAlpha = computeFinalAlphaValue();

    for ( int i = 0 ; i < 4 ; i++ )
    {
        matrixBoundaryConditions( 4, i + 6 ) =
                initialValueAlpha * compositeElevationFunction_->getComponentFunctionSecondDerivative( i, initialAzimuthalAngle_ );
        matrixBoundaryConditions( 5, i + 6 ) =
                finalValueAlpha * compositeElevationFunction_->getComponentFunctionSecondDerivative( i, finalAzimuthalAngle_ );
        matrixBoundaryConditions( 6, i + 6 ) = compositeElevationFunction_->getComponentFunctionCurrentValue( i, initialAzimuthalAngle_ );
        matrixBoundaryConditions( 7, i + 6 ) = compositeElevationFunction_->getComponentFunctionCurrentValue( i, finalAzimuthalAngle_ );
        matrixBoundaryConditions( 8, i + 6 ) = compositeElevationFunction_->getComponentFunctionFirstDerivative( i, initialAzimuthalAngle_ );
        matrixBoundaryConditions( 9, i + 6 ) = compositeElevationFunction_->getComponentFunctionFirstDerivative( i, finalAzimuthalAngle_ );
    }

    // Compute and return the inverse of the boundary conditions matrix.
    return matrixBoundaryConditions.inverse();
}

double SphericalShaping::computeInitialAlphaValue( )
{
    return - ( initialStateThetaParametrized_[ 3 ] * initialStateThetaParametrized_[ 5 ] / initialStateThetaParametrized_[ 0 ] )
            / ( std::pow( initialStateThetaParametrized_[ 5 ] / initialStateThetaParametrized_[ 0 ], 2.0 )
            + std::pow( std::cos( initialStateThetaParametrized_[ 2 ] ), 2.0 ) );
}

double SphericalShaping::computeFinalAlphaValue( )
{
    return - ( finalStateThetaParametrized_[ 3 ] *  finalStateThetaParametrized_[ 5 ] / finalStateThetaParametrized_[ 0 ] )
            / ( std::pow( finalStateThetaParametrized_[ 5 ] / finalStateThetaParametrized_[ 0 ], 2.0 )
            + std::pow( std::cos( finalStateThetaParametrized_[ 2 ] ), 2.0 ) );
}


void SphericalShaping::satisfyBoundaryConditions( )
{
    Eigen::VectorXd vectorBoundaryValues;
    vectorBoundaryValues.resize( 10 );

    vectorBoundaryValues[ 0 ] = 1.0 / initialStateThetaParametrized_[ 0 ];
    vectorBoundaryValues[ 1 ] = 1.0 / finalStateThetaParametrized_[ 0 ];
    vectorBoundaryValues[ 2 ] = - initialStateThetaParametrized_[ 3 ] / std::pow( initialStateThetaParametrized_[ 0 ], 2.0 );
    vectorBoundaryValues[ 3 ] = - finalStateThetaParametrized_[ 3 ] / std::pow( finalStateThetaParametrized_[ 0 ], 2.0 );
    vectorBoundaryValues[ 4 ] = computeInitialValueBoundaryConstant()
            - 2.0 * std::pow( initialStateThetaParametrized_[ 3 ], 2.0 ) / initialStateThetaParametrized_[ 0 ];
    vectorBoundaryValues[ 5 ] = computeFinalValueBoundaryConstant()
            - 2.0 * std::pow( finalStateThetaParametrized_[ 3 ], 2.0 ) / finalStateThetaParametrized_[ 0 ];
    vectorBoundaryValues[ 6 ] = initialStateThetaParametrized_[ 2 ];
    vectorBoundaryValues[ 7 ] = finalStateThetaParametrized_[ 2 ];
    vectorBoundaryValues[ 8 ] = initialStateThetaParametrized_[ 5 ] / initialStateThetaParametrized_[ 0 ];
    vectorBoundaryValues[ 9 ] = finalStateThetaParametrized_[ 5 ] / finalStateThetaParametrized_[ 0 ];

    Eigen::VectorXd vectorSecondComponentContribution;
    vectorSecondComponentContribution.resize( 10.0 );

    vectorSecondComponentContribution[ 0 ] = compositeRadialFunction_->getComponentFunctionCurrentValue( 2, initialAzimuthalAngle_ );
    vectorSecondComponentContribution[ 1 ] = compositeRadialFunction_->getComponentFunctionCurrentValue( 2, finalAzimuthalAngle_ );
    vectorSecondComponentContribution[ 2 ] = compositeRadialFunction_->getComponentFunctionFirstDerivative( 2, initialAzimuthalAngle_ );
    vectorSecondComponentContribution[ 3 ] = compositeRadialFunction_->getComponentFunctionFirstDerivative( 2, finalAzimuthalAngle_ );
    vectorSecondComponentContribution[ 4 ] = - std::pow( initialStateThetaParametrized_[ 0 ], 2 )
              * compositeRadialFunction_->getComponentFunctionSecondDerivative( 2, initialAzimuthalAngle_ );
    vectorSecondComponentContribution[ 5 ] = - std::pow( finalStateThetaParametrized_[ 0 ], 2 )
              * compositeRadialFunction_->getComponentFunctionSecondDerivative( 2, finalAzimuthalAngle_ );
    vectorSecondComponentContribution[ 6 ] = 0.0;
    vectorSecondComponentContribution[ 7 ] = 0.0;
    vectorSecondComponentContribution[ 8 ] = 0.0;
    vectorSecondComponentContribution[ 9 ] = 0.0;

    vectorSecondComponentContribution *= initialValueCoefficientRadialInversePolynomial_;

    Eigen::MatrixXd inverseMatrixBoundaryConditions = computeInverseMatrixBoundaryConditions();

    Eigen::MatrixXd compositeFunctionCoefficients = inverseMatrixBoundaryConditions * ( vectorBoundaryValues - vectorSecondComponentContribution );

    for ( int i = 0 ; i < 6 ; i++ )
    {
        if ( i < 2 )
        {
            freeCoefficientsRadialFunction_( i ) = compositeFunctionCoefficients( i );
        }
        else
        {
            freeCoefficientsRadialFunction_( i + 1 ) = compositeFunctionCoefficients( i );
        }
    }
    freeCoefficientsRadialFunction_( 2 ) = initialValueCoefficientRadialInversePolynomial_;

    for ( int i = 0 ; i < 4 ; i++ )
    {
        freeCoefficientsElevationFunction_( i ) = compositeFunctionCoefficients( i + 6 );
    }

    compositeRadialFunction_->resetCompositeFunctionCoefficients( freeCoefficientsRadialFunction_ );
    compositeElevationFunction_->resetCompositeFunctionCoefficients( freeCoefficientsElevationFunction_ );

}

double SphericalShaping::computeInitialValueBoundaryConstant( )
{
    double radialDistance = initialStateThetaParametrized_[ 0 ];
    double elevationAngle = initialStateThetaParametrized_[ 2 ];
    double derivativeRadialDistance = initialStateThetaParametrized_[ 3 ];
    double derivativeElevationAngle = initialStateThetaParametrized_[ 5 ] / initialStateThetaParametrized_[ 0 ];
    double derivativeOfTimeWrtAzimuthAngle = ( initialStateThetaParametrized_[ 0 ] * std::cos( initialStateThetaParametrized_[ 2 ] ) )
            / initialStateSphericalCoordinates_[ 4 ];

    return - centralBodyGravitationalParameter_ * std::pow( derivativeOfTimeWrtAzimuthAngle, 2.0 ) / std::pow( radialDistance, 2.0 )
            + 2.0 * std::pow( derivativeRadialDistance, 2.0 ) / radialDistance
            + radialDistance * ( std::pow( derivativeElevationAngle, 2.0 ) + std::pow( std::cos( elevationAngle ), 2.0 ) )
            - derivativeRadialDistance * derivativeElevationAngle * ( std::sin( elevationAngle ) * std::cos( elevationAngle ) )
            / ( std::pow( derivativeElevationAngle, 2.0 ) + std::pow( std::cos( elevationAngle ), 2.0 ) );
}

double SphericalShaping::computeFinalValueBoundaryConstant( )
{
    double radialDistance = finalStateThetaParametrized_[ 0 ];
    double elevationAngle = finalStateThetaParametrized_[ 2 ];
    double derivativeRadialDistance = finalStateThetaParametrized_[ 3 ];
    double derivativeElevationAngle = finalStateThetaParametrized_[ 5 ] / finalStateThetaParametrized_[ 0 ];
    double derivativeOfTimeWrtAzimuthAngle = ( finalStateThetaParametrized_[ 0 ] * std::cos( finalStateThetaParametrized_[ 2 ] ) )
            / finalStateSphericalCoordinates_[ 4 ];

    return - centralBodyGravitationalParameter_ * std::pow( derivativeOfTimeWrtAzimuthAngle, 2.0 ) / std::pow( radialDistance, 2.0 )
            + 2.0 * std::pow( derivativeRadialDistance, 2.0 ) / radialDistance
            + radialDistance * ( std::pow( derivativeElevationAngle, 2.0 ) + std::pow( std::cos( elevationAngle ), 2.0 ) )
            - derivativeRadialDistance * derivativeElevationAngle * ( std::sin( elevationAngle ) * std::cos( elevationAngle ) )
            / ( std::pow( derivativeElevationAngle, 2.0 ) + std::pow( std::cos( elevationAngle ), 2.0 ) );
}

double SphericalShaping::computeScalarFunctionTimeEquation( double currentAzimuthalAngle )
{
    double radialFunctionValue = compositeRadialFunction_->evaluateCompositeFunction( currentAzimuthalAngle );
    double firstDerivativeRadialFunction = compositeRadialFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthalAngle );
    double secondDerivativeRadialFunction = compositeRadialFunction_->evaluateCompositeFunctionSecondDerivative( currentAzimuthalAngle );

    double elevationFunctionValue = compositeElevationFunction_->evaluateCompositeFunction( currentAzimuthalAngle );
    double firstDerivativeElevationFunction = compositeElevationFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthalAngle );
    double secondDerivativeElevationFunction = compositeElevationFunction_->evaluateCompositeFunctionSecondDerivative( currentAzimuthalAngle );

    return - secondDerivativeRadialFunction + 2.0 * std::pow( firstDerivativeRadialFunction, 2.0 ) / radialFunctionValue
            + firstDerivativeRadialFunction * firstDerivativeElevationFunction
            * ( secondDerivativeElevationFunction - std::sin( elevationFunctionValue ) * std::cos( elevationFunctionValue ) )
            / ( std::pow( firstDerivativeElevationFunction, 2.0 ) + std::pow( std::cos( elevationFunctionValue ), 2.0 ) )
            + radialFunctionValue * ( std::pow( firstDerivativeElevationFunction, 2.0 ) + std::pow( std::cos( elevationFunctionValue ), 2.0 ) );
}

double SphericalShaping::computeDerivativeScalarFunctionTimeEquation( double currentAzimuthalAngle )
{
    double radialFunctionValue = compositeRadialFunction_->evaluateCompositeFunction( currentAzimuthalAngle );
    double firstDerivativeRadialFunction = compositeRadialFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthalAngle );
    double secondDerivativeRadialFunction = compositeRadialFunction_->evaluateCompositeFunctionSecondDerivative( currentAzimuthalAngle );
    double thirdDerivativeRadialFunction = compositeRadialFunction_->evaluateCompositeFunctionThirdDerivative( currentAzimuthalAngle );

    double elevationFunctionValue = compositeElevationFunction_->evaluateCompositeFunction( currentAzimuthalAngle );
    double firstDerivativeElevationFunction = compositeElevationFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthalAngle );
    double secondDerivativeElevationFunction = compositeElevationFunction_->evaluateCompositeFunctionSecondDerivative( currentAzimuthalAngle );
    double thirdDerivativeElevationFunction = compositeElevationFunction_->evaluateCompositeFunctionThirdDerivative( currentAzimuthalAngle );

    // Define constants F1, F2, F3 and F4 as proposed in... (ADD REFERENCE).
    double F1 = std::pow( firstDerivativeElevationFunction, 2.0 ) + std::pow( std::cos( elevationFunctionValue ), 2.0 );
    double F2 = secondDerivativeElevationFunction - std::sin( 2.0 * elevationFunctionValue ) / 2.0;
    double F3 = std::cos( 2.0 * elevationFunctionValue ) + 2.0 * std::pow( firstDerivativeElevationFunction, 2.0 ) + 1.0;
    double F4 = 2.0 * secondDerivativeElevationFunction - std::sin( 2.0 * elevationFunctionValue );

    return  F1 * firstDerivativeRadialFunction - thirdDerivativeRadialFunction
            - 2.0 * std::pow( firstDerivativeRadialFunction, 3.0 ) / std::pow( radialFunctionValue, 2.0 )
            + 4.0 * firstDerivativeRadialFunction * secondDerivativeRadialFunction / radialFunctionValue
            + F4 * firstDerivativeElevationFunction * radialFunctionValue
            + 2.0 * firstDerivativeElevationFunction * firstDerivativeRadialFunction
            * ( thirdDerivativeElevationFunction - firstDerivativeElevationFunction * std::cos( 2.0 * elevationFunctionValue ) ) / F3
            + F2 * firstDerivativeElevationFunction * secondDerivativeRadialFunction / F1
            + F2 * firstDerivativeRadialFunction * secondDerivativeElevationFunction / F1
            - 4.0 * F4 * F2 * std::pow( firstDerivativeElevationFunction, 2.0 ) * firstDerivativeRadialFunction / std::pow( F3, 2.0 );
}

double SphericalShaping::computeTimeOfFlight()
{
    // Compute step size.
    double stepSize = 2.0 * mathematical_constants::PI / 100.0;

    // Check that the trajectory is feasible, ie curved toward the central body.
    for ( int i = 0 ; i < std::ceil( ( finalAzimuthalAngle_ - initialAzimuthalAngle_ ) / stepSize ) ; i++ )
    {
        if ( computeScalarFunctionTimeEquation( initialAzimuthalAngle_ + i * stepSize ) < 0.0 )
        {
            std::cerr << "Error, trajectory not curved toward the central body, and thus not feasible." << "\n\n";
        }
    }

    // Define the derivative of the time function w.r.t the azimuthal angle theta.
    std::function< double( const double ) > derivativeTimeFunction = [ = ] ( const double currentAzimuthalAngle ){

        return std::sqrt( computeScalarFunctionTimeEquation( currentAzimuthalAngle ) *
                          std::pow( compositeRadialFunction_->evaluateCompositeFunction( currentAzimuthalAngle ), 2.0 )
                          / centralBodyGravitationalParameter_ );
    };

    // Create numerical quadrature from quadrature settings.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeTimeFunction, quadratureSettings_, finalAzimuthalAngle_ );

    return quadrature->getQuadrature( );
}

//! Compute current time from azimuth angle.
double SphericalShaping::computeCurrentTimeFromAzimuthAngle( const double currentAzimuthAngle )
{
//    // Compute step size.
//    double stepSize = 2.0 * mathematical_constants::PI / 100.0;

//    // Check that the trajectory is feasible, ie curved toward the central body.
//    for ( int i = 0 ; i < std::ceil( ( finalAzimuthalAngle_ - initialAzimuthalAngle_ ) / stepSize ) ; i++ )
//    {
//        if ( computeScalarFunctionTimeEquation( initialAzimuthalAngle_ + i * stepSize ) < 0.0 )
//        {
//            std::cerr << "Error, trajectory not curved toward the central body, and thus not feasible." << "\n\n";
//        }
//    }

    // Define the derivative of the time function w.r.t the azimuthal angle theta.
    std::function< double( const double ) > derivativeTimeFunction = [ = ] ( const double currentAzimuthalAngle ){

        return std::sqrt( computeScalarFunctionTimeEquation( currentAzimuthalAngle ) *
                          std::pow( compositeRadialFunction_->evaluateCompositeFunction( currentAzimuthalAngle ), 2.0 )
                          / centralBodyGravitationalParameter_ );
    };

    // Create numerical quadrature from quadrature settings.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeTimeFunction, quadratureSettings_, currentAzimuthAngle );

    return quadrature->getQuadrature( );
}

//! Iterate to match the required time of flight
void SphericalShaping::iterateToMatchRequiredTimeOfFlight( std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings,
                                                           const double lowerBound,
                                                           const double upperBound,
                                                           const double initialGuess )
{

    // Define the structure updating the time of flight from the free coefficient value, while still satisfying the boundary conditions.
    std::function< void ( const double ) > resetFreeCoefficientFunction = std::bind( &SphericalShaping::resetValueFreeCoefficient, this, std::placeholders::_1 );
    std::function< void( ) > satisfyBoundaryConditionsFunction = std::bind( &SphericalShaping::satisfyBoundaryConditions, this );
    std::function< double ( ) >  computeTOFfunction = std::bind( &SphericalShaping::computeTimeOfFlight, this );
    std::function< double ( ) > getRequiredTOFfunction = std::bind( &SphericalShaping::getRequiredTimeOfFlight, this );

    std::shared_ptr< basic_mathematics::Function< double, double > > timeOfFlightFunction =
            std::make_shared< SphericalShaping::TOFfunction >( resetFreeCoefficientFunction, satisfyBoundaryConditionsFunction, computeTOFfunction, getRequiredTOFfunction );

    // Create root finder from root finder settings.
    std::shared_ptr< root_finders::RootFinderCore< double > > rootFinder = root_finders::createRootFinder( rootFinderSettings, lowerBound, upperBound, initialGuess );

    // Iterate to find the free coefficient value that matches the required time of flight.
    double updatedFreeCoefficient = rootFinder->execute( timeOfFlightFunction, initialGuess );

    std::cout << "updated free coefficient: " << updatedFreeCoefficient << "\n\n";
    std::cout << "final difference in TOF: " << getRequiredTOFfunction() - computeTOFfunction() << "\n\n";
}

//! Compute current spherical position.
Eigen::Vector3d SphericalShaping::computeCurrentSphericalPosition( const double currentAzimuthalAngle )
{

    return ( Eigen::Vector3d() <<
             compositeRadialFunction_->evaluateCompositeFunction( currentAzimuthalAngle ),
             currentAzimuthalAngle,
             compositeElevationFunction_->evaluateCompositeFunction( currentAzimuthalAngle ) ).finished();

}

//! Compute current cartesian position.
Eigen::Vector6d SphericalShaping::computeCurrentCartesianPosition( const double currentAzimuthalAngle )
{
    Eigen::Vector3d sphericalPosition = computeCurrentSphericalPosition( currentAzimuthalAngle );
    return coordinate_conversions::convertSphericalToCartesianState( ( Eigen::Vector6d() << sphericalPosition[ 0 ],
                                                                sphericalPosition[ 1 ], sphericalPosition[ 2 ], 0.0, 0.0, 0.0 ).finished() );
}

//! Compute current derivative of the azimuthal angle.
double SphericalShaping::computeFirstDerivativeAzimuthalAngle( const double currentAzimuthalAngle )
{
    // Compute scalar function time equation.
    double scalarFunctionTimeEquation = computeScalarFunctionTimeEquation( currentAzimuthalAngle );

    // Compute current radial distance.
    double radialDistance = compositeRadialFunction_->evaluateCompositeFunction( currentAzimuthalAngle );

    // Compute and return the first derivative of the azimuth angle w.r.t. time.
    return std::sqrt( centralBodyGravitationalParameter_ / ( scalarFunctionTimeEquation * std::pow( radialDistance, 2.0 ) ) );
}

//! Compute second derivative of the azimuthal angle.
double SphericalShaping::computeSecondDerivativeAzimuthalAngle( const double currentAzimuthalAngle )
{
    // Compute first derivative azimuth angle w.r.t. time.
    double firstDerivativeAzimuthAngle = computeFirstDerivativeAzimuthalAngle( currentAzimuthalAngle );

    // Compute scalar function of the time equation, and its derivative w.r.t. azimuth angle.
    double scalarFunctionTimeEquation = computeScalarFunctionTimeEquation( currentAzimuthalAngle );
    double derivativeScalarFunctionTimeEquation = computeDerivativeScalarFunctionTimeEquation( currentAzimuthalAngle );

    // Compute radial distance, and its derivative w.r.t. azimuth angle.
    double radialDistance = compositeRadialFunction_->evaluateCompositeFunction( currentAzimuthalAngle );
    double firstDerivativeRadialDistance = compositeRadialFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthalAngle );

    return - std::pow( firstDerivativeAzimuthAngle, 2.0 )
            * ( derivativeScalarFunctionTimeEquation / ( 2.0 * scalarFunctionTimeEquation ) + firstDerivativeRadialDistance / radialDistance );
}

//! Compute current velocity in spherical coordinates.
Eigen::Vector3d SphericalShaping::computeCurrentSphericalVelocity(const double currentAzimuthalAngle )
{
    // Compute first derivative of the azimuth angle w.r.t. time.
    double derivativeAzimuthalAngle = computeFirstDerivativeAzimuthalAngle( currentAzimuthalAngle );

    // Compute and return current velocity vector in spherical coordinates.
    return derivativeAzimuthalAngle * computeCurrentVelocityParametrizedByAzimuthAngle( currentAzimuthalAngle );
}

//! Compute current velocity parametrized by azimuthal angle theta.
Eigen::Vector3d SphericalShaping::computeCurrentVelocityParametrizedByAzimuthAngle( const double currentAzimuthalAngle )
{

    // Retrieve current radial distance and elevation angle, as well as their derivatives w.r.t. azimuth angle.
    double radialDistance = compositeRadialFunction_->evaluateCompositeFunction( currentAzimuthalAngle );
    double derivativeRadialDistance = compositeRadialFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthalAngle );
    double elevationAngle = compositeElevationFunction_->evaluateCompositeFunction( currentAzimuthalAngle );
    double derivativeElevationAngle = compositeElevationFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthalAngle );

    // Compute and return velocity vector parametrized by azimuth angle.
    return ( Eigen::Vector3d() << derivativeRadialDistance, radialDistance * std::cos( elevationAngle ), radialDistance * derivativeElevationAngle ).finished();
}

//! Compute current acceleration in spherical coordinates.
Eigen::Vector3d SphericalShaping::computeCurrentSphericalAcceleration( const double currentAzimuthalAngle )
{
    // Compute first derivative of the azimuth angle w.r.t. time.
    double derivativeAzimuthalAngle = computeFirstDerivativeAzimuthalAngle( currentAzimuthalAngle );

    // Compute second derivative of the azimuth angle w.r.t. time.
    double secondDerivativeAzimuthalAngle = computeSecondDerivativeAzimuthalAngle( currentAzimuthalAngle );

    // Compute velocity vector parametrized by the azimuthal angle theta.
    Eigen::Vector3d velocityVectorThetaParametrized = computeCurrentVelocityParametrizedByAzimuthAngle( currentAzimuthalAngle );

    // Compute acceleration vector parametrized by the azimuthal angle theta.
    Eigen::Vector3d accelerationParametrizedByAzimuthAngle = computeCurrentAccelerationParametrizedByAzimuthAngle( currentAzimuthalAngle );

    // Compute and return acceleration vector parametrized w.r.t. time in spherical coordinates.
    return secondDerivativeAzimuthalAngle * velocityVectorThetaParametrized + std::pow( derivativeAzimuthalAngle, 2.0 ) * accelerationParametrizedByAzimuthAngle;
}

//! Compute current acceleration parametrized by azimuthal angle theta.
Eigen::Vector3d SphericalShaping::computeCurrentAccelerationParametrizedByAzimuthAngle( const double currentAzimuthalAngle )
{
    // Retrieve spherical coordinates and their derivatives w.r.t. to the azimuthal angle.
    double radialDistance = compositeRadialFunction_->evaluateCompositeFunction( currentAzimuthalAngle );
    double firstDerivativeRadialDistance = compositeRadialFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthalAngle );
    double secondDerivativeRadialDistance = compositeRadialFunction_->evaluateCompositeFunctionSecondDerivative( currentAzimuthalAngle );
    double elevationAngle = compositeElevationFunction_->evaluateCompositeFunction( currentAzimuthalAngle );
    double firstDerivativeElevationAngle = compositeElevationFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthalAngle );
    double secondDerivativeElevationAngle = compositeElevationFunction_->evaluateCompositeFunctionSecondDerivative( currentAzimuthalAngle );

    // Compute and return acceleration vector parametrized by the azimuthal angle theta.
    Eigen::Vector3d accelerationParametrizedByAzimuthAngle;
    accelerationParametrizedByAzimuthAngle[ 0 ] = secondDerivativeRadialDistance
            - radialDistance * ( std::pow( firstDerivativeElevationAngle, 2.0 ) + std::pow( std::cos( elevationAngle ), 2.0 ) );
    accelerationParametrizedByAzimuthAngle[ 1 ] = 2.0 * firstDerivativeRadialDistance * std::cos( elevationAngle )
            - 2.0 * radialDistance * firstDerivativeElevationAngle * std::sin( elevationAngle );
    accelerationParametrizedByAzimuthAngle[ 2 ] = 2.0 * firstDerivativeRadialDistance * firstDerivativeElevationAngle
            + radialDistance * ( secondDerivativeElevationAngle + std::sin( elevationAngle ) * std::cos( elevationAngle ) );

    return accelerationParametrizedByAzimuthAngle;

}

//! Compute control acceleration vector in spherical coordinates.
Eigen::Vector3d SphericalShaping::computeSphericalControlAccelerationVector( const double currentAzimuthalAngle )
{
    // Compute current radial distance.
    double radialDistance = compositeRadialFunction_->evaluateCompositeFunction( currentAzimuthalAngle );

    // Compute first and second derivatives of the azimuth angle w.r.t. time.
    double firstDerivativeAzimuthAngleWrtTime = computeFirstDerivativeAzimuthalAngle( currentAzimuthalAngle );
    double secondDerivativeAzimuthAngleWrtTime = computeSecondDerivativeAzimuthalAngle( currentAzimuthalAngle );

    // Compute velocity vector parametrized by azimuth angle theta.
    Eigen::Vector3d velocityParametrizedByAzimuthAngle = computeCurrentVelocityParametrizedByAzimuthAngle( currentAzimuthalAngle );

    // Compute acceleration vector parametrized by azimuth angle theta.
    Eigen::Vector3d accelerationParametrizedByAzimuthAngle = computeCurrentAccelerationParametrizedByAzimuthAngle( currentAzimuthalAngle );


    // Compute and return the current thrust acceleration vector in spherical coordinates.
    return std::pow( firstDerivativeAzimuthAngleWrtTime, 2.0 ) * accelerationParametrizedByAzimuthAngle
            + secondDerivativeAzimuthAngleWrtTime * velocityParametrizedByAzimuthAngle
            + centralBodyGravitationalParameter_ / std::pow( radialDistance, 3.0 ) * ( Eigen::Vector3d() << radialDistance, 0.0, 0.0 ).finished();

}

//! Compute current control acceleration in cartesian coordinates.
Eigen::Vector3d SphericalShaping::computeCartesianControlAcceleration( const double currentAzimuthAngle )
{
    Eigen::Vector3d cartesianAcceleration;

    Eigen::Vector6d sphericalStateToBeConverted;
    sphericalStateToBeConverted.segment( 0, 3 ) = computeCurrentSphericalPosition( currentAzimuthAngle );
    sphericalStateToBeConverted.segment( 3, 3 ) = computeSphericalControlAccelerationVector( currentAzimuthAngle );

    cartesianAcceleration = coordinate_conversions::convertSphericalToCartesianState( sphericalStateToBeConverted ).segment( 3, 3 );

    return cartesianAcceleration;
}

//! Compute magnitude cartesian acceleration.
double  SphericalShaping::computeMagnitudeCartesianAcceleration( double currentAzimuthAngle )
{
    return computeCartesianControlAcceleration( currentAzimuthAngle ).norm();
}

//! Compute direction cartesian acceleration.
Eigen::Vector3d SphericalShaping::computeDirectionCartesianAcceleration( double currentAzimuthAngle )
{
    return computeCartesianControlAcceleration( currentAzimuthAngle ).normalized();
}

//! Compute final deltaV.
double SphericalShaping::computeDeltav( )
{
    // Define the derivative of the deltaV, ie thrust acceleration function, as a function of the azimuth angle.
    std::function< double( const double ) > derivativeFunctionDeltaV = [ = ] ( const double currentAzimuthAngle ){

        double thrustAcceleration = computeSphericalControlAccelerationVector( currentAzimuthAngle ).norm()
                * std::sqrt( computeScalarFunctionTimeEquation( currentAzimuthAngle )
                             * std::pow( compositeRadialFunction_->evaluateCompositeFunction( currentAzimuthAngle ), 2.0 ) / centralBodyGravitationalParameter_ );

        return thrustAcceleration;

    };

    // Define numerical quadrature from quadratrure settings.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeFunctionDeltaV, quadratureSettings_, finalAzimuthalAngle_ );

    return quadrature->getQuadrature( );
}

//! Compute current spherical state.
Eigen::Vector6d SphericalShaping::computeCurrentSphericalState( const double currentAzimuthalAngle )
{
    Eigen::Vector6d currentSphericalState;
    currentSphericalState.segment( 0, 3 ) = computeCurrentSphericalPosition( currentAzimuthalAngle );
    currentSphericalState.segment( 3, 3 ) = computeCurrentSphericalVelocity( currentAzimuthalAngle );

    return currentSphericalState;
}

//! Compute current cartesian state.
Eigen::Vector6d SphericalShaping::computeCurrentCartesianState( const double currentAzimuthalAngle )
{
    return coordinate_conversions::convertSphericalToCartesianState( computeCurrentSphericalState( currentAzimuthalAngle ) );
}


//! Get low-thrust acceleration model from shaping method.
std::shared_ptr< propulsion::ThrustAcceleration > SphericalShaping::getLowThrustAccelerationModel(
        simulation_setup::NamedBodyMap& bodyMap,
        const std::string& bodyToPropagate,
        std::function< double( const double ) > specificImpulseFunction,
        std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > interpolatorPolarAngleFromTime )
{

    std::shared_ptr< simulation_setup::Body > vehicle = bodyMap[ bodyToPropagate ];

    // Define thrust magnitude function from the shaped trajectory.
    std::function< double( const double ) > thrustMagnitudeFunction = [ = ]( const double currentTime )
    {

        // Compute current azimuth angle.
        double currentAzimuthAngle = interpolatorPolarAngleFromTime->interpolate( currentTime );

        // Compute current acceleration.
        double currentAcceleration = computeMagnitudeCartesianAcceleration( currentAzimuthAngle ) * physical_constants::ASTRONOMICAL_UNIT
                / std::pow( physical_constants::JULIAN_YEAR, 2.0 );

        // Compute current mass of the vehicle.
        double currentMass = vehicle->getBodyMass();

        // Compute and return magnitude of the low-thrust force.
        return currentAcceleration * currentMass;
    };

    // Define thrust magnitude settings from thrust magnitude function.
    std::shared_ptr< simulation_setup::FromFunctionThrustMagnitudeSettings > thrustMagnitudeSettings
            = std::make_shared< simulation_setup::FromFunctionThrustMagnitudeSettings >(
                thrustMagnitudeFunction, specificImpulseFunction );


    // Define thrust direction function from the shaped trajectory.
    std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction = [ = ]( const double currentTime )
    {
        // Compute current azimuth angle.
        double currentAzimuthAngle = interpolatorPolarAngleFromTime->interpolate( currentTime );

        // Compute current acceleration.
        Eigen::Vector3d currentAccelerationDirection = computeDirectionCartesianAcceleration( currentAzimuthAngle );

        // Return direction of the low-thrust acceleration.
        return currentAccelerationDirection;
    };

    // Define thrust direction settings from the direction of thrust acceleration retrieved from the shaping method.
    std::shared_ptr< simulation_setup::CustomThrustDirectionSettings > thrustDirectionSettings =
            std::make_shared< simulation_setup::CustomThrustDirectionSettings >( thrustDirectionFunction );

    // Define thrust acceleration settings.
    std::shared_ptr< simulation_setup::ThrustAccelerationSettings > thrustAccelerationSettings =
            std::make_shared< simulation_setup::ThrustAccelerationSettings >(
                thrustDirectionSettings, thrustMagnitudeSettings );

    // Create low thrust acceleration model.
    std::shared_ptr< propulsion::ThrustAcceleration > lowThrustAccelerationModel = createThrustAcceleratioModel( thrustAccelerationSettings, bodyMap, bodyToPropagate );

    return lowThrustAccelerationModel;
}


void SphericalShaping::computeShapingTrajectoryAndFullPropagation(simulation_setup::NamedBodyMap& bodyMap,
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

    // Retrieve initial step size.
    double initialStepSize = integratorSettings->initialTimeStep_;

    // Vector of azimuth angles at which the time should be computed.
    Eigen::VectorXd azimuthAnglesToComputeAssociatedEpochs =
            Eigen::VectorXd::LinSpaced( std::ceil( computeTimeOfFlight() * physical_constants::JULIAN_YEAR / initialStepSize ),
                                        initialAzimuthalAngle_, finalAzimuthalAngle_ );

    std::map< double, double > dataToInterpolate;
    for ( int i = 0 ; i < azimuthAnglesToComputeAssociatedEpochs.size() ; i++ )
    {
        dataToInterpolate[ computeCurrentTimeFromAzimuthAngle( azimuthAnglesToComputeAssociatedEpochs[ i ] ) * physical_constants::JULIAN_YEAR ]
                = azimuthAnglesToComputeAssociatedEpochs[ i ];
    }

    // Create interpolator.
    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
            std::make_shared< interpolators::LagrangeInterpolatorSettings >( 10 );

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, double > > interpolator =
            interpolators::createOneDimensionalInterpolator( dataToInterpolate, interpolatorSettings );

    // Compute halved time of flight.
    double halvedTimeOfFlight = computeTimeOfFlight() / 2.0;

    // Compute azimuth angle at half of the time of flight.
    double azimuthAngleAtHalvedTimeOfFlight = interpolator->interpolate( halvedTimeOfFlight * physical_constants::JULIAN_YEAR );

    // Compute state at half of the time of flight.
    Eigen::Vector6d initialStateAtHalvedTimeOfFlight = computeCurrentCartesianState( azimuthAngleAtHalvedTimeOfFlight );
    initialStateAtHalvedTimeOfFlight.segment( 0, 3 ) *= physical_constants::ASTRONOMICAL_UNIT;
    initialStateAtHalvedTimeOfFlight.segment( 3, 3 ) *= physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_YEAR;

    // Create low thrust acceleration model.
    std::shared_ptr< propulsion::ThrustAcceleration > lowThrustAccelerationModel =
            getLowThrustAccelerationModel( bodyMap, bodyToPropagate, specificImpulseFunction, interpolator );

    accelerationMap[ bodyToPropagate ][ bodyToPropagate ].push_back( lowThrustAccelerationModel );


    std::vector< std::string > centralBodies;
    centralBodies.push_back( centralBody );

    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( bodyToPropagate );

    // Define forward propagator settings variables.
    integratorSettings->initialTime_ = halvedTimeOfFlight * physical_constants::JULIAN_YEAR;

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
        double currentAzimuthAngle = interpolator->interpolate( itr->first );

        Eigen::Vector6d currentNormalisedState = computeCurrentCartesianState( currentAzimuthAngle );
        Eigen::Vector6d currentState;
        currentState.segment( 0, 3 ) = currentNormalisedState.segment( 0, 3 ) * physical_constants::ASTRONOMICAL_UNIT;
        currentState.segment( 3, 3 ) = currentNormalisedState.segment( 3, 3 ) * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_YEAR;
        shapingMethodResults[ itr->first ] = currentState;
        fullPropagationResults[ itr->first ] = itr->second;
        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryForwardPropagation[ itr->first ];
    }


    // Define backward propagator settings variables.
    integratorSettings->initialTimeStep_ = - integratorSettings->initialTimeStep_;
    integratorSettings->initialTime_ = halvedTimeOfFlight * physical_constants::JULIAN_YEAR;

    // Perform the backward propagation.
    propagators::SingleArcDynamicsSimulator< > dynamicsSimulatorIntegrationBackwards( bodyMap, integratorSettings, propagatorSettingsBackwardPropagation );
    std::map< double, Eigen::VectorXd > stateHistoryFullProblemBackwardPropagation = dynamicsSimulatorIntegrationBackwards.getEquationsOfMotionNumericalSolution( );
    std::map< double, Eigen::VectorXd > dependentVariableHistoryBackwardPropagation = dynamicsSimulatorIntegrationBackwards.getDependentVariableHistory( );

    // Compute and save full propagation and shaping method results along the backward propagation direction
    for( std::map< double, Eigen::VectorXd >::iterator itr = stateHistoryFullProblemBackwardPropagation.begin( );
         itr != stateHistoryFullProblemBackwardPropagation.end( ); itr++ )
    {
        double currentAzimuthAngle = interpolator->interpolate( itr->first );

        Eigen::Vector6d currentNormalisedState = computeCurrentCartesianState( currentAzimuthAngle );
        Eigen::Vector6d currentState;
        currentState.segment( 0, 3 ) = currentNormalisedState.segment( 0, 3 ) * physical_constants::ASTRONOMICAL_UNIT;
        currentState.segment( 3, 3 ) = currentNormalisedState.segment( 3, 3 ) * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_YEAR;
        shapingMethodResults[ itr->first ] = currentState;
        fullPropagationResults[ itr->first ] = itr->second;
        dependentVariablesHistory[ itr->first ] = dependentVariableHistoryBackwardPropagation[ itr->first ];
    }

    // Reset initial integrator settings.
    integratorSettings->initialTimeStep_ = - integratorSettings->initialTimeStep_;

}


} // namespace shape_based_methods
} // namespace tudat
