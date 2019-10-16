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
                                    int numberOfRevolutions,
                                    simulation_setup::NamedBodyMap& bodyMap,
                                    const std::string bodyToPropagate,
                                    const std::string centralBody,
                                    double initialValueFreeCoefficient,
                                    std::shared_ptr< root_finders::RootFinderSettings >& rootFinderSettings,
                                    const double lowerBoundFreeCoefficient,
                                    const double upperBoundFreeCoefficient,
                                    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings ):
    ShapeBasedMethodLeg( initialState, finalState, requiredTimeOfFlight, bodyMap, bodyToPropagate, centralBody, integratorSettings ),
    numberOfRevolutions_( numberOfRevolutions ),
    initialValueFreeCoefficient_( initialValueFreeCoefficient ),
    rootFinderSettings_( rootFinderSettings ),
    lowerBoundFreeCoefficient_( lowerBoundFreeCoefficient ),
    upperBoundFreeCoefficient_( upperBoundFreeCoefficient ),
    integratorSettings_( integratorSettings )
{
    // Retrieve gravitational parameter of the central body.
    double centralBodyGravitationalParameter = bodyMap_[ centralBody_ ]->getGravityFieldModel()->getGravitationalParameter( );

    // Normalize the initial state.
    stateAtDeparture_.segment( 0, 3 ) = initialState.segment( 0, 3 ) / physical_constants::ASTRONOMICAL_UNIT;
    stateAtDeparture_.segment( 3, 3 ) = initialState.segment( 3, 3 ) * physical_constants::JULIAN_YEAR / physical_constants::ASTRONOMICAL_UNIT;

    // Normalize the final state.
    stateAtArrival_.segment( 0, 3 ) = finalState.segment( 0, 3 ) / physical_constants::ASTRONOMICAL_UNIT;
    stateAtArrival_.segment( 3, 3 ) = finalState.segment( 3, 3 ) * physical_constants::JULIAN_YEAR / physical_constants::ASTRONOMICAL_UNIT;

    // Normalize the required time of flight.
    requiredTimeOfFlight_ = requiredTimeOfFlight / physical_constants::JULIAN_YEAR;

    // Normalize the gravitational parameter of the central body.
    centralBodyGravitationalParameter_ = centralBodyGravitationalParameter * std::pow( physical_constants::JULIAN_YEAR, 2.0 )
            / std::pow( physical_constants::ASTRONOMICAL_UNIT, 3.0 );


    // Compute initial state in spherical coordinates.
    initialStateSphericalCoordinates_ = coordinate_conversions::convertCartesianToSphericalState( stateAtDeparture_ );

    // Compute final state in spherical coordinates.
    finalStateSphericalCoordinates_ = coordinate_conversions::convertCartesianToSphericalState( stateAtArrival_ );

    if ( initialStateSphericalCoordinates_( 1 ) < 0.0 )
    {
        initialStateSphericalCoordinates_( 1 ) += 2.0 * mathematical_constants::PI;
    }
    if ( finalStateSphericalCoordinates_( 1 ) < 0.0 )
    {
        finalStateSphericalCoordinates_( 1 ) += 2.0 * mathematical_constants::PI;
    }

    // Retrieve the initial value of the azimuth angle.
    initialAzimuthAngle_ = initialStateSphericalCoordinates_[ 1 ];

    // Compute final value of the azimuth angle.
    if ( ( finalStateSphericalCoordinates_( 1 ) - initialStateSphericalCoordinates_( 1 ) ) < 0.0 )
    {
        finalAzimuthAngle_ = finalStateSphericalCoordinates_( 1 ) + 2.0 * mathematical_constants::PI * ( numberOfRevolutions_ + 1.0 );
    }
    else
    {
        finalAzimuthAngle_ = finalStateSphericalCoordinates_( 1 ) + 2.0 * mathematical_constants::PI * ( numberOfRevolutions_ );
    }


    // Compute initial and final values of the derivative of the azimuth angle w.r.t. time.
    double initialDerivativeAzimuthAngle = initialStateSphericalCoordinates_[ 4 ]
            / ( initialStateSphericalCoordinates_[ 0 ] * std::cos( initialStateSphericalCoordinates_[ 2 ] ) );
    double finalDerivativeAzimuthAngle = finalStateSphericalCoordinates_[ 4 ]
            / ( finalStateSphericalCoordinates_[ 0 ] * std::cos( finalStateSphericalCoordinates_[ 2 ] ) );

    // Compute initial state parametrized by azimuth angle theta.
    initialStateParametrizedByAzimuthAngle_ = ( Eigen::Vector6d() << initialStateSphericalCoordinates_[ 0 ],
            initialStateSphericalCoordinates_[ 1 ],
            initialStateSphericalCoordinates_[ 2 ],
            initialStateSphericalCoordinates_[ 3 ] / initialDerivativeAzimuthAngle,
            initialStateSphericalCoordinates_[ 4 ] / initialDerivativeAzimuthAngle,
            initialStateSphericalCoordinates_[ 5 ] / initialDerivativeAzimuthAngle ).finished();

    // Compute final state parametrized by azimuth angle theta.
    finalStateParametrizedByAzimuthAngle_ = ( Eigen::Vector6d() << finalStateSphericalCoordinates_[ 0 ],
            finalStateSphericalCoordinates_[ 1 ],
            finalStateSphericalCoordinates_[ 2 ],
            finalStateSphericalCoordinates_[ 3 ] / finalDerivativeAzimuthAngle,
            finalStateSphericalCoordinates_[ 4 ] / finalDerivativeAzimuthAngle,
            finalStateSphericalCoordinates_[ 5 ] / finalDerivativeAzimuthAngle ).finished();


    // Initialise coefficients for radial distance and elevation angle functions.
    coefficientsRadialDistanceFunction_.resize( 7 );
    for ( int i = 0 ; i < 7 ; i++ )
    {
        coefficientsRadialDistanceFunction_[ i ] = 1.0;
    }
    coefficientsElevationAngleFunction_.resize( 4 );
    for ( int i = 0 ; i < 4 ; i++ )
    {
        coefficientsElevationAngleFunction_[ i ] = 1.0;
    }


    // Define coefficients for radial distance and elevation angle composite functions.
    radialDistanceCompositeFunction_ = std::make_shared< CompositeRadialFunctionSphericalShaping >( coefficientsRadialDistanceFunction_ );
    elevationAngleCompositeFunction_ = std::make_shared< CompositeElevationFunctionSphericalShaping >( coefficientsElevationAngleFunction_ );


    // Define settings for numerical quadrature, to be used to compute time of flight and final deltaV.
    quadratureSettings_ = std::make_shared< numerical_quadrature::GaussianQuadratureSettings < double > >( initialAzimuthAngle_, 16 );

    // Iterate on the free coefficient value until the time of flight matches its required value.
    iterateToMatchRequiredTimeOfFlight( rootFinderSettings_, lowerBoundFreeCoefficient_, upperBoundFreeCoefficient_, initialValueFreeCoefficient_ );


    // Retrieve initial step size.
    double initialStepSize = integratorSettings->initialTimeStep_;

    // Vector of azimuth angles at which the time should be computed.
    Eigen::VectorXd azimuthAnglesToComputeAssociatedEpochs =
            Eigen::VectorXd::LinSpaced( std::ceil( computeNormalizedTimeOfFlight() * physical_constants::JULIAN_YEAR / initialStepSize ),
                                        initialAzimuthAngle_, finalAzimuthAngle_ );

    std::map< double, double > dataToInterpolate;
    for ( int i = 0 ; i < azimuthAnglesToComputeAssociatedEpochs.size() ; i++ )
    {
        dataToInterpolate[ computeCurrentTimeFromAzimuthAngle( azimuthAnglesToComputeAssociatedEpochs[ i ] ) * physical_constants::JULIAN_YEAR ]
                = azimuthAnglesToComputeAssociatedEpochs[ i ];
    }

    // Create interpolator.
    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings =
            std::make_shared< interpolators::LagrangeInterpolatorSettings >( 10 );

    interpolator_ = interpolators::createOneDimensionalInterpolator( dataToInterpolate, interpolatorSettings );

}


//! Convert time to independent variable.
double SphericalShaping::convertTimeToIndependentVariable( const double time )
{
    return  interpolator_->interpolate( time ); ;
}


//! Convert independent variable to time.
double SphericalShaping::convertIndependentVariableToTime( const double independentVariable )
{
    // Define the derivative of the time function w.r.t the azimuth angle theta.
    std::function< double( const double ) > derivativeTimeFunction = [ = ] ( const double currentAzimuthAngle ){

        double scalarFunctionTimeEquation = computeScalarFunctionTimeEquation( currentAzimuthAngle );

        // Check that the trajectory is feasible, ie curved toward the central body.
        if ( scalarFunctionTimeEquation < 0.0 )
        {
            throw std::runtime_error ( "Error, trajectory not curved toward the central body, and thus not feasible." );
        }
        else
        {
            return std::sqrt( computeScalarFunctionTimeEquation( currentAzimuthAngle ) *
                              std::pow( radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle ), 2.0 )
                              / centralBodyGravitationalParameter_ );
        }
    };

    // Create numerical quadrature from quadrature settings.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeTimeFunction, quadratureSettings_, independentVariable );

    double time = quadrature->getQuadrature( );

    return time;
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
        matrixBoundaryConditions( 0, i ) = radialDistanceCompositeFunction_->getComponentFunctionCurrentValue( index, initialAzimuthAngle_ );
        matrixBoundaryConditions( 1, i ) = radialDistanceCompositeFunction_->getComponentFunctionCurrentValue( index, finalAzimuthAngle_ );
        matrixBoundaryConditions( 2, i ) = radialDistanceCompositeFunction_->getComponentFunctionFirstDerivative( index, initialAzimuthAngle_ );
        matrixBoundaryConditions( 3, i ) = radialDistanceCompositeFunction_->getComponentFunctionFirstDerivative( index, finalAzimuthAngle_ );
        matrixBoundaryConditions( 4, i ) = - std::pow( initialStateSphericalCoordinates_[ 0 ], 2.0 )
                * radialDistanceCompositeFunction_->getComponentFunctionSecondDerivative( index, initialAzimuthAngle_ );
        matrixBoundaryConditions( 5, i ) = - std::pow( finalStateSphericalCoordinates_[ 0 ], 2.0 )
                * radialDistanceCompositeFunction_->getComponentFunctionSecondDerivative( index, finalAzimuthAngle_ );
    }

    // Compute value of variable alpha at initial time.
    double initialValueAlpha = computeInitialAlphaValue();

    // Compute value of variable alpha at final time.
    double finalValueAlpha = computeFinalAlphaValue();

    for ( int i = 0 ; i < 4 ; i++ )
    {
        matrixBoundaryConditions( 4, i + 6 ) =
                initialValueAlpha * elevationAngleCompositeFunction_->getComponentFunctionSecondDerivative( i, initialAzimuthAngle_ );
        matrixBoundaryConditions( 5, i + 6 ) =
                finalValueAlpha * elevationAngleCompositeFunction_->getComponentFunctionSecondDerivative( i, finalAzimuthAngle_ );
        matrixBoundaryConditions( 6, i + 6 ) = elevationAngleCompositeFunction_->getComponentFunctionCurrentValue( i, initialAzimuthAngle_ );
        matrixBoundaryConditions( 7, i + 6 ) = elevationAngleCompositeFunction_->getComponentFunctionCurrentValue( i, finalAzimuthAngle_ );
        matrixBoundaryConditions( 8, i + 6 ) = elevationAngleCompositeFunction_->getComponentFunctionFirstDerivative( i, initialAzimuthAngle_ );
        matrixBoundaryConditions( 9, i + 6 ) = elevationAngleCompositeFunction_->getComponentFunctionFirstDerivative( i, finalAzimuthAngle_ );
    }

    // Compute and return the inverse of the boundary conditions matrix.
    return matrixBoundaryConditions.inverse();
}

double SphericalShaping::computeInitialAlphaValue( )
{
    return - ( initialStateParametrizedByAzimuthAngle_[ 3 ] * initialStateParametrizedByAzimuthAngle_[ 5 ] / initialStateParametrizedByAzimuthAngle_[ 0 ] )
            / ( std::pow( initialStateParametrizedByAzimuthAngle_[ 5 ] / initialStateParametrizedByAzimuthAngle_[ 0 ], 2.0 )
            + std::pow( std::cos( initialStateParametrizedByAzimuthAngle_[ 2 ] ), 2.0 ) );
}

double SphericalShaping::computeFinalAlphaValue( )
{
    return - ( finalStateParametrizedByAzimuthAngle_[ 3 ] *  finalStateParametrizedByAzimuthAngle_[ 5 ] / finalStateParametrizedByAzimuthAngle_[ 0 ] )
            / ( std::pow( finalStateParametrizedByAzimuthAngle_[ 5 ] / finalStateParametrizedByAzimuthAngle_[ 0 ], 2.0 )
            + std::pow( std::cos( finalStateParametrizedByAzimuthAngle_[ 2 ] ), 2.0 ) );
}

double SphericalShaping::computeInitialValueBoundariesConstant( )
{
    double radialDistance = initialStateParametrizedByAzimuthAngle_[ 0 ];
    double elevationAngle = initialStateParametrizedByAzimuthAngle_[ 2 ];
    double derivativeRadialDistance = initialStateParametrizedByAzimuthAngle_[ 3 ];
    double derivativeElevationAngle = initialStateParametrizedByAzimuthAngle_[ 5 ] / initialStateParametrizedByAzimuthAngle_[ 0 ];
    double derivativeOfTimeWrtAzimuthAngle = ( initialStateParametrizedByAzimuthAngle_[ 0 ] * std::cos( initialStateParametrizedByAzimuthAngle_[ 2 ] ) )
            / initialStateSphericalCoordinates_[ 4 ];

    return - centralBodyGravitationalParameter_ * std::pow( derivativeOfTimeWrtAzimuthAngle, 2.0 ) / std::pow( radialDistance, 2.0 )
            + 2.0 * std::pow( derivativeRadialDistance, 2.0 ) / radialDistance
            + radialDistance * ( std::pow( derivativeElevationAngle, 2.0 ) + std::pow( std::cos( elevationAngle ), 2.0 ) )
            - derivativeRadialDistance * derivativeElevationAngle * ( std::sin( elevationAngle ) * std::cos( elevationAngle ) )
            / ( std::pow( derivativeElevationAngle, 2.0 ) + std::pow( std::cos( elevationAngle ), 2.0 ) );
}

double SphericalShaping::computeFinalValueBoundariesConstant( )
{
    double radialDistance = finalStateParametrizedByAzimuthAngle_[ 0 ];
    double elevationAngle = finalStateParametrizedByAzimuthAngle_[ 2 ];
    double derivativeRadialDistance = finalStateParametrizedByAzimuthAngle_[ 3 ];
    double derivativeElevationAngle = finalStateParametrizedByAzimuthAngle_[ 5 ] / finalStateParametrizedByAzimuthAngle_[ 0 ];
    double derivativeOfTimeWrtAzimuthAngle = ( finalStateParametrizedByAzimuthAngle_[ 0 ] * std::cos( finalStateParametrizedByAzimuthAngle_[ 2 ] ) )
            / finalStateSphericalCoordinates_[ 4 ];

    return - centralBodyGravitationalParameter_ * std::pow( derivativeOfTimeWrtAzimuthAngle, 2.0 ) / std::pow( radialDistance, 2.0 )
            + 2.0 * std::pow( derivativeRadialDistance, 2.0 ) / radialDistance
            + radialDistance * ( std::pow( derivativeElevationAngle, 2.0 ) + std::pow( std::cos( elevationAngle ), 2.0 ) )
            - derivativeRadialDistance * derivativeElevationAngle * ( std::sin( elevationAngle ) * std::cos( elevationAngle ) )
            / ( std::pow( derivativeElevationAngle, 2.0 ) + std::pow( std::cos( elevationAngle ), 2.0 ) );
}

void SphericalShaping::satisfyBoundaryConditions( )
{
    Eigen::VectorXd vectorBoundaryValues;
    vectorBoundaryValues.resize( 10 );

    vectorBoundaryValues[ 0 ] = 1.0 / initialStateParametrizedByAzimuthAngle_[ 0 ];
    vectorBoundaryValues[ 1 ] = 1.0 / finalStateParametrizedByAzimuthAngle_[ 0 ];
    vectorBoundaryValues[ 2 ] = - initialStateParametrizedByAzimuthAngle_[ 3 ] / std::pow( initialStateParametrizedByAzimuthAngle_[ 0 ], 2.0 );
    vectorBoundaryValues[ 3 ] = - finalStateParametrizedByAzimuthAngle_[ 3 ] / std::pow( finalStateParametrizedByAzimuthAngle_[ 0 ], 2.0 );
    vectorBoundaryValues[ 4 ] = computeInitialValueBoundariesConstant()
            - 2.0 * std::pow( initialStateParametrizedByAzimuthAngle_[ 3 ], 2.0 ) / initialStateParametrizedByAzimuthAngle_[ 0 ];
    vectorBoundaryValues[ 5 ] = computeFinalValueBoundariesConstant()
            - 2.0 * std::pow( finalStateParametrizedByAzimuthAngle_[ 3 ], 2.0 ) / finalStateParametrizedByAzimuthAngle_[ 0 ];
    vectorBoundaryValues[ 6 ] = initialStateParametrizedByAzimuthAngle_[ 2 ];
    vectorBoundaryValues[ 7 ] = finalStateParametrizedByAzimuthAngle_[ 2 ];
    vectorBoundaryValues[ 8 ] = initialStateParametrizedByAzimuthAngle_[ 5 ] / initialStateParametrizedByAzimuthAngle_[ 0 ];
    vectorBoundaryValues[ 9 ] = finalStateParametrizedByAzimuthAngle_[ 5 ] / finalStateParametrizedByAzimuthAngle_[ 0 ];

    Eigen::VectorXd vectorSecondComponentContribution;
    vectorSecondComponentContribution.resize( 10.0 );

    vectorSecondComponentContribution[ 0 ] = radialDistanceCompositeFunction_->getComponentFunctionCurrentValue( 2, initialAzimuthAngle_ );
    vectorSecondComponentContribution[ 1 ] = radialDistanceCompositeFunction_->getComponentFunctionCurrentValue( 2, finalAzimuthAngle_ );
    vectorSecondComponentContribution[ 2 ] = radialDistanceCompositeFunction_->getComponentFunctionFirstDerivative( 2, initialAzimuthAngle_ );
    vectorSecondComponentContribution[ 3 ] = radialDistanceCompositeFunction_->getComponentFunctionFirstDerivative( 2, finalAzimuthAngle_ );
    vectorSecondComponentContribution[ 4 ] = - std::pow( initialStateParametrizedByAzimuthAngle_[ 0 ], 2 )
              * radialDistanceCompositeFunction_->getComponentFunctionSecondDerivative( 2, initialAzimuthAngle_ );
    vectorSecondComponentContribution[ 5 ] = - std::pow( finalStateParametrizedByAzimuthAngle_[ 0 ], 2 )
              * radialDistanceCompositeFunction_->getComponentFunctionSecondDerivative( 2, finalAzimuthAngle_ );
    vectorSecondComponentContribution[ 6 ] = 0.0;
    vectorSecondComponentContribution[ 7 ] = 0.0;
    vectorSecondComponentContribution[ 8 ] = 0.0;
    vectorSecondComponentContribution[ 9 ] = 0.0;

    vectorSecondComponentContribution *= initialValueFreeCoefficient_;

    Eigen::MatrixXd inverseMatrixBoundaryConditions_ = computeInverseMatrixBoundaryConditions();

    Eigen::MatrixXd compositeFunctionCoefficients = inverseMatrixBoundaryConditions_ * ( vectorBoundaryValues - vectorSecondComponentContribution );

    for ( int i = 0 ; i < 6 ; i++ )
    {
        if ( i < 2 )
        {
            coefficientsRadialDistanceFunction_( i ) = compositeFunctionCoefficients( i );
        }
        else
        {
            coefficientsRadialDistanceFunction_( i + 1 ) = compositeFunctionCoefficients( i );
        }
    }
    coefficientsRadialDistanceFunction_( 2 ) = initialValueFreeCoefficient_;

    for ( int i = 0 ; i < 4 ; i++ )
    {
        coefficientsElevationAngleFunction_( i ) = compositeFunctionCoefficients( i + 6 );
    }

    radialDistanceCompositeFunction_->resetCompositeFunctionCoefficients( coefficientsRadialDistanceFunction_ );
    elevationAngleCompositeFunction_->resetCompositeFunctionCoefficients( coefficientsElevationAngleFunction_ );

}



double SphericalShaping::computeScalarFunctionTimeEquation( double currentAzimuthAngle )
{
    double radialFunctionValue = radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );
    double firstDerivativeRadialFunction = radialDistanceCompositeFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthAngle );
    double secondDerivativeRadialFunction = radialDistanceCompositeFunction_->evaluateCompositeFunctionSecondDerivative( currentAzimuthAngle );

    double elevationFunctionValue = elevationAngleCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );
    double firstDerivativeElevationFunction = elevationAngleCompositeFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthAngle );
    double secondDerivativeElevationFunction = elevationAngleCompositeFunction_->evaluateCompositeFunctionSecondDerivative( currentAzimuthAngle );

    return - secondDerivativeRadialFunction + 2.0 * std::pow( firstDerivativeRadialFunction, 2.0 ) / radialFunctionValue
            + firstDerivativeRadialFunction * firstDerivativeElevationFunction
            * ( secondDerivativeElevationFunction - std::sin( elevationFunctionValue ) * std::cos( elevationFunctionValue ) )
            / ( std::pow( firstDerivativeElevationFunction, 2.0 ) + std::pow( std::cos( elevationFunctionValue ), 2.0 ) )
            + radialFunctionValue * ( std::pow( firstDerivativeElevationFunction, 2.0 ) + std::pow( std::cos( elevationFunctionValue ), 2.0 ) );
}

double SphericalShaping::computeDerivativeScalarFunctionTimeEquation( double currentAzimuthAngle )
{
    double radialFunctionValue = radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );
    double firstDerivativeRadialFunction = radialDistanceCompositeFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthAngle );
    double secondDerivativeRadialFunction = radialDistanceCompositeFunction_->evaluateCompositeFunctionSecondDerivative( currentAzimuthAngle );
    double thirdDerivativeRadialFunction = radialDistanceCompositeFunction_->evaluateCompositeFunctionThirdDerivative( currentAzimuthAngle );

    double elevationFunctionValue = elevationAngleCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );
    double firstDerivativeElevationFunction = elevationAngleCompositeFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthAngle );
    double secondDerivativeElevationFunction = elevationAngleCompositeFunction_->evaluateCompositeFunctionSecondDerivative( currentAzimuthAngle );
    double thirdDerivativeElevationFunction = elevationAngleCompositeFunction_->evaluateCompositeFunctionThirdDerivative( currentAzimuthAngle );

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

double SphericalShaping::computeNormalizedTimeOfFlight()
{

    // Define the derivative of the time function w.r.t the azimuth angle theta.
    std::function< double( const double ) > derivativeTimeFunction = [ = ] ( const double currentAzimuthAngle ){

        double scalarFunctionTimeEquation = computeScalarFunctionTimeEquation( currentAzimuthAngle );

        // Check that the trajectory is feasible, ie curved toward the central body.
        if ( scalarFunctionTimeEquation < 0.0 )
        {
            throw std::runtime_error ( "Error, trajectory not curved toward the central body, and thus not feasible." );
        }
        else
        {
            return std::sqrt( computeScalarFunctionTimeEquation( currentAzimuthAngle ) *
                              std::pow( radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle ), 2.0 )
                              / centralBodyGravitationalParameter_ );
        }
    };

    // Create numerical quadrature from quadrature settings.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeTimeFunction, quadratureSettings_, finalAzimuthAngle_ );

    return quadrature->getQuadrature( );
}

//! Compute current time from azimuth angle.
double SphericalShaping::computeCurrentTimeFromAzimuthAngle( const double currentAzimuthAngle )
{
    // Define the derivative of the time function w.r.t the azimuth angle theta.
    std::function< double( const double ) > derivativeTimeFunction = [ = ] ( const double currentAzimuthAngle ){

        double scalarFunctionTimeEquation = computeScalarFunctionTimeEquation( currentAzimuthAngle );

        // Check that the trajectory is feasible, ie curved toward the central body.
        if ( scalarFunctionTimeEquation < 0.0 )
        {
            throw std::runtime_error ( "Error, trajectory not curved toward the central body, and thus not feasible." );
        }
        else
        {
            return std::sqrt( computeScalarFunctionTimeEquation( currentAzimuthAngle ) *
                              std::pow( radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle ), 2.0 )
                              / centralBodyGravitationalParameter_ );
        }
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
    std::function< double ( ) >  computeTOFfunction = std::bind( &SphericalShaping::computeNormalizedTimeOfFlight, this );
    std::function< double ( ) > getRequiredTOFfunction = std::bind( &SphericalShaping::getNormalizedRequiredTimeOfFlight, this );

    std::shared_ptr< basic_mathematics::Function< double, double > > timeOfFlightFunction =
            std::make_shared< SphericalShaping::TimeOfFlightFunction >( resetFreeCoefficientFunction, satisfyBoundaryConditionsFunction, computeTOFfunction, getRequiredTOFfunction );

    // Create root finder from root finder settings.
    std::shared_ptr< root_finders::RootFinderCore< double > > rootFinder = root_finders::createRootFinder( rootFinderSettings, lowerBound, upperBound, initialGuess );

    // Iterate to find the free coefficient value that matches the required time of flight.
    double updatedFreeCoefficient = rootFinder->execute( timeOfFlightFunction, initialGuess );

}


//! Compute current derivative of the azimuth angle.
double SphericalShaping::computeFirstDerivativeAzimuthAngleWrtTime( const double currentAzimuthAngle )
{
    // Compute scalar function time equation.
    double scalarFunctionTimeEquation = computeScalarFunctionTimeEquation( currentAzimuthAngle );

    // Compute current radial distance.
    double radialDistance = radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );

    // Compute and return the first derivative of the azimuth angle w.r.t. time.
    return std::sqrt( centralBodyGravitationalParameter_ / ( scalarFunctionTimeEquation * std::pow( radialDistance, 2.0 ) ) );
}

//! Compute second derivative of the azimuth angle.
double SphericalShaping::computeSecondDerivativeAzimuthAngleWrtTime( const double currentAzimuthAngle )
{
    // Compute first derivative azimuth angle w.r.t. time.
    double firstDerivativeAzimuthAngle = computeFirstDerivativeAzimuthAngleWrtTime( currentAzimuthAngle );

    // Compute scalar function of the time equation, and its derivative w.r.t. azimuth angle.
    double scalarFunctionTimeEquation = computeScalarFunctionTimeEquation( currentAzimuthAngle );
    double derivativeScalarFunctionTimeEquation = computeDerivativeScalarFunctionTimeEquation( currentAzimuthAngle );

    // Compute radial distance, and its derivative w.r.t. azimuth angle.
    double radialDistance = radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );
    double firstDerivativeRadialDistance = radialDistanceCompositeFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthAngle );

    return - std::pow( firstDerivativeAzimuthAngle, 2.0 )
            * ( derivativeScalarFunctionTimeEquation / ( 2.0 * scalarFunctionTimeEquation ) + firstDerivativeRadialDistance / radialDistance );
}


//! Compute current state vector in spherical coordinates.
Eigen::Vector6d SphericalShaping::computeStateVectorInSphericalCoordinates( const double currentAzimuthAngle )
{
    Eigen::Vector6d currentSphericalState;
    currentSphericalState.segment( 0, 3 ) = ( Eigen::Vector3d() << radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle ),
                                              currentAzimuthAngle, elevationAngleCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle ) ).finished();

    // Compute first derivative of the azimuth angle w.r.t. time.
    double derivativeAzimuthAngle = computeFirstDerivativeAzimuthAngleWrtTime( currentAzimuthAngle );

    // Compute and return current velocity vector in spherical coordinates.
    currentSphericalState.segment( 3, 3 ) = derivativeAzimuthAngle * computeVelocityVectorParametrizedByAzimuthAngle( currentAzimuthAngle );

    return currentSphericalState;
}


//! Compute current cartesian state.
Eigen::Vector6d SphericalShaping::computeCurrentStateVector( const double currentAzimuthAngle )
{
    Eigen::Vector6d normalizedStateVector =  coordinate_conversions::convertSphericalToCartesianState( computeStateVectorInSphericalCoordinates( currentAzimuthAngle ) );

    Eigen::Vector6d dimensionalStateVector;
    dimensionalStateVector.segment( 0, 3 ) = normalizedStateVector.segment( 0, 3 )
            * physical_constants::ASTRONOMICAL_UNIT;
    dimensionalStateVector.segment( 3, 3 ) = normalizedStateVector.segment( 3, 3 )
            * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_YEAR;

    return dimensionalStateVector;
}


//! Compute current velocity in spherical coordinates parametrized by azimuth angle theta.
Eigen::Vector3d SphericalShaping::computeVelocityVectorParametrizedByAzimuthAngle( const double currentAzimuthAngle )
{

    // Retrieve current radial distance and elevation angle, as well as their derivatives w.r.t. azimuth angle.
    double radialDistance = radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );
    double derivativeRadialDistance = radialDistanceCompositeFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthAngle );
    double elevationAngle = elevationAngleCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );
    double derivativeElevationAngle = elevationAngleCompositeFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthAngle );

    // Compute and return velocity vector parametrized by azimuth angle.
    return ( Eigen::Vector3d() << derivativeRadialDistance,
             radialDistance * std::cos( elevationAngle ),
             radialDistance * derivativeElevationAngle ).finished();
}


//! Compute current acceleration in spherical coordinates parametrized by azimuth angle theta.
Eigen::Vector3d SphericalShaping::computeThrustAccelerationVectorParametrizedByAzimuthAngle( const double currentAzimuthAngle )
{
    // Retrieve spherical coordinates and their derivatives w.r.t. to the azimuth angle.
    double radialDistance = radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );
    double firstDerivativeRadialDistance = radialDistanceCompositeFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthAngle );
    double secondDerivativeRadialDistance = radialDistanceCompositeFunction_->evaluateCompositeFunctionSecondDerivative( currentAzimuthAngle );
    double elevationAngle = elevationAngleCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );
    double firstDerivativeElevationAngle = elevationAngleCompositeFunction_->evaluateCompositeFunctionFirstDerivative( currentAzimuthAngle );
    double secondDerivativeElevationAngle = elevationAngleCompositeFunction_->evaluateCompositeFunctionSecondDerivative( currentAzimuthAngle );

    // Compute and return acceleration vector parametrized by the azimuth angle theta.
    Eigen::Vector3d accelerationParametrizedByAzimuthAngle;
    accelerationParametrizedByAzimuthAngle[ 0 ] = secondDerivativeRadialDistance
            - radialDistance * ( std::pow( firstDerivativeElevationAngle, 2.0 ) + std::pow( std::cos( elevationAngle ), 2.0 ) );
    accelerationParametrizedByAzimuthAngle[ 1 ] = 2.0 * firstDerivativeRadialDistance * std::cos( elevationAngle )
            - 2.0 * radialDistance * firstDerivativeElevationAngle * std::sin( elevationAngle );
    accelerationParametrizedByAzimuthAngle[ 2 ] = 2.0 * firstDerivativeRadialDistance * firstDerivativeElevationAngle
            + radialDistance * ( secondDerivativeElevationAngle + std::sin( elevationAngle ) * std::cos( elevationAngle ) );

    return accelerationParametrizedByAzimuthAngle;

}

//! Compute thrust acceleration vector in spherical coordinates.
Eigen::Vector3d SphericalShaping::computeThrustAccelerationVectorInSphericalCoordinates( const double currentAzimuthAngle )
{
    // Compute current radial distance.
    double radialDistance = radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle );

    // Compute first and second derivatives of the azimuth angle w.r.t. time.
    double firstDerivativeAzimuthAngleWrtTime = computeFirstDerivativeAzimuthAngleWrtTime( currentAzimuthAngle );
    double secondDerivativeAzimuthAngleWrtTime = computeSecondDerivativeAzimuthAngleWrtTime( currentAzimuthAngle );

    // Compute velocity vector parametrized by azimuth angle theta.
    Eigen::Vector3d velocityParametrizedByAzimuthAngle = computeVelocityVectorParametrizedByAzimuthAngle( currentAzimuthAngle );

    // Compute acceleration vector parametrized by azimuth angle theta.
    Eigen::Vector3d accelerationParametrizedByAzimuthAngle = computeThrustAccelerationVectorParametrizedByAzimuthAngle( currentAzimuthAngle );


    // Compute and return the current thrust acceleration vector in spherical coordinates.
    return std::pow( firstDerivativeAzimuthAngleWrtTime, 2.0 ) * accelerationParametrizedByAzimuthAngle
            + secondDerivativeAzimuthAngleWrtTime * velocityParametrizedByAzimuthAngle
            + centralBodyGravitationalParameter_ / std::pow( radialDistance, 3.0 ) * ( Eigen::Vector3d() << radialDistance, 0.0, 0.0 ).finished();

}

//! Compute current thrust acceleration in cartesian coordinates.
Eigen::Vector3d SphericalShaping::computeThrustAccelerationVector( const double currentAzimuthAngle )
{
    Eigen::Vector6d sphericalStateToBeConverted;
    sphericalStateToBeConverted.segment( 0, 3 ) = ( Eigen::Vector3d() << radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle ),
                                                    currentAzimuthAngle, elevationAngleCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle ) ).finished();
    sphericalStateToBeConverted.segment( 3, 3 ) = computeThrustAccelerationVectorInSphericalCoordinates( currentAzimuthAngle );

    Eigen::Vector3d normalizedThrustAccelerationVector = coordinate_conversions::convertSphericalToCartesianState( sphericalStateToBeConverted ).segment( 3, 3 );

    return normalizedThrustAccelerationVector * physical_constants::ASTRONOMICAL_UNIT
            / std::pow( physical_constants::JULIAN_YEAR, 2.0 );
}

//! Compute magnitude cartesian acceleration.
double  SphericalShaping::computeCurrentThrustAccelerationMagnitude(
        double currentAzimuthAngle, std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    return computeThrustAccelerationVector( currentAzimuthAngle ).norm( );
}

//! Compute direction thrust acceleration in cartesian coordinates.
Eigen::Vector3d SphericalShaping::computeCurrentThrustAccelerationDirection(
        double currentAzimuthAngle, std::function< double ( const double ) > specificImpulseFunction,
        std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings )
{
    return computeThrustAccelerationVector( currentAzimuthAngle ).normalized( );
}

//! Compute final deltaV.
double SphericalShaping::computeDeltaV( )
{
    // Define the derivative of the deltaV, ie thrust acceleration function, as a function of the azimuth angle.
    std::function< double( const double ) > derivativeFunctionDeltaV = [ = ] ( const double currentAzimuthAngle ){

        double thrustAcceleration = computeThrustAccelerationVectorInSphericalCoordinates( currentAzimuthAngle ).norm()
                * std::sqrt( computeScalarFunctionTimeEquation( currentAzimuthAngle )
                             * std::pow( radialDistanceCompositeFunction_->evaluateCompositeFunction( currentAzimuthAngle ), 2.0 )
                             / centralBodyGravitationalParameter_ );

        return thrustAcceleration;

    };

    // Define numerical quadrature from quadratrure settings.
    std::shared_ptr< numerical_quadrature::NumericalQuadrature< double, double > > quadrature =
            numerical_quadrature::createQuadrature( derivativeFunctionDeltaV, quadratureSettings_, finalAzimuthAngle_ );

    // Return dimensional deltaV
    return quadrature->getQuadrature( ) * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_YEAR;
}


} // namespace shape_based_methods
} // namespace tudat
