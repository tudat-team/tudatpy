/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/numericalAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Dummy function used for update, performs no calculations.
void emptyFunction( ){ }

//! Dummy function used for update, performs no calculations.
void emptyTimeFunction( const double time ){ }

//! Function to numerical compute the partial derivative of an acceleration w.r.t. a body state.
Eigen::Matrix3d calculateAccelerationWrtStatePartials(
        boost::function< void( Eigen::Vector6d ) > setBodyState,
        boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        Eigen::Vector6d originalState,
        Eigen::Vector3d statePerturbation,
        int startIndex,
        boost::function< void( ) > updateFunction,
        const double evaluationTime )
{
    Eigen::Matrix3d upAccelerations = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d downAccelerations = Eigen::Matrix3d::Zero( );

    Eigen::Vector6d perturbedState = originalState;

    accelerationModel->resetTime( TUDAT_NAN );

    // Calculate perturbed accelerations for up-perturbed state entries.
    for( int i = 0; i < 3; i++ )
    {
        perturbedState( i + startIndex ) += statePerturbation( i );
        setBodyState( perturbedState );
        updateFunction( );
        upAccelerations.block( 0, i, 3, 1 ) = basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                    accelerationModel, evaluationTime );
        accelerationModel->resetTime( TUDAT_NAN );
        perturbedState = originalState;
    }

    // Calculate perturbed accelerations for down-perturbed state entries.
    for( int i = 0; i < 3; i++ )
    {
        perturbedState( i + startIndex ) -= statePerturbation( i );
        setBodyState( perturbedState );
        updateFunction( );
        downAccelerations.block( 0, i, 3, 1 ) = basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                    accelerationModel, evaluationTime );
        accelerationModel->resetTime( TUDAT_NAN );
        perturbedState = originalState;
    }


    // Reset state/environment to original state.
    setBodyState( perturbedState );
    updateFunction( );

    basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                        accelerationModel, evaluationTime );

    // Numerically compute partial derivatives.
    Eigen::Matrix3d accelerationPartials = upAccelerations - downAccelerations;
    for( int i = 0; i < 3; i++ )
    {
        accelerationPartials.block( 0, i, 3, 1 ) /= ( 2.0 * statePerturbation( i ) );
    }

    return accelerationPartials;
}

//! Function to numerical compute the partial derivative of an acceleration w.r.t. a double parameter
Eigen::Vector3d calculateAccelerationWrtParameterPartials(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        double parameterPerturbation,
        boost::function< void( ) > updateDependentVariables,
        const double currentTime,
        boost::function< void( const double ) > timeDependentUpdateDependentVariables )
{
    // Store uperturbed value.
    double unperturbedParameterValue = parameter->getParameterValue( );

    // Calculate up-perturbation
    parameter->setParameterValue(
                unperturbedParameterValue + parameterPerturbation );
    updateDependentVariables( );
    timeDependentUpdateDependentVariables( currentTime );
    accelerationModel->resetTime( TUDAT_NAN );

    Eigen::Vector3d upPerturbedAcceleration = basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                accelerationModel, currentTime );
    accelerationModel->resetTime( TUDAT_NAN );
    // Calculate down-perturbation.
    parameter->setParameterValue(
                unperturbedParameterValue - parameterPerturbation );
    updateDependentVariables( );
    timeDependentUpdateDependentVariables( currentTime );

    Eigen::Vector3d downPerturbedAcceleration = basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                accelerationModel, currentTime );

    accelerationModel->resetTime( TUDAT_NAN );

    // Reset to original value.
    parameter->setParameterValue(
                unperturbedParameterValue ) ;
    updateDependentVariables( );
    timeDependentUpdateDependentVariables( currentTime );

    accelerationModel->resetTime( TUDAT_NAN );
    basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                    accelerationModel, currentTime );

    // Calculate partial using central difference.
    return ( upPerturbedAcceleration - downPerturbedAcceleration ) / ( 2.0 * parameterPerturbation );

}

//! Function to numerical compute the partial derivative of an acceleration w.r.t. a vector parameter
Eigen::Matrix< double, 3, Eigen::Dynamic > calculateAccelerationWrtParameterPartials(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        Eigen::VectorXd parameterPerturbation,
        boost::function< void( ) > updateDependentVariables,
        const double currentTime,
        boost::function< void( const double ) > timeDependentUpdateDependentVariables )
{
    // Store uperturbed value.

    Eigen::VectorXd unperturbedParameterValue = parameter->getParameterValue( );


    if( unperturbedParameterValue.size( ) != parameterPerturbation.size( ) )
    {
        throw std::runtime_error( "Error when calculating numerical parameter partial of acceleration, parameter and perturbations are not the same size" );
    }

    Eigen::Matrix< double, 3, Eigen::Dynamic > partialMatrix = Eigen::MatrixXd::Zero( 3, unperturbedParameterValue.size( ) );

    accelerationModel->resetTime( TUDAT_NAN );

    Eigen::VectorXd perturbedParameterValue;
    for( int i = 0; i < unperturbedParameterValue.size( ); i++ )
    {
        perturbedParameterValue = unperturbedParameterValue;
        perturbedParameterValue( i ) += parameterPerturbation( i );

        // Calculate up-perturbation
        parameter->setParameterValue( perturbedParameterValue );
        updateDependentVariables( );
        timeDependentUpdateDependentVariables( currentTime );
        Eigen::Vector3d upPerturbedAcceleration = basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                    accelerationModel, currentTime );
        accelerationModel->resetTime( TUDAT_NAN );

        // Calculate down-perturbation.
        perturbedParameterValue = unperturbedParameterValue;
        perturbedParameterValue( i ) -= parameterPerturbation( i );
        parameter->setParameterValue( perturbedParameterValue );
        updateDependentVariables( );
        timeDependentUpdateDependentVariables( currentTime );
        Eigen::Vector3d downPerturbedAcceleration = basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                    accelerationModel, currentTime );
        accelerationModel->resetTime( TUDAT_NAN );

        // Compute partial entry.
        partialMatrix.block( 0, i, 3, 1 ) =
                ( upPerturbedAcceleration - downPerturbedAcceleration ) / ( 2.0 * parameterPerturbation( i ) );

    }

    // Reset to original value.
    parameter->setParameterValue(
                unperturbedParameterValue ) ;
    updateDependentVariables( );
    timeDependentUpdateDependentVariables( currentTime );
    accelerationModel->resetTime( TUDAT_NAN );

    basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                        accelerationModel, currentTime );

    // Calculate partial using central difference.
    return partialMatrix;
}


}

}
