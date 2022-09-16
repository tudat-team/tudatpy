/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/numericalAccelerationPartial.h"

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
        std::function< void( Eigen::Vector6d ) > setBodyState,
        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        Eigen::Vector6d originalState,
        Eigen::Vector3d statePerturbation,
        int startIndex,
        std::function< void( ) > updateFunction,
        const double evaluationTime )
{
    Eigen::Matrix3d upAccelerations = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d downAccelerations = Eigen::Matrix3d::Zero( );

    Eigen::Vector6d perturbedState = originalState;

    accelerationModel->resetCurrentTime( );

    // Calculate perturbed accelerations for up-perturbed state entries.
    for( int i = 0; i < 3; i++ )
    {
        perturbedState( i + startIndex ) += statePerturbation( i );
        setBodyState( perturbedState );
        updateFunction( );
        upAccelerations.block( 0, i, 3, 1 ) = basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                    accelerationModel, evaluationTime );
        accelerationModel->resetCurrentTime( );
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
        accelerationModel->resetCurrentTime( );
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

Eigen::Vector3d calculateAccelerationWrtMassPartials(
        std::function< void( double ) > setBodyMass,
        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        double originalMass,
        double massPerturbation,
        std::function< void( ) > updateFunction,
        const double evaluationTime )
{
    Eigen::Vector3d upAccelerations = Eigen::Vector3d::Zero( );
    Eigen::Vector3d downAccelerations = Eigen::Vector3d::Zero( );

    double perturbedMass = originalMass;

    accelerationModel->resetCurrentTime( );

    // Calculate perturbed accelerations for up-perturbed state entries.
    perturbedMass += massPerturbation;
    setBodyMass( perturbedMass );
    updateFunction( );
    upAccelerations = basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                accelerationModel, evaluationTime );
    accelerationModel->resetCurrentTime( );
    perturbedMass = originalMass;

    // Calculate perturbed accelerations for down-perturbed state entries.
    perturbedMass -= massPerturbation;
    setBodyMass( perturbedMass );
    updateFunction( );
    downAccelerations = basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                accelerationModel, evaluationTime );
    accelerationModel->resetCurrentTime( );
    perturbedMass = originalMass;


    // Reset state/environment to original state.
    setBodyMass( perturbedMass );
    updateFunction( );

    basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                accelerationModel, evaluationTime );

    // Numerically compute partial derivatives.
    return ( upAccelerations - downAccelerations ) / ( 2.0 * massPerturbation );
}

//! Function to numerical compute the partial derivative of a torque w.r.t. a body translational state.
Eigen::MatrixXd calculateTorqueWrtTranslationalStatePartials(
        std::function< void( Eigen::Vector6d ) > setBodyState,
        std::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel,
        Eigen::Vector6d originalState,
        Eigen::Vector3d statePerturbation,
        int startIndex,
        std::function< void( ) > updateFunction,
        const double evaluationTime )
{
    Eigen::Matrix3d upTorques = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d downTorques = Eigen::Matrix3d::Zero( );

    Eigen::Vector6d perturbedState = originalState;

    torqueModel->resetCurrentTime( );

    // Calculate perturbed torques for up-perturbed state entries.
    for( int i = 0; i < 3; i++ )
    {
        perturbedState( i + startIndex ) += statePerturbation( i );
        setBodyState( perturbedState );
        updateFunction( );
        upTorques.block( 0, i, 3, 1 ) = basic_astrodynamics::updateAndGetTorque(
                    torqueModel, evaluationTime );
        torqueModel->resetCurrentTime( );
        perturbedState = originalState;
    }

    // Calculate perturbed torques for down-perturbed state entries.
    for( int i = 0; i < 3; i++ )
    {
        perturbedState( i + startIndex ) -= statePerturbation( i );
        setBodyState( perturbedState );
        updateFunction( );
        downTorques.block( 0, i, 3, 1 ) = basic_astrodynamics::updateAndGetTorque(
                    torqueModel, evaluationTime );
        torqueModel->resetCurrentTime( );
        perturbedState = originalState;
    }


    // Reset state/environment to original state.
    setBodyState( perturbedState );
    updateFunction( );

    basic_astrodynamics::updateAndGetTorque( torqueModel, evaluationTime );

    // Numerically compute partial derivatives.
    Eigen::Matrix3d torquePartials = upTorques - downTorques;
    for( int i = 0; i < 3; i++ )
    {
        torquePartials.block( 0, i, 3, 1 ) /= ( 2.0 * statePerturbation( i ) );
    }

    return torquePartials;
}

//! Function to numerical compute the partial derivative of a torque w.r.t. a body rotational state.
Eigen::MatrixXd calculateTorqueWrtRotationalStatePartials(
        std::function< void( Eigen::Vector7d ) > setBodyRotationalState,
        std::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel,
        Eigen::Vector7d originalRotationalState,
        Eigen::VectorXd statePerturbations,
        int startIndex,
        int numberOfEntries,
        std::function< void( ) > updateFunction,
        const double evaluationTime )
{
    Eigen::MatrixXd upTorques = Eigen::MatrixXd::Zero( 3, numberOfEntries );
    Eigen::MatrixXd downTorques = Eigen::MatrixXd::Zero( 3, numberOfEntries );

    Eigen::Vector7d perturbedState = originalRotationalState;

    torqueModel->resetCurrentTime( );

    // Calculate perturbed accelerations for up-perturbed state entries.
    for( int i = 0; i < numberOfEntries; i++ )
    {
        perturbedState( i + startIndex ) += statePerturbations( i );
        setBodyRotationalState( perturbedState );
        updateFunction( );
        upTorques.block( 0, i, 3, 1 ) = basic_astrodynamics::updateAndGetTorque(
                    torqueModel, evaluationTime );
        torqueModel->resetCurrentTime( );
        perturbedState = originalRotationalState;
    }

    // Calculate perturbed accelerations for down-perturbed state entries.
    for( int i = 0; i < numberOfEntries; i++ )
    {
        perturbedState( i + startIndex ) -= statePerturbations( i );
        setBodyRotationalState( perturbedState );
        updateFunction( );
        downTorques.block( 0, i, 3, 1 ) = basic_astrodynamics::updateAndGetTorque(
                    torqueModel, evaluationTime );
        torqueModel->resetCurrentTime( );
        perturbedState = originalRotationalState;
    }


    // Reset state/environment to original state.
    setBodyRotationalState( perturbedState );
    updateFunction( );

    basic_astrodynamics::updateAndGetTorque(
                torqueModel, evaluationTime );

    // Numerically compute partial derivatives.
    Eigen::MatrixXd accelerationPartials = upTorques - downTorques;
    for( int i = 0; i < numberOfEntries; i++ )
    {
        accelerationPartials.block( 0, i, 3, 1 ) /= ( 2.0 * statePerturbations( i ) );
    }

    return accelerationPartials;
}

//! Function to numerical compute the partial derivative of a acceleration w.r.t. a body rotational quaternion.
Eigen::MatrixXd calculateAccelerationDeviationDueToOrientationChange(
        const std::function< void( Eigen::Vector7d ) > setBodyRotationalState,
        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        const Eigen::Vector7d& originalRotationalState,
        const Eigen::Vector4d& commandedQuaternionPerturbation,
        std::vector< Eigen::Vector4d >& appliedQuaternionPerturbation,
        std::function< void( ) > updateFunction,
        const double evaluationTime )
{
    appliedQuaternionPerturbation.resize( 4 );

    setBodyRotationalState( originalRotationalState );
    Eigen::Vector3d nominalAcceleration = basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                accelerationModel, evaluationTime );
    Eigen::MatrixXd deviations = Eigen::MatrixXd::Zero( 3, 3 );

    Eigen::Vector7d perturbedState = originalRotationalState;

    accelerationModel->resetCurrentTime( );

    // Calculate perturbed accelerations for up-perturbed state entries.
    for( int i = 1; i < 4; i++ )
    {
        perturbedState( i ) += commandedQuaternionPerturbation( i );
        perturbedState( 0 ) = ( originalRotationalState( 0 ) > 0 ? 1.0 : -1.0 ) * std::sqrt( 1.0 - std::pow( perturbedState.segment( 1, 3 ).norm( ), 2 ) );

        appliedQuaternionPerturbation[ i ] = perturbedState.segment( 0, 4 ).normalized( ) -
                originalRotationalState.segment( 0, 4 );
        setBodyRotationalState( perturbedState );
        updateFunction( );
        deviations.block( 0, i - 1, 3, 1 ) = basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                    accelerationModel, evaluationTime ) -
                nominalAcceleration;
        accelerationModel->resetCurrentTime( );
        perturbedState = originalRotationalState;

    }

    setBodyRotationalState( perturbedState );
    updateFunction( );

    basic_astrodynamics::updateAndGetAcceleration( accelerationModel, evaluationTime );

    return deviations;
}

//! Function to numerical compute the partial derivative of a torque w.r.t. a body rotational quaternion.
Eigen::MatrixXd calculateTorqueDeviationDueToOrientationChange(
        const std::function< void( Eigen::Vector7d ) > setBodyRotationalState,
        const std::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel,
        const Eigen::Vector7d& originalRotationalState,
        const Eigen::Vector4d& commandedQuaternionPerturbation,
        std::vector< Eigen::Vector4d >& appliedQuaternionPerturbation,
        std::function< void( ) > updateFunction,
        const double evaluationTime )
{
    appliedQuaternionPerturbation.resize( 4 );
    Eigen::Vector3d nominalTorque = basic_astrodynamics::updateAndGetTorque( torqueModel, evaluationTime );
    Eigen::MatrixXd deviations = Eigen::MatrixXd::Zero( 3, 3 );

    Eigen::Vector7d perturbedState = originalRotationalState;

    torqueModel->resetCurrentTime( );

    // Calculate perturbed accelerations for up-perturbed state entries.
    for( int i = 1; i < 4; i++ )
    {
        perturbedState( i ) += commandedQuaternionPerturbation( i );
        perturbedState( 0 ) = std::sqrt( 1.0 - std::pow( perturbedState.segment( 1, 3 ).norm( ), 2 ) );

        appliedQuaternionPerturbation[ i ] = perturbedState.segment( 0, 4 ).normalized( ) -
                originalRotationalState.segment( 0, 4 );

        setBodyRotationalState( perturbedState );
        updateFunction( );
        deviations.block( 0, i - 1, 3, 1 ) = basic_astrodynamics::updateAndGetTorque( torqueModel, evaluationTime ) -
                nominalTorque;
        torqueModel->resetCurrentTime( );
        perturbedState = originalRotationalState;

    }

    setBodyRotationalState( perturbedState );
    updateFunction( );

    basic_astrodynamics::updateAndGetTorque( torqueModel, evaluationTime );

    return deviations;
}

//! Function to numerical compute the partial derivative of an acceleration w.r.t. a double parameter
Eigen::Vector3d calculateAccelerationWrtParameterPartials(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        double parameterPerturbation,
        std::function< void( ) > updateDependentVariables,
        const double currentTime,
        std::function< void( const double ) > timeDependentUpdateDependentVariables )
{
    // Store uperturbed value.
    double unperturbedParameterValue = parameter->getParameterValue( );

    // Calculate up-perturbation
    parameter->setParameterValue(
                unperturbedParameterValue + parameterPerturbation );
    updateDependentVariables( );
    timeDependentUpdateDependentVariables( currentTime );
    accelerationModel->resetCurrentTime( );

    Eigen::Vector3d upPerturbedAcceleration = basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                accelerationModel, currentTime );
    accelerationModel->resetCurrentTime( );
    // Calculate down-perturbation.
    parameter->setParameterValue(
                unperturbedParameterValue - parameterPerturbation );
    updateDependentVariables( );
    timeDependentUpdateDependentVariables( currentTime );

    Eigen::Vector3d downPerturbedAcceleration = basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                accelerationModel, currentTime );

    accelerationModel->resetCurrentTime( );

    // Reset to original value.
    parameter->setParameterValue(
                unperturbedParameterValue ) ;
    updateDependentVariables( );
    timeDependentUpdateDependentVariables( currentTime );

    accelerationModel->resetCurrentTime( );
    basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                accelerationModel, currentTime );

    // Calculate partial using central difference.
    return ( upPerturbedAcceleration - downPerturbedAcceleration ) / ( 2.0 * parameterPerturbation );

}

//! Function to numerical compute the partial derivative of an torque w.r.t. a double parameter
Eigen::Vector3d calculateTorqueWrtParameterPartials(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        std::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel,
        double parameterPerturbation,
        std::function< void( ) > updateDependentVariables,
        const double currentTime,
        std::function< void( const double ) > timeDependentUpdateDependentVariables )
{
    // Store uperturbed value.
    double unperturbedParameterValue = parameter->getParameterValue( );

    // Calculate up-perturbation
    parameter->setParameterValue(
                unperturbedParameterValue + parameterPerturbation );
    updateDependentVariables( );
    timeDependentUpdateDependentVariables( currentTime );
    torqueModel->resetCurrentTime( );

    Eigen::Vector3d upPerturbedTorque = basic_astrodynamics::updateAndGetTorque(
                torqueModel, currentTime );
    torqueModel->resetCurrentTime( );

    // Calculate down-perturbation.
    parameter->setParameterValue(
                unperturbedParameterValue - parameterPerturbation );
    updateDependentVariables( );
    timeDependentUpdateDependentVariables( currentTime );

    Eigen::Vector3d downPerturbedTorque = basic_astrodynamics::updateAndGetTorque(
                torqueModel, currentTime );

    torqueModel->resetCurrentTime( );

    // Reset to original value.
    parameter->setParameterValue(
                unperturbedParameterValue ) ;
    updateDependentVariables( );
    timeDependentUpdateDependentVariables( currentTime );

    torqueModel->resetCurrentTime( );
    basic_astrodynamics::updateAndGetTorque(
                torqueModel, currentTime );

    // Calculate partial using central difference.
    return ( upPerturbedTorque - downPerturbedTorque ) / ( 2.0 * parameterPerturbation );

}


//! Function to numerical compute the partial derivative of an acceleration w.r.t. a double parameter
double calculateMassRateWrtParameterPartials(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
        std::shared_ptr< basic_astrodynamics::MassRateModel > massRateModel,
        double parameterPerturbation,
        std::function< void( ) > updateDependentVariables,
        const double currentTime,
        std::function< void( const double ) > timeDependentUpdateDependentVariables )
{
    // Store uperturbed value.
    double unperturbedParameterValue = parameter->getParameterValue( );

    // Calculate up-perturbation
    parameter->setParameterValue(
                unperturbedParameterValue + parameterPerturbation );
    updateDependentVariables( );
    timeDependentUpdateDependentVariables( currentTime );
    massRateModel->resetCurrentTime( );

    massRateModel->updateMembers( currentTime );
    double upPerturbedMassRate = massRateModel->getMassRate( );
    massRateModel->resetCurrentTime( );

    // Calculate down-perturbation.
    parameter->setParameterValue(
                unperturbedParameterValue - parameterPerturbation );
    updateDependentVariables( );
    timeDependentUpdateDependentVariables( currentTime );

    massRateModel->updateMembers( currentTime );
    double downPerturbedMassRate = massRateModel->getMassRate( );
    massRateModel->resetCurrentTime( );

    // Reset to original value.
    parameter->setParameterValue(
                unperturbedParameterValue ) ;
    updateDependentVariables( );
    timeDependentUpdateDependentVariables( currentTime );

    massRateModel->resetCurrentTime( );
    massRateModel->updateMembers( currentTime );

    // Calculate partial using central difference.
    return ( upPerturbedMassRate - downPerturbedMassRate ) / ( 2.0 * parameterPerturbation );

}


//! Function to numerical compute the partial derivative of an acceleration w.r.t. a vector parameter
Eigen::Matrix< double, 3, Eigen::Dynamic > calculateAccelerationWrtParameterPartials(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        Eigen::VectorXd parameterPerturbation,
        std::function< void( ) > updateDependentVariables,
        const double currentTime,
        std::function< void( const double ) > timeDependentUpdateDependentVariables )
{
    // Store uperturbed value.

    Eigen::VectorXd unperturbedParameterValue = parameter->getParameterValue( );


    if( unperturbedParameterValue.size( ) != parameterPerturbation.size( ) )
    {
        throw std::runtime_error( "Error when calculating numerical parameter partial of acceleration, parameter and perturbations are not the same size" );
    }

    Eigen::Matrix< double, 3, Eigen::Dynamic > partialMatrix = Eigen::MatrixXd::Zero( 3, unperturbedParameterValue.size( ) );

    accelerationModel->resetCurrentTime( );

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
        accelerationModel->resetCurrentTime( );

        // Calculate down-perturbation.
        perturbedParameterValue = unperturbedParameterValue;
        perturbedParameterValue( i ) -= parameterPerturbation( i );
        parameter->setParameterValue( perturbedParameterValue );
        updateDependentVariables( );
        timeDependentUpdateDependentVariables( currentTime );
        Eigen::Vector3d downPerturbedAcceleration = basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                    accelerationModel, currentTime );
        accelerationModel->resetCurrentTime( );

        // Compute partial entry.
        partialMatrix.block( 0, i, 3, 1 ) =
                ( upPerturbedAcceleration - downPerturbedAcceleration ) / ( 2.0 * parameterPerturbation( i ) );

    }

    // Reset to original value.
    parameter->setParameterValue(
                unperturbedParameterValue ) ;
    updateDependentVariables( );
    timeDependentUpdateDependentVariables( currentTime );
    accelerationModel->resetCurrentTime( );

    basic_astrodynamics::updateAndGetAcceleration< Eigen::Vector3d >(
                accelerationModel, currentTime );

    // Calculate partial using central difference.
    return partialMatrix;
}

//! Function to numerical compute the partial derivative of an torque w.r.t. a vector parameter
Eigen::Matrix< double, 3, Eigen::Dynamic > calculateTorqueWrtParameterPartials(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
        std::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel,
        Eigen::VectorXd parameterPerturbation,
        std::function< void( ) > updateDependentVariables,
        const double currentTime,
        std::function< void( const double ) > timeDependentUpdateDependentVariables )
{
    // Store uperturbed value.

    Eigen::VectorXd unperturbedParameterValue = parameter->getParameterValue( );


    if( unperturbedParameterValue.size( ) != parameterPerturbation.size( ) )
    {
        throw std::runtime_error( "Error when calculating numerical parameter partial of torque, parameter and perturbations are not the same size" );
    }

    Eigen::Matrix< double, 3, Eigen::Dynamic > partialMatrix = Eigen::MatrixXd::Zero( 3, unperturbedParameterValue.size( ) );

    torqueModel->resetCurrentTime( );

    Eigen::VectorXd perturbedParameterValue;
    for( int i = 0; i < unperturbedParameterValue.size( ); i++ )
    {
        perturbedParameterValue = unperturbedParameterValue;
        perturbedParameterValue( i ) += parameterPerturbation( i );

        // Calculate up-perturbation
        parameter->setParameterValue( perturbedParameterValue );
        updateDependentVariables( );
        timeDependentUpdateDependentVariables( currentTime );
        Eigen::Vector3d upPerturbedTorque = basic_astrodynamics::updateAndGetTorque(
                    torqueModel, currentTime );
        torqueModel->resetCurrentTime( );

        // Calculate down-perturbation.
        perturbedParameterValue = unperturbedParameterValue;
        perturbedParameterValue( i ) -= parameterPerturbation( i );
        parameter->setParameterValue( perturbedParameterValue );
        updateDependentVariables( );
        timeDependentUpdateDependentVariables( currentTime );
        Eigen::Vector3d downPerturbedTorque = basic_astrodynamics::updateAndGetTorque(
                    torqueModel, currentTime );
        torqueModel->resetCurrentTime( );

        // Compute partial entry.
        partialMatrix.block( 0, i, 3, 1 ) =
                ( upPerturbedTorque - downPerturbedTorque ) / ( 2.0 * parameterPerturbation( i ) );

    }

    // Reset to original value.
    parameter->setParameterValue(
                unperturbedParameterValue ) ;
    updateDependentVariables( );
    timeDependentUpdateDependentVariables( currentTime );
    torqueModel->resetCurrentTime( );

    basic_astrodynamics::updateAndGetTorque(
                torqueModel, currentTime );

    // Calculate partial using central difference.
    return partialMatrix;
}


}

}
