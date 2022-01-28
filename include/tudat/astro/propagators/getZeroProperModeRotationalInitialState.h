/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GETZEROPROPERMODEROTATIONALINITIALSTATE_H
#define TUDAT_GETZEROPROPERMODEROTATIONALINITIALSTATE_H

#include <Eigen/Core>

#include "tudat/astro/basic_astro/dissipativeTorqueModel.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/io/basicInputOutput.h"

namespace tudat
{

namespace propagators
{

//! Function to obtain the synthetic dissipation matrix from the damping time
/*!
 * Function to obtain the synthetic dissipation matrix from the damping time
 * \param dampingTime Damping time
 * \param inertiaTensor Inertia Tensor
 * \return Dissipation Matrix
 */
inline Eigen::Matrix3d getDissipationMatrix(
        const double dampingTime,
        const Eigen::Matrix3d& inertiaTensor )
{
    return inertiaTensor / dampingTime;
}

//! Function to integrate forward in time with synthetic dissipation, and subsequently backwards without dissipation
/*!
 *  Function to integrate forward in time with synthetic dissipation, and subsequently backwards without dissipation. This
 *  function is used to find the initial rotational state with a damped free mode.
 * \brief integrateForwardWithDissipationAndBackwardsWithout
 * \param dynamicsSimulator Object to numerically propagate the dynamics
 * \param dissipativeTorque Synthetic dissipative torque model
 * \param propagatedStates Pair of numerically propagated states (first: forward propagation; second: backwards propagation).
 * Returned by reference.
 * \param dependentVariables Pair of numerically computed dependent variables (first: forward propagation; second:
 * backwards propagation). Returned by reference.
 */
template< typename StateScalarType = double, typename TimeType = double >
void integrateForwardWithDissipationAndBackwardsWithout(
        const std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > dynamicsSimulator,
        const std::shared_ptr< basic_astrodynamics::DissipativeTorqueModel > dissipativeTorque,
        std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1> >,
        std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >& propagatedStates,
        std::pair< std::map< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1> >,
        std::map< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1 > > >& dependentVariables )
{
    using namespace tudat::numerical_integrators;

    std::shared_ptr< SingleArcPropagatorSettings< TimeType > > propagatorSettings = dynamicsSimulator->getPropagatorSettings( );
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1> initialState = propagatorSettings->getInitialStates( );
    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings =
            dynamicsSimulator->getIntegratorSettings( );

    // Integrate forward with dissipation and retrieve results
    dynamicsSimulator->integrateEquationsOfMotion( initialState );
    std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > forwardIntegrated =
            dynamicsSimulator->getEquationsOfMotionNumericalSolution( );
    std::map< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1 > > forwardIntegratedDependent =
            dynamicsSimulator->getDependentVariableHistory( );

    // Save original integration/propagation settings
    TimeType originalStartTime = integratorSettings->initialTime_;
    TimeType originalEndTime = std::dynamic_pointer_cast< PropagationTimeTerminationSettings >(
                propagatorSettings->getTerminationSettings( ) )->terminationTime_;
    double originalTimeStep = integratorSettings->initialTimeStep_;

    // Turn off dissipation
    dissipativeTorque->setDampingMatrixFunction( Eigen::Matrix3d::Zero( ) );

    // Reset propagation/integration settings for backwards propagation
    auto outputMapIterator = forwardIntegrated.rbegin( );
    integratorSettings->initialTime_ = outputMapIterator->first;
    dynamicsSimulator->resetInitialPropagationTime( outputMapIterator->first );
    propagatorSettings->resetTerminationSettings(
                std::make_shared< PropagationTimeTerminationSettings >( originalStartTime ) );
    integratorSettings->initialTimeStep_ = -originalTimeStep;

    // Integrate backward without dissipation and retrieve results
    dynamicsSimulator->integrateEquationsOfMotion(
                outputMapIterator->second );
    std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1> > backwardIntegrated =
            dynamicsSimulator->getEquationsOfMotionNumericalSolution( );
    std::map< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1 > > backwardIntegratedDependent =
            dynamicsSimulator->getDependentVariableHistory( );

    // Reset original integration/propagation settings
    integratorSettings->initialTime_ = originalStartTime;
    dynamicsSimulator->resetInitialPropagationTime( originalStartTime );
    propagatorSettings->resetTerminationSettings(
                std::make_shared< PropagationTimeTerminationSettings >( originalEndTime ) );
    integratorSettings->initialTimeStep_ = originalTimeStep;

    // Save states and dependent variables for both forward and backward integration
    propagatedStates = std::make_pair( forwardIntegrated, backwardIntegrated );
    dependentVariables  = std::make_pair( forwardIntegratedDependent, backwardIntegratedDependent );
}

//! Function to determine the initial rotational state for which the free mode is fully damped
/*!
 *  Function to determine the initial rotational state for which the free mode is fully damped, using the approach of
 *  Rambaux et al. (2012). Using this method, the dynamics is propagated forward in time with an artificial damping, and backwards
 *  in time without this damping. The process is repeated for an ever increasing value of the dissipation time. Using this
 *  function, the numerical state and dependent variables at each of the iterations is returned (by reference) by the function.
 *  Additionally, the results can be written to files after each iteration.
 *  \param bodies Map of bodies in the propagation, with keys the names of the bodies.
 *  \param integratorSettings Settings for the numerical integration scheme
 *  \param propagatorSettings Settings for the propagator, must include rotational dynamics of only a single body, but may also
 *  include translational dynamics.
 *  \param bodyMeanRotationRate Mean rotation rate of body for which the initial state is to be determined, used in the
 *  computation of the dissipation matrix.
 *  \param dissipationTimes List of characteristic times for dissipation matrix (in order of which they will be used)
 *  \param propagatedStates Numerical states for each iteration
 *  \param dependentVariables Dependent variables for each iteration
 *  \param propagateNominal Boolean designating whether the undamped dynamics is also to be propagated (saved in index 0)
 *  \param writeToFileInLoop Boolean designating whether the results are written to files after each iteration. NOTE: setting this
 *  variable to true will result in propagatedStates and dependentVariables being returned empty (the data being stored in
 *  files instead)
 *  \param baseFileName File name prefix used for files, if writeToFileInLoop is true
 *  \param outputFolder Directory to which files will be written, if writeToFileInLoop is true
 */
template< typename TimeType, typename StateScalarType >
Eigen::VectorXd getZeroProperModeRotationalState(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings,
        const double bodyMeanRotationRate,
        const std::vector< double > dissipationTimes,
        std::vector< std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >,
        std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >& propagatedStates,
        std::vector< std::pair< std::map< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1 > >,
        std::map< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1 > > > >& dependentVariables,
        const bool propagateNominal = true,
        const bool writeToFileInLoop = false,
        const std::string baseFileName = "",
        const std::string outputFolder = paths::get_resources_path() + "/outputFolder/"  )
{
    propagatedStates.resize( dissipationTimes.size( ) + 1 );
    dependentVariables.resize( dissipationTimes.size( ) + 1 );

    basic_astrodynamics::TorqueModelMap torqueModelMap;

    // Find rotational propagator settings
    std::shared_ptr< RotationalStatePropagatorSettings< StateScalarType > > rotationPropagationSettings_ =
            std::dynamic_pointer_cast< RotationalStatePropagatorSettings< StateScalarType > >(
                propagatorSettings );
    if( rotationPropagationSettings_ == nullptr )
    {
        // Retrieve rotational propagator settings from multi-type settings
        std::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
                std::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( propagatorSettings );
        if( multiTypePropagatorSettings != nullptr )
        {
            // Retrieve rotational propagator settings
            std::map< IntegratedStateType, std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > >
                    propagatorSettingsMap = multiTypePropagatorSettings->propagatorSettingsMap_;
            if( propagatorSettingsMap.count( rotational_state ) == 0 )
            {
                throw std::runtime_error( "Error when finding initial rotational state, no rotational dynamics in multi-type list" );
            }
            else if( propagatorSettingsMap.at( rotational_state ).size( ) != 1 )
            {
                throw std::runtime_error( "Error when finding initial rotational state, multiple rotational dynamics in multi-type list" );
            }
            else
            {
                rotationPropagationSettings_ = std::dynamic_pointer_cast< RotationalStatePropagatorSettings< StateScalarType > >(
                            propagatorSettingsMap.at( rotational_state ).at( 0 ) );
            }
        }
        else
        {
            throw std::runtime_error( "Error when finding initial rotational state, no propagator settings found" );
        }
    }

    // Get torque map
    torqueModelMap = rotationPropagationSettings_->getTorqueModelsMap( );
    if( torqueModelMap.size( ) != 1 )
    {
        std::cerr<<"Error when finding initial rotational state, "<<torqueModelMap.size( )<<" bodies are propagated."<<std::endl;
    }

    //Retrieve body inertia tensor and create synthetic dissipation model
    Eigen::Matrix3d inertiaTensor = bodies.at( torqueModelMap.begin( )->first )->getBodyInertiaTensor( );
    std::shared_ptr< basic_astrodynamics::DissipativeTorqueModel > dissipativeTorque =
            std::make_shared< basic_astrodynamics::DissipativeTorqueModel >(
                std::bind( &simulation_setup::Body::getCurrentAngularVelocityVectorInLocalFrame,
                           bodies.at( torqueModelMap.begin( )->first ) ),
                [ = ]( ){ return Eigen::Matrix3d::Zero( ); },
                bodyMeanRotationRate );
    torqueModelMap[ torqueModelMap.begin( )->first ][ torqueModelMap.begin( )->first ].push_back(
                dissipativeTorque );
    rotationPropagationSettings_->resetTorqueModelsMap( torqueModelMap );

    // Create object to propagate dynamics
    std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > dynamicsSimulator =
            std::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                bodies, integratorSettings, propagatorSettings, 0, 0, 0 );

    // Propagate nominal dynamics forward and backward (without dissipation) if required
    if( propagateNominal )
    {
        integrateForwardWithDissipationAndBackwardsWithout< StateScalarType, TimeType >(
                    dynamicsSimulator, dissipativeTorque, propagatedStates.at( 0 ), dependentVariables.at( 0 ) );

        // Write data to files if required
        int i = 0;
        if( writeToFileInLoop )
        {
            input_output::writeDataMapToTextFile(
                        dependentVariables.at( i ).second,
                        baseFileName + "_dependent_damped_" + std::to_string( i ) + ".dat", outputFolder );
            input_output::writeDataMapToTextFile(
                        dependentVariables.at( i ).first,
                        baseFileName + "_dependent_damped_forward_" + std::to_string( i ) + ".dat", outputFolder );
            input_output::writeDataMapToTextFile(
                        propagatedStates.at( i ).second,
                        baseFileName + "_states_damped_" + std::to_string( i ) + ".dat", outputFolder );
            input_output::writeDataMapToTextFile(
                        propagatedStates.at( i ).first,
                        baseFileName + "_states_damped_forward_" + std::to_string( i ) + ".dat", outputFolder );
            propagatedStates[ 0 ].first.clear( );
            dependentVariables[ 0 ].first.clear( );

            propagatedStates[ 0 ].second.clear( );
            dependentVariables[ 0 ].second.clear( );
        }

    }

    // Retrieve original initial state
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentInitialState =
            propagatorSettings->getInitialStates( );

    double newFinalTime;
    for( unsigned int i = 0; i < dissipationTimes.size( ); i++ )
    {
        std::cout<<"Getting zero proper mode, iteration "<<i<<std::endl;

        // Set damping for current iteration
        dissipativeTorque->setDampingMatrixFunction(
                    getDissipationMatrix( dissipationTimes.at( i ), inertiaTensor ) );

        // Reset final time (propagate for 10 times the dissipation time)
        newFinalTime = integratorSettings->initialTime_ + 10.0 * dissipationTimes.at( i );
        propagatorSettings->resetTerminationSettings(
                    std::make_shared< PropagationTimeTerminationSettings >( newFinalTime ) );

        // Propagate forward (with dissipation) and backward (without dissipation)
        integrateForwardWithDissipationAndBackwardsWithout< StateScalarType, TimeType >(
                    dynamicsSimulator, dissipativeTorque, propagatedStates.at( i + 1 ), dependentVariables.at( i + 1 ) );

        // Update initial state to current damped result
        currentInitialState = propagatedStates.at( i + 1 ).second.begin( )->second;
        propagatorSettings->resetInitialStates( currentInitialState );

        // Write data to files if required
        if( writeToFileInLoop )
        {
            input_output::writeDataMapToTextFile(
                        dependentVariables.at(i + 1).second,
                        baseFileName + "_dependent_damped_" + std::to_string(i + 1) + ".dat", outputFolder );
            input_output::writeDataMapToTextFile(
                        dependentVariables.at(i + 1).first,
                        baseFileName + "_dependent_damped_forward_" + std::to_string(i + 1) + ".dat", outputFolder );
            input_output::writeDataMapToTextFile(
                        propagatedStates.at(i + 1).second,
                        baseFileName + "_states_damped_" + std::to_string(i + 1) + ".dat", outputFolder );
            input_output::writeDataMapToTextFile(
                        propagatedStates.at(i + 1).first,
                        baseFileName + "_states_damped_forward_" + std::to_string(i + 1) + ".dat", outputFolder );
            propagatedStates[ i + 1 ].first.clear( );
            dependentVariables[ i + 1 ].first.clear( );

            propagatedStates[ i + 1 ].second.clear( );
            dependentVariables[ i + 1 ].second.clear( );
        }


    }

    // Return damped initial state that is computed
    return currentInitialState;
}

template< typename TimeType = double, typename StateScalarType = double >
std::tuple<
Eigen::VectorXd,
std::vector< std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >,
std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >,
std::vector< std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >,
std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > >
> getZeroProperModeRotationalState(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings,
        const double bodyMeanRotationRate,
        const std::vector< double > dissipationTimes,
        const bool propagateNominal = true )
{
    std::vector< std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >,
    std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > > propagatedStates;
    std::vector< std::pair< std::map< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1 > >,
    std::map< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1 > > > > dependentVariables;

    Eigen::VectorXd initialState =
            getZeroProperModeRotationalState(
                bodies, integratorSettings, propagatorSettings, bodyMeanRotationRate, dissipationTimes,
                propagatedStates, dependentVariables,
                propagateNominal, false);

    return std::make_tuple( initialState, propagatedStates, dependentVariables );
}

//! Function to determine the initial rotational state for which the free mode is fully damped
/*!
 *  Function to determine the initial rotational state for which the free mode is fully damped, using the approach of
 *  Rambaux et al. (2012). Using this method, the dynamics is propagated forward in time with an artificial damping, and backwards
 *  in time without this damping. The process is repeated for an ever increasing value of the dissipation time
 *  \param bodies Map of bodies in the propagation, with keys the names of the bodies.
 *  \param integratorSettings Settings for the numerical integration scheme
 *  \param propagatorSettings Settings for the propagator, must include rotational dynamics of only a single body, but may also
 *  include translational dynamics.
 *  \param bodyMeanRotationRate Mean rotation rate of body for which the initial state is to be determined, used in the
 *  computation of the dissipation matrix.
 *  \param dissipationTimes List of characteristic times for dissipation matrix (in order of which they will be used)
 *
 */
template< typename TimeType = double, typename StateScalarType = double >
Eigen::VectorXd getZeroProperModeRotationalState(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > propagatorSettings,
        const double bodyMeanRotationRate,
        const std::vector< double > dissipationTimes )
        //const bool propagateNominal = true )
{

    std::vector< std::pair< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >,
            std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > > propagatedStates;
    std::vector< std::pair< std::map< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1 > >,
            std::map< TimeType, Eigen::Matrix< double, Eigen::Dynamic, 1 > > > > dependentVariables;
    return getZeroProperModeRotationalState( bodies, integratorSettings, propagatorSettings, bodyMeanRotationRate, dissipationTimes,
                                             propagatedStates, dependentVariables );//, propagateNominal );
}

}

}

#endif // TUDAT_GETZEROPROPERMODEROTATIONALINITIALSTATE_H
