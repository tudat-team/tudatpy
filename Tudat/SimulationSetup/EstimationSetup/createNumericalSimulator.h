/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef TUDAT_CREATENUMERICALSIMULATOR_H
#define TUDAT_CREATENUMERICALSIMULATOR_H

#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"
#include "Tudat/SimulationSetup/EstimationSetup/variationalEquationsSolver.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create single-arc dynamics simulator object
/*!
 *  Function to create single-arc dynamics simulator object
 *  \param bodyMap Map of bodies (with names) of all bodies in integration.
 *  \param integratorSettings Settings for numerical integrator.
 *  \param propagatorSettings Settings for propagator.
 *  \param areEquationsOfMotionToBeIntegrated Boolean to denote whether equations of motion should be integrated
 *  immediately at the end of the contructor or not (default true).
 *  \param clearNumericalSolutions Boolean to determine whether to clear the raw numerical solution member variables
 *  after propagation and resetting ephemerides (default false).
 *  \param setIntegratedResult Boolean to determine whether to automatically use the integrated results to set
 *  ephemerides (default false).
 *  \return Single-arc dynamics simulator object
 */
template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType > > createSingleArcDynamicsSimulator(
        const  simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated = true,
        const bool clearNumericalSolutions = true,
        const bool setIntegratedResult = true )
{
    return std::make_shared< propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                bodyMap,
                integratorSettings, propagatorSettings, areEquationsOfMotionToBeIntegrated, clearNumericalSolutions,
                setIntegratedResult );
}

//! Function to create single-arc variational equations solver object
/*!
 *  Function to create single-arc variational equations solver object
 *  \param bodyMap Map of bodies (with names) of all bodies in integration.
 *  \param integratorSettings Settings for numerical integrator of combined propagation of variational equations
 *  and equations of motion.
 *  \param propagatorSettings Settings for propagation of equations of motion.
 *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current
 *  settings and values.
 *  \param integrateDynamicalAndVariationalEquationsConcurrently Boolean defining whether variational and dynamical
 *  equations are to be propagated concurrently (if true) or sequentially (of false)
 *  \param variationalOnlyIntegratorSettings Settings for numerical integrator when integrating only variational
 *  equations.
 *  \param clearNumericalSolution Boolean to determine whether to clear the raw numerical solution member variables
 *  (default true) after propagation and resetting of state transition interface.
 *  \param integrateEquationsOnCreation Boolean to denote whether equations should be integrated immediately at the
 *  end of this contructor.
 *  \return Single-arc variational equations solver object
 */
template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType > >
createSingleArcVariationalEquationsSolver(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  StateScalarType > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently = 1,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings
        = std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
        const bool clearNumericalSolution = 1,
        const bool integrateEquationsOnCreation = 1 )
{
    return std::make_shared< propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                bodyMap, integratorSettings, propagatorSettings, parametersToEstimate,
                integrateDynamicalAndVariationalEquationsConcurrently, variationalOnlyIntegratorSettings,
                clearNumericalSolution, integrateEquationsOnCreation );
}

//! Function to create multi-arc variational equations solver object
/*!
 *  Function to create multi-arc variational equations solver object
 *  \param bodyMap Map of bodies (with names) of all bodies in integration.
 *  \param integratorSettings Settings for numerical integrator.
 *  \param propagatorSettings Settings for propagator.
 *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current settings and values.
 *  \param integrateDynamicalAndVariationalEquationsConcurrently Boolean defining whether variational and dynamical
 *  equations are to be propagated concurrently (if true) or sequentially (of false)
 *  \param variationalOnlyIntegratorSettings Settings for numerical integrator when integrating only variational
 *  equations.
 *  \param clearNumericalSolution Boolean to determine whether to clear the raw numerical solution member variables
 *  (default true) after propagation and resetting of state transition interface.
 *  \param integrateEquationsOnCreation Boolean to denote whether equations should be integrated immediately at the
 *  end of this contructor.
 *  \return Multi-arc variational equations solver object
 */
template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< propagators::MultiArcVariationalEquationsSolver< StateScalarType, TimeType > >
createMultiArcVariationalEquationsSolver(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  StateScalarType > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently = 1,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings
        = std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
        const bool clearNumericalSolution = 1,
        const bool integrateEquationsOnCreation = 1 )
{
    std::vector< double > arcStartTimes = estimatable_parameters::getMultiArcStateEstimationArcStartTimes(
                parametersToEstimate, true );
    return std::make_shared< propagators::MultiArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, arcStartTimes,
                integrateDynamicalAndVariationalEquationsConcurrently, variationalOnlyIntegratorSettings,
                clearNumericalSolution, integrateEquationsOnCreation );
}


//! Function to create hybrid-arc variational equations solver object
/*!
 *  Function to create hybrid-arc variational equations solver object
 *  \param bodyMap Map of bodies (with names) of all bodies in integration.
 *  \param integratorSettings Settings for numerical integrator.
 *  \param propagatorSettings Settings for propagator.
 *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current settings and values.
 *  \param integrateDynamicalAndVariationalEquationsConcurrently Boolean defining whether variational and dynamical
 *  equations are to be propagated concurrently (if true) or sequentially (of false)
 *  \param variationalOnlyIntegratorSettings Settings for numerical integrator when integrating only variational
 *  equations.
 *  \param clearNumericalSolution Boolean to determine whether to clear the raw numerical solution member variables
 *  (default true) after propagation and resetting of state transition interface.
 *  \param integrateEquationsOnCreation Boolean to denote whether equations should be integrated immediately at the
 *  end of this contructor.
 *  \return Hybrid-arc variational equations solver object
 */
template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< propagators::HybridArcVariationalEquationsSolver< StateScalarType, TimeType > >
createHybridArcVariationalEquationsSolver(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  StateScalarType > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently = 1,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings
        = std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
        const bool clearNumericalSolution = 1,
        const bool integrateEquationsOnCreation = 1 )
{
    std::vector< double > arcStartTimes = estimatable_parameters::getMultiArcStateEstimationArcStartTimes(
                parametersToEstimate, false );
    return std::make_shared< propagators::HybridArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, arcStartTimes,
                integrateDynamicalAndVariationalEquationsConcurrently,
                clearNumericalSolution, integrateEquationsOnCreation );
}

//! Function to create variational equations solver object
/*!
*  Function to create variational equations solver object
*  \param bodyMap Map of bodies (with names) of all bodies in integration.
*  \param integratorSettings Settings for numerical integrator.
*  \param propagatorSettings Settings for propagator.
*  \param parametersToEstimate Object containing all parameters that are to be estimated and their current settings and values.
*  \param integrateDynamicalAndVariationalEquationsConcurrently Boolean defining whether variational and dynamical
*  equations are to be propagated concurrently (if true) or sequentially (of false)
*  \param variationalOnlyIntegratorSettings Settings for numerical integrator when integrating only variational
*  equations.
*  \param clearNumericalSolution Boolean to determine whether to clear the raw numerical solution member variables
*  (default true) after propagation and resetting of state transition interface.
*  \param integrateEquationsOnCreation Boolean to denote whether equations should be integrated immediately at the
*  end of this contructor.
*  \return Variational equations solver object
*/
template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< propagators::VariationalEquationsSolver< StateScalarType, TimeType > >
createVariationalEquationsSolver(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  StateScalarType > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently = 1,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings
        = std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
        const bool clearNumericalSolution = 1,
        const bool integrateEquationsOnCreation = 1 )
{
    if( std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != nullptr )
    {
        return createSingleArcVariationalEquationsSolver(
                    bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, integrateDynamicalAndVariationalEquationsConcurrently,
                    variationalOnlyIntegratorSettings, clearNumericalSolution, integrateEquationsOnCreation );
    }
    else if( std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != nullptr )
    {
        return createMultiArcVariationalEquationsSolver(
                    bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, integrateDynamicalAndVariationalEquationsConcurrently,
                    variationalOnlyIntegratorSettings, clearNumericalSolution, integrateEquationsOnCreation );
    }
    else if( std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != nullptr )
    {
        return createHybridArcVariationalEquationsSolver(
                    bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, integrateDynamicalAndVariationalEquationsConcurrently,
                    variationalOnlyIntegratorSettings, clearNumericalSolution, integrateEquationsOnCreation );
    }
    else
    {
        throw std::runtime_error( "Error when creating variational equations solver, dynamics type not recognized" );
    }
}


//! Function to crate interface for state transition/sensitivity matrix results
/*!
 *  Function to crate interface for state transition/sensitivity matrix results
 *  \param propagatorSettings Settings used for propagation
 *  \param dynamicalStateSize Size of dynamical state that is estimation
 *  \param totalParameterSize Total size of estimated state vector
 */
template< typename StateScalarType = double >
std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > createStateTransitionAndSensitivityMatrixInterface(
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const int dynamicalStateSize,
        const int totalParameterSize )
{
    if( std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != nullptr )
    {
        return  std::make_shared<
                propagators::SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >( ),
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >( ),
                    dynamicalStateSize, totalParameterSize );
    }
    else if( std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != nullptr )
    {
        return  std::make_shared<
                propagators::MultiArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >( ),
                    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >( ),
                    std::vector< double >( ),
                    dynamicalStateSize, totalParameterSize );
    }
    else if( std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != nullptr )
    {
        throw std::runtime_error( "Error, cannot yet create empty hybrid arc state transition interface" );
        return nullptr;
    }
    else
    {
        throw std::runtime_error( "Error when creating state transition interface, dynamics type not recognized" );
    }
}


extern template std::shared_ptr< propagators::SingleArcVariationalEquationsSolver< double, double > >
createSingleArcVariationalEquationsSolver< double, double >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

extern template std::shared_ptr< propagators::SingleArcDynamicsSimulator< double, double > >
createSingleArcDynamicsSimulator< double, double >(
        const  simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated,
        const bool clearNumericalSolutions,
        const bool setIntegratedResult );

extern template std::shared_ptr< propagators::MultiArcVariationalEquationsSolver< double, double > >
createMultiArcVariationalEquationsSolver< double, double >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

extern template std::shared_ptr< propagators::HybridArcVariationalEquationsSolver< double, double > >
createHybridArcVariationalEquationsSolver< double, double >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );


extern template std::shared_ptr< propagators::VariationalEquationsSolver< double, double > >
createVariationalEquationsSolver< double, double >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )


extern template std::shared_ptr< propagators::SingleArcDynamicsSimulator< double, Time > >
createSingleArcDynamicsSimulator< double, Time >(
        const  simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated,
        const bool clearNumericalSolutions,
        const bool setIntegratedResult );

extern template std::shared_ptr< propagators::SingleArcDynamicsSimulator< long double, double > >
createSingleArcDynamicsSimulator< long double, double >(
        const  simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated,
        const bool clearNumericalSolutions,
        const bool setIntegratedResult );

extern template std::shared_ptr< propagators::SingleArcDynamicsSimulator< long double, Time > >
createSingleArcDynamicsSimulator< long double, Time >(
        const  simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated,
        const bool clearNumericalSolutions,
        const bool setIntegratedResult );

extern template std::shared_ptr< propagators::SingleArcVariationalEquationsSolver< double, Time > >
createSingleArcVariationalEquationsSolver< double, Time >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

extern template std::shared_ptr< propagators::SingleArcVariationalEquationsSolver< long double, Time > >
createSingleArcVariationalEquationsSolver< long double, Time >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  long double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

extern template std::shared_ptr< propagators::SingleArcVariationalEquationsSolver< long double, double > >
createSingleArcVariationalEquationsSolver< long double, double >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  long double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );



extern template std::shared_ptr< propagators::HybridArcVariationalEquationsSolver< double, Time > >
createHybridArcVariationalEquationsSolver< double, Time >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

extern template std::shared_ptr< propagators::HybridArcVariationalEquationsSolver< long double, Time > >
createHybridArcVariationalEquationsSolver< long double, Time >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  long double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

extern template std::shared_ptr< propagators::HybridArcVariationalEquationsSolver< long double, double > >
createHybridArcVariationalEquationsSolver< long double, double >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  long double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );



extern template std::shared_ptr< propagators::MultiArcVariationalEquationsSolver< double, Time > >
createMultiArcVariationalEquationsSolver< double, Time >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

extern template std::shared_ptr< propagators::MultiArcVariationalEquationsSolver< long double, Time > >
createMultiArcVariationalEquationsSolver< long double, Time >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  long double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

extern template std::shared_ptr< propagators::MultiArcVariationalEquationsSolver< long double, double > >
createMultiArcVariationalEquationsSolver< long double, double >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  long double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );




extern template std::shared_ptr< propagators::VariationalEquationsSolver< double, Time > >
createVariationalEquationsSolver< double, Time >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

extern template std::shared_ptr< propagators::VariationalEquationsSolver< long double, Time > >
createVariationalEquationsSolver< long double, Time >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< Time > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  long double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

extern template std::shared_ptr< propagators::VariationalEquationsSolver< long double, double > >
createVariationalEquationsSolver< long double, double >(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< long double > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  long double > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings,
        const bool clearNumericalSolution,
        const bool integrateEquationsOnCreation );

#endif

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATENUMERICALSIMULATOR_H
