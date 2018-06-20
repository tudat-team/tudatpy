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
#include "Tudat/SimulationSetup/PropagationSetup/variationalEquationsSolver.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"

namespace tudat
{

namespace simulation_setup
{

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

template< typename StateScalarType = double, typename TimeType = double, typename ParameterType = double >
std::shared_ptr< propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType > >
createSingleArcVariationalEquationsSolver(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
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

template< typename StateScalarType = double, typename TimeType = double, typename ParameterType = double >
std::shared_ptr< propagators::HybridArcVariationalEquationsSolver< StateScalarType, TimeType > >
createHybridArcVariationalEquationsSolver(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
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

template< typename StateScalarType = double, typename TimeType = double, typename ParameterType = double >
std::shared_ptr< propagators::MultiArcVariationalEquationsSolver< StateScalarType, TimeType > >
createMultiArcVariationalEquationsSolver(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
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


template< typename StateScalarType = double, typename TimeType = double, typename ParameterType = double >
std::shared_ptr< propagators::VariationalEquationsSolver< StateScalarType, TimeType > >
createVariationalEquationsSolver(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const bool integrateDynamicalAndVariationalEquationsConcurrently = 1,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings
        = std::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
        const bool clearNumericalSolution = 1,
        const bool integrateEquationsOnCreation = 1 )
{
    if( boost::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != NULL )
    {
        return createSingleArcVariationalEquationsSolver(
                    bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, integrateDynamicalAndVariationalEquationsConcurrently,
                    variationalOnlyIntegratorSettings, clearNumericalSolution, integrateEquationsOnCreation );
    }
    else if( boost::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != NULL )
    {
        return createMultiArcVariationalEquationsSolver(
                    bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, integrateDynamicalAndVariationalEquationsConcurrently,
                    variationalOnlyIntegratorSettings, clearNumericalSolution, integrateEquationsOnCreation );
    }
    else if( boost::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != NULL )
    {
        return createHybridArcVariationalEquationsSolver(
                    bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, integrateDynamicalAndVariationalEquationsConcurrently,
                    variationalOnlyIntegratorSettings, clearNumericalSolution, integrateEquationsOnCreation );
    }
    else
    {
        throw std::runtime_error( "Error when creating variational equations solved, dynamics type not recognized" );
    }
}

template< typename StateScalarType = double >
std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > createStateTransitionAndSensitivityMatrixInterface(
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const int dynamicalStateSize,
        const int totalParameterSize )
{
    if( boost::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != NULL )
    {
        return  std::make_shared<
                propagators::SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >( ),
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >( ),
                    dynamicalStateSize, totalParameterSize );
    }
    else if( boost::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != NULL )
    {
        return  std::make_shared<
                propagators::MultiArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >( ),
                    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >( ),
                    std::vector< double >( ),
                    dynamicalStateSize, totalParameterSize );
    }
    else if( boost::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != NULL )
    {
        throw std::runtime_error( "Error, cannot yet create empty hybrid arc state transition interface" );
        return NULL;
    }
    else
    {
        throw std::runtime_error( "Error when creating state transition interface, dynamics type not recognized" );
    }
}

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATENUMERICALSIMULATOR_H
