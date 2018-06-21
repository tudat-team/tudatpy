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

extern template std::shared_ptr< propagators::SingleArcDynamicsSimulator< double, double > >
createSingleArcDynamicsSimulator< double, double >(
        const  simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
        const std::shared_ptr< propagators::PropagatorSettings< double > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated,
        const bool clearNumericalSolutions,
        const bool setIntegratedResult );

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
    if( std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != NULL )
    {
        return createSingleArcVariationalEquationsSolver(
                    bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, integrateDynamicalAndVariationalEquationsConcurrently,
                    variationalOnlyIntegratorSettings, clearNumericalSolution, integrateEquationsOnCreation );
    }
    else if( std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != NULL )
    {
        return createMultiArcVariationalEquationsSolver(
                    bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, integrateDynamicalAndVariationalEquationsConcurrently,
                    variationalOnlyIntegratorSettings, clearNumericalSolution, integrateEquationsOnCreation );
    }
    else if( std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != NULL )
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

template< typename StateScalarType = double >
std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > createStateTransitionAndSensitivityMatrixInterface(
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const int dynamicalStateSize,
        const int totalParameterSize )
{
    if( std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != NULL )
    {
        return  std::make_shared<
                propagators::SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >( ),
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >( ),
                    dynamicalStateSize, totalParameterSize );
    }
    else if( std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != NULL )
    {
        return  std::make_shared<
                propagators::MultiArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >( ),
                    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >( ),
                    std::vector< double >( ),
                    dynamicalStateSize, totalParameterSize );
    }
    else if( std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != NULL )
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
