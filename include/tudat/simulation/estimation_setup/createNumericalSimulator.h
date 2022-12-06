/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create single-arc dynamics simulator object
/*!
 *  Function to create single-arc dynamics simulator object
 *  \param bodies Map of bodies (with names) of all bodies in integration.
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
        const  simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated = true,
        const bool clearNumericalSolutions = true,
        const bool setIntegratedResult = true )
{
    return std::make_shared< propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                bodies,
                propagatorSettings, areEquationsOfMotionToBeIntegrated, clearNumericalSolutions,
                setIntegratedResult );
}

template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< propagators::DynamicsSimulator< StateScalarType, TimeType > > createDynamicsSimulator(
        const  simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated = true )
{
    if( std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) != nullptr )
    {
        return std::make_shared< propagators::SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                    bodies, std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ),
                    areEquationsOfMotionToBeIntegrated );
    }
    else if( std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) != nullptr )
    {
        return std::make_shared< propagators::MultiArcDynamicsSimulator< StateScalarType, TimeType > >(
                    bodies, std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ),
                    areEquationsOfMotionToBeIntegrated );
    }
    else if( std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != nullptr )
    {
        return std::make_shared< propagators::HybridArcDynamicsSimulator< StateScalarType, TimeType > >(
                    bodies, std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ),
                    areEquationsOfMotionToBeIntegrated );
    }
    else
    {
        throw std::runtime_error( "Error when creating variational equations solver, dynamics type not recognized" );
    }
}

//! Function to create variational equations solver object
/*!
*  Function to create variational equations solver object
*  \param bodies Map of bodies (with names) of all bodies in integration.
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
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet<  StateScalarType > > parametersToEstimate,
        const bool integrateEquationsOnCreation = 1 )
{
    if( std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) != nullptr )
    {
        return std::make_shared< propagators::SingleArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                    bodies, std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< StateScalarType, TimeType > >(
                        propagatorSettings ), parametersToEstimate, integrateEquationsOnCreation );
    }
    else if( std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) != nullptr )
    {
        return std::make_shared< propagators::MultiArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                    bodies, std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< StateScalarType, TimeType > >(
                        propagatorSettings ), parametersToEstimate, integrateEquationsOnCreation );
    }
    else if( std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) != nullptr )
    {
        return std::make_shared< propagators::HybridArcVariationalEquationsSolver< StateScalarType, TimeType > >(
                    bodies, std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< StateScalarType, TimeType > >(
                        propagatorSettings ), parametersToEstimate, integrateEquationsOnCreation );
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
template< typename StateScalarType, typename TimeType >
std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > createStateTransitionAndSensitivityMatrixInterface(
        const std::shared_ptr< propagators::PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
        const int dynamicalStateSize,
        const int totalParameterSize )
{
    if( std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) != nullptr )
    {
        return  std::make_shared<
                propagators::SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >( ),
                    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >( ),
                    dynamicalStateSize, totalParameterSize, std::vector< std::pair< int, int > >( ) );
    }
    else if( std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) != nullptr )
    {
        return  std::make_shared<
                propagators::MultiArcCombinedStateTransitionAndSensitivityMatrixInterface< StateScalarType > >(
                    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >( ),
                    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > >( ),
                    std::vector< double >( ),
                    std::vector< double >( ),
                    parametersToEstimate,
                    dynamicalStateSize, totalParameterSize,
                    std::vector< std::vector< std::pair< int, int > > >( ));
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


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATENUMERICALSIMULATOR_H
