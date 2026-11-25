/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_HYBRIDARCDYNAMICSSIMULATORIMPLEMENTATION_H
#define TUDAT_HYBRIDARCDYNAMICSSIMULATORIMPLEMENTATION_H

#include "tudat/simulation/propagation_setup/hybridArcDynamicsSimulator.h"

namespace tudat
{

namespace propagators
{

template< typename StateScalarType, typename TimeType >
HybridArcDynamicsSimulator< StateScalarType, TimeType >::HybridArcDynamicsSimulator(
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > hybridPropagatorSettings,
        const bool areEquationsOfMotionToBeIntegrated,
        const bool addSingleArcBodiesToMultiArcDynamics ):
    DynamicsSimulator< StateScalarType, TimeType >( bodies, hybridPropagatorSettings ),
    hybridPropagatorSettings_( hybridPropagatorSettings )
{
    if( hybridPropagatorSettings_ == nullptr )
    {
        throw std::runtime_error( "Error when making HybridArcDynamicsSimulator, propagator settings are incompatible" );
    }

    singleArcDynamicsSize_ = hybridPropagatorSettings_->getSingleArcPropagatorSettings( )->getPropagatedStateSize( );
    multiArcDynamicsSize_ = hybridPropagatorSettings_->getMultiArcPropagatorSettings( )->getPropagatedStateSize( );

    if( !addSingleArcBodiesToMultiArcDynamics )
    {
        singleArcDynamicsSimulator_ = std::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                bodies,
                hybridPropagatorSettings_->getSingleArcPropagatorSettings( ),
                false,
                PredefinedSingleArcStateDerivativeModels< StateScalarType, TimeType >( ) );
        multiArcDynamicsSimulator_ = std::make_shared< MultiArcDynamicsSimulator< StateScalarType, TimeType > >(
                bodies, hybridPropagatorSettings_->getMultiArcPropagatorSettings( ), false );
        propagationResults_ = std::make_shared< HybridArcResults >( singleArcDynamicsSimulator_->getSingleArcPropagationResults( ),
                                                                    multiArcDynamicsSimulator_->getMultiArcPropagationResults( ) );
    }
    else
    {
        throw std::runtime_error( "Cannot yet add single-arc bodies to multi-arc propagation" );
    }

    if( areEquationsOfMotionToBeIntegrated )
    {
        integrateEquationsOfMotion( hybridPropagatorSettings_->getInitialStates( ) );
    }
}

template< typename StateScalarType, typename TimeType >
void HybridArcDynamicsSimulator< StateScalarType, TimeType >::integrateEquationsOfMotion(
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialGlobalStates )
{
    printPrePropagationMessages( );
    singleArcDynamicsSimulator_->integrateEquationsOfMotion( initialGlobalStates.block( 0, 0, singleArcDynamicsSize_, 1 ) );
    multiArcDynamicsSimulator_->integrateEquationsOfMotion(
            initialGlobalStates.block( singleArcDynamicsSize_, 0, multiArcDynamicsSize_, 1 ) );
    printPostPropagationMessages( );
}

template< typename StateScalarType, typename TimeType >
void HybridArcDynamicsSimulator< StateScalarType, TimeType >::integrate(
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& concatenatedInitialStates )
{
    integrateEquationsOfMotion( concatenatedInitialStates );
}

template< typename StateScalarType, typename TimeType >
void HybridArcDynamicsSimulator< StateScalarType, TimeType >::processNumericalEquationsOfMotionSolution( )
{
    singleArcDynamicsSimulator_->processNumericalEquationsOfMotionSolution( );
    multiArcDynamicsSimulator_->processNumericalEquationsOfMotionSolution( );
    if( hybridPropagatorSettings_->getOutputSettings( )->getUpdateDependentVariableInterpolator( ) )
    {
        propagationResults_->updateDependentVariableInterface( );
    }
}

template< typename StateScalarType, typename TimeType >
bool HybridArcDynamicsSimulator< StateScalarType, TimeType >::integrationCompletedSuccessfully( ) const
{
    return propagationResults_->integrationCompletedSuccessfully( );
}

template< typename StateScalarType, typename TimeType >
std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >
HybridArcDynamicsSimulator< StateScalarType, TimeType >::getEquationsOfMotionNumericalSolutionBase( )
{
    return propagationResults_->getConcatenatedEquationsOfMotionResults( );
}

template< typename StateScalarType, typename TimeType >
std::vector< std::map< TimeType, Eigen::VectorXd > >
HybridArcDynamicsSimulator< StateScalarType, TimeType >::getDependentVariableNumericalSolutionBase( )
{
    return propagationResults_->getConcatenatedDependentVariableResults( );
}

template< typename StateScalarType, typename TimeType >
std::vector< std::map< TimeType, double > >
HybridArcDynamicsSimulator< StateScalarType, TimeType >::getCumulativeComputationTimeHistoryBase( )
{
    return propagationResults_->getConcatenatedCumulativeComputationTimeHistory( );
}

template< typename StateScalarType, typename TimeType >
void HybridArcDynamicsSimulator< StateScalarType, TimeType >::printPrePropagationMessages( )
{
    if( hybridPropagatorSettings_->getOutputSettings( )->printAnyOutput( ) )
    {
        std::cout << hybridPropagatorSettings_->getOutputSettings( )->getPropagationStartHeader( ) << std::endl << std::endl;
    }
}

template< typename StateScalarType, typename TimeType >
void HybridArcDynamicsSimulator< StateScalarType, TimeType >::printPostPropagationMessages( )
{
    if( hybridPropagatorSettings_->getOutputSettings( )->printAnyOutput( ) )
    {
        std::cout << hybridPropagatorSettings_->getOutputSettings( )->getPropagationEndHeader( ) << std::endl << std::endl;
    }
}

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_HYBRIDARCDYNAMICSSIMULATORIMPLEMENTATION_H
