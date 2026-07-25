/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_HYBRIDARCDYNAMICSSIMULATOR_H
#define TUDAT_HYBRIDARCDYNAMICSSIMULATOR_H

#include "tudat/simulation/propagation_setup/multiArcDynamicsSimulator.h"

namespace tudat
{

namespace propagators
{
template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > validateDeprecatedHybridArcSettings(
        const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
        const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const std::vector< double > arcStartTimes,
        const bool clearNumericalSolutions = true,
        const bool setIntegratedResult = true,
        const bool updateDependentVariableInterpolator = false )
{
    return validateDeprecatedHybridArcSettings( integratorSettings,
                                                utilities::deepcopyPointer( integratorSettings ),
                                                propagatorSettings,
                                                arcStartTimes,
                                                clearNumericalSolutions,
                                                setIntegratedResult,
                                                updateDependentVariableInterpolator );
}

//! Class for performing full numerical integration of a dynamical system, with a compbination of single and multi-arc propagations
/*!
 *  Class for performing full numerical integration of a dynamical system, with a compbination of single and multi-arc
 *  propagations. it is assumed that the single-arc propagations are not influence by the multi-arc propagations: first the
 *  single-arc propagation is done, followed by the multi-arc one (using the single-arc propagaton result for the environment)
 *  In this class, the governing equations are set once,  but can be re-integrated for different initial conditions using
 * the same instance of the class.
 */
template< typename StateScalarType = double, typename TimeType = double >
class HybridArcDynamicsSimulator : public DynamicsSimulator< StateScalarType, TimeType >
{
public:
    //! Using statemebts
    using DynamicsSimulator< StateScalarType, TimeType >::bodies_;

    typedef HybridArcSimulationResults< SingleArcSimulationResults, StateScalarType, TimeType > HybridArcResults;

    HybridArcDynamicsSimulator( const simulation_setup::SystemOfBodies& bodies,
                                const std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > hybridPropagatorSettings,
                                const bool areEquationsOfMotionToBeIntegrated = true,
                                const bool addSingleArcBodiesToMultiArcDynamics = false ):
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

    //! Constructor of multi-arc simulator for same integration settings per arc.
    /*!
     *  Constructor of multi-arc simulator for same integration settings per arc.
     *  \param bodies Map of bodies (with names) of all bodies in integration.
     *  \param integratorSettings Integrator settings for numerical integrator, used for all arcs.
     *  \param propagatorSettings Propagator settings for dynamics (must be of multi arc type)
     *  \param arcStartTimes Times at which the separate arcs start, for the multi-arc case
     *  \param areEquationsOfMotionToBeIntegrated Boolean to denote whether equations of motion should be integrated at
     *  the end of the contructor or not.//extern template class MultiArcDynamicsSimulator< double, double >;


     *  \param clearNumericalSolutions Boolean to determine whether to clear the raw numerical solution member variables
     *  after propagation and resetting ephemerides (default true).
     *  \param setIntegratedResult Boolean to determine whether to automatically use the integrated results to set
     *  ephemerides (default true).
     *  \param addSingleArcBodiesToMultiArcDynamics Boolean denoting whether to add single arc bodies to multi-arc
     *  dynamics (default false).
     */
    HybridArcDynamicsSimulator( const simulation_setup::SystemOfBodies& bodies,
                                const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
                                const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
                                const std::vector< double > arcStartTimes,
                                const bool areEquationsOfMotionToBeIntegrated = true,
                                const bool clearNumericalSolutions = true,
                                const bool setIntegratedResult = true,
                                const bool addSingleArcBodiesToMultiArcDynamics = false,
                                const bool updateDependentVariableInterpolator = false ):
        HybridArcDynamicsSimulator( bodies,
                                    validateDeprecatedHybridArcSettings( integratorSettings,
                                                                         propagatorSettings,
                                                                         arcStartTimes,
                                                                         clearNumericalSolutions,
                                                                         setIntegratedResult,
                                                                         updateDependentVariableInterpolator ),
                                    areEquationsOfMotionToBeIntegrated,
                                    addSingleArcBodiesToMultiArcDynamics )
    {}

    //    HybridArcDynamicsSimulator(
    //            const simulation_setup::SystemOfBodies& bodies,
    //            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > singleArcIntegratorSettings,
    //            const std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > multiArcIntegratorSettings,
    //            const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
    //            const std::vector< double > arcStartTimes,
    //            const bool areEquationsOfMotionToBeIntegrated = true,
    //            const bool clearNumericalSolutions = true,
    //            const bool setIntegratedResult = true,
    //            const bool addSingleArcBodiesToMultiArcDynamics = false,
    //            const bool updateDependentVariableInterpolator = false ):
    //        HybridArcDynamicsSimulator( bodies, validateDeprecatedHybridArcSettings(
    //                                        singleArcIntegratorSettings, multiArcIntegratorSettings, propagatorSettings, arcStartTimes,
    //                                        clearNumericalSolutions, setIntegratedResult ),
    //                                    areEquationsOfMotionToBeIntegrated, addSingleArcBodiesToMultiArcDynamics,
    //                                    updateDependentVariableInterpolator ){ }
    //! Destructor
    ~HybridArcDynamicsSimulator( ) {}

    //! This function numerically (re-)integrates the equations of motion, using concatenated states for single and multi-arcs
    /*!
     *  This function numerically (re-)integrates the equations of motion, using the settings set through the constructor
     *  and a new initial state vector provided here. The raw results are set in the equationsOfMotionNumericalSolution_
     *  \param initialGlobalStates Initial state vector that is to be used for numerical integration. Note that this state
     *  should be in the correct frame (i.e. corresponding to centralBodies in propagatorSettings_), but not in the propagator-
     *  specific form (i.e Encke, Gauss, etc. for translational dynamics). The states for all arcs must be concatenated in
     *  order into a single Eigen Vector, starting with the single-arc states, followed by the mulit-arc states
     */
    void integrateEquationsOfMotion( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialGlobalStates )
    {
        printPrePropagationMessages( );
        singleArcDynamicsSimulator_->integrateEquationsOfMotion( initialGlobalStates.block( 0, 0, singleArcDynamicsSize_, 1 ) );
        multiArcDynamicsSimulator_->integrateEquationsOfMotion(
                initialGlobalStates.block( singleArcDynamicsSize_, 0, multiArcDynamicsSize_, 1 ) );
        printPostPropagationMessages( );
    }

    void integrate( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& concatenatedInitialStates )
    {
        integrateEquationsOfMotion( concatenatedInitialStates );
    }

    //! This function updates the environment with the numerical solution of the propagation
    /*!
     *  This function updates the environment with the numerical solution of the propagation
     *  (no additional functionality in hybrid arc). Function may be used to process manually updtaed propagation results in
     *  single and/or multi-arc model
     */
    void processNumericalEquationsOfMotionSolution( )
    {
        singleArcDynamicsSimulator_->processNumericalEquationsOfMotionSolution( );
        multiArcDynamicsSimulator_->processNumericalEquationsOfMotionSolution( );
        if( hybridPropagatorSettings_->getOutputSettings( )->getUpdateDependentVariableInterpolator( ) )
        {
            propagationResults_->updateDependentVariableInterface( );
        }
    }

    //! Function to retrieve the single-arc dynamics simulator
    /*!
     * Function to retrieve the single-arc dynamics simulator
     * \return Single-arc dynamics simulator
     */
    std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > getSingleArcDynamicsSimulator( )
    {
        return singleArcDynamicsSimulator_;
    }

    //! Function to retrieve the multi-arc dynamics simulator
    /*!
     * Function to retrieve the multi-arc dynamics simulator
     * \return Multi-arc dynamics simulator
     */
    std::shared_ptr< MultiArcDynamicsSimulator< StateScalarType, TimeType > > getMultiArcDynamicsSimulator( )
    {
        return multiArcDynamicsSimulator_;
    }

    //! Get whether the integration was completed successfully.
    /*!
     * Get whether the integration was completed successfully.
     * \return Whether the integration was completed successfully by reaching the termination condition.
     */
    virtual bool integrationCompletedSuccessfully( ) const
    {
        return propagationResults_->integrationCompletedSuccessfully( );
    }

    //! Function to return the numerical solution to the equations of motion (base class interface).
    /*!
     *  Function to return the numerical solution to the equations of motion for last numerical integration. First vector entry
     *  contains single-arc results. Each subsequent vector entry contains one of the multi-arcs. Key of map denotes time,
     *  values are full propagated state vectors.
     *  \return List of maps of history of numerically integrated states.
     */
    std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > > getEquationsOfMotionNumericalSolutionBase( )
    {
        return propagationResults_->getConcatenatedEquationsOfMotionResults( );
    }

    //! Function to return the numerical solution to the dependent variables (base class interface).
    /*!
     *  Function to return the numerical solution to the dependent variables for last numerical integration. First vector entry
     *  contains single-arc results. Each subsequent vector entry contains one of the multi-arcs. Key of map denotes time,
     *  values are dependent variables vectors
     *  \return List of maps of dependent variable history
     */
    std::vector< std::map< TimeType, Eigen::VectorXd > > getDependentVariableNumericalSolutionBase( )
    {
        return propagationResults_->getConcatenatedDependentVariableResults( );
    }

    //! Function to return the map of cumulative computation time history that was saved during numerical propagation.
    /*!
     *  Function to return the map of cumulative computation time history that was saved during numerical propagation.  First vector
     *  entry contains single-arc results. Each subsequent vector entry contains one of the multi-arcs. Key of map denotes time,
     *  values are computation times.
     *  \return Vector is size 1, with entry: map of cumulative computation time history that was saved during numerical propagation.
     */
    std::vector< std::map< TimeType, double > > getCumulativeComputationTimeHistoryBase( )
    {
        return propagationResults_->getConcatenatedCumulativeComputationTimeHistory( );
    }

    void printPrePropagationMessages( )
    {
        if( hybridPropagatorSettings_->getOutputSettings( )->printAnyOutput( ) )
        {
            std::cout << hybridPropagatorSettings_->getOutputSettings( )->getPropagationStartHeader( ) << std::endl << std::endl;
        }
    }

    void printPostPropagationMessages( )
    {
        if( hybridPropagatorSettings_->getOutputSettings( )->printAnyOutput( ) )
        {
            std::cout << hybridPropagatorSettings_->getOutputSettings( )->getPropagationEndHeader( ) << std::endl << std::endl;
        }
    }

    std::shared_ptr< SimulationResults< StateScalarType, TimeType > > getPropagationResults( )
    {
        return propagationResults_;
    }

    std::shared_ptr< HybridArcResults > getHybridArcPropagationResults( )
    {
        return propagationResults_;
    }

protected:
    std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > hybridPropagatorSettings_;

    //! Object used to propagate single-arc dynamics
    std::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > singleArcDynamicsSimulator_;

    //! Object used to propagate multi-arc dynamics
    std::shared_ptr< MultiArcDynamicsSimulator< StateScalarType, TimeType > > multiArcDynamicsSimulator_;

    //! Size of single-arc (initial) state vector
    int singleArcDynamicsSize_;

    //! Size of multi-arc concatenated initial state vector
    int multiArcDynamicsSize_;

    std::shared_ptr< HybridArcResults > propagationResults_;
};

template< typename StateScalarType = double, typename TimeType = double >
std::shared_ptr< PropagatorSettings< StateScalarType > > validateDeprecatePropagatorSettings(
        const std::vector< std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > >& integratorSettings,
        const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings )
{
    if( std::dynamic_pointer_cast< propagators::SingleArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) !=
        nullptr )
    {
        if( integratorSettings.size( ) == 0 )
        {
            throw std::runtime_error(
                    "Error when validating deprecated propagator settings; did not find integrator settings for single-arc propagation" );
            ;
        }
        return validateDeprecatedSingleArcSettings( integratorSettings.at( 0 ), propagatorSettings );
    }
    else if( std::dynamic_pointer_cast< propagators::MultiArcPropagatorSettings< StateScalarType, TimeType > >( propagatorSettings ) !=
             nullptr )
    {
        return validateDeprecatedMultiArcSettings( integratorSettings, propagatorSettings );
    }
    else if( std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< StateScalarType > >( propagatorSettings ) != nullptr )
    {
        std::shared_ptr< propagators::HybridArcPropagatorSettings< StateScalarType > > hybridArcSettings =
                std::dynamic_pointer_cast< propagators::HybridArcPropagatorSettings< StateScalarType > >( propagatorSettings );

        if( integratorSettings.size( ) == 0 )
        {
            throw std::runtime_error(
                    "Error when validating deprecated propagator settings; did not find integrator settings for hybrid-arc propagation" );
            ;
        }
        validateDeprecatedSingleArcSettings< StateScalarType, TimeType >( integratorSettings.at( 0 ),
                                                                          hybridArcSettings->getSingleArcPropagatorSettings( ) );
        validateDeprecatedMultiArcSettings< StateScalarType, TimeType >( { integratorSettings.begin( ) + 1, integratorSettings.end( ) },
                                                                         hybridArcSettings->getMultiArcPropagatorSettings( ) );
        return hybridArcSettings;
    }
    else
    {
        throw std::runtime_error( "Error when validating deprecated propagator settings" );
        return nullptr;
    }
}

#if TUDAT_BUILD_EXPLICIT_INSTANTIATIONS
extern template class HybridArcDynamicsSimulator< double, double >;
#endif

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_HYBRIDARCDYNAMICSSIMULATOR_H
