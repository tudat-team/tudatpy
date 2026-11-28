/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_HYBRIDARCVARIATIONALEQUATIONSSOLVER_H
#define TUDAT_HYBRIDARCVARIATIONALEQUATIONSSOLVER_H

#include "tudat/simulation/estimation_setup/singleArcVariationalEquationsSolver.h"
#include "tudat/simulation/estimation_setup/multiArcVariationalEquationsSolver.h"
#include "tudat/simulation/propagation_setup/hybridArcDynamicsSimulator.h"
#include "tudat/simulation/estimation_setup/getAssociatedMultiArcParameter.h"
#include "tudat/simulation/estimation_setup/setInitialStateVectorFromParameterSet.h"

namespace tudat
{

namespace propagators
{
//! Class to manage and execute the numerical integration of variational equations of a dynamical system in a combination
//! of single and multiple arcs
/*!
 *  Class to manage and execute the numerical integration of variational equations of a dynamical system, in addition
 *  to the dynamics itself, in a combination  of single and multiple arcs. In this class, the governing equations are set once,
 *  but can be re-integrated for different initial conditions using the same instance of the class.
 */
template< typename StateScalarType = double, typename TimeType = double >
class HybridArcVariationalEquationsSolver : public VariationalEquationsSolver< StateScalarType, TimeType >
{
public:
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > MatrixType;
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > VectorType;
    typedef HybridArcSimulationResults< SingleArcVariationalSimulationResults, StateScalarType, TimeType > HybridArcResults;

    using VariationalEquationsSolver< StateScalarType, TimeType >::parametersToEstimate_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::bodies_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::stateTransitionMatrixSize_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::parameterVectorSize_;
    using VariationalEquationsSolver< StateScalarType, TimeType >::stateTransitionInterface_;

    HybridArcVariationalEquationsSolver(
            const simulation_setup::SystemOfBodies& bodies,
            const std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings,
            const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate,
            const bool integrateEquationsOnCreation = false );

    void initializeHybridArcVariationalEquationsSolver( const simulation_setup::SystemOfBodies& bodies,
                                                        const std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
                                                        const bool integrateEquationsOnCreation );

    //! Destructor
    ~HybridArcVariationalEquationsSolver( ) {}

    //! Function to integrate variational equations and equations of motion.
    /*!
     *  Function to integrate variational equations and equations of motion. At the end of this function,
     *  the stateTransitionInterface_ is reset with the new state transition and sensitivity matrices. If dynamical
     *  solution is to be processed, the environment is also updated to the new solution.
     *  \param initialStateEstimate Initial statez of the equations of motion that is to be used (in same order as in
     *  parametersToEstimate_). The initial states of single and multi-arcs propagations are concatenated into a single
     *  vector.
     *  \param integrateEquationsConcurrently Variable determining whether the equations of motion are to be
     *  propagated concurrently with variational equations of motion (if true), or before variational equations (if false).
     */
    void integrateVariationalAndDynamicalEquations( const VectorType& initialStateEstimate, const bool integrateEquationsConcurrently );

    //! Function to integrate equations of motion only.
    /*!
     *  Function to integrate equations of motion only.  If dynamical
     *  solution is to be processed, the environment is also updated to the new solution.
     *  \param initialStateEstimate Initial state of the equations of motion that is to be used (in same order as in
     *  parametersToEstimate_). The initial states of single and multi-arcs propagations are concatenated into a single
     *  vector.
     */
    void integrateDynamicalEquationsOfMotionOnly( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& initialStateEstimate );

    //! Function to reset parameter estimate and re-integrate equations of motion and, if desired, variational equations.
    /*!
     *  Function to reset parameter estimate and re-integrate equations of motion and, if desired, variational equations
     *  using the new physical parameters/body initial states.
     *  \param newParameterEstimate New estimate of parameters that are to be estimated, in same order as defined
     *  in parametersToEstimate_ member.
     *  \param areVariationalEquationsToBeIntegrated Boolean defining whether the variational equations are to be
     *  reintegrated with the new parameter values.
     */
    void resetParameterEstimate( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > newParameterEstimate,
                                 const bool areVariationalEquationsToBeIntegrated = true );

    //! Function to retrieve propagator settings used for equations of motion
    /*!
     * Function to retrieve propagator settings used for equations of motion
     * \return Propagator settings used for equations of motion
     */
    std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > getPropagatorSettings( )
    {
        return propagatorSettings_;
    }

    //! Function to retrieve the dynamics simulator object (as base-class pointer)
    /*!
     * Function to retrieve the dynamics simulator object (as base-class pointer). This function is not yet implemented
     * in hybric-arc model, as no single DynamicsSimulator model is used. Calling this function throws an error
     * \return Dynamics simulator object (as base-class pointer)
     */
    virtual std::shared_ptr< DynamicsSimulator< StateScalarType, TimeType > > getDynamicsSimulatorBase( );

    std::shared_ptr< MultiArcVariationalEquationsSolver< StateScalarType, TimeType > > getMultiArcSolver( )
    {
        return multiArcSolver_;
    }

    //! Object to solve single-arc variational equations.
    std::shared_ptr< SingleArcVariationalEquationsSolver< StateScalarType, TimeType > > getSingleArcSolver( )
    {
        return singleArcSolver_;
    }

    std::shared_ptr< HybridArcResults > getHybridArcVariationalPropagationResults( )
    {
        return variationalPropagationResults_;
    }

    std::shared_ptr< SimulationResults< StateScalarType, TimeType > > getVariationalPropagationResults( );

protected:
    //! Function to set and process the arc start times of the multi-arc propagation
    /*!
     * Function to set and process the arc start times of the multi-arc propagation
     * \param arcStartTimes Arc start times of the multi-arc propagation
     */
    void setExtendedMultiArcParameters( const std::vector< double >& arcStartTimes );

    //! Function to reset the initial multi-arc states
    /*!
     * Function to reset the initial multi-arc states
     * \param manualMultiArcStates New multi-arc states
     */
    void resetMultiArcInitialStates( const VectorType& manualMultiArcStates );

    //! Function that removes the single-arc body data from propagation results before processing data
    /*!
     * Function that removes the single-arc body data from propagation results before processing data
     * \param numericalMultiArcSolution Full numerical solution of single and multi-arc bodies.
     */
    void removeSingleArcBodiesFromMultiArcSolultion( std::vector< std::map< TimeType, VectorType > >& numericalMultiArcSolution );

    //! Update original propagator settings
    void copyExtendedMultiArcInitialStatesToOriginalSettins( );

    //! Object to solve multi-arc variational equations (multi-arc bodies only).
    std::shared_ptr< MultiArcVariationalEquationsSolver< StateScalarType, TimeType > > originalMultiArcSolver_;

    //! Object to solve single-arc variational equations.
    std::shared_ptr< SingleArcVariationalEquationsSolver< StateScalarType, TimeType > > singleArcSolver_;

    //! Object to solve multi-arc variational equations (single- and multi-arc bodies).
    std::shared_ptr< MultiArcVariationalEquationsSolver< StateScalarType, TimeType > > multiArcSolver_;

    //! Propagator settings, with single- and multi-arc separately
    std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > originalPopagatorSettings_;

    //! Propagator settings, with single-arc bodies added to multi-arc list
    std::shared_ptr< HybridArcPropagatorSettings< StateScalarType, TimeType > > propagatorSettings_;

    //! Settings to be used for integrator
    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > singleArcIntegratorSettings_;

    std::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > multiArcIntegratorSettings_;

    //! Size of estimated single-arc dynamical parameters
    int singleArcDynamicsSize_;

    //! Vector containing, for each arc, the size of original estimated multi-arc dynamical parameters
    std::vector< int > originalMultiArcDynamicsSingleArcSize_;

    //! Total size of original estimated multi-arc dynamical parameters
    int originalMultiArcDynamicsSize_;

    //! Size of single arc of extended estimated multi-arc dynamical parameters
    int multiArcDynamicsSize_;

    //! Vector containing, for each arc, the total size of extended estimated multi-arc dynamical parameters
    std::vector< int > multiArcDynamicsSingleArcSize_;

    //! Estimated parameter set with single-arc dynamical parameters only
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > singleArcParametersToEstimate_;

    //! Estimated parameter set with original multi-arc dynamical parameters only
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > originalMultiArcParametersToEstimate_;

    //! Estimated parameter set with extended multi-arc dynamical parameters only
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > multiArcParametersToEstimate_;

    std::shared_ptr< HybridArcResults > variationalPropagationResults_;

    //! Times at which arcs for multi-arc solution start
    std::vector< double > arcStartTimes_;

    //! Function that retrieves the single-arc bodies' initial states as a function of time
    /*!
     *  Function that retrieves the single-arc bodies' initial states as a function of time, is used to update the multi-arc
     *  initial states after the single-arc propagation
     */
    std::function< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >( const double ) > initialStatesFromSingleArcPropagation_;

    double singleArcInitialTime_;
};

#if TUDAT_BUILD_TESTS_WITH_EXPLICIT_INSTANTIATIONS
extern template class HybridArcVariationalEquationsSolver< double, double >;
#endif

}  // namespace propagators

}  // namespace tudat

#include "tudat/simulation/estimation_setup/hybridArcVariationalEquationsSolverImplementation.h"

#endif  // TUDAT_HYBRIDARCVARIATIONALEQUATIONSSOLVER_H
