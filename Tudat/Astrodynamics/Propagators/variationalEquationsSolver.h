#ifndef TUDAT_VARIATIONALEQUATIONSSOLVER_H
#define TUDAT_VARIATIONALEQUATIONSSOLVER_H

#include <boost/make_shared.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>

#include "Tudat/Basics/utilities.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Mathematics/Interpolators/interpolator.h"

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/Propagators/stateTransitionMatrixInterface.h"
#include "Tudat/Astrodynamics/Propagators/dynamicsSimulator.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/SimulationSetup/createStateDerivativePartials.h"

namespace tudat
{

namespace propagators
{


//! Base class to manage and execute the numerical integration of equations of motion and variational equations.
/*!
 *  Base class to manage and execute the numerical integration of equations of motion and variational equations.
 *  Governing equations are set once, but can be re-integrated for different initial conditions using the same
 *  instance of the class. Derived classes define the specific kind of integration that is performed
 *  (single-arc/multi-arc; dynamics/variational equations, etc.)
 */
template< typename StateScalarType = double, typename TimeType = double, typename ParameterType = double >
class VariationalEquationsSolver
{
public:

    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > MatrixType;
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > VectorType;

    //! Constructor
    /*!
     *  Constructor, sets up object for automatic evaluation and numerical integration of variational equations and
     *  equations of motion.
     *  \param bodyMap Map of bodies (with names) of all bodies in integration.
     *  \param integratorSettings Settings for numerical integrator of combined propagation of variational equations
     *  and equations of motion.
     *  \param propagatorSettings Settings for propagation of equations of motion.
     *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current
     *  settings and values.
     *  \param variationalOnlyIntegratorSettings Settings for numerical integrator when integrating only variational
     *  equations.
     *  \param clearNumericalSolution Boolean to determine whether to clear the raw numerical solution member variables
     *  (default true) after propagation and resetting of state transition interface.
     */
    VariationalEquationsSolver(
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings=
            boost::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
            const bool clearNumericalSolution = 1 ):
        parametersToEstimate_( parametersToEstimate ),
        bodyMap_( bodyMap ),
        propagatorSettings_( propagatorSettings ), integratorSettings_( integratorSettings ),
        variationalOnlyIntegratorSettings_( variationalOnlyIntegratorSettings ),
        stateTransitionMatrixSize_( parametersToEstimate_->getInitialDynamicalStateParameterSize( ) ),
        parameterVectorSize_( parametersToEstimate_->getParameterSetSize( ) ),
        clearNumericalSolution_( clearNumericalSolution )
    { }

    //! Destructor
    virtual ~VariationalEquationsSolver( ){ }

    //! Pure virtual function to integrate variational equations and equations of motion.
    /*!
     *  Pure virtual function to integrate variational equations and equations of motion, to be implemented in derived
     *  class
     *  \param initialStateEstimate Initial state of the equations of motion that is to be used.
     *  \param integrateEquationsConcurrently Variable determining whether the equations of motion are to be
     *  propagated concurrently with variational equations of motion (if true), or before variational equations (if false).
     */
    virtual void integrateVariationalAndDynamicalEquations(
            const VectorType& initialStateEstimate, const bool integrateEquationsConcurrently ) = 0;

    //! Pure virtual function to integrate equations of motion only.
    /*!
     *  Pure virtual function to integrate equations of motion only, to be implemented in derived
     *  class
     *  \param initialStateEstimate Initial state of the equations of motion that is to be used.
     */
    virtual void integrateDynamicalEquationsOfMotionOnly(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialStateEstimate ) = 0;


    //! Function to get the list of objects representing the parameters that are to be integrated.
    /*!
     *  Function to get the list of objects representing the parameters that are to be integrated.
     *  \return List of objects representing the parameters that are to be integrated.
     */
    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > getParametersToEstimate( )
    {
        return parametersToEstimate_;
    }

    //! Function to reset parameter estimate and re-integrate equations of motion and, if desired, variational equations.
    /*!
     *  Function to reset parameter estimate and re-integrate equations of motion and, if desired, variational equations
     *  using the new physical parameters/body initial states.
     *  \param newParameterEstimate New estimate of parameters that are to be estimated, in same order as defined
     *  in parametersToEstimate_ member.
     *  \param areVariationalEquationsToBeIntegrated Boolean defining whether the variational equations are to be
     *  reintegrated with the new parameter values.
     */
    void resetParameterEstimate( const Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > newParameterEstimate,
                                 const bool areVariationalEquationsToBeIntegrated = true )

    {
        // Reset values of parameters.
        parametersToEstimate_->template resetParameterValues< ParameterType >( newParameterEstimate );
        propagatorSettings_->resetInitialStates(
                    estimatable_parameters::getInitialStateVectorOfBodiesToEstimate( parametersToEstimate_ ) );

        dynamicsStateDerivative_->template updateStateDerivativeModelSettings(
                    propagatorSettings_->getInitialStates( ), 0 );

        // Check if re-integration of variational equations is requested
        if( areVariationalEquationsToBeIntegrated )
        {
            // Integrate variational and state equations.
            this->integrateVariationalAndDynamicalEquations( propagatorSettings_->getInitialStates( ), 1 );
        }
        else
        {
            this->integrateDynamicalEquationsOfMotionOnly( propagatorSettings_->getInitialStates( ) );
        }
    }

    //! Function to get the state transition matric interface object.
    /*!
     *  Function to get the state transition matric interface object.
     *  \return The state transition matric interface object.
     */
    boost::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > getStateTransitionMatrixInterface( )
    {
        return stateTransitionInterface_;
    }


protected:


    //! Create initial matrix of numerical soluation to variational + dynamical equations.
    /*!
     *  Create initial matrix of numerical soluation to variational + dynamical equations. The structure of the matrix is
     *  [Phi;S;y], with Phi the state transition matrix, S the sensitivity matrix y the state vector.
     *  \param initialStateEstimate vector of initial state (position/velocity) of bodies to be integrated numerically.
     *  order determined by order of bodiesToIntegrate_.
     *  \return Initial matrix of numerical soluation to variation + state equations.
     */
    MatrixType createInitialConditions( const VectorType initialStateEstimate )
    {
        // Initialize initial conditions to zeros.
        MatrixType varSystemInitialState = MatrixType( stateTransitionMatrixSize_,
                                                       parameterVectorSize_ + 1 ).setZero( );

        // Set initial state transition matrix to identity
        varSystemInitialState.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ).setIdentity( );

        // Set initial body states to current estimate of initial body states.
        varSystemInitialState.block( 0, parameterVectorSize_,
                                     stateTransitionMatrixSize_, 1 ) = initialStateEstimate;

        return varSystemInitialState;
    }

    //! Create initial matrix of numerical soluation to variational equations
    /*!
     *  Create initial matrix of numerical soluation to variational equations, with structure [Phi;S]. Initial state
     *  transition matrix Phi is identity matrix. Initial sensitivity matrix S is all zeros.
     *  \return Initial matrix solution to variational equations.
     */
    Eigen::MatrixXd createInitialVariationalEquationsSolution( )
    {
        // Initialize initial conditions to zeros.
        Eigen::MatrixXd varSystemInitialState = Eigen::MatrixXd::Zero(
                    stateTransitionMatrixSize_, parameterVectorSize_ );

        // Set initial state transition matrix to identity
        varSystemInitialState.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ).setIdentity( );

        return varSystemInitialState;
    }

    //! Object containing all parameters that are to be estimated and their current  settings and values.
    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate_ ;

    //! Map of bodies (with names) of all bodies in integration.
    simulation_setup::NamedBodyMap bodyMap_;

    //! Settings for propagation of equations of motion.
    boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings_;

    //! Settings for numerical integrator of combined propagation of variational equations and equations of motion.
    boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings_;

    //! Settings for numerical integrator when integrating only variational equations.
    boost::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings_;

    //! Size (rows and columns are equal) of state transition matrix.
    int stateTransitionMatrixSize_;

    //! Number of rows in sensitivity matrix
    int parameterVectorSize_;

    //! Boolean to determine whether to clear the raw numerical solution member variables after propagation
    /*!
     *  Boolean to determine whether to clear the raw numerical solution member variables after propagation
     *  and resetting of state transition interface.
     */
    bool clearNumericalSolution_;

    //! Object used for interpolating numerical results of state transition and sensitivity matrix.
    boost::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface_;  

    //! Object used to compute the full state derivative in equations of motion and variational equations.
    /*!
     *  Object used to compute the full state derivative in equations of motion and variational equations,
     *  including relevant updates of environment from current state and time. Object may be used for
     *  either full or separate propagation of equations.
     */
    boost::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > dynamicsStateDerivative_;

};

//! Function to separate the time histories of the sensitivity and state transition matrices from a full numerical solution.
/*!
 *  Function to separate the time histories of the sensitivity and state transition matrices from a full numerical solution,
 *  in which the solution is represented as a single matrix block per time value.
 *  NOTE: numericalIntegrationResult contents are deleted by this function (all information is conserved in
 *  variationalEquationsSolution.
 *  \param numericalIntegrationResult Full time history from which separate matrix histories are to be retrieved.
 *  \param variationalEquationsSolution Vector of two matrix histories (returned by reference). First vector entry
 *  is state transition matrix history, second entry is sensitivity matrix history.
 *  \param stateTransitionStartIndices First row and column (first and second) of state transition matrix in entries of
 *  numericalIntegrationResult.
 *  \param sensitivityStartIndices First row and column (first and second) of sensitivity matrix in entries of
 *  numericalIntegrationResult.
 *  \param stateTransitionMatrixSize Size (rows and columns are equal) of state transition matrix.
 *  \param parameterSetSize Number of rows in sensitivity matrix
 */
template< typename TimeType, typename StateScalarType >
void setVariationalEquationsSolution(
        std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >&
        numericalIntegrationResult,
        std::vector< std::map< double, Eigen::MatrixXd > >& variationalEquationsSolution,
        const std::pair< int, int > stateTransitionStartIndices,
        const std::pair< int, int > sensitivityStartIndices,
        const int stateTransitionMatrixSize,
        const int parameterSetSize )
{
    variationalEquationsSolution.clear( );
    variationalEquationsSolution.resize( 2 );

    for( typename std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > >::iterator
         integrationIterator = numericalIntegrationResult.begin( );
         integrationIterator != numericalIntegrationResult.end( ); )
    {
        // Set result for state transition matrix in each time step.
        variationalEquationsSolution[ 0 ][ integrationIterator->first ] =
                ( integrationIterator->second.block( stateTransitionStartIndices.first, stateTransitionStartIndices.second,
                                                     stateTransitionMatrixSize,
                                                     stateTransitionMatrixSize ) ).template cast< double >( );

        // Set result for sensitivity matrix in each time step.
        variationalEquationsSolution[ 1 ][ integrationIterator->first ] =
                ( integrationIterator->second.block( sensitivityStartIndices.first, sensitivityStartIndices.second,
                                                     stateTransitionMatrixSize,
                                                     parameterSetSize -
                                                     stateTransitionMatrixSize ) ).template cast< double >( );
        numericalIntegrationResult.erase( integrationIterator++ );
    }
}

//! Function to create interpolators for state transition and sensitivity matrices from numerical results.
/*!
 * Function to create interpolators for state transition and sensitivity matrices from numerical results.
 * \param stateTransitionMatrixInterpolator Interpolator object for state transition matrix (returned by reference).
 * \param sensitivityMatrixInterpolator Interpolator object for sensitivity matrix (returned by reference).
 * \param variationalEquationsSolution Vector of two matrix histories. First vector entry
 *  is state transition matrix history, second entry is sensitivity matrix history.
 * \param clearRawSolution Boolean denoting whether to clear entries of variationalEquationsSolution after creation
 * of interpolators.
 */
void createStateTransitionAndSensitivityMatrixInterpolator(
        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >&
        stateTransitionMatrixInterpolator,
        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >&
        sensitivityMatrixInterpolator,
        std::vector< std::map< double, Eigen::MatrixXd > >& variationalEquationsSolution,
        const bool clearRawSolution = 1 );

//! Function to check the consistency between propagation settings of equations of motion, and estimated parameters.
/*!
 *  Function to check the consistency between propagation settings of equations of motion, and estimated parameters.
 *  In particular, it is presently required that the set of propagated states is equal to the set of estimated states.
 *  \param propagatorSettings Settings for propagation of equations of motion.
 *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current
 *  settings and values.
 */
template< typename StateScalarType = double, typename TimeType = double, typename ParameterType = double >
bool checkPropagatorSettingsAndParameterEstimationConsistency(
        const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate )
{
    bool isInputConsistent = 1;

    // Check type of dynamics
    switch( propagatorSettings->stateType_ )
    {
    case transational_state:
    {
        boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > translationalPropagatorSettings =
                boost::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );

        // Retrieve estimated and propagated translational states, and check equality.
        std::vector< std::string > propagatedBodies = translationalPropagatorSettings->bodiesToIntegrate_;
        std::vector< std::string > estimatedBodies = estimatable_parameters::getListOfBodiesWithTranslationalStateToEstimate(
                    parametersToEstimate );
        if( propagatedBodies.size( ) != estimatedBodies.size( ) )
        {
            std::cerr<<"Error, propagated and estimated body vector sizes are inconsistent"<<
                       propagatedBodies.size( )<<" "<<estimatedBodies.size( )<<std::endl;
            isInputConsistent = 0;
        }
        else
        {
            for( unsigned int i = 0; i < propagatedBodies.size( ); i++ )
            {
                if( propagatedBodies.at( i ) != estimatedBodies.at( i ) )
                {
                    std::cerr<<"Error, propagated and estimated body vectors inconsistent at index "<<i<<": "<<
                               propagatedBodies.at( i )<<" "<<estimatedBodies.at( i )<<std::endl;
                    isInputConsistent = 0;
                }
            }

        }
        break;
    }
    default:
        std::cerr<<"Error, cannot yet check consistency of propagator settings for type: "<<
                   propagatorSettings->stateType_<<std::endl;
    }
    return isInputConsistent;
}

//! Class to manage and execute the numerical integration of variational equations of a dynamical system in a single arc.
/*!
 *  Class to manage and execute the numerical integration of variational equations of a dynamical system, in addition
 *  to the dynamics itself, in a single arc: i.e. the governing equations a single initial time, and are propagated once
 *  for the full prescribed time interval. This is in contrast to multi-arc dynamics, where the time interval is cut into
 *  pieces. In this class, the governing equations are set once, but can be re-integrated for
 *  different initial conditions using the same instance of the class.
 */
template< typename StateScalarType = double, typename TimeType = double, typename ParameterType = double >
class SingleArcVariationalEquationsSolver: public VariationalEquationsSolver< StateScalarType, TimeType, ParameterType >
{
public:

    //! Local typedefs for vector and matrix of given scalar type
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > MatrixType;
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > VectorType;

    //! Base class using statements
    using VariationalEquationsSolver< StateScalarType, TimeType, ParameterType >::parametersToEstimate_;
    using VariationalEquationsSolver< StateScalarType, TimeType, ParameterType >::bodyMap_;
    using VariationalEquationsSolver< StateScalarType, TimeType, ParameterType >::dynamicsStateDerivative_;
    using VariationalEquationsSolver< StateScalarType, TimeType, ParameterType >::propagatorSettings_;
    using VariationalEquationsSolver< StateScalarType, TimeType, ParameterType >::integratorSettings_;
    using VariationalEquationsSolver< StateScalarType, TimeType, ParameterType >::stateTransitionMatrixSize_;
    using VariationalEquationsSolver< StateScalarType, TimeType, ParameterType >::parameterVectorSize_;
    using VariationalEquationsSolver< StateScalarType, TimeType, ParameterType >::variationalOnlyIntegratorSettings_;
    using VariationalEquationsSolver< StateScalarType, TimeType, ParameterType >::stateTransitionInterface_;

    //! Constructor
    /*!
     *  Constructor, sets up object for automatic evaluation and numerical integration of variational equations and
     *  equations of motion.
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
     */
    SingleArcVariationalEquationsSolver(
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const bool integrateDynamicalAndVariationalEquationsConcurrently = 1,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings
            = boost::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
            const bool clearNumericalSolution = 1,
            const bool integrateEquationsOnCreation = 1 ):
        VariationalEquationsSolver< StateScalarType, TimeType, ParameterType >(
            bodyMap, integratorSettings, propagatorSettings, parametersToEstimate,
            variationalOnlyIntegratorSettings, clearNumericalSolution )
    {
        // Check input consistency
        if( !checkPropagatorSettingsAndParameterEstimationConsistency< StateScalarType, TimeType, ParameterType >(
                    propagatorSettings, parametersToEstimate ) )
        {
            throw std::runtime_error(
                        "Error when making single arc variational equations solver, estimated and propagated bodies are inconsistent" );
        }
        else
        {
            // Create simulation object for dynamics only.
            dynamicsSimulator_ =  boost::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                        bodyMap, integratorSettings, propagatorSettings, false, clearNumericalSolution );
            dynamicsStateDerivative_ = dynamicsSimulator_->getDynamicsStateDerivative( );

            // Create state derivative partials
            std::map< IntegratedStateType, orbit_determination::partial_derivatives::StateDerivativePartialsMap >
                    stateDerivativePartials =
                    orbit_determination::partial_derivatives::createStateDerivativePartials
                    < StateScalarType, TimeType, ParameterType >(
                        dynamicsStateDerivative_->getStateDerivativeModels( ), bodyMap, parametersToEstimate );

            // Create variational equations objects.
            variationalEquationsObject_ = boost::make_shared< VariationalEquations >(
                        stateDerivativePartials, parametersToEstimate_,
                        dynamicsStateDerivative_->getStateTypeStartIndices( ) );
            dynamicsStateDerivative_->addVariationalEquations( variationalEquationsObject_ );

            // Resize solution of variational equations to 2 (state transition and sensitivity matrices)
            variationalEquationsSolution_.resize( 2 );

            // Integrate variational equations from initial state estimate.
            if( integrateEquationsOnCreation )
            {
                if( integrateDynamicalAndVariationalEquationsConcurrently )
                {
                    integrateVariationalAndDynamicalEquations( propagatorSettings->getInitialStates( ), 1 );
                }
                else
                {
                    integrateVariationalAndDynamicalEquations( propagatorSettings->getInitialStates( ), 0 );
                }
            }
        }
    }

    //! Destructor
    ~SingleArcVariationalEquationsSolver( ){ }

    //! Function to integrate equations of motion only.
    /*!
     *  Function to integrate equations of motion only (in single arc).  If dynamical
     *  solution is to be processed, the environment is also updtaed to teh new solution.
     *  \param initialStateEstimate Initial state of the equations of motion that is to be used (in same order as in
     *  parametersToEstimate_)
     */
    void integrateDynamicalEquationsOfMotionOnly(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialStateEstimate )
    {
        dynamicsStateDerivative_->setPropagationSettings( std::vector< IntegratedStateType >( ), 1, 0 );
        dynamicsSimulator_->integrateEquationsOfMotion( initialStateEstimate );
    }

    //! Function to integrate variational equations and equations of motion.
    /*!
     *  Function to integrate variational equations and equations of motion (in single arc). At the end of this function,
     *  the stateTransitionInterface_ is reset with the new state transition and sensitivity matrices. If dynamical
     *  solution is to be processed, the environment is also updtaed to the new solution.
     *  \param initialStateEstimate Initial state of the equations of motion that is to be used (in same order as in
     *  parametersToEstimate_).
     *  \param integrateEquationsConcurrently Variable determining whether the equations of motion are to be
     *  propagated concurrently with variational equations of motion (if true), or before variational equations (if false).
     */
    void integrateVariationalAndDynamicalEquations(
            const VectorType& initialStateEstimate, const bool integrateEquationsConcurrently )
    {
        variationalEquationsSolution_[ 0 ].clear( );
        variationalEquationsSolution_[ 1 ].clear( );


        if( integrateEquationsConcurrently )
        {

            // Create initial conditions from new estimate.
            MatrixType initialVariationalState = this->createInitialConditions(
                        dynamicsStateDerivative_->convertFromOutputSolution(
                            initialStateEstimate, integratorSettings_->initialTime_ ) );


            // Integrate variational and state equations.
            dynamicsStateDerivative_->setPropagationSettings( std::vector< IntegratedStateType >( ), 1, 1 );
            std::map< TimeType, Eigen::VectorXd > dependentVariableHistory;
            std::map< TimeType, MatrixType > rawNumericalSolution;
            integrateEquations< MatrixType, TimeType >(
                        dynamicsSimulator_->getStateDerivativeFunction( ), rawNumericalSolution,
                        initialVariationalState, integratorSettings_,
                        boost::bind( &PropagationTerminationCondition::checkStopCondition,
                                     dynamicsSimulator_->getPropagationTerminationCondition( ), _1 ),
                        dependentVariableHistory );

            std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > equationsOfMotionNumericalSolution;
            utilities::createVectorBlockMatrixHistory(
                        rawNumericalSolution, equationsOfMotionNumericalSolution,
                        std::make_pair( 0, parameterVectorSize_ ), stateTransitionMatrixSize_ );

            equationsOfMotionNumericalSolution = convertNumericalStateSolutionsToOutputSolutions(
                        equationsOfMotionNumericalSolution, dynamicsStateDerivative_ );
            dynamicsSimulator_->manuallySetAndProcessRawNumericalEquationsOfMotionSolution(
                        equationsOfMotionNumericalSolution );

            // Reset solution for state transition and sensitivity matrices.
            setVariationalEquationsSolution< TimeType, StateScalarType >(
                        rawNumericalSolution, variationalEquationsSolution_,
                        std::make_pair( 0, 0 ), std::make_pair( 0, stateTransitionMatrixSize_ ),
                        stateTransitionMatrixSize_, parameterVectorSize_ );
        }
        else
        {

            dynamicsStateDerivative_->setPropagationSettings( std::vector< IntegratedStateType >( ), 1, 0 );
            dynamicsSimulator_->integrateEquationsOfMotion( initialStateEstimate );

            // Integrate variational equations.
            dynamicsStateDerivative_->setPropagationSettings( boost::assign::list_of( transational_state ), 0, 1 );
            Eigen::MatrixXd initialVariationalState = this->createInitialVariationalEquationsSolution( );
            std::map< double, Eigen::MatrixXd > rawNumericalSolution;
            std::map< TimeType, Eigen::VectorXd > dependentVariableHistory;

            integrateEquations< Eigen::MatrixXd, double >(
                        dynamicsSimulator_->getDoubleStateDerivativeFunction( ), rawNumericalSolution, initialVariationalState,
                        variationalOnlyIntegratorSettings_,
                        boost::bind( &PropagationTerminationCondition::checkStopCondition,
                                     dynamicsSimulator_->getPropagationTerminationCondition( ), _1 ),
                        dependentVariableHistory );

            setVariationalEquationsSolution< double, double >(
                        rawNumericalSolution, variationalEquationsSolution_, std::make_pair( 0, 0 ),
                        std::make_pair( 0, stateTransitionMatrixSize_ ),
                        stateTransitionMatrixSize_, parameterVectorSize_ );

        }

        // Reset solution for state transition and sensitivity matrices.
        resetVariationalEquationsInterpolators( );

    }

    //! Function to return the numerical solution history of numerically integrated variational equations.
    /*!
     *  Function to return the numerical solution history of numerically integrated variational equations.
     *  \return Vector of mapa of state transition matrix history (first vector entry)
     *  and sensitivity matrix history (second vector entry)
     */
    std::vector< std::map< double, Eigen::MatrixXd > >& getNumericalVariationalEquationsSolution( )
    {
        return variationalEquationsSolution_;
    }

    //! Function to return object used for numerically propagating and managing the solution of the equations of motion.
    /*!
     * Function to return object used for numerically propagating and managing the solution of the equations of motion.
     * \return Object used for numerically propagating and managing the solution of the equations of motion.
     */
    boost::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > getDynamicsSimulator( )
    {
        return dynamicsSimulator_;
    }

protected:

private:


    //! Reset solutions of variational equations.
    /*!
     *  Reset solutions of variational equations (stateTransitionMatrixInterpolator_ and sensitivityMatrixInterpolator_),
     *  i.e. use numerical integration results to create new look-up tables
     *  and interpolators of state transition and sensitivity matrix through the createInterpolatorsForVariationalSolution
     *  function
     */
    void resetVariationalEquationsInterpolators( )
    {
        using namespace interpolators;
        using namespace utilities;

        // Create interpolators.
        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
                stateTransitionMatrixInterpolator;
        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >
                sensitivityMatrixInterpolator;
        createStateTransitionAndSensitivityMatrixInterpolator(
                    stateTransitionMatrixInterpolator, sensitivityMatrixInterpolator, variationalEquationsSolution_,
                    this->clearNumericalSolution_ );

        // Create (if non-existent) or reset state transition matrix interface
        if( stateTransitionInterface_ == NULL )
        {
            stateTransitionInterface_ = boost::make_shared< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                        stateTransitionMatrixInterpolator, sensitivityMatrixInterpolator,
                        propagatorSettings_->getStateSize( ), parameterVectorSize_ );
        }
        else
        {
            boost::dynamic_pointer_cast< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                        stateTransitionInterface_ )->updateMatrixInterpolators(
                        stateTransitionMatrixInterpolator, sensitivityMatrixInterpolator );
        }
    }

    //! Object used for numerically propagating and managing the solution of the equations of motion.
    boost::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > dynamicsSimulator_;

    //!  Object that is used to evaluate the variational equations at the given state and time.
    boost::shared_ptr< VariationalEquations > variationalEquationsObject_;

    //! Map of history of numerically integrated variational equations.
    /*!
     *  Map of history of numerically integrated variational equations. Key of map denotes time, values are
     *  state transition matrix Phi (first vector entry) and sensitivity matrix S (second vector entry)
     */
    std::vector< std::map< double, Eigen::MatrixXd > > variationalEquationsSolution_;

};

}

}




#endif // TUDAT_VARIATIONALEQUATIONSSOLVER_H
