#ifndef DYNAMICSSIMULATOR3_H
#define DYNAMICSSIMULATOR3_H

#include <boost/make_shared.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple_io.hpp>

#include "Tudat/Basics/utilities.h"

#include "Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Mathematics/NumericalIntegrators/euler.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Mathematics/Interpolators/interpolator.h"
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h"
#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/Propagators/variationalEquations.h"
#include "Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/dynamicsStateDerivativeModel.h"
#include "Tudat/Astrodynamics/Propagators/stateTransitionMatrixInterface.h"
#include "Tudat/Astrodynamics/Propagators/dynamicsSimulator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/SimulationSetup/createStateDerivativePartials.h"

namespace tudat
{

namespace propagators
{


template< typename StateScalarType = double, typename ParameterType = double >
class VariationalEquationsSolverBase
{
public:
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > VectorType;

    VariationalEquationsSolverBase(
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate ):
        parametersToEstimate_( parametersToEstimate ){ }

    virtual ~VariationalEquationsSolverBase( ){ }

    virtual void resetParameterEstimate( const Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > newParameterEstimate,
                                         const bool areVariationalEquationsToBeIntegrated = true ) = 0;

    virtual void integrateVariationalAndDynamicalEquations(
            const VectorType& initialStateEstimate, const bool integrateEquationsConcurrently ) = 0;

    virtual void integrateDynamicalEquationsOfMotionOnly(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialCartesianStates ) = 0;


    virtual boost::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > getStateTransitionMatrixInterface( ) = 0;

    //! Function to get the list of objects representing the parameters that are to be integrated.
    /*!
     *  Function to get the list of objects representing the parameters that are to be integrated.
     *  \return List of objects representing the parameters that are to be integrated.
     */
    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > getParametersToEstimate( )
    {
        return parametersToEstimate_;
    }

protected:
    //! Object containing all parameters that are to be estimated.
    /*!
     *  Object containing all parameters that are to be estimated and their current settings and values.
     */
    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate_ ;
};

//! Class to manage integration of equations of motion and variational equations.
/*!
 *  Class to manage integration of equations of motion and variational equations. Derives from DynamicsSimulator that deals solely with
 *  equations of motion
 */
template< typename StateScalarType = double, typename TimeType = double, typename ParameterType = double >
class VariationalEquationsSolver : public VariationalEquationsSolverBase< StateScalarType, ParameterType >
{
public:

    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > MatrixType;
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > VectorType;

    using VariationalEquationsSolverBase< StateScalarType, ParameterType >::parametersToEstimate_;

    //! Constructor
    /*!
     *  Constructor, sets up object for automatic evaluation and numerical integration of variational equations and equations of motion.
     *  \param bodyMap Map of bodies (with names) of all bodies in integration.
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagatorSettings Settings for propagator.
     *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current settings and values.
     *  \param bodiesToEstime List of bodies to estimate
     */
    VariationalEquationsSolver(
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings =
            boost::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
            const bool clearNumericalSolution = 1 ):
        VariationalEquationsSolverBase< StateScalarType, ParameterType >( parametersToEstimate ),
        bodyMap_( bodyMap ),
        propagatorSettings_( propagatorSettings ), integratorSettings_( integratorSettings ),
        variationalOnlyIntegratorSettings_( variationalOnlyIntegratorSettings ),
        stateTransitionMatrixSize_( parametersToEstimate->getInitialDynamicalStateParameterSize( ) ),
        parameterVectorSize_( parametersToEstimate_->getParameterSetSize( ) ),
        clearNumericalSolution_( clearNumericalSolution )
    { }

    //! Destructor
    /*!
     *  Destructor
     */
    virtual ~VariationalEquationsSolver( ){ }

    //! Function to reset parameter estimate and re-integrate equations of motion and, if desired, variational equations.
    /*!
     *  Function to reset parameter estimate and re-integrate equations of motion and, if desired, variational equations
     *  using the new physical parameters/body initial states.
     *  \param newParameterEstimate New estimate of parameters that are to be estimated. First entries are initial
     *  states of bodies to be integrated ( order determined by bodiesToIntegrate_ ),
     *  followed by physical parameters that are to be estimated ( order determined by parametersToEstimate_).
     */
    void resetParameterEstimate( const Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > newParameterEstimate,
                                 const bool areVariationalEquationsToBeIntegrated = true )

    {
        // Reset values of parameters.
        parametersToEstimate_->template resetParameterValues< ParameterType >( newParameterEstimate );
        propagatorSettings_->resetInitialStates( estimatable_parameters::getInitialStateVectorOfBodiesToEstimate( parametersToEstimate_ ) );

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

    boost::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > getStateTransitionMatrixInterface( )
    {
        return stateTransitionInterface_;
    }


protected:


    //! Create initial matrix of numerical soluation to variation + state equations.
    /*!
     *  Create initial matrix of numerical soluation to variation + state equations. The structure of the matrix is:
     * [Phi(t,t0)_{nb*6,nb*6};S(t)_{6*nb,np};[y0_{6,1}...ynb_{6,1} ;0_{(nb-1)*6,nb}] ]. Subscripts denote the size of
     *  the componets, Phi the state transition matrix, S the sensitivity matrix yi the state of body i, nb the number of bodies
     *  and np the number of parameters to be estimated.
     *  \param initialStateEstimate vector of initial state (position/velocity) of bodies to be integrated numerically.
     *  order determined by order of bodiesToIntegrate_.
     *  \return Iinitial matrix of numerical soluation to variation + state equations.
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

    Eigen::MatrixXd createInitialVariationalEquationsSolution( )
    {
        // Initialize initial conditions to zeros.
        Eigen::MatrixXd varSystemInitialState = Eigen::MatrixXd::Zero(
                    stateTransitionMatrixSize_, parameterVectorSize_ );

        // Set initial state transition matrix to identity
        varSystemInitialState.block( 0, 0, stateTransitionMatrixSize_, stateTransitionMatrixSize_ ).setIdentity( );

        return varSystemInitialState;
    }

    boost::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface_;

    boost::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings_;

    int stateTransitionMatrixSize_;

    int parameterVectorSize_;


    boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings_;

    boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings_;

    boost::shared_ptr< DynamicsStateDerivativeModel< TimeType, StateScalarType > > dynamicsStateDerivative_;

    simulation_setup::NamedBodyMap bodyMap_;

    bool clearNumericalSolution_;
};

template< typename MapTimeType, typename MapStateScalarType >
void setVariationalEquationsSolution(
        std::map< MapTimeType, Eigen::Matrix< MapStateScalarType, Eigen::Dynamic, Eigen::Dynamic > >& numericalIntegrationResult,
        std::vector< std::map< double, Eigen::MatrixXd > >& variationalEquationsSolution,
        const std::pair< int, int > stateTransitionStartIndices,
        const std::pair< int, int > sensitivityStartIndices,
        const int stateTransitionMatrixSize,
        const int parameterSetSize )
{
    for( typename std::map< MapTimeType, Eigen::Matrix< MapStateScalarType, Eigen::Dynamic, Eigen::Dynamic > >::iterator
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


void createStateTransitionAndSensitivityMatrixInterpolator(
        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >& stateTransitionMatrixInterpolator,
        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > >& sensitivityMatrixInterpolator,
        std::vector< std::map< double, Eigen::MatrixXd > >& variationalEquationsSolution,
        const bool clearRawSolution = 1 );

template< typename StateScalarType = double, typename TimeType = double, typename ParameterType = double >
bool checkPropagatorSettingsAndParameterEstimationConsistency(
        const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate )
{
    bool isInputConsistent = 1;
    switch( propagatorSettings->stateType_ )
    {
    case transational_state:
    {
        boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > translationalPropagatorSettings =
                boost::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );
        std::vector< std::string > propagatedBodies = translationalPropagatorSettings->bodiesToIntegrate_;
        std::vector< std::string > estimatedBodies = estimatable_parameters::getListOfBodiesToEstimate(
                    parametersToEstimate );
        if( propagatedBodies.size( ) != estimatedBodies.size( ) )
        {
            std::cerr<<"Error, propagated and estimated body vector sizes are inconsistent"<<propagatedBodies.size( )<<" "<<estimatedBodies.size( )<<std::endl;
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
        std::cerr<<"Error, cannot yet check consistency of propagator settings for type: "<<propagatorSettings->stateType_<<std::endl;
    }
    return isInputConsistent;
}

template< typename StateScalarType = double, typename TimeType = double, typename ParameterType = double >
class SingleArcVariationalEquationsSolver: public VariationalEquationsSolver< StateScalarType, TimeType, ParameterType >
{
public:

    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > MatrixType;
    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > VectorType;

    using VariationalEquationsSolverBase< StateScalarType, ParameterType >::parametersToEstimate_;

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
     *  Constructor, sets up object for automatic evaluation and numerical integration of variational equations and equations of motion.
     *  \param bodyMap Map of bodies (with names) of all bodies in integration.
     *  \param integratorSettings Settings for numerical integrator.
     *  \param propagatorSettings Settings for propagator.
     *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current settings and values.
     *  \param bodiesToEstime List of bodies to estimate
     */
    SingleArcVariationalEquationsSolver(
            const simulation_setup::NamedBodyMap& bodyMap,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< TimeType > > integratorSettings,
            const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const bool integrateDynamicalAndVariationalEquationsConcurrently = 1,
            const boost::shared_ptr< numerical_integrators::IntegratorSettings< double > > variationalOnlyIntegratorSettings =
            boost::shared_ptr< numerical_integrators::IntegratorSettings< double > >( ),
            const bool clearNumericalSolution = 1,
            const bool integrateEquationsOnCreation = 1 ):
        VariationalEquationsSolver< StateScalarType, TimeType, ParameterType >(
            bodyMap, integratorSettings, propagatorSettings, parametersToEstimate, variationalOnlyIntegratorSettings, clearNumericalSolution )
    {
        if( !checkPropagatorSettingsAndParameterEstimationConsistency< StateScalarType, TimeType, ParameterType >(
                    propagatorSettings, parametersToEstimate ) )
        {
            std::cerr<<"Error when making single arc variational equations solver, estimated and propagated bodies are inconsistent"<<std::endl;
        }
        else
        {
            dynamicsSimulator_ =  boost::make_shared< SingleArcDynamicsSimulator< StateScalarType, TimeType > >(
                        bodyMap, integratorSettings, propagatorSettings, false, clearNumericalSolution );
            dynamicsStateDerivative_ = dynamicsSimulator_->getDynamicsStateDerivative( );

            // Create variational equations objects.
            std::map< IntegratedStateType, orbit_determination::partial_derivatives::StateDerivativePartialsMap > stateDerivativePartials =
                    orbit_determination::partial_derivatives::createStateDerivativePartials< StateScalarType, TimeType, ParameterType >(
                        dynamicsStateDerivative_->getStateDerivativeModels( ), bodyMap, parametersToEstimate );

            variationalEquationsObject_ = boost::make_shared< VariationalEquations >(
                        stateDerivativePartials, parametersToEstimate_, dynamicsStateDerivative_->getStateTypeStartIndices( ) );

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
    /*!
     *  Destructor
     */
    ~SingleArcVariationalEquationsSolver( ){ }

    void integrateDynamicalEquationsOfMotionOnly(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& initialStateEstimate )
    {
        std::cout<<"resetting only dynamics"<<std::endl;

        dynamicsStateDerivative_->setPropagationSettings( std::vector< IntegratedStateType >( ), 1, 0 );
        dynamicsSimulator_->integrateEquationsOfMotion( initialStateEstimate );
    }

    void integrateVariationalAndDynamicalEquations(
            const VectorType& initialStateEstimate, const bool integrateEquationsConcurrently )
    {
        variationalEquationsSolution_[ 0 ].clear( );
        variationalEquationsSolution_[ 1 ].clear( );


        if( integrateEquationsConcurrently )
        {

            // Create initial conditions from new estimate.
            MatrixType initialVariationalState = this->createInitialConditions(
                        dynamicsStateDerivative_->convertFromOutputSolution( initialStateEstimate, integratorSettings_->initialTime_ ) );

            // Integrate variational and state equations.
            dynamicsStateDerivative_->setPropagationSettings( std::vector< IntegratedStateType >( ), 1, 1 );
            std::map< TimeType, MatrixType > rawNumericalSolution =
                    integrateEquations< MatrixType, TimeType >(
                        dynamicsSimulator_->getStateDerivativeFunction( ), initialVariationalState, integratorSettings_ );

            std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > equationsOfMotionNumericalSolution;
            utilities::createVectorBlockMatrixHistory(
                        rawNumericalSolution, equationsOfMotionNumericalSolution,
                        std::make_pair( 0, parameterVectorSize_ ), stateTransitionMatrixSize_ );

            equationsOfMotionNumericalSolution = convertNumericalStateSolutionsToOutputSolutions(
                        equationsOfMotionNumericalSolution, dynamicsStateDerivative_ );
            dynamicsSimulator_->manuallySetAndProcessRawNumericalEquationsOfMotionSolution( equationsOfMotionNumericalSolution );

            // Reset solution for state transition and sensitivity matrices.
            setVariationalEquationsSolution< TimeType, StateScalarType >(
                        rawNumericalSolution, variationalEquationsSolution_, std::make_pair( 0, 0 ), std::make_pair( 0, stateTransitionMatrixSize_ ),
                        stateTransitionMatrixSize_, parameterVectorSize_ );
        }
        else
        {

            dynamicsStateDerivative_->setPropagationSettings( std::vector< IntegratedStateType >( ), 1, 0 );
            dynamicsSimulator_->integrateEquationsOfMotion( initialStateEstimate );

            // Integrate variational equations.
            dynamicsStateDerivative_->setPropagationSettings( boost::assign::list_of( transational_state ), 0, 1 );
            Eigen::MatrixXd initialVariationalState = this->createInitialVariationalEquationsSolution( );
            std::map< double, Eigen::MatrixXd > rawNumericalSolution = integrateEquations< Eigen::MatrixXd, double >(
                        dynamicsSimulator_->getDoubleStateDerivativeFunction( ), initialVariationalState, variationalOnlyIntegratorSettings_ );
            setVariationalEquationsSolution< double, double >(
                        rawNumericalSolution, variationalEquationsSolution_, std::make_pair( 0, 0 ), std::make_pair( 0, stateTransitionMatrixSize_ ),
                        stateTransitionMatrixSize_, parameterVectorSize_ );

        }

        // Reset solution for state transition and sensitivity matrices.
        resetVariationalEquationsInterpolators( );

    }

    //! Function to return the numerical solution history of numerically integrated variational equations.
    /*!
     *  Function to return the numerical solution history of numerically integrated variational equations. Key of map denotes time, values are concatenated
     *  matrices state transition Phi and sensitivity S: [Phi;S].
     *  \return Map of state history of numerically integrated bodies.
     */
    std::vector< std::map< double, Eigen::MatrixXd > >& getNumericalVariationalEquationsSolution( )
    {
        return variationalEquationsSolution_;
    }

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
     *  and interpolators of state transition and sensitivity matrix through the createInterpolatorsForVariationalSolution( )
     *  function
     */
    void resetVariationalEquationsInterpolators( )
    {
        using namespace interpolators;
        using namespace utilities;

        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > stateTransitionMatrixInterpolator;
        boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::MatrixXd > > sensitivityMatrixInterpolator;

        createStateTransitionAndSensitivityMatrixInterpolator(
                    stateTransitionMatrixInterpolator, sensitivityMatrixInterpolator, variationalEquationsSolution_, this->clearNumericalSolution_ );

        if( stateTransitionInterface_ == NULL )
        {
            stateTransitionInterface_ = boost::make_shared< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                        stateTransitionMatrixInterpolator, sensitivityMatrixInterpolator, propagatorSettings_->getStateSize( ),
                        parameterVectorSize_ );
        }
        else
        {
            boost::dynamic_pointer_cast< SingleArcCombinedStateTransitionAndSensitivityMatrixInterface >(
                        stateTransitionInterface_ )->updateMatrixInterpolators(
                        stateTransitionMatrixInterpolator, sensitivityMatrixInterpolator );
        }
    }

    boost::shared_ptr< SingleArcDynamicsSimulator< StateScalarType, TimeType > > dynamicsSimulator_;

    //!  Function to return the evaluated variational equations
    /*!
     *   Function to return the evaluated variational equations and update environment to current integration step.
     */
    boost::shared_ptr< VariationalEquations > variationalEquationsObject_;

    //! Map of history of numerically integrated variational equations.
    /*!
     *  Map of history of numerically integrated variational equations. Key of map denotes time, values are concatenated
     *  matrices state transition Phi and sensitivity S: [Phi;S].
     */
    std::vector< std::map< double, Eigen::MatrixXd > > variationalEquationsSolution_;

};

}

}




#endif // DYNAMICSMANAGER_H
