#ifndef TUDAT_VARIATIONALEQUATIONS_H
#define TUDAT_VARIATIONALEQUATIONS_H

#include <map>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"
#include "Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h"
#include "Tudat/Astrodynamics/Propagators/setNumericallyIntegratedStates.h"

namespace tudat
{

namespace propagators
{

//! Class from which the variational equations can be evaluated.
class VariationalEquations
{
    
public:
    
    //! Constructor of variation equations class.
    /*!
     * Constructor of variation equations class. Since the vehicle state must be integrated along with
     * the variational equations, an object calculating the state derivative is required.
     * \param stateDerivativePartialList List partials of state derivative models from which the variational equations
     * are set up. The key is the type of dynamics for which partials are taken, the values are StateDerivativePartialsMap
     * (see StateDerivativePartialsMap definition for details)
     *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current settings and
     *  values.
     * \param stateTypeStartIndices Start index (value) in vector of propagated state for each type of state (key)
     */
    template< typename ParameterType >
    VariationalEquations(
            const std::map< IntegratedStateType, orbit_determination::partial_derivatives::StateDerivativePartialsMap >
            stateDerivativePartialList,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const std::map< IntegratedStateType, int >& stateTypeStartIndices ):
        stateDerivativePartialList_( stateDerivativePartialList ), stateTypeStartIndices_( stateTypeStartIndices ),
        estimatedUnintegratedBodies_( std::vector< std::string >( ) )
    {
        dynamicalStatesToEstimate_ = estimatable_parameters::getListOfInitialDynamicalStateParametersEstimate< ParameterType >(
                    parametersToEstimate );
        
        // Get size of dynamical state to estimate
        numberOfParameterValues_ = estimatable_parameters::getSingleArcParameterSetSize( parametersToEstimate );
        totalDynamicalStateSize_ = 0;        
        for( std::map< IntegratedStateType, orbit_determination::partial_derivatives::StateDerivativePartialsMap >::iterator partialTypeIterator =
             stateDerivativePartialList_.begin( ); partialTypeIterator != stateDerivativePartialList_.end( ); partialTypeIterator++ )
        {
            
            if( dynamicalStatesToEstimate_.count( partialTypeIterator->first ) == 0 )
            {
                std::cerr<<"Error when making variational equations object, found no state to estimate of type "<<partialTypeIterator->first <<std::endl;
            }
            else if( dynamicalStatesToEstimate_.at( partialTypeIterator->first ).size( ) != partialTypeIterator->second.size( ) )
            {
                std::cerr<<"Error when making variational equations object, input partial list size is inconsistent"<<std::endl;
            }
            
            totalDynamicalStateSize_ += getSingleIntegrationSize( partialTypeIterator->first ) * partialTypeIterator->second.size( );
        }
        
        // Initialize matrices.
        currentMatrixDerivative_ = Eigen::MatrixXd::Zero( totalDynamicalStateSize_, numberOfParameterValues_ );
        currentLongMatrixDerivative_ = Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >::Zero( totalDynamicalStateSize_, numberOfParameterValues_ );
        variationalMatrix_ = Eigen::MatrixXd::Zero( totalDynamicalStateSize_, totalDynamicalStateSize_ );
        variationalParameterMatrix_ = Eigen::MatrixXd::Zero( totalDynamicalStateSize_, numberOfParameterValues_ - totalDynamicalStateSize_ );

        // Set parameter partial functions.
        setStatePartialFunctionList( );
        setTranslationalStatePartialFrameScalingFunctions( parametersToEstimate );
        setParameterPartialFunctionList( parametersToEstimate );
    }
    
    //! Calculates matrix containing partial derivatives of state derivatives w.r.t. body state.
    /*!
     *  Calculates matrix containing partial derivatives of state derivatives w.r.t. body state, i.e.
     *  first matrix in right hand side of Eq. (7.45) in (Montenbruck & Gill, 2000).
     *  \return Matrix containing partial derivatives of state derivative w.r.t. body state
     */
    void setBodyStatePartialMatrix( );

    //! Function to compute the contribution of the derivatives w.r.t. current states in the variational equations
    /*!
     *  Function to compute the contribution of the derivatives w.r.t. current states in the variational equations,
     *  e.g. first term in Eq. (7.45) in (Montenbruck & Gill, 2000).
     *  \param stateTransitionAndSensitivityMatrices Current combined state transition and sensitivity matric
     *  \param currentMatrixDerivative Matrix block which is to return (by reference) the given contribution to the
     *  variational equations.
     */
    template< typename StateScalarType >
    void getBodyInitialStatePartialMatrix(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& stateTransitionAndSensitivityMatrices,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > currentMatrixDerivative )
    {
        setBodyStatePartialMatrix( );

        // Add partials of body positions and velocities.
        currentMatrixDerivative.block( 0, 0, totalDynamicalStateSize_, numberOfParameterValues_ ) =
                ( variationalMatrix_ * stateTransitionAndSensitivityMatrices ).template cast< StateScalarType >( );
    }

    //! Calculates matrix containing partial derivatives of state derivatives w.r.t. parameters.
    /*!
     *  Calculates matrix containing partial derivatives of state derivatives  w.r.t. parameters, i.e.
     *  second matrix in rhs of Eq. (7.45) in (Montenbruck & Gill, 2000).
     *  \param currentMatrixDerivative Matrix block containing partial derivatives of accelerarion w.r.t. parameters
     *  (returned by reference).
     */
    template< typename StateScalarType >
    void getParameterPartialMatrix(
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > currentMatrixDerivative )
    {
        // Initialize matrix to zeros
        variationalParameterMatrix_.setZero( );


        // Iterate over all bodies undergoing accelerations for which initial condition is to be estimated.
        for( std::map< IntegratedStateType, std::vector< std::multimap< std::pair< int, int >,
             boost::function< void( Eigen::Block< Eigen::MatrixXd > ) > > > >::iterator typeIterator =
             parameterPartialList_.begin( ); typeIterator != parameterPartialList_.end( ); typeIterator++ )
        {
            int startIndex = stateTypeStartIndices_.at( typeIterator->first );
            int currentStateSize = getSingleIntegrationSize( typeIterator->first );
            int entriesToSkipPerEntry = currentStateSize - currentStateSize / getSingleIntegrationDifferentialEquationOrder( typeIterator->first );

            // Iterate over all bodies being estimated.
            for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
            {
                // Iterate over all parameter partial functions determined by setParameterPartialFunctionList( )
                for( functionIterator = typeIterator->second[ i ].begin( ); functionIterator != typeIterator->second[ i ].end( );
                     functionIterator++ )
                {
                    functionIterator->second(
                                variationalParameterMatrix_.block(
                                startIndex + entriesToSkipPerEntry + currentStateSize * i, functionIterator->first.first - totalDynamicalStateSize_,
                                currentStateSize - entriesToSkipPerEntry, functionIterator->first.second ) );
                }
            }

        }

        currentMatrixDerivative.block( 0, totalDynamicalStateSize_, totalDynamicalStateSize_, totalDynamicalStateSize_ - numberOfParameterValues_ ) +=
                variationalParameterMatrix_.template cast< StateScalarType >( );
    }
    
    //! Evaluates the complete variational equations.
    /*!
     *  Evaluates the complete variational equations at a given time and state transition matrix, sensitivity matrix and
     *  state (accessed indirectly). This function evaluates the complete Eq. (7.45) from (Montenbruck & Gill, 2000).
     *  \param time Current time
     *  \param stateTransitionAndSensitivityMatrices Combined state transition and sensitivity matrix.
     *  \param currentMatrixDerivative Variation equations result (returned by reference).
     */
    template< typename StateScalarType >
    void evaluateVariationalEquations(
            const double time, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& stateTransitionAndSensitivityMatrices,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > currentMatrixDerivative )
    {
        // Compute and add state partials.
        getBodyInitialStatePartialMatrix< StateScalarType >( stateTransitionAndSensitivityMatrices, currentMatrixDerivative );

        if( numberOfParameterValues_ > totalDynamicalStateSize_ )
        {
            // Add partials of parameters.
            getParameterPartialMatrix< StateScalarType >( currentMatrixDerivative );
        }
    }
    
    //! This function updates all state derivative models to the current time and state.
    /*!
     *  This function updates all state derivative models to the current time and state.
     *  \param currentTime Time to  which the system is to be updated.
     */
    void updatePartials( const double currentTime );
    
    //! Returns the number of parameter values.
    /*!
     *  Returns the number of parameter values (i.e. number of columns in state transition matrix).
     *  \return Number of parameter values.
     */
    double getNumberOfParameterValues( )
    {
        return numberOfParameterValues_;
    }
    
protected:
    
private:
    
    void setStatePartialFunctionList( );
        
    template< typename CurrentParameterType >
    void addParameterPartialToList(
            const std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< CurrentParameterType > > >& parameterList,
            const boost::shared_ptr< orbit_determination::partial_derivatives::StateDerivativePartial > partialObject,
            std::multimap< std::pair< int, int >, boost::function< void( Eigen::Block< Eigen::MatrixXd > ) > >& functionListOfBody,
            const int totalParameterVectorIndicesToSubtract = 0 )
    {
        using namespace orbit_determination::partial_derivatives;

        // Iterate over all parameters.
        for( typename std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< CurrentParameterType > > >::const_iterator
             parameterIterator = parameterList.begin( ); parameterIterator != parameterList.end( ); parameterIterator++ )
        {
            // Add current parameter to list of partials to be computed for current acceleration (if dependency exists)
            int functionToEvaluate =
                    partialObject->setParameterPartialUpdateFunction( parameterIterator->second );
            
            // If function is non-NULL, add to list
            if( functionToEvaluate != 0 )
            {
                // Make pair of indices for generating parameter partial matrix:
                //first is start column in matrix, second is number of entries (1 for double parameter)
                std::pair< int, int > indexPair = std::make_pair(
                            parameterIterator->first - totalParameterVectorIndicesToSubtract, functionToEvaluate );
                
                // Add to list.
                functionListOfBody.insert(
                            std::pair< std::pair< int, int >, boost::function< void( Eigen::Block< Eigen::MatrixXd > ) > >
                            ( indexPair, boost::bind(
                                  static_cast< void ( StateDerivativePartial::* )
                                  ( const boost::shared_ptr< estimatable_parameters::EstimatableParameter< CurrentParameterType > >,
                                    Eigen::Block< Eigen::MatrixXd > )>
                                  ( &StateDerivativePartial::getCurrentParameterPartial ), partialObject, parameterIterator->second, _1  ) ) );
            }
        }
    }
    
    //! This function creates the list of partial derivatives of the state w.r.t. parameter values.
    /*!
     *  This function creates the list of partial derivatives of the state w.r.t. parameter values. The function is called once by the constructor
     *  and the resulting functions are set as memebr variables. This prevents having to check whether an acceleration model depends on every parameter
     *  during every time step.
     *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current settings and values.
     */
    template< typename ParameterType >
    void setParameterPartialFunctionList(
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate )
    {
        // Get double parameters.
        std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > doubleParametersToEstimate =
                parametersToEstimate->getDoubleParameters( );
        std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > vectorParametersToEstimate =
                parametersToEstimate->getVectorParameters( );
        
        int totalParameterVectorIndicesToSubtract = parametersToEstimate->getInitialDynamicalStateParameterSize( ) -
                estimatable_parameters::getSingleArcInitialDynamicalStateParameterSetSize( parametersToEstimate );

        for( std::map< propagators::IntegratedStateType, orbit_determination::partial_derivatives::StateDerivativePartialsMap >::iterator
             stateDerivativeTypeIterator = stateDerivativePartialList_.begin( ); stateDerivativeTypeIterator != stateDerivativePartialList_.end( );
             stateDerivativeTypeIterator++ )
        {
            
            // Initialize vector of lists to correct size.
            parameterPartialList_[ stateDerivativeTypeIterator->first ].resize( stateDerivativeTypeIterator->second.size( ) );
            
            // Iterate over all bodies of which initial position is being estimated.
            for( unsigned int i = 0; i < stateDerivativeTypeIterator->second.size( ); i++ )
            {
                // Initialize list of parameter partial functions for single body.
                std::multimap< std::pair< int, int >, boost::function< void( Eigen::Block< Eigen::MatrixXd > ) > > functionListOfBody;
                
                // Iterate over all accelerations due to this body on current body.
                for( unsigned int j = 0; j < stateDerivativeTypeIterator->second.at( i ).size( ); j++ )
                {
                    addParameterPartialToList< double >(
                                doubleParametersToEstimate, stateDerivativeTypeIterator->second.at( i ).at( j ),
                                functionListOfBody, totalParameterVectorIndicesToSubtract );
                    addParameterPartialToList< Eigen::VectorXd >(
                                vectorParametersToEstimate, stateDerivativeTypeIterator->second.at( i ).at( j ),
                                functionListOfBody, totalParameterVectorIndicesToSubtract );
                }
                
                
                // Add generated parameter partial list of current body.
                parameterPartialList_[ stateDerivativeTypeIterator->first ][ i ] = functionListOfBody;
            }
        }
    }

    template< typename ParameterType >
    void setTranslationalStatePartialFrameScalingFunctions(
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate )
    {
        std::vector< boost::shared_ptr< estimatable_parameters::EstimatableParameter<
                Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
                parametersToEstimate->getEstimatedInitialStateParameters( );

        std::vector< std::string > propagatedBodies;
        std::vector< std::string > centralBodies;

        for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
        {
            if( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::initial_body_state )
            {
                propagatedBodies.push_back(
                            initialDynamicalParameters.at( i )->getParameterName( ).second.first );
                centralBodies.push_back( boost::dynamic_pointer_cast< estimatable_parameters::InitialTranslationalStateParameter< ParameterType > >(
                                             initialDynamicalParameters.at( i ) )->getCentralBody( ) );
            }
        }

        std::vector< std::string > updateOrder = determineEphemerisUpdateorder(
                    propagatedBodies, centralBodies, centralBodies );

        for( int i = updateOrder.size( ) - 1; i >= 0 ; i-- )
        {
            int currentBodyIndex = std::distance( propagatedBodies.begin( ),
                                                  std::find( propagatedBodies.begin( ), propagatedBodies.end( ), updateOrder.at( i ) ) );
            for( unsigned int j = 0; j < propagatedBodies.size( ); j++ )
            {
                if( centralBodies.at( currentBodyIndex ) == propagatedBodies.at( j ) )
                {

                    statePartialAdditionIndices_.push_back( std::make_pair( stateTypeStartIndices_[ propagators::transational_state ] +
                                                            currentBodyIndex * propagators::getSingleIntegrationSize( propagators::transational_state ),
                                                            stateTypeStartIndices_[ propagators::transational_state ] +
                            j * propagators::getSingleIntegrationSize( propagators::transational_state ) ) );
                }
            }
        }
    }

    
    //! Map listing all (named) acceleration models exerted by the (named) bodies.
    /*!
     *  Map listing all (named) acceleration models exerted by the (named) bodies.
     */
    std::map< propagators::IntegratedStateType, orbit_determination::partial_derivatives::StateDerivativePartialsMap > stateDerivativePartialList_;
    
    std::map< IntegratedStateType, int > stateTypeStartIndices_;
    
    
    std::map< IntegratedStateType, std::vector< std::multimap< std::pair< int, int >, boost::function< void( Eigen::Block< Eigen::MatrixXd > ) > > > > statePartialList_;
    
    std::multimap< std::pair< int, int >, boost::function< void( Eigen::Block< Eigen::MatrixXd > ) > >::iterator statePartialIterator_;
    
    std::vector< std::pair< int, int > > statePartialAdditionIndices_;

    
    //! Vector of multimaps of parameter partial function, with index information in variational equations as key
    /*!
     *  Vector of multimaps of parameter partial function. Functions are set after checking dependencies of each accelration partial.
     *  The vector indices coincide with the bodiesToEstimate indices. The key pair of the multimap indicates the start column on the sensitivity matrix
     *  part of the variational equations and the number of columns that the return of the function occupies (number of rows = 3 from size of acceleration)
     */
    std::map< IntegratedStateType, std::vector< std::multimap< std::pair< int, int >,
    boost::function< void( Eigen::Block< Eigen::MatrixXd > ) > > > > parameterPartialList_;
    
    //! Pre-declared iterator over all parameter partial functions
    /*!
     *  Pre-declared iterator over all parameter partial functions. Declared to prevent large number of iterator crations and destructions (performance)
     */
    std::multimap< std::pair< int, int >, boost::function< void( Eigen::Block< Eigen::MatrixXd > ) > >::iterator functionIterator;

    std::map< propagators::IntegratedStateType, orbit_determination::partial_derivatives::StateDerivativePartialsMap >::iterator stateDerivativeTypeIterator_;
    
    

    
    std::map< propagators::IntegratedStateType, std::vector< std::pair< std::string, std::string > > > dynamicalStatesToEstimate_;
    
    
    
    std::vector< std::string > estimatedUnintegratedBodies_;
    
    
    //! Number of parameter values in estimation
    /*!
     *  Number of parameter values in estimation (i.e. number of columns in sensitivity matrix)
     */
    int numberOfParameterValues_;
    
    
    int totalDynamicalStateSize_;

    Eigen::MatrixXd currentMatrixDerivative_;

    Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic > currentLongMatrixDerivative_;

    Eigen::MatrixXd variationalMatrix_;

    Eigen::MatrixXd variationalParameterMatrix_;
};


}

}

#endif // TUDAT_VARIATIONALEQUATIONS_H
