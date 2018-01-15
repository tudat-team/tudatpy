/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDAT_VARIATIONALEQUATIONS_H
#define TUDAT_VARIATIONALEQUATIONS_H

#include <map>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"

#include "Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"

namespace tudat
{

namespace propagators
{

//! Class from which the variational equations can be evaluated.
/*!
 *  Class from which the variational equations can be evaluated. The time derivative of the state transition  and
 *  sensitivity matrices are computed from a set of state derivative partials objects, at the current time and state.
 *  This class performs all required bookkeeping to update, evaluate and combine these state derivative partials into
 *  the variational equations. The VariationalEquationsSolver object is used to manage and execute the full numerical
 *  integration of these variational equations and equations of motion.
 */
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
     * \param parametersToEstimate Object containing all parameters that are to be estimated and their current settings and
     * values.
     * \param stateTypeStartIndices Start index (value) in vector of propagated state for each type of state (key)
     */
    template< typename ParameterType >
    VariationalEquations(
            const std::map< IntegratedStateType, orbit_determination::StateDerivativePartialsMap >
            stateDerivativePartialList,
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
            const std::map< IntegratedStateType, int >& stateTypeStartIndices ):
        stateDerivativePartialList_( stateDerivativePartialList ), stateTypeStartIndices_( stateTypeStartIndices )
    {        
        dynamicalStatesToEstimate_ =
                estimatable_parameters::getListOfInitialDynamicalStateParametersEstimate< ParameterType >(
                    parametersToEstimate );
        
        // Get size of dynamical state to estimate
        numberOfParameterValues_ = estimatable_parameters::getSingleArcParameterSetSize( parametersToEstimate );
        totalDynamicalStateSize_ = 0;        
        for( std::map< IntegratedStateType, orbit_determination::StateDerivativePartialsMap >::iterator
             partialTypeIterator = stateDerivativePartialList_.begin( );
             partialTypeIterator != stateDerivativePartialList_.end( ); partialTypeIterator++ )
        {
            
            if( dynamicalStatesToEstimate_.count( partialTypeIterator->first ) == 0 )
            {
                std::string errorMessage = "Error when making variational equations object, found no state to estimate of type " +
                        std::to_string( partialTypeIterator->first );
                throw std::runtime_error( errorMessage );
            }
            else if( dynamicalStatesToEstimate_.at( partialTypeIterator->first ).size( ) !=
                     partialTypeIterator->second.size( ) )
            {
                throw std::runtime_error( "Error when making variational equations object, input partial list size is inconsistent" );
            }
            
            totalDynamicalStateSize_ +=
                    getSingleIntegrationSize( partialTypeIterator->first ) * partialTypeIterator->second.size( );
        }
        
        // Initialize matrices.
        variationalMatrix_ = Eigen::MatrixXd::Zero( totalDynamicalStateSize_, totalDynamicalStateSize_ );
        variationalParameterMatrix_ =
                Eigen::MatrixXd::Zero( totalDynamicalStateSize_, numberOfParameterValues_ - totalDynamicalStateSize_ );

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
                ( variationalMatrix_.template cast< StateScalarType >( ) * stateTransitionAndSensitivityMatrices );
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
            int entriesToSkipPerEntry = currentStateSize -
                    currentStateSize / getSingleIntegrationDifferentialEquationOrder( typeIterator->first );

            // Iterate over all bodies being estimated.
            for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
            {
                // Iterate over all parameter partial functions determined by setParameterPartialFunctionList( )
                for( functionIterator = typeIterator->second[ i ].begin( );
                     functionIterator != typeIterator->second[ i ].end( );
                     functionIterator++ )
                {
                    functionIterator->second(
                                variationalParameterMatrix_.block(
                                    startIndex + entriesToSkipPerEntry + currentStateSize * i,
                                    functionIterator->first.first - totalDynamicalStateSize_,
                                    currentStateSize - entriesToSkipPerEntry,
                                    functionIterator->first.second ) );
                }
            }

        }

        currentMatrixDerivative.block( 0, totalDynamicalStateSize_, totalDynamicalStateSize_,
                                       numberOfParameterValues_ - totalDynamicalStateSize_ ) +=
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
            const double time, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >&
            stateTransitionAndSensitivityMatrices,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > currentMatrixDerivative )
    {
        // Compute and add state partials.
        getBodyInitialStatePartialMatrix< StateScalarType >( stateTransitionAndSensitivityMatrices,currentMatrixDerivative );

        if( numberOfParameterValues_ > totalDynamicalStateSize_ )
        {
            // Add partials of parameters.
            getParameterPartialMatrix< StateScalarType >( currentMatrixDerivative );
        }
    }

    //! Function to clear reference/cached values of state derivative partials.
    /*!
     * Function to clear reference/cached values of state derivative partials, to ensure that they are all recalculated.
     */
    void clearPartials( );
    
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
    
    //! Function (called by constructor) to set up the statePartialList_ member from the state derivative partials
    /*!
     * Function (called by constructor) to set up the functions to evaluate the partial derivatives of the state derivatives
     * w.r.t. a current state (stored in the statePartialList_ member) from the state derivative partials.
     */
    void setStatePartialFunctionList( );
        
    //! Function to add parameter partial functions for single state derivative model, and set of parameter objects.
    /*!
     *  Function to add parameter partial functions for single state derivative model, and set of parameter objects.
     *  Partial derivative functions that are not-empty are added to the functionListOfBody input (returned by reference).
     *  A list of parameters of a single type (double or vector) are handled a single function call.
     *  \param parameterList Map of parameters for which partial functions are to checked and created. Map keys are
     *  start entry of parameter in total parameter vector.
     *  \param partialObject State derivative partial object from which partial functions are to be retrieved.
     *  \param functionListOfBody Multimap of partial derivative functions to which entries are to be added by this function.
     *  Map key is start index and size of given parameter in sensitivity matrix. Map value is partial function.
     *  \param totalParameterVectorIndicesToSubtract Number of entries by which to shift start index in sensitivity
     *  matrix from entry in parameter vector (used for multi-arc estimation).
     */
    template< typename CurrentParameterType >
    void addParameterPartialToList(
            const std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< CurrentParameterType > > >&
            parameterList,
            const boost::shared_ptr< orbit_determination::StateDerivativePartial > partialObject,
            std::multimap< std::pair< int, int >, boost::function< void( Eigen::Block< Eigen::MatrixXd > ) > >&
            functionListOfBody,
            const int totalParameterVectorIndicesToSubtract = 0 )
    {
        using namespace acceleration_partials;

        // Iterate over all parameters.
        for( typename std::map< int,
             boost::shared_ptr< estimatable_parameters::EstimatableParameter< CurrentParameterType > > >::const_iterator
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
                                  static_cast< void ( orbit_determination::StateDerivativePartial::* )
                                  ( const boost::shared_ptr<
                                    estimatable_parameters::EstimatableParameter< CurrentParameterType > >,
                                    Eigen::Block< Eigen::MatrixXd > )>
                                  ( &orbit_determination::StateDerivativePartial::getCurrentParameterPartial ),
                                  partialObject, parameterIterator->second, _1  ) ) );
            }
        }
    }
    
    //! This function creates the list of partial derivatives of the state w.r.t. parameter values.
    /*!
     *  This function creates the list of partial derivatives of the state w.r.t. parameter values.
     *  The function is called once by the constructor and the resulting functions are set as member variables.
     *  This prevents having to check whether an acceleration model depends on every parameter during every time step.
     *  \param parametersToEstimate Object containing all parameters that are to be estimated and their current settings and
     *  values.
     */
    template< typename ParameterType >
    void setParameterPartialFunctionList(
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > >
            parametersToEstimate )
    {
        // Get double parameters.
        std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >
                doubleParametersToEstimate =
                parametersToEstimate->getDoubleParameters( );

        // Get vector parameters.
        std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >
                vectorParametersToEstimate =
                parametersToEstimate->getVectorParameters( );
        
        int totalParameterVectorIndicesToSubtract = parametersToEstimate->getInitialDynamicalStateParameterSize( ) -
                estimatable_parameters::getSingleArcInitialDynamicalStateParameterSetSize( parametersToEstimate );

        for( std::map< propagators::IntegratedStateType,
             orbit_determination::StateDerivativePartialsMap >::iterator
             stateDerivativeTypeIterator = stateDerivativePartialList_.begin( );
             stateDerivativeTypeIterator != stateDerivativePartialList_.end( );
             stateDerivativeTypeIterator++ )
        {
            
            // Initialize vector of lists to correct size.
            parameterPartialList_[ stateDerivativeTypeIterator->first ].resize(
                        stateDerivativeTypeIterator->second.size( ) );
            
            // Iterate over all bodies of which initial position is being estimated.
            for( unsigned int i = 0; i < stateDerivativeTypeIterator->second.size( ); i++ )
            {
                // Initialize list of parameter partial functions for single body.
                std::multimap< std::pair< int, int >, boost::function< void( Eigen::Block< Eigen::MatrixXd > ) > >
                        functionListOfBody;
                
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

    //! Function called by constructor to handle estimation of hierarchical translational dynamics
    /*!
     *  Function called by constructor to handle estimation of hierarchical translational dynamics, i.e. where
     *  the state of body A is estimated w.r.t. body B, and body B is itself estimated w.r.t. to some third body (or inertial
     *  point) C.
     *  \param parametersToEstimate Total list of parameters to estimate.
     */
    template< typename ParameterType >
    void setTranslationalStatePartialFrameScalingFunctions(
            const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > >
            parametersToEstimate )
    {
        std::vector< boost::shared_ptr< estimatable_parameters::EstimatableParameter<
                Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
                parametersToEstimate->getEstimatedInitialStateParameters( );

        std::vector< std::string > propagatedBodies;
        std::vector< std::string > centralBodies;

        // Retrieve propagated bodies and central bodies of estimation.
        for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
        {
            if( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::initial_body_state )
            {
                propagatedBodies.push_back(
                            initialDynamicalParameters.at( i )->getParameterName( ).second.first );
                centralBodies.push_back( boost::dynamic_pointer_cast
                                         < estimatable_parameters::InitialTranslationalStateParameter< ParameterType > >(
                                             initialDynamicalParameters.at( i ) )->getCentralBody( ) );
            }
            else if( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::arc_wise_initial_body_state )
            {
                propagatedBodies.push_back(
                            initialDynamicalParameters.at( i )->getParameterName( ).second.first );
                centralBodies.push_back( boost::dynamic_pointer_cast< estimatable_parameters::ArcWiseInitialTranslationalStateParameter< ParameterType > >(
                                             initialDynamicalParameters.at( i ) )->getCentralBody( ) );
            }
            else
            {
                throw std::runtime_error( "Error when settinf up variational equations, did not recognize initial state type" );
            }
        }

        // Get order in which ephemerides were to be updated.
        std::vector< std::string > updateOrder = determineEphemerisUpdateorder(
                    propagatedBodies, centralBodies, centralBodies );

        // Iterate over central bodies and propagated bodies and check for dependencies
        for( int i = updateOrder.size( ) - 1; i >= 0 ; i-- )
        {
            int currentBodyIndex = std::distance(
                        propagatedBodies.begin( ),
                        std::find( propagatedBodies.begin( ), propagatedBodies.end( ), updateOrder.at( i ) ) );
            for( unsigned int j = 0; j < propagatedBodies.size( ); j++ )
            {
                if( centralBodies.at( currentBodyIndex ) == propagatedBodies.at( j ) )
                {

                    statePartialAdditionIndices_.push_back(
                                std::make_pair( stateTypeStartIndices_[ propagators::transational_state ] +
                                currentBodyIndex * propagators::getSingleIntegrationSize( propagators::transational_state ),
                                stateTypeStartIndices_[ propagators::transational_state ] +
                            j * propagators::getSingleIntegrationSize( propagators::transational_state ) ) );
                }
            }
        }
    }

    
    //! Map with list of StateDerivativePartialsMaps, with state type as key.
    /*!
     *  List partials of state derivative models from which the variational equations
     *  are set up. The key is the type of dynamics for which partials are taken, the values are StateDerivativePartialsMap
     *  (see StateDerivativePartialsMap definition for details)
     */
    std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >
    stateDerivativePartialList_;
    
    //! Map of start entry in sensitivity matrix of each type of estimated dynamics.
    std::map< IntegratedStateType, int > stateTypeStartIndices_;
    
    //! List of all functions returning current partial derivative w.r.t. a current dynamical state
    /*!
     *  List of all functions returning current partial derivative w.r.t. a current dynamical state (map key). The
     *  vector entries correspond to the entries in the outer vector of StateDerivativePartialsMaps in
     *  stateDerivativePartialList_. The multimaps inside the vector provide the functions (as values) adding the
     *  partials to a given matrix block and the start column and number of columns in matrix partial (as keys).
     */
    std::map< IntegratedStateType,
    std::vector< std::multimap< std::pair< int, int >, boost::function< void( Eigen::Block< Eigen::MatrixXd > ) > > > >
    statePartialList_;
    
    //! Pre-defined iterator for efficiency.
    std::multimap< std::pair< int, int >, boost::function< void( Eigen::Block< Eigen::MatrixXd > ) > >::iterator
    statePartialIterator_;
    
    //! Vector of pair providing indices of column blocks of variational equations to add to other column blocks
    /*!
     * Vector of pair providing indices of column blocks of variational equations to add to other column blocks,
     * which is needed by the hierarchical estimation fo dynamics. The second pair entry is the start of the column
     * block (of size 3) to where it should be copied. The first entry denotes from where it should be copied.
     * \sa setTranslationalStatePartialFrameScalingFunctions
     */
    std::vector< std::pair< int, int > > statePartialAdditionIndices_;

    
    //! List of all functions returning current partial derivative w.r.t. a parameter
    /*!
     *  List of all functions returning current partial derivative w.r.t. a parameter.
     *  Map key denotes associated dynamics type w.r.t which partial is taken. The
     *  vector entries correspond to the entries in the outer vector of StateDerivativePartialsMaps in
     *  stateDerivativePartialList_. The multimaps inside the vector provide the functions (as values) adding the
     *  partials to a given matrix block and the start column and number of columns in matrix partial (as keys).
     */
    std::map< IntegratedStateType, std::vector< std::multimap< std::pair< int, int >,
    boost::function< void( Eigen::Block< Eigen::MatrixXd > ) > > > > parameterPartialList_;
    
    //! Pre-declared iterator over all parameter partial functions.
    std::multimap< std::pair< int, int >, boost::function< void( Eigen::Block< Eigen::MatrixXd > ) > >
    ::iterator functionIterator;

    //! Pre-declared iterator over all state types
    std::map< propagators::IntegratedStateType, orbit_determination::StateDerivativePartialsMap >
    ::iterator stateDerivativeTypeIterator_;
    
    //! List of identifiers for points/bodies for which initial dynamical state is to be estimated.
    std::map< propagators::IntegratedStateType, std::vector< std::pair< std::string, std::string > > >
    dynamicalStatesToEstimate_;

    
    //! Number of parameter values in estimation (i.e. number of columns in sensitivity matrix)
    int numberOfParameterValues_;
    
    //! Total size of (single-arc) state vector of dynamics that is to be estimated.
    int totalDynamicalStateSize_;

    //! Total matrix of partial derivatives of state derivatives w.r.t. current states.
    Eigen::MatrixXd variationalMatrix_;

    //! Total matrix of partial derivatives of state derivatives w.r.t. parameter vectors.
    Eigen::MatrixXd variationalParameterMatrix_;
};


} // namespace propagators

} // namespace tudat

#endif // TUDAT_VARIATIONALEQUATIONS_H
