/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_STATEDERIVATIVEPARTIAL_H
#define TUDAT_STATEDERIVATIVEPARTIAL_H

#include <string>
#include <map>
#include <Eigen/Core>

#include <boost/bind/bind.hpp>
using namespace boost::placeholders;


#include "tudat/astro/basic_astro/accelerationModel.h"

#include "tudat/astro/basic_astro/accelerationModelTypes.h"
#include "tudat/astro/propagators/singleStateTypeDerivative.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"


namespace tudat
{

namespace orbit_determination
{

//! Base class for computing the partial derivatives of a state derivative model
/*!
 * Base class for computing the partial derivatives of a state derivative model (i.e. acceleration model for
 * translational dynamics, torque model for rottional dynamics, etc.).
 * Typically, two levels of derived classes are required: one for the type of dynamics, and one for the type of model
 * (e.g. one level for acceleration model, and one level for central gravitational, spherical harmonic model, etc.).
 */
class StateDerivativePartial
{

public:

    //! Constructor.
    /*!
     * Constructor
     * \param integratedStateType Type of dynamics for which partials are to be computed
     * \param integrationReferencePoint Reference point (i.e. propagated body and point) for which the dynamics is
     * propagated. First entry denotes the full body, second entry the reference point on the body that is propagated
     * (empty for translational, rotational dynamics).
     */
    StateDerivativePartial( const propagators::IntegratedStateType integratedStateType,
                            const std::pair< std::string, std::string >& integrationReferencePoint ):
        integratedStateType_( integratedStateType ), integrationReferencePoint_( integrationReferencePoint )
    {
        currentTime_ = TUDAT_NAN;
        accelerationSize_ = propagators::getGeneralizedAccelerationSize( integratedStateType );
    }

    //! Destructor.
    virtual ~StateDerivativePartial( ) { }

    //! Pure virtual function to retrieve the function that returns the partial derivative w.r.t. a propagated state.
    /*!
     * Pure virtual function to retrieve the function that returns the partial derivative w.r.t. a propagated state.
     * \param stateReferencePoint Reference point (id) for propagated state (i.e. body name for translational dynamics).
     * \param integratedStateType Type of propagated state.
     * \return Pair with function, returning partial derivative, and number of columns in partial vector,
     */
    virtual std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int >
    getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ) = 0;

    //! Pure virtual function to check whether a partial w.r.t. some integrated state is non-zero.
    /*!
     * Pure virtual function to check whether a partial w.r.t. some integrated state is non-zero.
     * \param stateReferencePoint Reference point (id) for propagated state (i.e. body name for translational dynamics).
     * \param integratedStateType Type of propagated state.
     * \return True if dependency exists, false otherwise.
     */
    virtual bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ) = 0;

    //! Function to directly compute the partial of the state derivative w.r.t. a double parameter.
    /*!
     * Function to directly compute the partial of the state derivative w.r.t. a double parameter. NOTE: This function is
     * incldued for testing purposes, and is not to be used during propagation (highly inefficient).
     * \param parameter Parameter w.r.t. which partial is to be computed
     * \return Partial of state derivative w.r.t. given parameter.
     */
    Eigen::MatrixXd wrtParameter(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        // Initialize partial
        Eigen::MatrixXd partial = Eigen::MatrixXd( accelerationSize_, 1 );

        // Get partial computation function.
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunction =
                getParameterPartialFunction( parameter );

        // If parameter dependency exists, compute it, otherwise set partial to zero.
        if( partialFunction.second > 0 )
        {
             partialFunction.first( partial );
        }
        else
        {
            partial = Eigen::MatrixXd::Zero(
                        propagators::getGeneralizedAccelerationSize( integratedStateType_ ), 1 );
        }

        return partial;
    }

    //! Function to retrieve the function that computes (by reference) a given double parameter partial.
    /*!
     *  Function to retrieve the function that computes (by reference) a given double parameter partial. NOTE: this function
     *  is implemented in this base class with default no dependency. If any double parameter dependencioes exists, this
     *  function should be overriden in derived class.
     *  \param parameter Parameter w.r.t. which partial is to be computed
     */
    virtual std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    //! Function to directly compute the partial of the state derivative w.r.t. a vector parameter.
    /*!
     * Function to directly compute the partial of the state derivative w.r.t. a vector parameter. NOTE: This function is
     * incldued for testing purposes, and is not to be used during propagation (highly inefficient).
     * \param parameter Parameter w.r.t. which partial is to be computed
     * \return Partial of state derivative w.r.t. given parameter.
     */
    Eigen::MatrixXd wrtParameter(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        // Initialize partial
        Eigen::MatrixXd partial = Eigen::MatrixXd( accelerationSize_, parameter->getParameterSize( ) );

        // Get partial computation function.
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunction =
                getParameterPartialFunction( parameter );

        // If parameter dependency exists, compute it, otherwise set partial to zero.
        if( partialFunction.second > 0 )
        {
             partialFunction.first( partial );
        }
        else
        {
            partial = Eigen::MatrixXd::Zero( propagators::getSingleIntegrationSize( integratedStateType_ ),
                                             parameter->getParameterSize( ) );
        }
        return partial;
    }

    //! Function to retrieve the function that computes (by reference) a given vector parameter partial.
    /*!
     *  Function to retrieve the function that computes (by reference) a given vector parameter partial. NOTE: this function
     *  is implemented in this base class with default no dependency. If any double parameter dependencioes exists, this
     *  function should be overriden in derived class.
     *  \param parameter Parameter w.r.t. which partial is to be computed
     */
    virtual std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    //! Pure virtual for updating the common parts of the partial object to current state and time.
    /*!
     *  Pure virtual for updating the common parts of the partial object to current state and time. Note that this
     *  function is distinct from the updateParameterPartials, which only updates the specific required partial
     *  derivatives w.r.t. required parameters. This function computes variables that are always required when
     *  the derived class partial is used.
     *  \param currentTime Time to which partials are to be updated.
     */
    virtual void update( const double currentTime ) = 0;

    //! Function to get the type of state for which partials are to be computed.
    /*!
     * Function to get the type of state for which partials are to be computed.
     * \return Type of state for which partials are to be computed.
     */
    propagators::IntegratedStateType getIntegratedStateType( )
    {
        return integratedStateType_;
    }

    //! Function to get the identifier for body/reference point for which propagation is performed.
    /*!
     * Function to get the identifier for body/reference point for which propagation is performed.
     * \return Identifier for body/reference point for which propagation is performed.
     */
    std::pair< std::string, std::string > getIntegrationReferencePoint( )
    {
        return integrationReferencePoint_;
    }

    //! Function to reset the  object to the current time
    /*!
     * Function to reset the  object to the current time, recomputing partials to current state.
     *  \param currentTime Time to which partials are to be updated.
     */
    void resetCurrentTime( )
    {
        resetCurrentParameterValues( );
        currentTime_ = TUDAT_NAN;

        // Perform updates of member objects if needed.
        resetCurrentTimeOfMemberObjects( );
    }

    //! Function to retrieve a partial w.r.t. a double parameter
    /*!
     * Function to retrieve a partial w.r.t. a double parameter. An error is thrown if there is no dependency w.r.t.
     * the requested parameter. A warning is printed if the dependency exists, but has not yet been computed for the
     * current time step.
     * \param parameter Partial w.r.t. which a parameter is to be computed
     * \param partialMatrix Partial of state derivative w.r.t. given parameter (return by 'reference', amking use of
     * Eigen::Block architecture).
     */
    void getCurrentParameterPartial(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
            Eigen::Block< Eigen::MatrixXd > partialMatrix )
    {
        // Check if dependecy is computed
        if( currentDoubleParameterPartials_.count( parameter ) == 0 )
        {
            // Check if dependecy exists at all
            if( parameterDoublePartialFunctions_.count( parameter ) == 0 )
            {
                throw std::runtime_error(
                            "Parameter of type " +
                            std::to_string( parameter->getParameterName( ).first ) + ", " +
                            std::string( parameter->getParameterName( ).second.first ) + ", " +
                            std::string( parameter->getParameterName( ).second.second ) + ", " +
                            " not found in list of existing partials" );
            }
            else
            {
                std::cerr << "Warning, double partial should already be calculated" << std::endl;
                parameterDoublePartialFunctions_.at( parameter )( currentDoubleParameterPartials_[ parameter ] );
            }
        }
        partialMatrix += currentDoubleParameterPartials_[ parameter ];
    }

    //! Function to retrieve a partial w.r.t. a double parameter
    /*!
     * Function to retrieve a partial w.r.t. a double parameter. An error is thrown if there is no dependency w.r.t.
     * the requested parameter. A warning is printed if the dependency exists, but has not yet been computed for the
     * current time step.
     * \param parameter Partial w.r.t. which a parameter is to be computed
     * \param parameterPartial Partial of state derivative w.r.t. given parameter (return by reference)
     */
    void getCurrentDoubleParameterPartial(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
            Eigen::MatrixXd& parameterPartial )
    {
        parameterPartial.block( 0, 0, accelerationSize_, 1 ) .setZero( );
        getCurrentParameterPartial( parameter, parameterPartial.block( 0, 0, accelerationSize_, 1 ) );
    }

    //! Function to retrieve a partial w.r.t. a vector parameter
    /*!
     * Function to retrieve a partial w.r.t. a vector parameter. An error is thrown if there is no dependency w.r.t.
     * the requested parameter. A warning is printed if the dependency exists, but has not yet been computed for the
     * current time step.
     * \param parameter Partial w.r.t. which a parameter is to be computed
     * \param partialMatrix Partial of state derivative w.r.t. given parameter (return by 'reference', amking use of
     * Eigen::Block architecture).
     */
    void getCurrentParameterPartial(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
            Eigen::Block< Eigen::MatrixXd > partialMatrix )
    {
        if( currentVectorParameterPartials_.count( parameter ) == 0 )
        {
            if( parameterVectorPartialFunctions_.count( parameter ) == 0 )
            {
                std::string errorMessage = "Parameter of type " +
                        std::to_string( parameter->getParameterName( ).first ) + ", " +
                        std::string( parameter->getParameterName( ).second.first ) + ", " +
                        std::string( parameter->getParameterName( ).second.second ) + ", " +
                        " not found in list of existing partials";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                std::cerr << "Warning, vector partial should already be calculated" << std::endl;
                parameterVectorPartialFunctions_.at( parameter )( currentVectorParameterPartials_[ parameter ] );
            }
        }
        partialMatrix += currentVectorParameterPartials_[ parameter ];
    }

    //! Function to retrieve a partial w.r.t. a vector parameter
    /*!
     * Function to retrieve a partial w.r.t. a vector parameter. An error is thrown if there is no dependency w.r.t.
     * the requested parameter. A warning is printed if the dependency exists, but has not yet been computed for the
     * current time step.
     * \param parameter Partial w.r.t. which a parameter is to be computed
     * \param parameterPartial Partial of state derivative w.r.t. given parameter (return by reference)
     */
    void getCurrentVectorParameterPartial(
            const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
            Eigen::MatrixXd& parameterPartial )
    {
        parameterPartial.block( 0, 0, accelerationSize_, parameter->getParameterSize( ) ).setZero( );
        getCurrentParameterPartial( parameter, parameterPartial.block(
                                        0, 0, accelerationSize_, parameter->getParameterSize( ) ) );
    }

    //! Function to set a dependency of this partial object w.r.t. a given double parameter.
    /*!
     * Function to set a dependency of this partial object w.r.t. a given double parameter. If a dependency exists, the given
     * partial is recomputed on every call of updateParameterPartials.
     * \param parameter Partial w.r.t. which dependency is to be checked and set.
     * \return Size (number of columns) of parameter partial. Zero if no dependency, 1 otherwise.
     */
    virtual int setParameterPartialUpdateFunction(
                std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        // Get partial function.
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > parameterPartialFunction =
                getParameterPartialFunction( parameter );

        // If partial function found, add function to list of computations to be performed when calling
        // updateParameterPartials.
        if( parameterPartialFunction.second > 0 && parameterDoublePartialFunctions_.count( parameter ) == 0 )
        {
            if( parameterDoublePartialFunctions_.count( parameter ) == 0 )
            {
                parameterDoublePartialFunctions_[ parameter ] = parameterPartialFunction.first;
                isCurrentDoubleParameterPartialSet_[ parameter ] = 0;
                currentDoubleParameterPartials_[ parameter ] = Eigen::MatrixXd( accelerationSize_, 1 );
            }
        }

        return parameterPartialFunction.second;
    }

    //! Function to set a dependency of this partial object w.r.t. a given vector parameter.
    /*!
     * Function to set a dependency of this partial object w.r.t. a given vector parameter. If a dependency exists, the given
     * partial is recomputed on every call of updateParameterPartials.
     * \param parameter Partial w.r.t. which dependency is to be checked and set.
     * \return Size (number of columns) of parameter partial. Zero if no dependency, size of parameter otherwise.
     */
    virtual int setParameterPartialUpdateFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        // Get partial function.
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > parameterPartialFunction =
                getParameterPartialFunction( parameter );

        // If partial function found, add function to list of computations to be performed when calling
        // updateParameterPartials.
        if( parameterPartialFunction.second > 0 && currentVectorParameterPartials_.count( parameter ) == 0 )
        {
            if( parameterVectorPartialFunctions_.count( parameter ) == 0 )
            {
                parameterVectorPartialFunctions_[ parameter ] = parameterPartialFunction.first;
                isCurrentVectorParameterPartialSet_[ parameter ] = 0;
                currentVectorParameterPartials_[ parameter ] =
                        Eigen::MatrixXd( accelerationSize_, parameter->getParameterSize( ) );
            }
        }

        return parameterPartialFunction.second;
    }

    //! Function to update the values partial derivatives to current state and time.
    /*!
     * Function to update the values partial derivatives to current state and time.
     */
    void updateParameterPartials( )
    {
        // Update member object partials.
        updateParameterPartialsOfMemberObjects( );

        // Update double parameter partials
        for( parameterDoublePartialFunctionIterator_ = parameterDoublePartialFunctions_.begin( );
             parameterDoublePartialFunctionIterator_ != parameterDoublePartialFunctions_.end( );
             parameterDoublePartialFunctionIterator_++ )
        {
            if( isCurrentDoubleParameterPartialSet_.at( parameterDoublePartialFunctionIterator_->first ) == 0 )
            {
                parameterDoublePartialFunctionIterator_->second( currentDoubleParameterPartials_[ parameterDoublePartialFunctionIterator_->first ]  );
                isCurrentDoubleParameterPartialSet_[ parameterDoublePartialFunctionIterator_->first ] = 1;
            }
        }

        // Update vector parameter partials
        for( parameterVectorPartialFunctionIterator_ = parameterVectorPartialFunctions_.begin( );
             parameterVectorPartialFunctionIterator_ != parameterVectorPartialFunctions_.end( );
             parameterVectorPartialFunctionIterator_++ )
        {
            if( isCurrentVectorParameterPartialSet_.at( parameterVectorPartialFunctionIterator_->first ) == 0 )
            {
                parameterVectorPartialFunctionIterator_->second( currentVectorParameterPartials_[ parameterVectorPartialFunctionIterator_->first ]  );
                isCurrentVectorParameterPartialSet_[  parameterVectorPartialFunctionIterator_->first ] = 1;
            }
        }
    }

protected:

    //! Function to compute parameter partials of member objects.
    /*!
     *  Function to  compute parameter partials of member objects. By default (implemented here) no computations are
     *  performed. For certain derived classed (i.e. ThirdBodyGravityPartial), there are member StateDerivativePartial
     *  objects that need to be updated when calling updateParameterPartials, for which this function should be redefined.
     */
    virtual void updateParameterPartialsOfMemberObjects( )
    {

    }

    //! Function to reset the member object to the current time
    /*!
     *  Function to reset the member object to the current time. By default (implemented here) no computations are performed.
     *  For certain derived classed (i.e. ThirdBodyGravityPartial), there are member StateDerivativePartial objects that
     *  need to be updated when calling resetCurrentTime, for which this function should be redefined.
     */
    virtual void resetCurrentTimeOfMemberObjects( )
    {

    }

    //! Function to define all current parameter partials as 'not computed'
    void resetCurrentParameterValues( )
    {
        for( isCurrentDoubleParameterPartialSetIterator_ = isCurrentDoubleParameterPartialSet_.begin( );
             isCurrentDoubleParameterPartialSetIterator_ !=  isCurrentDoubleParameterPartialSet_.end( );
             isCurrentDoubleParameterPartialSetIterator_++ )
        {
            isCurrentDoubleParameterPartialSet_[ isCurrentDoubleParameterPartialSetIterator_->first ] = 0;
        }

        for( isCurrentVectorParameterPartialSetIterator_ = isCurrentVectorParameterPartialSet_.begin( );
             isCurrentVectorParameterPartialSetIterator_ !=  isCurrentVectorParameterPartialSet_.end( );
             isCurrentVectorParameterPartialSetIterator_++ )
        {
            isCurrentVectorParameterPartialSet_[ isCurrentVectorParameterPartialSetIterator_->first ] = 0;
        }
    }

    //! Type of state for which partials are to be computed.
    propagators::IntegratedStateType integratedStateType_;

    //! Size of the single order state derivative model (i.e. 3 for translational dynamics).
    int accelerationSize_;

    //! Identifier for body/reference point for which propagation is performed.
    /*!
     *  Identifier for body/reference point for which propagation is performed. First entry represents the body being
     *  propagated. Second entry is empty for propagation of entire bodies (translational, rotational), but must be
     *  set for propagation of local dynamics (i.e. observer proper time).
     */
    std::pair< std::string, std::string > integrationReferencePoint_;


    //! List of booleans defining whether the partial w.r.t. the current double parameter (key) has been computed.
    std::map< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > >, bool >
    isCurrentDoubleParameterPartialSet_;

    //! Iterator for list defining whether the partial w.r.t. the current double parameter (key) has been computed.
    std::map< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > >, bool >::iterator
    isCurrentDoubleParameterPartialSetIterator_;



    //! List of current values of partials w.r.t. double parameter values (emptied at beginning of every time step).
    std::map< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > >, Eigen::MatrixXd >
    currentDoubleParameterPartials_;

    //! List of functions to compute (return by reference) values of partials w.r.t. doule parameter partials
    std::map< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > >,
    std::function< void( Eigen::MatrixXd& ) > > parameterDoublePartialFunctions_;

    //! Iterator over list of functions to compute values of partials w.r.t. double parameter partials
    std::map< std::shared_ptr< estimatable_parameters::EstimatableParameter< double > >,
    std::function< void( Eigen::MatrixXd& ) > >::iterator parameterDoublePartialFunctionIterator_;



    //! List of booleans defining whether the partial w.r.t. the current vector parameter (key) has been computed.
    std::map< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >, bool >
    isCurrentVectorParameterPartialSet_;

    //! Iterator for list defining whether the partial w.r.t. the current vector parameter (key) has been computed.
    std::map< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >, bool >::iterator
    isCurrentVectorParameterPartialSetIterator_;



    //! List of current values of partials w.r.t. double parameter values (emptied at beginning of every time step).
    std::map< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >,
    Eigen::MatrixXd > currentVectorParameterPartials_;

    //! List of functions to compute (return by reference) values of partials w.r.t. vector parameter partials
    std::map< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >,
    std::function< void( Eigen::MatrixXd&  ) > > parameterVectorPartialFunctions_;

    //! Iterator over list of functions to compute values of partials w.r.t. vector parameter partials
    std::map< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >,
    std::function< void( Eigen::MatrixXd& ) > >::iterator parameterVectorPartialFunctionIterator_;


    double currentTime_;

};

//! Typedef for double vector of StateDerivativePartial objects.
/*!
 *  Typedef for double vector of StateDerivativePartial objects. First (outer) vector is typically the
 *  bodies undergoing 'acceleration (and being estimated), the second (inner) vector is the list of partials
 *  being exerted on a single body.
 */
typedef std::vector< std::vector< std::shared_ptr< orbit_determination::StateDerivativePartial > > >
StateDerivativePartialsMap;

//! Function to evaluate the negative value of a parameter partial.
/*!
 *  Function to evaluate the negative value of a parameter partial.
 *  \param parameterPartialFunction Function to compute the regular paramater partial (by reference).
 *  \param partial Negative value of partial computed by parameterPartialFunction (returned by reference).
 */
void evaluateNegativeParameterPartialFunction(
        const std::function< void( Eigen::MatrixXd& ) > parameterPartialFunction,
        Eigen::MatrixXd& partial );

//! Function to evaluate the subtraction of two parameter partials.
/*!
 *  Function to evaluate the subtraction of two parameter partials.
 *  \param firstParameterPartialFunction Function to compute the first paramater partial (by reference), from which
 *  the value computed by parameterPartialFunctionToSubtract is subtracted.
 *  \param parameterPartialFunctionToSubtract Function to compute the paramater partial (by reference) that is to be
 *  subtracted.
 *  \param partial Value of partial returned by parameterPartialFunctionToSubtract, subtracted from value returned by
 *  firstParameterPartialFunction (returned by reference).
 */
void evaluateSubtractedParameterPartialFunction(
        const std::function< void( Eigen::MatrixXd& ) > firstParameterPartialFunction,
        const std::function< void( Eigen::MatrixXd& ) > parameterPartialFunctionToSubtract,
        Eigen::MatrixXd& partial );

//! Function to evaluate the addition of two parameter partials.
/*!
 *  Function to evaluate the addition of two parameter partials.
 *  \param firstParameterPartialFunction Function to compute the first paramater partial (by reference) that is to be
 *  added.
 *  \param parameterPartialFunctionToAdd Function to compute the second paramater partial (by reference) that is to be
 *  added.
 *  \param partial Value of partial returned by firstParameterPartialFunction, added to value returned by
 *  parameterPartialFunctionToAdd (returned by reference).
 */
void evaluateAddedParameterPartialFunction(
        const std::function< void( Eigen::MatrixXd& ) > firstParameterPartialFunction,
        const std::function< void( Eigen::MatrixXd& ) > parameterPartialFunctionToAdd,
        Eigen::MatrixXd& partial );

//! Create a parameter partial function obtained from the subtraction of two such function results.
/*!
 * Create a parameter partial function, as returned by the StateDerivativePartial::getParameterPartialFunction function
 * The partial created here is obtained from the subtraction of two such function results. The two input variables
 * may be both empty, both define a function, or only one of them may define a function.
 * \param partialFunctionOfAccelerationToAdd Function and associated parameter size (first and second of pair) that
 * are to be added to the total partial.
 * \param partialFunctionOfAccelerationToSubtract Function and associated parameter size (first and second of pair) that
 * are to be subtracted from the total partial.
 * \return Function and parameter size obtained from 'subtracting' partialFunctionOfAccelerationToSubtract from
 * partialFunctionOfAccelerationToAdd
 */
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > createMergedParameterPartialFunction(
        const std::pair< std::function< void( Eigen::MatrixXd& ) >, int >& partialFunctionOfAccelerationToAdd,
        const std::pair< std::function< void( Eigen::MatrixXd& ) >, int >& partialFunctionOfAccelerationToSubtract );


//! Function to create a parameter partial evaluation function, obtained by adding or subtracting a given partial
//! w.r.t. a double parameter from 2 state derivative partial models.
/*!
 * Function to create a parameter partial evaluation function, obtained by adding or subtracting a given partial
 * w.r.t. a double parameter from 2 state derivative partial models. Note that this function creates a merged
 * function from two getCurrentDoubleParameterPartial functions of the two input partial objects. The automatic
 * computation when updating the partial (done by call to setParameterPartialUpdateFunction) is not handled by
 * this function.
 * \param firstPartial First object for computing partial derivatives.
 * \param secondPartial Second object for computing partial derivatives.
 * \param parameterObject Parameter w.r.t. which a partial is to be taken.
 * \param firstPartialSize Size of partial from firstPartial object (only 0 and 1 are valid).
 * \param secondPartialSize Size of partial from secondPartial object (only 0 and 1 are valid).
 * \param subtractPartials Boolean denoting whether the second parameter is to be subtracted or added to the
 * total partial.
 * \return Function computing and returning (by reference) the combined partial according to the required settings.
 */
std::function< void( Eigen::MatrixXd& ) > getCombinedCurrentDoubleParameterFunction(
        const std::shared_ptr< StateDerivativePartial > firstPartial,
        const std::shared_ptr< StateDerivativePartial > secondPartial,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameterObject,
        const int firstPartialSize, const int secondPartialSize,
        const bool subtractPartials = 0 );

//! Function to create a parameter partial evaluation function, obtained by adding or subtracting a given partial
//! w.r.t. a vector parameter from 2 state derivative partial models.
/*!
 * Function to create a parameter partial evaluation function, obtained by adding or subtracting a given partial
 * w.r.t. a vector parameter from 2 state derivative partial models. Note that this function creates a merged
 * function from two getCurrentVectorParameterPartial functions of the two input partial objects. The automatic
 * computation when updating the partial (done by call to setParameterPartialUpdateFunction) is not handled by
 * this function.
 * \param firstPartial First object for computing partial derivatives.
 * \param secondPartial Second object for computing partial derivatives.
 * \param parameterObject Parameter w.r.t. which a partial is to be taken.
 * \param firstPartialSize Size of partial from firstPartial object
 * \param secondPartialSize Size of partial from secondPartial object
 * \param subtractPartials Boolean denoting whether the second parameter is to be subtracted or added to the
 * total partial.
 * \return Function computing and returning (by reference) the combined partial according to the required settings.
 */
std::function< void( Eigen::MatrixXd& ) > getCombinedCurrentVectorParameterFunction(
        const std::shared_ptr< StateDerivativePartial > firstPartial,
        const std::shared_ptr< StateDerivativePartial > secondPartial,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameterObject,
        const int firstPartialSize, const int secondPartialSize,
        const bool subtractPartials = 0 );

} // namespace orbit_determination

} // namespace tudat

#endif // TUDAT_STATEDERIVATIVEPARTIAL_H
