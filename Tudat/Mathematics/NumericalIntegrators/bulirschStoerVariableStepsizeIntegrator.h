/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120316    K. Kumar          File created.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of
 *          Scientific Computing. Cambridge University Press, February 2002.
 *
 */

#ifndef TUDAT_BULIRSCH_STOER_VARIABLE_STEP_SIZE_INTEGRATOR_H
#define TUDAT_BULIRSCH_STOER_VARIABLE_STEP_SIZE_INTEGRATOR_H

#include <cmath>

#include <boost/assign/std/vector.hpp>

#include <Eigen/Core>

#include <Tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>


namespace tudat
{
namespace numerical_integrators
{

enum ExtrapolationMethodStepSequences
{
    bulirsch_stoer_sequence,
    deufelhard_sequence
};


std::vector< unsigned int > getBulirschStoerStepSequence(
        const ExtrapolationMethodStepSequences& extrapolationMethodStepSequenceType = bulirsch_stoer_sequence,
        const unsigned int lengthOfSequence = 12 );

//! Class that implements the Bulirsch-Stoer variable stepsize integrator.
/*!
 * Class that implements the Bulirsch-Stoer variable step size integrator.
 * \tparam StateType The type of the state. This type should be an Eigen::Matrix derived type.
 * \tparam StateDerivativeType The type of the state derivative. This type should be an
 *          Eigen::Matrix derived type.
 * \tparam IndependentVariableType The type of the independent variable. This type should be
 *          either a float or double.
 * \sa NumericalIntegrator.
 */
template < typename IndependentVariableType = double, typename StateType = Eigen::VectorXd,
           typename StateDerivativeType = Eigen::VectorXd >
class BulirschStoerVariableStepSizeIntegrator :
        public NumericalIntegrator< IndependentVariableType, StateType, StateDerivativeType >
{
public:

    //! Typedef of the base class.
    /*!
     * Typedef of the base class with all template parameters filled in.
     */
    typedef NumericalIntegrator< IndependentVariableType, StateType, StateDerivativeType > Base;

    //! Typedef to the state derivative function.
    /*!
     * Typedef to the state derivative function inherited from the base class.
     * \sa NumericalIntegrator::StateDerivativeFunction.
     */
    typedef typename Base::StateDerivativeFunction StateDerivativeFunction;

    //! Default constructor.
    /*!
     * Default constructor, taking sequence, a state derivative function, initial conditions,
     * minimum step size and relative error tolerance per item in the state vector as argument.
     * \param sequence Rational function sequence used by algorithm.
     * \param stateDerivativeFunction State derivative function.
     * \param intervalStart The start of the integration interval.
     * \param initialState The initial state.
     * \param minimumStepSize The minimum step size to take. If this constraint is violated, a
     *          flag will be set that can be retrieved with isMinimumStepSizeViolated( ).
     * \param relativeErrorTolerance The relative error tolerance, for each individual state
     *          vector element.
     * \sa NumericalIntegrator::NumericalIntegrator.
     */
    BulirschStoerVariableStepSizeIntegrator(
            const std::vector< unsigned int >& sequence,
            const StateDerivativeFunction& stateDerivativeFunction,
            const IndependentVariableType intervalStart,  const StateType& initialState,
            const IndependentVariableType minimumStepSize,
            const StateType& relativeErrorTolerance ) :
        Base( stateDerivativeFunction ), currentIndependentVariable_( intervalStart ),
        currentState_( initialState ), lastIndependentVariable_( intervalStart ),
        sequence_( sequence ), minimumStepSize_( minimumStepSize ),
        relativeErrorTolerance_( relativeErrorTolerance ), isMinimumStepSizeViolated_( false )
    {
        maximumNumberOfAttempts_ = sequence_.size( );
        subSteps_.resize( maximumNumberOfAttempts_ );

        integratedStates_.resize( maximumNumberOfAttempts_ );
        for( unsigned int i = 0; i < maximumNumberOfAttempts_; i++ )
        {
            integratedStates_[ i ].resize( maximumNumberOfAttempts_ );
        }
    }

    //! Default constructor.
    /*!
     * Default constructor, taking coefficients, a state derivative function, initial conditions,
     * minimum step size and relative error tolerance for all items in the state vector as argument.
     * \param coefficients Coefficients to use with this integrator.
     * \param stateDerivativeFunction State derivative function.
     * \param intervalStart The start of the integration interval.
     * \param initialState The initial state.
     * \param minimumStepSize The minimum step size to take. If this constraint is violated, a
     *          flag will be set that can be retrieved with isMinimumStepSizeViolated( ).
     * \param relativeErrorTolerance The relative error tolerance, equal for all individual state
     *          vector elements.
     * \sa NumericalIntegrator::NumericalIntegrator.
     */
    BulirschStoerVariableStepSizeIntegrator(
            const std::vector< unsigned int >& sequence,
            const StateDerivativeFunction& stateDerivativeFunction,
            const IndependentVariableType intervalStart, const StateType& initialState,
            const IndependentVariableType minimumStepSize = 1.0e-15,
            const typename StateType::Scalar relativeErrorTolerance = 1.0e-12,
            const typename StateType::Scalar absoluteErrorTolerance = 1.0e-12) :
        Base( stateDerivativeFunction ), currentIndependentVariable_( intervalStart ),
        currentState_( initialState ), lastIndependentVariable_( intervalStart ),
        sequence_( sequence ), minimumStepSize_( minimumStepSize ),
        relativeErrorTolerance_( StateType::Constant( initialState.rows( ), initialState.cols( ),
                                                      relativeErrorTolerance ) ),
        absoluteErrorTolerance_( StateType::Constant( initialState.rows( ), initialState.cols( ),
                                                      absoluteErrorTolerance ) ),
        isMinimumStepSizeViolated_( false )
    {
        maximumNumberOfAttempts_ = sequence_.size( );
        subSteps_.resize( maximumNumberOfAttempts_ );

        integratedStates_.resize( maximumNumberOfAttempts_ );
        for( unsigned int i = 0; i < maximumNumberOfAttempts_; i++ )
        {
            integratedStates_[ i ].resize( maximumNumberOfAttempts_ );
        }
    }

    //! Get step size of the next step.
    /*!
     * Returns the step size of the next step.
     * \return Step size to be used for the next step.
     */
    virtual IndependentVariableType getNextStepSize( ) const { return stepSize_; }

    //! Get current state.
    /*!
     * Returns the current state of the integrator.
     * \return Current integrated state.
     */
    virtual StateType getCurrentState( ) const { return currentState_; }

    //! Returns the current independent variable.
    /*!
     * Returns the current value of the independent variable of the integrator.
     * \return Current independent variable.
     */
    virtual IndependentVariableType getCurrentIndependentVariable( ) const
    {
        return currentIndependentVariable_;
    }

    //! Perform a single integration step.
    /*!
     * Perform a single integration step and compute a new step size.
     * \param stepSize The step size to take. If the time step is too large to satisfy the error
     *          constraints, the step is redone until the error constraint is satisfied.
     * \return The state at the end of the interval.
     */
    virtual StateType performIntegrationStep( const IndependentVariableType stepSize );

    //! Rollback internal state to the last state.
    /*!
     * Performs rollback of the internal state to the last state. This function can only be called
     * once after calling integrateTo( ) or performIntegrationStep( ) unless specified otherwise by
     * implementations, and can not be called before any of these functions have been called. Will
     * return true if the rollback was successful, and false otherwise.
     * \return True if the rollback was successful.
     */
    virtual bool rollbackToPreviousState( )
    {
        if ( currentIndependentVariable_ == lastIndependentVariable_ )
        {
            return false;
        }

        currentIndependentVariable_ = lastIndependentVariable_;
        currentState_ = lastState_;
        return true;
    }

    //! Check if minimum step size constraint was violated.
    /*!
     * Returns true if the minimum step size constraint has been violated since this integrator
     * was constructed.
     * \return True if the minimum step size constraint was violated.
     */
    bool isMinimumStepSizeViolated( ) const { return isMinimumStepSizeViolated_; }

protected:

    //! Last used step size.
    /*!
     * Last used step size, passed to either integrateTo( ) or performIntegrationStep( ).
     */
    IndependentVariableType stepSize_;

    //! Current independent variable.
    /*!
     * Current independent variable as computed by performIntegrationStep().
     */
    IndependentVariableType currentIndependentVariable_;

    //! Current state.
    /*!
     * Current state as computed by performIntegrationStep( ).
     */
    StateType currentState_;

    //! Last independent variable.
    /*!
     * Last independent variable value as computed by performIntegrationStep().
     */
    IndependentVariableType lastIndependentVariable_;

    //! Last state.
    /*!
     * Last state as computed by performIntegrationStep( ).
     */
    StateType lastState_;

    //! Sequence for the integrator.
    /*!
     * Rational function sequence for the integrator.
     */
    std::vector< unsigned int > sequence_;

    //! Minimum step size.
    /*!
     * Minimum step size.
     */
    IndependentVariableType minimumStepSize_;

    //! Relative error tolerance.
    /*!
     * Relative error tolerance per element in the state.
     */
    StateType relativeErrorTolerance_;

    StateType absoluteErrorTolerance_;

    //! Flag to indicate whether the minimum step size constraint has been violated.
    /*!
     * Flag to indicate whether the minimum step size constraint has been violated.
     */
    bool isMinimumStepSizeViolated_;

    //! Execute mid-point method.
    /*!
     * Executes mid-point method, given a known state and state derivative.
     * \param stateAtFirstPoint State at first point.
     * \param stateAtCenterPoint State at center point.
     * \param independentVariableAtFirstPoint Independent variable at first point.
     * \param subStepSize Sub step size between successive states used by mid-point method.
     * \param stateAtLastPoint State at last point.
     */
    StateType executeMidPointMethod( StateType stateAtFirstPoint, StateType stateAtCenterPoint,
                                     const IndependentVariableType independentVariableAtFirstPoint,
                                     const IndependentVariableType subStepSize );

    std::vector< std::vector< StateType > > integratedStates_;

    unsigned int maximumNumberOfAttempts_;

    std::vector< double > subSteps_;


private:
};

//! Perform a single integration step.
template < typename IndependentVariableType, typename StateType, typename StateDerivativeType >
StateType
BulirschStoerVariableStepSizeIntegrator< IndependentVariableType, StateType, StateDerivativeType >
::performIntegrationStep( const IndependentVariableType stepSize )
{

   //  std::cout<<stepSize<<" "<<this->currentIndependentVariable_<<std::endl;

    //    std::cout<<"Start "<<stepSize<<" "<<currentIndependentVariable_<<std::endl;
    StateType stateAtFirstPoint_;
    StateType stateAtCenterPoint_;
    StateType stateAtLastPoint_;


    // Compute sub steps to take.
    for ( unsigned int p = 0; p <  subSteps_.size( ); p++ )
    {
        subSteps_.at( p ) = stepSize / static_cast< double >(
                    sequence_.at( p ) );
    }

    for ( unsigned int i = 0; i < sequence_.size( ); i++ )
    {
        // Compute Euler step and set as state at center point for use with mid-point method.
        stateAtCenterPoint_ = currentState_ + subSteps_.at( i )
                * this->stateDerivativeFunction_( currentIndependentVariable_, currentState_ );

        // Apply modified mid-point rule.
        stateAtFirstPoint_ = currentState_;
        IndependentVariableType independentVariableAtFirstPoint_ = currentIndependentVariable_;
        for ( unsigned int j = 0; j < sequence_.at( i ) - 1; j++ )
        {
            stateAtLastPoint_ = executeMidPointMethod(
                        stateAtFirstPoint_, stateAtCenterPoint_,
                        independentVariableAtFirstPoint_, subSteps_.at( i ) );

            if ( j < sequence_.at( i ) - 2 )
            {
                stateAtFirstPoint_ = stateAtCenterPoint_;
                stateAtCenterPoint_ = stateAtLastPoint_;
                independentVariableAtFirstPoint_ += subSteps_.at( i );
            }
        }

        // Apply end-point correction.
        integratedStates_[ i ][ 0 ]
                = 0.5 * ( stateAtLastPoint_ + stateAtCenterPoint_+ subSteps_.at( i ) * this->stateDerivativeFunction_(
                              currentIndependentVariable_ + stepSize, stateAtLastPoint_ ) );

        //        std::cout<<"Integrated state: "<<integratedStates_[ i ][ 0 ].transpose( )<<std::endl;

        for ( unsigned int k = 1; k < i + 1; k++ )
        {
            integratedStates_[ i ][ k ] =
                    integratedStates_[ i ][ k - 1 ] + 1.0 /
                    ( pow( subSteps_.at( i - k ), 2.0 ) / std::pow( subSteps_.at( i ), 2.0 ) - 1.0 )
                    * ( integratedStates_[ i ][ k - 1 ] - integratedStates_[ i - 1 ][ k - 1 ] );
        }



        //        if ( i == sequence_.size( )  )
        //        {
        //            std::cout << "Problem in BS integrator. " << std::endl;
        //        }

//        if( i > 1 )
//        {
//            std::cout<<"i: "<<i<<" "<<( integratedStates_.at( i ).at( i ) - integratedStates_.at( i ).at( i - 1 ) ).array( ).abs( ).maxCoeff( )<<std::endl;
//        }

        if ( ( i > 1 ) &&
             ( ( integratedStates_.at( i ).at( i ) - integratedStates_.at( i ).at( i - 1 ) ).array( ).abs( ).maxCoeff( )
               < relativeErrorTolerance_.array( ).maxCoeff( ) ) )
        {
            // Accept the current step.
            lastIndependentVariable_ = currentIndependentVariable_;
            lastState_ = currentState_;
            currentIndependentVariable_ += stepSize;
            currentState_ = integratedStates_[ i ][ i ];
            stepSize_ = stepSize;

//            std::cout<<"Accepting: "<<i<<" "<<std::endl;
//            std::cout<<"Tolerances : "<<relativeErrorTolerance_.array( ).maxCoeff( )<<" "<<
//                       lastIndependentVariable_<<" "<<( integratedStates_.at( i ).at( i ) - integratedStates_.at( i ).at( i - 1 ) ).array( ).abs( ).maxCoeff( )<<std::endl<<
//                       integratedStates_.at( i ).at( i ).transpose( )<<std::endl;
            //sleep( 1 );
            break;
        }
        else if( i == ( maximumNumberOfAttempts_ - 1 ) )
        {
            double maximumErrorValue = ( integratedStates_.at( i ).at( i ) - integratedStates_.at( i ).at( i - 1 ) ).array( ).abs( ).maxCoeff( );
            //std::cout<<"Reducing step: "<<stepSize<<" ";
            this->stepSize_ = 0.7 * stepSize * std::pow( relativeErrorTolerance_.array( ).maxCoeff( ) / maximumErrorValue,
                                                         1.0 / ( 2.0 * i - 1.0 ) );

//            std::cout<<"Resetting: "<<stepSize<<" "<<this->stepSize_<<" "<<maximumErrorValue<<std::endl;
//            sleep( 1 );

            return performIntegrationStep( this->stepSize_ );
        }
    }

    //    std::cout<<"End"<<std::endl<<std::endl;


    return currentState_;
}

//! Execute mid-point method.
template < typename IndependentVariableType, typename StateType, typename StateDerivativeType >
StateType
BulirschStoerVariableStepSizeIntegrator< IndependentVariableType, StateType, StateDerivativeType >
::executeMidPointMethod( StateType stateAtFirstPoint, StateType stateAtCenterPoint,
                         const IndependentVariableType independentVariableAtFirstPoint,
                         const IndependentVariableType subStepSize )
{
    return stateAtFirstPoint + 2.0 * subStepSize
            * this->stateDerivativeFunction_( independentVariableAtFirstPoint + subStepSize,
                                              stateAtCenterPoint );
}

//! Typedef of variable-step size Bulirsch-Stoer integrator (state/state derivative = VectorXd,
//! independent variable = double).
/*!
 * Typedef of a variable-step size Bulirsch-Stoer integrator with VectorXds as state and state
 * derivative and double as independent variable.
 */
typedef BulirschStoerVariableStepSizeIntegrator< > BulirschStoerVariableStepSizeIntegratorXd;

} // namespace numerical_integrators
} // namespace tudat

#endif // TUDAT_BULIRSCH_STOER_VARIABLE_STEP_SIZE_INTEGRATOR_H
