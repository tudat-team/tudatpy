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

//! Types of sequences available for extrapolation integration methods
enum ExtrapolationMethodStepSequences
{
    bulirsch_stoer_sequence,
    deufelhard_sequence
};

//! Function to retrieve the sequence of number of steps to used for Bulirsch-Stoer integration
/*!
 * Function to retrieve the sequence of number of steps to used for Bulirsch-Stoer integration
 * \param extrapolationMethodStepSequenceType Type of sequence that is to be retrieved
 * \param lengthOfSequence Length of the sequence that is to be retrieved (default 12)
 * \return Function to retrieve the sequence of number of steps to used for Bulirsch-Stoer integration
 */
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
           typename StateDerivativeType = Eigen::VectorXd, typename TimeStepType = double >
class BulirschStoerVariableStepSizeIntegrator :
        public NumericalIntegrator< IndependentVariableType, StateType, StateDerivativeType, TimeStepType >
{
public:

    //! Typedef of the base class.
    /*!
     * Typedef of the base class with all template parameters filled in.
     */
    typedef NumericalIntegrator< IndependentVariableType, StateType, StateDerivativeType, TimeStepType > Base;

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
     * \param maximumStepSize The maximum step size to take.
     * \param relativeErrorTolerance The relative error tolerance, for each individual state
     *          vector element.
     * \param absoluteErrorTolerance The absolute error tolerance, for each individual state
     *          vector element.
     * \param safetyFactorForNextStepSize Safety factor used to scale prediction of next step size.
     * \param maximumFactorIncreaseForNextStepSize Maximum factor increase for next step size.
     * \param minimumFactorDecreaseForNextStepSize Maximum factor decrease for next step size.
     * \sa NumericalIntegrator::NumericalIntegrator.
     */
    BulirschStoerVariableStepSizeIntegrator(
            const std::vector< unsigned int >& sequence,
            const StateDerivativeFunction& stateDerivativeFunction,
            const IndependentVariableType intervalStart,  const StateType& initialState,
            const TimeStepType minimumStepSize,
            const TimeStepType maximumStepSize,
            const StateType& relativeErrorTolerance,
            const StateType& absoluteErrorTolerance,
            const TimeStepType safetyFactorForNextStepSize = 0.6,
            const TimeStepType maximumFactorIncreaseForNextStepSize = 4.0,
            const TimeStepType minimumFactorDecreaseForNextStepSize = 0.1 ):
        Base( stateDerivativeFunction ), currentIndependentVariable_( intervalStart ),
        currentState_( initialState ), lastIndependentVariable_( intervalStart ),
        sequence_( sequence ), minimumStepSize_( minimumStepSize ), maximumStepSize_( maximumStepSize ),
        relativeErrorTolerance_( relativeErrorTolerance ),
        absoluteErrorTolerance_( absoluteErrorTolerance ),
        safetyFactorForNextStepSize_( safetyFactorForNextStepSize ),
        maximumFactorIncreaseForNextStepSize_( maximumFactorIncreaseForNextStepSize ),
        minimumFactorDecreaseForNextStepSize_( minimumFactorDecreaseForNextStepSize ),
        isMinimumStepSizeViolated_( false )
    {
        maximumStepIndex_ = sequence_.size( ) - 1;
        subSteps_.resize( maximumStepIndex_ + 1 );

        integratedStates_.resize( maximumStepIndex_ + 1  );
        for( unsigned int i = 0; i < maximumStepIndex_ + 1 ; i++ )
        {
            integratedStates_[ i ].resize( maximumStepIndex_ + 1  );
        }
    }

    //! Default constructor.
    /*!
     * Default constructor, taking coefficients, a state derivative function, initial conditions,
     * minimum step size and relative error tolerance for all items in the state vector as argument.
     * \param sequence Rational function sequence used by algorithm.
     * \param stateDerivativeFunction State derivative function.
     * \param intervalStart The start of the integration interval.
     * \param initialState The initial state.
     * \param minimumStepSize The minimum step size to take. If this constraint is violated, a
     *          flag will be set that can be retrieved with isMinimumStepSizeViolated( ).
     * \param maximumStepSize The maximum step size to take.
     * \param relativeErrorTolerance The relative error tolerance, equal for all individual state
     *          vector elements.
     * \param absoluteErrorTolerance The absolute error tolerance, for each individual state
     *          vector element.
     * \param safetyFactorForNextStepSize Safety factor used to scale prediction of next step size.
     * \param maximumFactorIncreaseForNextStepSize Maximum factor increase for next step size.
     * \param minimumFactorDecreaseForNextStepSize Maximum factor decrease for next step size.
     * \sa NumericalIntegrator::NumericalIntegrator.
     */
    BulirschStoerVariableStepSizeIntegrator(
            const std::vector< unsigned int >& sequence,
            const StateDerivativeFunction& stateDerivativeFunction,
            const IndependentVariableType intervalStart, const StateType& initialState,
            const TimeStepType minimumStepSize,
            const TimeStepType maximumStepSize,
            const typename StateType::Scalar relativeErrorTolerance = 1.0e-12,
            const typename StateType::Scalar absoluteErrorTolerance = 1.0e-12,
            const TimeStepType safetyFactorForNextStepSize = 0.75,
            const TimeStepType maximumFactorIncreaseForNextStepSize = 4.0,
            const TimeStepType minimumFactorDecreaseForNextStepSize = 0.1 ):
        Base( stateDerivativeFunction ), currentIndependentVariable_( intervalStart ),
        currentState_( initialState ), lastIndependentVariable_( intervalStart ),
        sequence_( sequence ), minimumStepSize_( minimumStepSize ),  maximumStepSize_( maximumStepSize ),
        relativeErrorTolerance_( StateType::Constant( initialState.rows( ), initialState.cols( ),
                                                      relativeErrorTolerance ) ),
        absoluteErrorTolerance_( StateType::Constant( initialState.rows( ), initialState.cols( ),
                                                      absoluteErrorTolerance ) ),
        safetyFactorForNextStepSize_( safetyFactorForNextStepSize ),
        maximumFactorIncreaseForNextStepSize_( maximumFactorIncreaseForNextStepSize ),
        minimumFactorDecreaseForNextStepSize_( minimumFactorDecreaseForNextStepSize ),
        isMinimumStepSizeViolated_( false )
    {
        maximumStepIndex_ = sequence_.size( ) - 1;
        subSteps_.resize( maximumStepIndex_ + 1 );

        integratedStates_.resize( maximumStepIndex_ + 1  );
        for( unsigned int i = 0; i < maximumStepIndex_ + 1 ; i++ )
        {
            integratedStates_[ i ].resize( maximumStepIndex_ + 1  );
        }
    }

    //! Get step size of the next step.
    /*!
     * Returns the step size of the next step.
     * \return Step size to be used for the next step.
     */
    virtual TimeStepType getNextStepSize( ) const { return stepSize_; }

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
    virtual StateType performIntegrationStep( const TimeStepType stepSize )
    {
        StateType stateAtFirstPoint_;
        StateType stateAtCenterPoint_;
        StateType stateAtLastPoint_;

        bool stepSuccessful = 0;

        // Compute sub steps to take.
        for ( unsigned int p = 0; p <  subSteps_.size( ); p++ )
        {
            subSteps_.at( p ) = stepSize / static_cast< double >(
                        sequence_.at( p ) );
        }

        double errorScaleTerm = TUDAT_NAN;
        for ( unsigned int i = 0; i <= maximumStepIndex_; i++ )
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

            for ( unsigned int k = 1; k < i + 1; k++ )
            {
                integratedStates_[ i ][ k ] =
                        integratedStates_[ i ][ k - 1 ] + 1.0 /
                        ( pow( subSteps_.at( i - k ), 2.0 ) / std::pow( subSteps_.at( i ), 2.0 ) - 1.0 )
                        * ( integratedStates_[ i ][ k - 1 ] - integratedStates_[ i - 1 ][ k - 1 ] );
            }

            if( i == maximumStepIndex_ )
            {
                const StateType errorTolerance_ =
                        ( integratedStates_.at( i ).at( i ).array( ).abs( ) * relativeErrorTolerance_.array( ) ).matrix( )
                        + absoluteErrorTolerance_;

                double maximumAllowableErrorValue = errorTolerance_.array( ).maxCoeff( );
                double maximumErrorValue = ( integratedStates_.at( i ).at( i ) - integratedStates_.at( i ).at( i - 1 ) ).array( ).abs( ).maxCoeff( );

                errorScaleTerm = safetyFactorForNextStepSize_ * std::pow( maximumAllowableErrorValue / maximumErrorValue,
                                           ( 1.0 / static_cast< double >( 2 * i - 1 ) ) );

                if( maximumErrorValue < maximumAllowableErrorValue )
                {
                    // Accept the current step.
                    lastIndependentVariable_ = currentIndependentVariable_;
                    lastState_ = currentState_;
                    currentIndependentVariable_ += stepSize;
                    currentState_ = integratedStates_[ i ][ i ];
                    stepSize_ = stepSize;
                    stepSuccessful = true;
                }
                else
                {
                    stepSuccessful = false;
                }
            }
        }
        if( !stepSuccessful )
        {
            if( safetyFactorForNextStepSize_ * errorScaleTerm < minimumFactorDecreaseForNextStepSize_ )
            {
                this->stepSize_ = stepSize * minimumFactorDecreaseForNextStepSize_ * errorScaleTerm;
            }
            else
            {
                this->stepSize_ = stepSize * safetyFactorForNextStepSize_ * errorScaleTerm;
            }

            if( std::fabs( stepSize_ ) < std::fabs( minimumStepSize_ ) )
            {
                isMinimumStepSizeViolated_ = true;
                throw std::runtime_error( "Error in BS integrator, minimum step size exceeded" );
            }
            performIntegrationStep( stepSize_ );
        }
        else
        {
            if( errorScaleTerm > maximumFactorIncreaseForNextStepSize_ )
            {
                this->stepSize_ = stepSize * maximumFactorIncreaseForNextStepSize_;
                if(  std::fabs( this->stepSize_ ) >=  std::fabs( maximumStepSize_ ) )
                {
                   this->stepSize_ = stepSize / ( std::fabs( stepSize ) ) * maximumStepSize_ ;
                }
            }
            else
            {
                this->stepSize_ = stepSize * errorScaleTerm;
            }
        }


        return currentState_;
    }

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

    IndependentVariableType getPreviousIndependentVariable( )
    {
        return lastIndependentVariable_;
    }

    StateType getPreviousState( )
    {
        return lastState_;
    }

private:

    //! Last used step size.
    /*!
     * Last used step size, passed to either integrateTo( ) or performIntegrationStep( ).
     */
    TimeStepType stepSize_;

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
    TimeStepType minimumStepSize_;

    //! Minimum step size.
    /*!
     * Maximum step size.
     */
    TimeStepType maximumStepSize_;

    //! Relative error tolerance.
    /*!
     * Relative error tolerance per element in the state.
     */
    StateType relativeErrorTolerance_;

    //! Absolute error tolerance.
    /*!
     *  Absolute error tolerance per element in the state.
     */
    StateType absoluteErrorTolerance_;

    //! Safety factor for next step size.
    /*!
     * Safety factor used to scale prediction of next step size. This is usually picked between
     * 0.8 and 0.9 (Burden and Faires, 2001).
     */
    TimeStepType safetyFactorForNextStepSize_;

    //! Maximum factor increase for next step size.
    /*!
     * The maximum factor by which the next step size can increase compared to the current value.
     * The need for this maximum stems from a need to ensure that the step size changes do not
     * alias with the dynamics of the model being integrated.
     */
    TimeStepType maximumFactorIncreaseForNextStepSize_;

    //! Minimum factor decrease for next step size.
    /*!
     * The minimum factor by which the next step size can decrease compared to the current value.
     * The need for this minimum stems from a need to ensure that the step size changes do not
     * alias with the dynamics of the model being integrated.
     */
    TimeStepType minimumFactorDecreaseForNextStepSize_;

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
     * \return Result of midpoint method
     */
    StateType executeMidPointMethod( StateType stateAtFirstPoint, StateType stateAtCenterPoint,
                                     const IndependentVariableType independentVariableAtFirstPoint,
                                     const IndependentVariableType subStepSize )
    {
        return stateAtFirstPoint + 2.0 * subStepSize
                * this->stateDerivativeFunction_( independentVariableAtFirstPoint + subStepSize,
                                                  stateAtCenterPoint );
    }

    std::vector< std::vector< StateType > > integratedStates_;

    unsigned int maximumStepIndex_;

    std::vector< double > subSteps_;

};

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
