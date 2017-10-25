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
        relativeErrorTolerance_( relativeErrorTolerance ), isMinimumStepSizeViolated_( false ),
        stepSizeChanged_( true )
    {
        maximumStepIndex_ = sequence_.size( ) - 1;
        subSteps_.resize( maximumStepIndex_ );

        integratedStates_.resize( maximumStepIndex_ );
        for( unsigned int i = 0; i < maximumStepIndex_; i++ )
        {
            integratedStates_[ i ].resize( maximumStepIndex_ );
        }

        initializeWorkParameters( );
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
            const typename StateType::Scalar absoluteErrorTolerance = 1.0e-12 ):
        Base( stateDerivativeFunction ), currentIndependentVariable_( intervalStart ),
        currentState_( initialState ), lastIndependentVariable_( intervalStart ),
        sequence_( sequence ), minimumStepSize_( minimumStepSize ),
        relativeErrorTolerance_( StateType::Constant( initialState.rows( ), initialState.cols( ),
                                                      relativeErrorTolerance ) ),
        absoluteErrorTolerance_( StateType::Constant( initialState.rows( ), initialState.cols( ),
                                                      absoluteErrorTolerance ) ),
        isMinimumStepSizeViolated_( false ),
        stepSizeChanged_( true )
    {
        maximumStepIndex_ = sequence_.size( ) - 1;
        subSteps_.resize( maximumStepIndex_ + 1 );

        integratedStates_.resize( maximumStepIndex_ + 1  );
        for( unsigned int i = 0; i < maximumStepIndex_ + 1 ; i++ )
        {
            integratedStates_[ i ].resize( maximumStepIndex_ + 1  );
        }

        initializeWorkParameters( );
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

    void initializeWorkParameters( )
    {

        aTerms_.resize( maximumStepIndex_ + 1  );
        aTerms_[ 0 ] = sequence_.at( 0 ) + 1;

        alphaTerms_.resize( maximumStepIndex_ + 1  );
        for( unsigned int i = 0; i <= maximumStepIndex_; i++ )
        {
            alphaTerms_[ i ].resize( maximumStepIndex_ + 1  );

            if( i > 0 )
            {
                aTerms_[ i ] = sequence_.at( i ) + aTerms_[ i - 1 ];
            }
        }

        for( unsigned int i = 0; i < maximumStepIndex_ ; i++ )
        {
            for( unsigned int j = 0; j < maximumStepIndex_ ; j++ )
            {
                alphaTerms_[ i ][ j ] = std::pow(
                            relativeErrorTolerance_.array( ).maxCoeff( ), static_cast< double >( aTerms_[ i + 1 ] - aTerms_[ j + 1 ] ) /
                         static_cast< double >( ( 2 * i  - 1 ) * ( aTerms_[ j + 1 ] - aTerms_[ 0 ] + 1 )  ) );
                std::cout<<"Computing alpha: "<<i<<" "<<j<<" "<<alphaTerms_[ i ][ j ]<<std::endl;
            }
        }
    }

    std::vector< std::vector< StateType > > integratedStates_;

    unsigned int maximumStepIndex_;

    unsigned int optimalOrder_;

    std::vector< double > subSteps_;

    bool stepSizeChanged_;

    std::vector< std::vector< double > > alphaTerms_;

    std::vector< int > aTerms_;

private:
};

//! Perform a single integration step.
template < typename IndependentVariableType, typename StateType, typename StateDerivativeType >
StateType
BulirschStoerVariableStepSizeIntegrator< IndependentVariableType, StateType, StateDerivativeType >
::performIntegrationStep( const IndependentVariableType stepSize )
{
    //std::cout<<"Starting: "<<currentIndependentVariable_<<" "<<stepSize<<std::endl;


    //  std::cout<<stepSize<<" "<<this->currentIndependentVariable_<<std::endl;

    //    std::cout<<"Start "<<stepSize<<" "<<currentIndependentVariable_<<std::endl;
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

    if( stepSizeChanged_ )
    {
        optimalOrder_ = maximumStepIndex_;
    }

    double stepSizeChangeFactor = TUDAT_NAN;

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

        if ( ( i > 1 ) )//&& ( ( i >= optimalOrder_ - 1 ) || stepSizeChanged_ ) )
        {
            const StateType errorTolerance_ =
                    ( integratedStates_.at( i ).at( i ).array( ).abs( ) * relativeErrorTolerance_.array( ) ).matrix( )
                    + absoluteErrorTolerance_;

            double maximumErrorValue = errorTolerance_.array( ).maxCoeff( );
            double errorValue = ( integratedStates_.at( i ).at( i ) - integratedStates_.at( i ).at( i - 1 ) ).array( ).abs( ).maxCoeff( );

            double errorScaleTerm = std::pow( errorValue / maximumErrorValue,
                                              static_cast< double >( 1 / ( 2 * i - 1 ) ) );

            //std::cout<<"Scaling errors: "<<i<<" "<<errorValue<<" "<<maximumErrorValue<<" "<<
            //           errorValue / maximumErrorValue<<" "<<errorScaleTerm<<std::endl;

            if( errorValue < maximumErrorValue )
            {
                //std::cout<<"In A"<<std::endl;
                // Accept the current step.
                lastIndependentVariable_ = currentIndependentVariable_;
                lastState_ = currentState_;
                currentIndependentVariable_ += stepSize;
                currentState_ = integratedStates_[ i ][ i ];
                stepSize_ = stepSize;
                stepSuccessful = true;

                if( i < optimalOrder_ )
                {
                    stepSizeChangeFactor = 4.0;
                }

                break;
            }

            if( ( i == maximumStepIndex_ ) || (  i == optimalOrder_ + 1 ) )
            {
                std::cout<<"In A"<<std::endl;
                stepSizeChangeFactor = 0.7 / errorScaleTerm;
                break;
            }
            else if( ( i == optimalOrder_ ) && ( alphaTerms_[ optimalOrder_ - 1 ][ optimalOrder_ ] < errorScaleTerm ) )
            {
                std::cout<<"In B"<<std::endl;
                stepSizeChangeFactor = 1.0 / errorScaleTerm;
                break;
            }
            else if( ( optimalOrder_ == maximumStepIndex_ ) && ( alphaTerms_[ i - 1 ][ maximumStepIndex_ - 1 ] < errorScaleTerm ) )
            {
                std::cout<<"In C "<<i - 1<<" "<<maximumStepIndex_ - 1<<" "<<alphaTerms_[ i - 1 ][ maximumStepIndex_ - 1 ]<<std::endl;
                stepSizeChangeFactor = alphaTerms_[ i - 1 ][ maximumStepIndex_ - 1 ]  * 0.7 / errorScaleTerm;
                break;
            }
            else if( alphaTerms_[ i - 1 ][ optimalOrder_ - 1 ] < errorScaleTerm )
            {
                std::cout<<"In D "<<i - 1<<" "<<optimalOrder_ - 1 <<" "<<alphaTerms_[ i - 1 ][ optimalOrder_ - 1 ] <<std::endl;
                stepSizeChangeFactor = alphaTerms_[ i - 1 ][ optimalOrder_ - 1 ] / errorScaleTerm;
                break;
            }
        }
    }
    if( !stepSuccessful )
    {
        if( ! ( stepSizeChangeFactor == stepSizeChangeFactor ) )
        {
            stepSizeChangeFactor = 0.5;
        }
        else if( stepSizeChangeFactor < 0.1 )
        {
            stepSizeChangeFactor = 0.1;
        }
        std::cout<<"Changing step size and rejecting: "<<stepSizeChangeFactor<<std::endl;
        stepSize_ = stepSize * stepSizeChangeFactor;
        stepSizeChanged_ = true;
        performIntegrationStep( stepSize_ );
    }
    else
    {
        if( stepSizeChangeFactor == stepSizeChangeFactor )
        {
            stepSize_ = stepSize * stepSizeChangeFactor;
            stepSizeChanged_ = true;
            std::cout<<"Changing step size and accepting: "<<stepSizeChangeFactor<<std::endl;
        }
        else
        {
            stepSizeChanged_ = false;
        }
    }

    //std::cout<<"Finished: "<<std::endl<<std::endl;

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
