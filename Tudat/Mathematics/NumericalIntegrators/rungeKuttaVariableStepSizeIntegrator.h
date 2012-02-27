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
 *      110826    K. Kumar          File created.
 *      110909    E.A.G. Heeren     Modified to be included in Tudat. Added if-statement to ensure
 *                                  a maximum stepsize change.
 *      110912    K. Kumar          Minor changes.
 *      111021    F.M. Engelen      Added interface for RKF45 and RKF56.
 *      111110    E.A.G Heeren      Minor changes.
 *      120204    B. Tong Minh      Adapted for new core integrators interface.
 *      120213    K. Kumar          Modified getCurrentInterval() to getIndependentVariable().
 *
 *    References
 *      Eagle, C.D. Celestial Computing with MATLAB, Science Software,
 *          http://www.cdeagle.com/ccmatlab/ccmatlab.pdf, last accessed: 28 August, 2011.
 *
 */

#ifndef TUDAT_RUNGE_KUTTA_VARIABLE_STEP_SIZE_INTEGRATOR_H
#define TUDAT_RUNGE_KUTTA_VARIABLE_STEP_SIZE_INTEGRATOR_H

#include <vector>
#include <TudatCore/Mathematics/NumericalIntegrators/numericalIntegrator.h>
#include "Tudat/Mathematics/NumericalIntegrators/rungeKuttaCoefficients.h"

namespace tudat
{
namespace mathematics
{
namespace numerical_integrators
{

//! Class that implements the Runge-Kutta variable stepsize integrator.
/*!
 * Class that implements the Runge-Kutta variable step size integrator.
 * \tparam StateType The type of the state. This type should be an Eigen::Matrix derived type.
 * \tparam StateDerivativeType The type of the state derivative. This type should be an
 *          Eigen::Matrix derived type.
 * \tparam IndependentVariableType The type of the independent variable. This type should be
 *          either a float or double.
 * \sa NumericalIntegrator.
 */
template < typename IndependentVariableType = double, typename StateType = Eigen::VectorXd,
           typename StateDerivativeType = Eigen::VectorXd >
class RungeKuttaVariableStepSizeIntegrator :
        public NumericalIntegrator<IndependentVariableType, StateType, StateDerivativeType>
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
     * Default constructor, taking coefficientsm a state derivative function, initial conditions,
     * minimum step size and relative error tolerance per item in the state vector as argument.
     * \param coefficients Coefficients to use with this integrator.
     * \param stateDerivativeFunction State derivative function.
     * \param intervalStart The start of the integration interval.
     * \param initialState The initial state.
     * \param minimumStepSize The minimum step size to take. If this constraint is violated, a
     *          flag will be set that can be retrieved with isMinimumStepSizeViolated( ).
     * \param relativeErrorTolerance The relative error tolerance, for each individual state
     *          vector element.
     * \sa NumericalIntegrator::NumericalIntegrator.
     */
    RungeKuttaVariableStepSizeIntegrator( const RungeKuttaCoefficients& coefficients,
                                          const StateDerivativeFunction& stateDerivativeFunction,
                                          const IndependentVariableType intervalStart,
                                          const StateType& initialState,
                                          const IndependentVariableType minimumStepSize,
                                          const StateType& relativeErrorTolerance ) :
        Base( stateDerivativeFunction ), currentIndependentVariable_( intervalStart ),
        currentState_( initialState ), lastIndependentVariable_( intervalStart ),
        coefficients_( coefficients ), minimumStepSize_( minimumStepSize ),
        relativeErrorTolerance_( relativeErrorTolerance ), isMinimumStepSizeViolated_( false ) { }

    //! Default constructor.
    /*!
     * Default constructor, taking coefficients a state derivative function, initial conditions,
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
    RungeKuttaVariableStepSizeIntegrator( const RungeKuttaCoefficients& coefficients,
                                          const StateDerivativeFunction& stateDerivativeFunction,
                                          const IndependentVariableType intervalStart,
                                          const StateType& initialState,
                                          const IndependentVariableType minimumStepSize = 1.0e-15,
                                          const typename StateType::Scalar
                                                relativeErrorTolerance = 1.0e-12 ) :
        Base( stateDerivativeFunction ), currentIndependentVariable_( intervalStart ),
        currentState_( initialState ), lastIndependentVariable_( intervalStart ),
        coefficients_( coefficients ), minimumStepSize_( minimumStepSize ),
        relativeErrorTolerance_( StateType::Constant( initialState.rows( ), initialState.cols( ),
                                                      relativeErrorTolerance ) ),
        isMinimumStepSizeViolated_( false ) { }

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
    bool const isMinimumStepSizeViolated( ) const { return isMinimumStepSizeViolated_; }

    //! Set minimum step size constraint.
    /*!
     * Sets the minimum step size constraint.
     * \param minimumStepSize The new minimum step size.
     */
    void setMinimumStepSize( const IndependentVariableType minimumStepSize )
    {
        minimumStepSize_ = minimumStepSize;
    }

    //! Set relative error tolerance per item in the state vector.
    /*!
     * Sets the relative error tolerance per item in the state vector.
     * \param relativeErrorTolerance The relative error tolerance, for each individual state
     *          vector element.
     */
    void setRelativeErrorTolerance( const StateType& relativeErrorTolerance )
    {
        relativeErrorTolerance_ = relativeErrorTolerance;
    }

    //! Set the relative error tolerance for all items in the state vector.
    /*!
     * Set the relative error tolerance for all items in the state vector.
     * \param relativeErrorTolerance The relative error tolerance, equal for all individual state
     *          vector elements.
     */
    void setRelativeErrorTolerance( const typename StateType::Scalar relativeErrorTolerance )
    {
        setRelativeErrorTolerance( StateType::Constant( currentState_.rows( ),
                                                        currentState_.cols( ),
                                                        relativeErrorTolerance ) );
    }

protected:

    //! Computes the next step size and validates the result.
    /*!
     * Computes the next step size based on a higher and lower order estimate, determines if the
     * error is within bounds and returns a new step size.
     * \param lowerOrderEstimate The integrated result with the lower order coefficients.
     * \param higherOrderEstimate The integrated result with the higher order coefficients.
     * \param stepSize The step size used to obtain these results.
     * \return True if teh error was within bounds, false otherwise.
     */
    virtual bool computeNextStepSizeAndValidateResult( const StateType& lowerOrderEstimate,
                                                       const StateType& higherOrderEstimate,
                                                       const IndependentVariableType stepSize );

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

    //! Coefficients for the integrator.
    /*!
     * Coefficients for the integrator, as defined by the Butcher tableau.
     */
    RungeKuttaCoefficients coefficients_;

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

    //! Flag to indicate whether the minimum step size constraint has been violated.
    /*!
     * Flag to indicate whether the minimum step size constraint has been violated.
     */
    bool isMinimumStepSizeViolated_;

};

//! Perform a single integration step.
template < typename IndependentVariableType, typename StateType, typename StateDerivativeType >
StateType
RungeKuttaVariableStepSizeIntegrator< IndependentVariableType, StateType, StateDerivativeType >
::performIntegrationStep( const IndependentVariableType stepSize )
{
    // Define and allocated vector for the number of stages
    std::vector<StateDerivativeType> stateDerivatives; // vector of k_i
    stateDerivatives.reserve( coefficients_.cCoefficients.rows( ) );

    StateType lowerOrderEstimate( currentState_ ), higherOrderEstimate( currentState_ );

    // Compute the k_i state derivatives per stage
    for ( int stage = 0; stage < coefficients_.cCoefficients.rows( ); stage++ )
    {
        // Compute the state to pass to the state derivative for this stage
        StateType state( currentState_ );
        for ( int column = 0; column < stage; column++ )
        {
            state += stepSize * coefficients_.aCoefficients( stage, column ) *
                    stateDerivatives[column];
        }

        // Compute the state derivative
        stateDerivatives.push_back(
                    stateDerivativeFunction_( currentIndependentVariable_ +
                                              coefficients_.cCoefficients( stage ) * stepSize,
                                              state ) );

        // Update the estimate
        lowerOrderEstimate += coefficients_.bCoefficients( 0, stage ) * stepSize *
                stateDerivatives[stage];
        higherOrderEstimate += coefficients_.bCoefficients( 1, stage ) * stepSize *
                stateDerivatives[stage];
    }

    // Determine if the error was within bounds and compute a new step size
    if ( computeNextStepSizeAndValidateResult( lowerOrderEstimate, higherOrderEstimate, stepSize ) )
    {
        // Accept the current step
        lastIndependentVariable_ = currentIndependentVariable_;
        lastState_ = currentState_;
        currentIndependentVariable_ += stepSize;
        currentState_ = higherOrderEstimate;

        return currentState_;
    }
    else
    {
        // Reject current step.
        return performIntegrationStep( stepSize_ );
    }
}

//! Compute the next step size and validate the result.
template < typename IndependentVariableType, typename StateType, typename StateDerivativeType >
bool
RungeKuttaVariableStepSizeIntegrator< IndependentVariableType, StateType, StateDerivativeType >
::computeNextStepSizeAndValidateResult(
        const StateType& lowerOrderEstimate, const StateType& higherOrderEstimate,
        const IndependentVariableType stepSize )
{
    // Calculate the error in the lower order estimate
    StateType truncationError_ = higherOrderEstimate - lowerOrderEstimate;
    StateType dampedAbsoluteTolerance_ = ( currentState_.array( ).abs( ) *
                                           relativeErrorTolerance_.array( ) ).matrix( )
            + relativeErrorTolerance_;

    StateType relativeTruncationError_ = truncationError_.array( ).abs( ) /
            dampedAbsoluteTolerance_.array( );

    typename StateType::Scalar maximumErrorInState_ = relativeTruncationError_.maxCoeff( );

    if ( std::fabs( maximumErrorInState_ ) <=
         std::numeric_limits< typename StateType::Scalar >::epsilon( ) )
    {
        // The estimates are equal, so there is no proper guess for a new step size. Continue with
        // the current step size.
        stepSize_ = stepSize;
    }
    else
    {
        // Compute the new step size.
        stepSize_ = 0.8 * stepSize * std::pow( 1.0 / maximumErrorInState_,
                                               1.0 / coefficients_.higherOrder );
    }

    // Check whether change in stepsize does not exceed bounds.
    // If the stepsize is reduced to less than 1/10th, set to 1/10th.
    // If the stepsize is increased to more than 4-times, set to 4-times.
    // These bounds are necessary to prevent the stepsize changes from aliasing
    // with the dynamics of the system of ODEs.
    // The bounds can be found in (Burden and Faires, 2001).
    if ( stepSize_ / stepSize < 0.1 )
    {
        stepSize_ = 0.1 * stepSize;
    }

    else if ( stepSize_ / stepSize > 4.0 )
    {
        stepSize_ = 4.0 * stepSize;
    }

    // Check if minimum stepsize condition is violated.
    if ( stepSize_ < minimumStepSize_ )
    {
        stepSize_ = minimumStepSize_;
        isMinimumStepSizeViolated_ = true;

        // Return true if the previous step size is already equal to the minimum step size to
        // allow the integrator to continue.
        return std::abs( stepSize - minimumStepSize_ ) <=
                std::numeric_limits< typename StateType::Scalar >::epsilon( );
    }

    // Check if computed error in state is too large and reject step if true.
    return maximumErrorInState_ <= 1.0;
}

//! Typedef of variable-step size Runge-Kutta integrator (state/state derivative = VectorXd,
//! independent variable = double).
/*!
 * Typedef of a variable-step size Runge-Kutta integrator with VectorXds as state and state
 * derivative and double as independent variable.
 */
typedef RungeKuttaVariableStepSizeIntegrator< > RungeKuttaVariableStepSizeIntegratorXd;

} // namespace numerical_integrators
} // namespace mathematics
} // namespace tudat

#endif // TUDAT_RUNGE_KUTTA_VARIABLE_STEP_SIZE_INTEGRATOR_H
