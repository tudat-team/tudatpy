/*! \file variableStepsizeRungeKuttaIntegrator.cpp
 *    Source file that defines the 7(8)th-order, variable stepsize, Runge-Kutta
 *    integrator.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author/Checker    : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author/Checker    : E.A.G. Heeren
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : E.A.G.Heeren@student.tudelft.nl
 *
 *    Author            : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Date created      : 26 August, 2011
 *    Last modified     : 10 November, 2011
 *
 *    References
 *      Burden, R.L., Faires, J.D. Numerical Analysis, 7th Edition, Books/Cole,
 *          2001.
 *      Eagle, C.D. Celestial Computing with MATLAB, Science Software,
 *          http://www.cdeagle.com/ccmatlab/ccmatlab.pdf, last accessed:
 *          28 August, 2011.
 *
 *    Notes
 *      This is a translation of the rkf78.m Matlab integrator.
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      111021    F.M. Engelen      Generalised and added coefficients for RKF45 and RKF56
 *      111110    E.A.G Heeren      Minor changes.
 *      111118    F.M. Engelen      Solved bug, not resetting previousStepsize_.
 */

// Include statements.
#include <cmath>
#include <iostream>
#include "Mathematics/basicMathematicsFunctions.h"
#include "Mathematics/NumericalIntegrators/rungeKuttaFehlbergVariableStepsize.h"

//! Tudat library namespace.
namespace tudat
{

//! Default constructor.
RungeKuttaFehlBergVariableStepsizeIntegrator::RungeKuttaFehlBergVariableStepsizeIntegrator(
        RungeKuttaFehlbergVariableStepsizeIntegratorType
        rungeKuttaFehlbergVariableStepsizeIntegratorType ) : currentTime_( - 0.0 ),
    dampedAbsoluteTolerance_( -0.0 ), minimumStepsize_( 1.0e-15 ), previousStepsize_ ( -0.0 ),
    relativeErrorTolerance_( 1.0e-12 ), relativeTruncationError_( -0.0 ), truncationError_( -0.0 )
{
    switch( rungeKuttaFehlbergVariableStepsizeIntegratorType )
    {
    case rungeKuttaFelhberg45VariableStepsize:

        numberOfStages_ = 6;
        integratorOrder_ = 5;
        break;

    case rungeKuttaFelhberg56VariableStepsize:

        numberOfStages_ = 8;
        integratorOrder_ = 6;
        break;

    case rungeKuttaFelhberg78VariableStepsize:

        numberOfStages_ = 13;
        integratorOrder_ = 8;
        break;

    default :
        std::cerr << "This is not an available Runge-Kutta-Fehlberg-type integrator.";
    }

    // Initialize coefficient vectors/matrices.
    aCoefficients_ = MatrixXd::Zero( numberOfStages_, numberOfStages_ );
    bCoefficients_ = MatrixXd::Zero( 2, numberOfStages_ );
    cCoefficients_ = VectorXd::Zero( numberOfStages_ );

    // Set Runge-Kutta-Fehlberg coefficients.
    setCoefficients_( rungeKuttaFehlbergVariableStepsizeIntegratorType );
}

//! Integrate.
void RungeKuttaFehlBergVariableStepsizeIntegrator::integrate( )
{
    // Using declarations.
    using std::fabs;
    using std::cerr;
    using std::endl;

    // Clear vector of current states.
    vectorOfCurrentStates_.clear( );

    // Set stepsize to value of initial stepsize.
    stepsize_ = initialStepsize_;

    // Set the previous stepsize to zero.
    previousStepsize_ = -0.0;

    // Compute size of integration interval.
    integrationInterval_ = integrationIntervalEnd_ - integrationIntervalStart_;

    // Add pointer to initial state to container vector.
    vectorOfCurrentStates_.push_back( *pointerToInitialState_ );

    // Set dimension of initial state.
    dimensionOfState_ = pointerToInitialState_->state.rows( );

    // Set dimension of final state to correspond with dimension of initial
    // state.
    finalState_.state.setZero( dimensionOfState_ );

    // Set dimension of state derivative to correpsond with dimension of
    // initial state.
    stateDerivative_.state.setZero( dimensionOfState_ );

    // Set current point in integration interval to start of
    // integration interval.
    integrationIntervalCurrentPoint_ = integrationIntervalStart_;
    currentTime_ = integrationIntervalCurrentPoint_;

    // Set dimensions of f-Matrix.
    fMatrix.setZero( dimensionOfState_, numberOfStages_ );

    // Loop until end of integration is reached.
    while ( true )
    {
        // Current state and current time.
        currentState_ = vectorOfCurrentStates_.at( 0 );
        currentTime_ = currentTime_ + previousStepsize_;

        integrationIntervalCurrentPoint_ = currentTime_;

        // Check for last step.
        if ( fabs( stepsize_ )
             > fabs( integrationIntervalEnd_ - integrationIntervalCurrentPoint_ ) )
        {
            stepsize_ = integrationIntervalEnd_ - integrationIntervalCurrentPoint_;
        }


        // Check for end of integration period and break from loop if end is reached.
        // Removed fabs statement, to prevent overrun (although this should not be possible)
        // Changed
        if (  integrationIntervalEnd_ - integrationIntervalCurrentPoint_
             < minimumStepsize_ )
        {
            break;
        }

        // Compute state derivative.
        computeStateDerivative_( integrationIntervalCurrentPoint_,
                                 &vectorOfCurrentStates_.at( 0 ),
                                 &stateDerivative_ );

        // Copy computed state derivative to first column of fMatrix.
        fMatrix.col( 0 ) = stateDerivative_.state;

        // Compute current pth-order estimate by looping over k-coefficients.
        for ( unsigned int k = 0; k < numberOfStages_; k++ )
        {
            // Loop over elements of state.
            for ( int i = 0; i < vectorOfCurrentStates_.at( 0 ).state.rows( ); i++ )
            {
                // Compute ith element of modified state.
                vectorOfCurrentStates_.at( 0 ).state[ i ] = currentState_.state[ i ]  + stepsize_
                        * ( aCoefficients_.row( k ).array( ) * fMatrix.row( i ).array( ) ).sum( );
            }

            // Compute modified current point in integration interval.
            integrationIntervalCurrentPoint_ = currentTime_ + cCoefficients_( k ) * stepsize_;

            // Compute state derivative using modified parameters.
            computeStateDerivative_( integrationIntervalCurrentPoint_,
                                     &vectorOfCurrentStates_.at( 0 ),
                                     &stateDerivative_ );

            // Copy computed state derivative to kth-column of fMatrix.
            fMatrix.col( k ) = stateDerivative_.state.col( 0 );
        }

        // Set error in state to relative error tolerance.
        errorInState_ = relativeErrorTolerance_;

        // Loop over state elements to compute relative error in pth-order estimate
        // using (p+1)th-order estimate.
        for ( int i = 0; i < vectorOfCurrentStates_.at( 0 ).state.rows( ); i++ )
        {
            // Compute ith element of state.
           vectorOfCurrentStates_.at( 0 ).state[ i ] = currentState_.state[ i ];

            double lowerOrderEstimate = 0.0;
            double higherOrderEstimate = 0.0;

            // Compute truncation error.
            for ( unsigned int k = 0; k < numberOfStages_; k++)
            {
                lowerOrderEstimate +=
                        stepsize_ * bCoefficients_( 0, k ) * fMatrix( i, k );

                higherOrderEstimate +=
                        stepsize_ * bCoefficients_( 1, k ) * fMatrix( i, k );
            }


            truncationError_ = fabs( higherOrderEstimate - lowerOrderEstimate );

            vectorOfCurrentStates_.at( 0 ).state[ i ] =
                    currentState_.state[ i ] + higherOrderEstimate;

            // Compute damped absolute tolerance.
            dampedAbsoluteTolerance_ = fabs( vectorOfCurrentStates_.at( 0).state[ i ] )
                                       * relativeErrorTolerance_ + relativeErrorTolerance_;

            // Compute relative truncation error.
            relativeTruncationError_ = truncationError_ / dampedAbsoluteTolerance_;

            // Check computed error.
            if ( relativeTruncationError_ > errorInState_ )
            {
                errorInState_ = relativeTruncationError_;
            }
        }

        // Set Previous stepsize to current stepsize.
        previousStepsize_ = stepsize_;

        // Compute new step size.

        stepsize_ = 0.8 * stepsize_ * pow( 1.0 / errorInState_, 1.0 / integratorOrder_ );


        // Check whether change in stepsize does not exceed bounds.
        // If the stepsize is reduced to less than 1/10th, set to 1/10th.
        // If the stepsize is increased to more than 4-times, set to 4-times.
        // These bounds are necessary to prevent the stepsize changes from aliasing
        // with the dynamics of the system of ODEs.
        // The bounds can be found in (Burden and Faires, 2001).
        if ( stepsize_ / previousStepsize_ < 0.1 )
        {
            stepsize_ = 0.1 * previousStepsize_;
        }

        else if ( stepsize_ / previousStepsize_ > 4.0 )
        {
            stepsize_ = 4.0 * previousStepsize_;
        }

        // Check if minimum stepsize condition is violated.
        if ( stepsize_ < minimumStepsize_ )
        {
            cerr << "Minimum stepsize condition is violated. Stepsize is: " << stepsize_ << endl;
            cerr << "RKF78 integrator failed!" << endl;
        }

        // Check if computed error in state is too large and reject step if true.

        if ( errorInState_ > 1.0 )
        {
            // Reject current step.
            previousStepsize_ = 0.0;
            vectorOfCurrentStates_.at( 0 ).state = currentState_.state;
        }

        // Else, accept step and store final state computed.
        else
        {
            // Set final state.
            finalState_ = vectorOfCurrentStates_.at( 0 );
        }
    }
}

//! Set coefficients of integration scheme.
void RungeKuttaFehlBergVariableStepsizeIntegrator::setCoefficients_(
        RungeKuttaFehlbergVariableStepsizeIntegratorType
        rungeKuttaFelhbergVariableStepsizeIntegratorCoefficientSet )
{
    switch( rungeKuttaFelhbergVariableStepsizeIntegratorCoefficientSet )
    {
    case rungeKuttaFelhberg45VariableStepsize:

        // Define a-coefficients for the Runge-Kutta-Fehlberg method of order 5
        // with an embedded 4th-order method for stepsize cont rol and a total of 6 stages.
        aCoefficients_( 1, 0 ) = 1.0 / 4.0;

        aCoefficients_( 2, 0 ) = 3.0 / 32.0;
        aCoefficients_( 2, 1 ) = 9.0 / 32.0;

        aCoefficients_( 3, 0 ) =  1932.0 / 2197.0;
        aCoefficients_( 3, 1 ) = -7200.0 / 2197.0;
        aCoefficients_( 3, 2 ) =  7296.0 / 2197.0;

        aCoefficients_( 4, 0 ) = 439.0 / 216.0;
        aCoefficients_( 4, 1 ) = -8.0;
        aCoefficients_( 4, 2 ) = 3680.0 / 513.0;
        aCoefficients_( 4, 3 ) = -845.0 / 4104.0;

        aCoefficients_( 5, 0 ) = -8.0 / 27.0;
        aCoefficients_( 5, 1 ) = 2.0;
        aCoefficients_( 5, 2 ) = -3544.0 / 2565.0;
        aCoefficients_( 5, 3 ) = 1859.0 / 4104.0;
        aCoefficients_( 5, 4 ) = -11.0 / 40.0;


        // Define c-coefficients for the Runge-Kutta-Fehlberg method of order 5
        // with an embedded 4th-order method for stepsize control and a total of 6 stages.
        cCoefficients_( 0 ) = 0.0;
        cCoefficients_( 1 ) = 1.0 / 4.0;
        cCoefficients_( 2 ) = 3.0 / 8.0;
        cCoefficients_( 3 ) = 12.0 / 13.0;
        cCoefficients_( 4 ) = 1.0;
        cCoefficients_( 5 ) = 1.0 / 2.0;


        // Define b-coefficients for the Runge-Kutta method of order 5
        // with an embedded 4th-order method for stepsize control and a total of 6 stages.
        bCoefficients_( 0, 0 ) = 25.0 / 216.0;
        bCoefficients_( 0, 2 ) = 1408.0 / 2565.0;
        bCoefficients_( 0, 3 ) = 2197.0 / 4104.0;
        bCoefficients_( 0, 4 ) = -1.0 / 5.0;

        bCoefficients_( 1, 0 ) = 16.0 / 135.0;
        bCoefficients_( 1, 2 ) = 6656.0 / 12825.0;
        bCoefficients_( 1, 3 ) = 28561.0 / 56430.0;
        bCoefficients_( 1, 4 ) = -9.0 / 50.0;
        bCoefficients_( 1, 5 ) = 2.0 / 55.0;

        break;

    case rungeKuttaFelhberg56VariableStepsize:

        // Define a-coefficients for the Runge-Kutta-Fehlberg method of order 6
        // with an embedded 5th-order method for stepsize control and a total of 8 stages.
        aCoefficients_( 1, 0 ) = 1.0 / 6.0;

        aCoefficients_( 2, 0 ) =  4.0 / 75.0;
        aCoefficients_( 2, 1 ) = 16.0 / 75.0;

        aCoefficients_( 3, 0 ) =  5.0 / 6.0;
        aCoefficients_( 3, 1 ) = -8.0 / 3.0;
        aCoefficients_( 3, 2 ) =  5.0 / 2.0;

        aCoefficients_( 4, 0 ) =  -8.0 / 5.0;
        aCoefficients_( 4, 1 ) = 144.0 / 25.0;
        aCoefficients_( 4, 2 ) =  -4.0;
        aCoefficients_( 4, 3 ) =  16.0 / 25.0;

        aCoefficients_( 5, 0 ) = 361.0 / 320.0;
        aCoefficients_( 5, 1 ) = -18.0 / 5.0;
        aCoefficients_( 5, 2 ) = 407.0 / 128.0;
        aCoefficients_( 5, 3 ) = -11.0 / 80.0;
        aCoefficients_( 5, 4 ) =  55.0 / 128.0;

        aCoefficients_( 6, 0 ) = -11.0 / 640.0;
        aCoefficients_( 6, 1 ) =   0.0;
        aCoefficients_( 6, 2 ) =  11.0 / 256.0;
        aCoefficients_( 6, 3 ) = -11.0 / 160.0;
        aCoefficients_( 6, 4 ) =  11.0 / 256.0;
        aCoefficients_( 6, 5 ) =   0.0;

        aCoefficients_( 7, 0 ) = 93.0 / 640.0;
        aCoefficients_( 7, 1 ) = -18.0 / 5.0;
        aCoefficients_( 7, 2 ) = 803.0 / 256.0;
        aCoefficients_( 7, 3 ) = -11.0 / 160.0;
        aCoefficients_( 7, 4 ) = 99.0 / 256.0;
        aCoefficients_( 7, 5 ) = 0.0;
        aCoefficients_( 7, 6 ) = 1.0;

        // Define c-coefficients for the Runge-Kutta-Fehlberg method of order 6
        // with an embedded 5th-order method for stepsize control and a total of 8 stages.
        cCoefficients_( 1 ) = 1.0 / 6.0;
        cCoefficients_( 2 ) = 4.0 / 15.0;
        cCoefficients_( 3 ) = 2.0 / 3.0;
        cCoefficients_( 4 ) = 4.0 / 5.0;
        cCoefficients_( 5 ) = 1.0;
        cCoefficients_( 7 ) = 1.0;

        // Define b-coefficients for the Runge-Kutta method of order 8
        // with an embedded 7th-order method for stepsize control and a total of 13 stages.
        bCoefficients_( 0, 0 ) = 31.0 / 384.0;
        bCoefficients_( 0, 2 ) = 1125.0 / 2816.0;
        bCoefficients_( 0, 3 ) = 9.0 / 32.0;
        bCoefficients_( 0, 4 ) = 125.0 / 768.0;
        bCoefficients_( 0, 5 ) = 5.0 / 66.0;

        bCoefficients_( 1, 0 ) = 7.0 / 1408.0;
        bCoefficients_( 1, 2 ) = 1125.0 / 2816.0;
        bCoefficients_( 1, 3 ) = 9.0 / 32.0;
        bCoefficients_( 1, 4 ) = 125.0 / 768.0;
        bCoefficients_( 1, 6 ) = 5.0 / 66.0;
        bCoefficients_( 1, 7 ) = 5.0 / 66.0;

        break;


    case rungeKuttaFelhberg78VariableStepsize:

        // Define a-coefficients for the Runge-Kutta-Fehlberg method of order 8
        // with an embedded 7th-order method for stepsize control and a total of 13 stages.
        aCoefficients_( 1, 0 ) = 2.0 / 27.0;

        aCoefficients_( 2, 0 ) = 1.0 / 36.0;
        aCoefficients_( 2, 1 ) = 1.0 / 12.0;

        aCoefficients_( 3, 0 ) = 1.0 / 24.0;
        aCoefficients_( 3, 1 ) = 0.0;
        aCoefficients_( 3, 2 ) = 1.0 / 8.0;

        aCoefficients_( 4, 0 ) = 5.0 / 12.0;
        aCoefficients_( 4, 1 ) = 0.0;
        aCoefficients_( 4, 2 ) = -25.0 / 16.0;
        aCoefficients_( 4, 3 ) = 25.0 / 16.0;

        aCoefficients_( 5, 0 ) = 1.0 / 20.0;
        aCoefficients_( 5, 1 ) = 0.0;
        aCoefficients_( 5, 2 ) = 0.0;
        aCoefficients_( 5, 3 ) = 1.0 / 4.0;
        aCoefficients_( 5, 4 ) = 1.0 / 5.0;

        aCoefficients_( 6, 0 ) = -25.0 / 108.0;
        aCoefficients_( 6, 1 ) = 0.0;
        aCoefficients_( 6, 2 ) = 0.0;
        aCoefficients_( 6, 3 ) = 125.0 / 108.0;
        aCoefficients_( 6, 4 ) = -65.0 / 27.0;
        aCoefficients_( 6, 5 ) = 125.0 / 54.0;

        aCoefficients_( 7, 0 ) = 31.0 / 300.0;
        aCoefficients_( 7, 1 ) = 0.0;
        aCoefficients_( 7, 2 ) = 0.0;
        aCoefficients_( 7, 3 ) = 0.0;
        aCoefficients_( 7, 4 ) = 61.0 / 225.0;
        aCoefficients_( 7, 5 ) = -2.0 / 9.0;
        aCoefficients_( 7, 6 ) = 13.0 / 900.0;

        aCoefficients_( 8, 0 ) = 2.0;
        aCoefficients_( 8, 1 ) = 0.0;
        aCoefficients_( 8, 2 ) = 0.0;
        aCoefficients_( 8, 3 ) = -53.0 / 6.0;
        aCoefficients_( 8, 4 ) = 704.0 / 45.0;
        aCoefficients_( 8, 5 ) = -107.0 / 9.0;
        aCoefficients_( 8, 6 ) = 67.0 / 90.0;
        aCoefficients_( 8, 7 ) = 3.0;

        aCoefficients_( 9, 0 ) = -91.0 / 108.0;
        aCoefficients_( 9, 1 ) = 0.0;
        aCoefficients_( 9, 2 ) = 0.0;
        aCoefficients_( 9, 3 ) = 23.0 / 108.0;
        aCoefficients_( 9, 4 ) = -976.0 / 135.0;
        aCoefficients_( 9, 5 ) = 311.0 / 54.0;
        aCoefficients_( 9, 6 ) = -19.0 / 60.0;
        aCoefficients_( 9, 7 ) = 17.0 / 6.0;
        aCoefficients_( 9, 8 ) = -1.0 / 12.0;

        aCoefficients_( 10, 0 ) = 2383.0 / 4100.0;
        aCoefficients_( 10, 1 ) = 0.0;
        aCoefficients_( 10, 2 ) = 0.0;
        aCoefficients_( 10, 3 ) = -341.0 / 164.0;
        aCoefficients_( 10, 4 ) = 4496.0 / 1025.0;
        aCoefficients_( 10, 5 ) = -301.0 / 82.0;
        aCoefficients_( 10, 6 ) = 2133.0 / 4100.0;
        aCoefficients_( 10, 7 ) = 45.0 / 82.0;
        aCoefficients_( 10, 8 ) = 45.0 / 164.0;
        aCoefficients_( 10, 9 ) = 18.0 / 41.0;

        aCoefficients_( 11, 0 ) = 3.0 / 205.0;
        aCoefficients_( 11, 1 ) = 0.0;
        aCoefficients_( 11, 2 ) = 0.0;
        aCoefficients_( 11, 3 ) = 0.0;
        aCoefficients_( 11, 4 ) = 0.0;
        aCoefficients_( 11, 5 ) = -6.0 / 41.0;
        aCoefficients_( 11, 6 ) = -3.0 / 205.0;
        aCoefficients_( 11, 7 ) = -3.0 / 41.0;
        aCoefficients_( 11, 8 ) = 3.0 / 41.0;
        aCoefficients_( 11, 9 ) = 6.0 / 41.0;
        aCoefficients_( 11, 10 ) = 0.0;

        aCoefficients_( 12, 0 ) = -1777.0 / 4100.0;
        aCoefficients_( 12, 1 ) = 0.0;
        aCoefficients_( 12, 2 ) = 0.0;
        aCoefficients_( 12, 3 ) = -341.0 / 164.0;
        aCoefficients_( 12, 4 ) = 4496.0 / 1025.0;
        aCoefficients_( 12, 5 ) = -289.0 / 82.0;
        aCoefficients_( 12, 6 ) = 2193.0 / 4100.0;
        aCoefficients_( 12, 7 ) = 51.0 / 82.0;
        aCoefficients_( 12, 8 ) = 33.0 / 164.0;
        aCoefficients_( 12, 9 ) = 12.0 / 41.0;
        aCoefficients_( 12, 10 ) = 0.0;
        aCoefficients_( 12, 11 ) = 1.0;

        // Define c-coefficients for the Runge-Kutta-Fehlberg method of order 8
        // with an embedded 7th-order method for stepsize control and a total of 13 stages.
        cCoefficients_( 0 ) = 0.0;
        cCoefficients_( 1 ) = 2.0 / 27.0;
        cCoefficients_( 2 ) = 1.0 / 9.0;
        cCoefficients_( 3 ) = 1.0 / 6.0;
        cCoefficients_( 4 ) = 5.0 / 12.0;
        cCoefficients_( 5 ) = 1.0 / 2.0;
        cCoefficients_( 6 ) = 5.0 / 6.0;
        cCoefficients_( 7 ) = 1.0 / 6.0;
        cCoefficients_( 8 ) = 2.0 / 3.0;
        cCoefficients_( 9 ) = 1.0 / 3.0;
        cCoefficients_( 10 ) = 1.0;
        cCoefficients_( 11 ) = 0.0;
        cCoefficients_( 12 ) = 1.0;

        // Define b-coefficients for the Runge-Kutta method of order 8
        // with an embedded 7th-order method for stepsize control and a total of 13 stages.
        bCoefficients_( 0, 0 ) = 41.0 / 840.0;
        bCoefficients_( 0, 5 ) = 34.0 / 105.0;
        bCoefficients_( 0, 6 ) = 9.0 / 35.0;
        bCoefficients_( 0, 7 ) = bCoefficients_( 0, 6 );
        bCoefficients_( 0, 8 ) = 9.0 / 280.0;
        bCoefficients_( 0, 9 ) = bCoefficients_( 0, 8 );
        bCoefficients_( 0, 10 ) = 41.0 / 840.0;

        bCoefficients_( 1, 5 ) = 34.0 / 105.0;
        bCoefficients_( 1, 6 ) = 9.0 / 35.0;
        bCoefficients_( 1, 7 ) = bCoefficients_( 1, 6 );
        bCoefficients_( 1, 8 ) = 9.0 / 280.0;
        bCoefficients_( 1, 9 ) = bCoefficients_( 1, 8 );
        bCoefficients_( 1, 11 ) = 41.0 / 840.0;
        bCoefficients_( 1, 12 ) = bCoefficients_( 1, 11 );

        break;

    default:

        std::cerr << "This is not an available coefficient set for the Runge-Kutta-Fehlberg "
                  << "integrator.";
    }
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
    RungeKuttaFehlBergVariableStepsizeIntegrator& rungeKuttaFehlBergVariableStepsizeIntegrator )
{
    // Using declarations.
    using std::endl;

    stream << "This is a VariableStepsizeIntegrator object" << endl;
    stream << "The initial state is set to: " << endl;
    stream << rungeKuttaFehlBergVariableStepsizeIntegrator.getInitialState( )->state << endl;
    stream << "The stepsize is set to: "
           << rungeKuttaFehlBergVariableStepsizeIntegrator.getStepsize( ) << endl;
    stream << "The start of the integration interval is set to: "
           << rungeKuttaFehlBergVariableStepsizeIntegrator.getIntegrationIntervalStart( ) << endl;
    stream << "The end of the integration interval is set to: "
           << rungeKuttaFehlBergVariableStepsizeIntegrator.getIntegrationIntervalEnd( ) << endl;
    stream << "The number of integration steps required is: "
           << rungeKuttaFehlBergVariableStepsizeIntegrator.getNumberOfIntegrationSteps( ) << endl;

    // Return stream.
    return stream;
}

}
// End of file.
