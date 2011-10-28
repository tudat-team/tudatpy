/*! \file rungeKuttaFehlberg78VariableStepsize.cpp
 *    Source file that defines the 7(8)th-order, variable stepsize, Runge-Kutta
 *    integrator.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : E.A.G. Heeren
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : E.A.G.Heeren@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl

 *    Checker           : E.A.G. Heeren
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : E.A.G.Heeren@student.tudelft.nl
 *
 *    Date created      : 26 August, 2011
 *    Last modified     : 12 September, 2011
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
 *      110909    E.A.G. Heeren     Modified to be included in Tudat. Added if-statement
 *                                  to ensure a maximum stepsize change.
 *      110912    K. Kumar          Minor changes.
 */

// Include statements.
#include <cmath>
#include <iostream>
#include "Mathematics/basicMathematicsFunctions.h"
#include "Mathematics/NumericalIntegrators/rungeKuttaFehlberg78VariableStepsize.h"

//! Tudat library namespace.
namespace tudat
{

//! Integrate.
void RungeKuttaFehlberg78VariableStepsize::integrate( )
{
    // Using declarations.
    using std::fabs;
    using std::cerr;
    using std::endl;

    // Clear vector of current states.
    vectorOfCurrentStates_.clear( );

    // Set stepsize to value of initial stepsize.
    stepsize_ = initialStepsize_;

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

    // Set dimensions of f-Matrix.
    fMatrix.setZero( dimensionOfState_, 13 );

    // Loop until end of integration is reached.
    while ( true )
    {
        // Current state and current time.
        currentState_ = vectorOfCurrentStates_.at( 0 );
        currentTime_ = integrationIntervalCurrentPoint_;

        // Check for last step.
        if ( fabs( stepsize_ )
             > fabs( integrationIntervalEnd_ - integrationIntervalCurrentPoint_ ) )
        {
            stepsize_ = integrationIntervalEnd_ - integrationIntervalCurrentPoint_;
        }

        // Check for end of integration period and break from loop if end is reached.
        if ( fabs( integrationIntervalCurrentPoint_ - integrationIntervalEnd_ )
             < mathematics::MACHINE_PRECISION_DOUBLES )
        {
            break;
        }

        // Compute state derivative.
        computeStateDerivative_( integrationIntervalCurrentPoint_, &vectorOfCurrentStates_.at( 0 ),
                                 &stateDerivative_ );

        // Copy computed state derivative to first column of fMatrix.
        fMatrix.col( 0 ) = stateDerivative_.state;

        // Compute current 7th-order estimate by looping over k-coefficients.
        for ( unsigned int k = 0; k < 13; k++ )
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
                                     &vectorOfCurrentStates_.at( 0 ), &stateDerivative_ );

            // Copy computed state derivative to kth-column of fMatrix.
            fMatrix.col( k ) = stateDerivative_.state.col( 0 );
        }

        // Set error in state to relative error tolerance.
        errorInState_ = relativeErrorTolerance_;

        // Loop over state elements to compute relative error in 7th-order estimate
        // using 8th-order estimate.
        for ( int i = 0; i < vectorOfCurrentStates_.at( 0 ).state.rows( ); i++ )
        {
            // Compute ith element of state.
            vectorOfCurrentStates_.at( 0 ).state[ i ] = currentState_.state[ i ] + stepsize_
                    * ( bCoefficients_.array( ) * fMatrix.row( i ).array( ) ).sum( );

            // Compute truncation error.
            truncationError_ = fabs( ( fMatrix( i, 0 ) + fMatrix( i, 10 )
                                       - fMatrix( i, 11 ) - fMatrix( i, 12 ) )
                                     * bCoefficients_( 11 ) * stepsize_ );

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
        stepsize_ = 0.8 * stepsize_ * pow( 1.0 / errorInState_, 1.0 / 8.0 );

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
            cerr << "Minimum stepsize condition is violated." << endl;
            cerr << "RKF78 integrator failed!" << endl;
        }

        // Check if computed error in state is too large and reject step if true.
        if ( errorInState_ > 1.0 )
        {
            // Reject current step.
            integrationIntervalCurrentPoint_ = currentTime_;
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
void RungeKuttaFehlberg78VariableStepsize::setCoefficients_( )
{
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

    bCoefficients_( 5 ) = 34.0 / 105.0;
    bCoefficients_( 6 ) = 9.0 / 35.0;
    bCoefficients_( 7 ) = bCoefficients_( 6 );
    bCoefficients_( 8 ) = 9.0 / 280.0;
    bCoefficients_( 9 ) = bCoefficients_( 8 );
    bCoefficients_( 11 ) = 41.0 / 840.0;
    bCoefficients_( 12 ) = bCoefficients_( 11 );
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
    RungeKuttaFehlberg78VariableStepsize& rungeKuttaFehlberg78VariableStepsize )
{
    // Using declarations.
    using std::endl;

    stream << "This is a RungeKuttaFehlberg78VariableStepsize object" << endl;
    stream << "The initial state is set to: " << endl;
    stream << rungeKuttaFehlberg78VariableStepsize.getInitialState( )->state << endl;
    stream << "The stepsize is set to: "
           << rungeKuttaFehlberg78VariableStepsize.getStepsize( ) << endl;
    stream << "The start of the integration interval is set to: "
           << rungeKuttaFehlberg78VariableStepsize.getIntegrationIntervalStart( ) << endl;
    stream << "The end of the integration interval is set to: "
           << rungeKuttaFehlberg78VariableStepsize.getIntegrationIntervalEnd( ) << endl;
    stream << "The number of integration steps required is: "
           << rungeKuttaFehlberg78VariableStepsize.getNumberOfIntegrationSteps( ) << endl;

    // Return stream.
    return stream;
}

}

// End of file.
