/*! \file rungeKutta4thOrderFixedStepsize.cpp
 *    Source file that defines the RungeKutta4thOrderFixedStepsize integrator
 *    implemented in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 6
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 29 September, 2010
 *    Last modified     : 7 February, 2011
 *
 *    References
 *
 *    Notes
 *      The integrate() function is in principle the same for all single-step
 *      integration methods, so it should be possible to inherit this function
 *      from the SingleStepIntegrationMethods class. The problem with this is
 *      however that SingletepIntegrationMethods is unable to make use of the
 *      computeNextValue() function necessary. Hence, to consolidate the code
 *      the architecture design has to be revisited in future.
 *
 *      The code is now not written such that it is possible to switch between
 *      order or stepsize type of a specific class of integrators easily. This
 *      must be incorporated in future, so that it should be equally easy to
 *      make use of a fixed-stepsize 4th-order or 5th-order Runge-Kutta
 *      integrator.
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      100929    K. Kumar          File created.
 *      100929    D. Dirkx          File checked.
 *      100929    K. Kumar          Minor comment modifications.
 *      110201    K. Kumar          Updated code to make use of State class.
 *      110203    J. Melman         Suggested two minor code updates, to make
 *                                  it more generally applicable.
 *      110204    K. Kumar          Added note about consolidating integrate();
 *                                  added a note about general RK integrator.
 *      110207    K. Kumar          Path changed; moved integrate() function
 *                                  to SingleStepIntegrationMethods.
 */

// Include statements.
#include "rungeKutta4thOrderFixedStepsize.h"

//! Default constructor.
RungeKutta4thOrderFixedStepsize::RungeKutta4thOrderFixedStepsize( )
{
}

//! Default destructor.
RungeKutta4thOrderFixedStepsize::~RungeKutta4thOrderFixedStepsize( )
{
}

//! Compute next state.
void RungeKutta4thOrderFixedStepsize::computeNextState_( const double& stepsize )
{
    // Compute next point in integration interval and set it to current point.
    integrationIntervalCurrentPoint_ += stepsize;

    // Compute k-Coefficients for 4th order, fixed stepsize,
    // Runge-Kutta integration scheme.
    // Compute k1.
    stateDerivative_ = *computeStateDerivative_(
            vectorOfCurrentStatePointers_.at( 0 ) );

    kCoefficients_[ 0 ] = stepsize * stateDerivative_.state;

    // Compute k2.
    modifiedInitialState_.state = vectorOfCurrentStatePointers_.at( 0 )->state
                                  + kCoefficients_[ 0 ] / 2.0;

    stateDerivative_ = *computeStateDerivative_( &modifiedInitialState_ );

    kCoefficients_[ 1 ] = stepsize * stateDerivative_.state;

    // Compute k3.
    modifiedInitialState_.state = vectorOfCurrentStatePointers_.at( 0 )->state
                                  + kCoefficients_[ 1 ] / 2.0;

    stateDerivative_ = *computeStateDerivative_( &modifiedInitialState_ );

    kCoefficients_[ 2 ] = stepsize * stateDerivative_.state;

    // Compute k4.
    modifiedInitialState_.state = vectorOfCurrentStatePointers_.at( 0 )->state
                                  + kCoefficients_[ 2 ];

    stateDerivative_ = *computeStateDerivative_( &modifiedInitialState_ );

    kCoefficients_[ 3 ] = stepsize * stateDerivative_.state;

    // Compute final state using 4th Order, fixed stepsize Runge-Kutta
    // algorithm.
    vectorOfCurrentStatePointers_.at( 0 ) ->state
            = vectorOfCurrentStatePointers_.at( 0 )->state
              + ( kCoefficients_[ 0 ] + 2.0 * kCoefficients_[ 1 ]
                  + 2.0 * kCoefficients_[ 2 ]
                  + kCoefficients_[ 3 ] ) / 6.0;
}

// End of file.
