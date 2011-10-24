/*! \file rungeKutta4thOrderFixedStepsize.cpp
 *    Source file that defines the RungeKutta4thOrderFixedStepsize integrator implemented in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 8
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 29 September, 2010
 *    Last modified     : 5 September, 2011
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
 *      100929    K. Kumar          File created.
 *      100929    D. Dirkx          File checked.
 *      100929    K. Kumar          Minor comment modifications.
 *      110201    K. Kumar          Updated code to make use of State class.
 *      110203    J. Melman         Suggested two minor code updates, to make it more generally
 *                                  applicable.
 *      110204    K. Kumar          Added note about consolidating integrate();
 *                                  added a note about general RK integrator.
 *      110207    K. Kumar          Path changed; moved integrate() function
 *                                  to SingleStepIntegrationMethods.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

// Include statements.
#include "Mathematics/NumericalIntegrators/rungeKutta4thOrderFixedStepsize.h"

//! Compute next state.
void RungeKutta4thOrderFixedStepsize::computeNextState_( const double& stepsize )
{
    // Clear k-Coefficients.
    kCoefficients_.clear( );

    // Set size of state derivative vector.
    stateDerivative_.state.setZero( vectorOfCurrentStates_.at( 0 ).state.rows( ) );

    // Compute k-Coefficients for 4th order, fixed stepsize,
    // Runge-Kutta integration scheme.
    // Compute k1.
    computeStateDerivative_( integrationIntervalCurrentPoint_, &vectorOfCurrentStates_.at( 0 ),
                             &stateDerivative_ );

    kCoefficients_.push_back( stepsize * stateDerivative_.state );

    // Compute k2.
    modifiedInitialState_.state = vectorOfCurrentStates_.at( 0 ).state
            + kCoefficients_.at( 0 ) / 2.0;

    modifiedIntegrationIntervalCurrentPoint_ = integrationIntervalCurrentPoint_ + stepsize / 2.0;

    computeStateDerivative_( modifiedIntegrationIntervalCurrentPoint_, &modifiedInitialState_,
                             &stateDerivative_ );

    kCoefficients_.push_back( stepsize * stateDerivative_.state );

    // Compute k3.
    modifiedInitialState_.state = vectorOfCurrentStates_.at( 0 ).state
            + kCoefficients_.at( 1 ) / 2.0;

    modifiedIntegrationIntervalCurrentPoint_ = integrationIntervalCurrentPoint_ + stepsize / 2.0;

    computeStateDerivative_( modifiedIntegrationIntervalCurrentPoint_, &modifiedInitialState_,
                             &stateDerivative_ );

    kCoefficients_.push_back( stepsize * stateDerivative_.state );

    // Compute k4.
    modifiedInitialState_.state = vectorOfCurrentStates_.at( 0 ).state + kCoefficients_.at( 2 );

    modifiedIntegrationIntervalCurrentPoint_ = integrationIntervalCurrentPoint_ + stepsize;

    computeStateDerivative_( modifiedIntegrationIntervalCurrentPoint_, &modifiedInitialState_,
                             &stateDerivative_ );

    kCoefficients_.push_back( stepsize * stateDerivative_.state );

    // Compute final state using 4th Order, fixed stepsize Runge-Kutta
    // algorithm.
    vectorOfCurrentStates_.at( 0 ).state = vectorOfCurrentStates_.at( 0 ).state
            + ( kCoefficients_.at( 0 ) + 2.0 * kCoefficients_.at( 1 )
                + 2.0 * kCoefficients_.at( 2 ) + kCoefficients_.at( 3 ) ) / 6.0;
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                          RungeKutta4thOrderFixedStepsize& rungeKutta4thOrderFixedStepsize )
{
    stream << "This is a RungeKutta4thOrderFixedStepsize object" << std::endl;
    stream << "The initial state is set to: " << std::endl;
    stream << rungeKutta4thOrderFixedStepsize.getInitialState( )->state << std::endl;
    stream << "The stepsize is set to: "
           << rungeKutta4thOrderFixedStepsize.getStepsize( ) << std::endl;
    stream << "The start of the integration interval is set to: "
           << rungeKutta4thOrderFixedStepsize.getIntegrationIntervalStart( ) << std::endl;
    stream << "The end of the integration interval is set to: "
           << rungeKutta4thOrderFixedStepsize.getIntegrationIntervalEnd( ) << std::endl;
    stream << "The number of integration steps required is: "
           << rungeKutta4thOrderFixedStepsize.getNumberOfIntegrationSteps( ) << std::endl;

    // Return stream.
    return stream;
}

// End of file.
