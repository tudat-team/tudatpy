/*! \file euler.cpp
 *    Source file that defines the Euler integrator implemented in Tudat.
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
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 30 July, 2010
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
 *      100909    K. Kumar          File header and footer added.
 *      100928    K. Kumar          Added fixed output interval option, added
 *                                  missing comments.
 *      100929    D. Dirkx          File checked, integration loop modified, as
 *                                  it was integrating one step extra.
 *      100929    K. Kumar          Minor comment modifications, modified
 *                                  integration loop to avoid problems with
 *                                  unsigned int values equally 0.
 *      110201    K. Kumar          Updated code to make use of State class.
 *      110203    J. Melman         File checked.
 *      110204    K. Kumar          Added note about consolidating integrate().
 *      110207    K. Kumar          Path changed; moved integrate() function
 *                                  to SingleStepIntegrationMethods.
 */

// Include statements.
#include "euler.h"

// Using declarations.
using std::endl;

//! Default constructor.
Euler::Euler( )
{
}

//! Default destructor.
Euler::~Euler( )
{
}

//! Compute next state.
void Euler::computeNextState_( const double& stepsize )
{
    // Compute next point in integration interval and set it to current point.
    integrationIntervalCurrentPoint_ += stepsize;

    // Compute state derivative given initial state.
    stateDerivative_ = *computeStateDerivative_(
            vectorOfCurrentStatePointers_.at( 0 ) );

    // Compute next state using Euler algorithm by updating current state
    // stored in vector container.
    vectorOfCurrentStatePointers_.at( 0 )->state
            = vectorOfCurrentStatePointers_.at( 0 )->state
              + stepsize * stateDerivative_.state;
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream, Euler& euler )
{
    stream << "This is an Euler object" << endl;
    stream << "The initial state is set to: " << endl;
    stream << euler.getInitialState( )->state << endl;
    stream << "The stepsize is set to: "
           << euler.getStepsize( ) << endl;
    stream << "The start of the integration interval is set to: "
           << euler.getIntegrationIntervalStart( ) << endl;
    stream << "The end of the integration interval is set to: "
           << euler.getIntegrationIntervalEnd( ) << endl;
    stream << "The number of integration steps required is: "
           << euler.getNumberOfIntegrationSteps( ) << endl;

    // Return stream.
    return stream;
}

// End of file.
