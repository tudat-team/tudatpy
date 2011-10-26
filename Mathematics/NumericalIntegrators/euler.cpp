/*! \file euler.cpp
 *    Source file that defines the Euler integrator implemented in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 9
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
 *    Date created      : 30 July, 2010
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
 *      100909    K. Kumar          File header and footer added.
 *      100928    K. Kumar          Added fixed output interval option, added missing comments.
 *      100929    D. Dirkx          File checked, integration loop modified, as it was integrating
 *                                  one step extra.
 *      100929    K. Kumar          Minor comment modifications, modified integration loop to avoid
 *                                  problems with unsigned int values equally 0.
 *      110201    K. Kumar          Updated code to make use of State class.
 *      110203    J. Melman         File checked.
 *      110204    K. Kumar          Added note about consolidating integrate().
 *      110207    K. Kumar          Path changed; moved integrate() function
 *                                  to SingleStepIntegrationMethods.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

// Include statements.
#include "Mathematics/NumericalIntegrators/euler.h"

//! Tudat library namespace.
namespace tudat
{

//! Compute next state.
void Euler::computeNextState_( const double& stepsize )
{
    // Compute state derivative initial given state and initial time.
    computeStateDerivative_( integrationIntervalCurrentPoint_, &vectorOfCurrentStates_.at( 0 ),
                             &stateDerivative_ );

    // Compute next state using Euler algorithm by updating current state
    // stored in vector container.
    vectorOfCurrentStates_.at( 0 ).state = vectorOfCurrentStates_.at( 0 ).state
            + stepsize * stateDerivative_.state;
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream, Euler& eulerIntegrator )
{
    stream << "This is an Euler object" << std::endl;
    stream << "The initial state is set to: " << std::endl;
    stream << eulerIntegrator.getInitialState( )->state << std::endl;
    stream << "The stepsize is set to: " << eulerIntegrator.getStepsize( ) << std::endl;
    stream << "The start of the integration interval is set to: "
           << eulerIntegrator.getIntegrationIntervalStart( ) << std::endl;
    stream << "The end of the integration interval is set to: "
           << eulerIntegrator.getIntegrationIntervalEnd( ) << std::endl;
    stream << "The number of integration steps required is: "
           << eulerIntegrator.getNumberOfIntegrationSteps( ) << std::endl;

    // Return stream.
    return stream;
}

}

// End of file.
