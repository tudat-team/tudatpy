/*! \file singleStepIntegrationMethods.h
 *    Header file that defines the base class for all single step integration
 *    methods included in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 1
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
 *    Date created      : 7 February, 2011
 *    Last modified     : 7 February, 2011
 *
 *    References
 *
 *    Notes
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
 *      110207    K. Kumar          First creation of file.
 */

// Include statements.
#include "singleStepIntegrationMethods.h"

//! Default constructor.
SingleStepIntegrationMethods::SingleStepIntegrationMethods( )
{
}

//! Default destructor.
SingleStepIntegrationMethods::~SingleStepIntegrationMethods( )
{
}

//! Integrate.
void SingleStepIntegrationMethods::integrate( )
{
    // Set current point in integration interval to start of
    // integration interval.
    integrationIntervalCurrentPoint_ = integrationIntervalStart_;

    // Loop over the number of integration steps to compute final
    // state at each integration step.
    for ( unsigned int i = 1; i < numberOfIntegrationSteps_; i++ )
    {
        // Compute final state.
        computeNextState_( stepsize_ );

        // Move current point in integration interval to next.
        integrationIntervalCurrentPoint_ += stepsize_;
    }

    // Compute stepsize of last step.
    lastStepStepsize_ = integrationIntervalEnd_
                        - integrationIntervalCurrentPoint_;

    // Compute final state.
    computeNextState_( lastStepStepsize_ );

    // Set final state.
    finalState_ = vectorOfCurrentStates_.at( 0 );
}

// End of file.
