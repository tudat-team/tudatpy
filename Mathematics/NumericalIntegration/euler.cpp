/*! \file euler.cpp
 *    Source file that defines the Euler integrator implemented in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegration/
 *    Version           : 4
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dikx@student.tudelft.nl
 *
 *    Date created      : 30 July, 2010
 *    Last modified     : 29 September, 2010
 *
 *    References
 *
 *    Notes
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
 *      YYMMDD    author              comment
 *      100909    K. Kumar            File header and footer added.
 *      100928    K. Kumar            Added fixed output interval option,
 *                                    added missing comments.
 *      100929    D. Dirkx            File checked, integration loop modified,
 *                                    as it was integrating one step extra.
 *      100929    K. Kumar            Minor comment modifications, modified
 *                                    integration loop to avoid problems with
 *                                    unsigned int values equally 0.
 */

// Include statments.
#include "euler.h"

//! Default constructor.
Euler::Euler( )
{
}

//! Default destructor.
Euler::~Euler( )
{
}

//! Perform Euler integration.
void Euler::integrate( )
{
    // Set current point in integration interval to start of
    // integration interval.
    integrationIntervalCurrentPoint_ = integrationIntervalStart_;

    // Loop over the number of integration steps to compute final
    // state vector at each integration step.
    for ( unsigned int i = 1; i < numberOfIntegrationSteps_; i++)
    {
        // Set initial state vector to final state vector before proceeding to
        // next step.
        finalStateVector_ = vectorOfCurrentStateVectors_.at( 0 );

        // Check if integration history should be saved.
        if ( isSaveIntegrationHistory_ )
        {
            // Store intermediate results.
            storeIntermediateIntegrationResult_( );
        }

        // Compute final state vector.
        computeNextStateVector( stepsize_ );
    }

    // Compute stepsize of last step.
    lastStepStepsize_ = integrationIntervalEnd_
                        - integrationIntervalCurrentPoint_;

    // Set initial state vector to final state vector before proceeding to
    // next step.
    finalStateVector_ = vectorOfCurrentStateVectors_.at( 0 );

    // Check if integration history should be saved.
    if ( isSaveIntegrationHistory_ )
    {
        // Store intermediate results.
        storeIntermediateIntegrationResult_( );
    }

    // Compute final state vector.
    computeNextStateVector( lastStepStepsize_ );

    // Set initial state vector to final state vector before proceeding to
    // next step.
    finalStateVector_ = vectorOfCurrentStateVectors_.at( 0 );

    // Check if integration history should be saved.
    if ( isSaveIntegrationHistory_ )
    {
        // Store intermediate results.
        storeIntermediateIntegrationResult_( );
    }
}

//! Compute next state vector.
void Euler::computeNextStateVector( const double& stepsize )
{
    // Compute next point in integration interval and set it to current point.
    integrationIntervalCurrentPoint_ += stepsize;

    // Compute state derivative given initial state vector.
    computeStateDerivatives_( vectorOfCurrentStateVectors_.at( 0 ),
                                       stateDerivativeVector_ );

    // Compute next state vector using Euler algorithm by updating current
    // state vector stored in vector container.
    // PLACEHOLDER
    vectorOfCurrentStateVectors_.at( 0 ) = vectorOfCurrentStateVectors_.at( 0 )
                                           + stepsize * stateDerivativeVector_;
}

// End of file.
