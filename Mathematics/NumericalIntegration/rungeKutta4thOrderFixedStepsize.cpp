/*! \file rungeKutta4thOrderFixedStepsize.cpp
 *    Source file that defines the RungeKutta4thOrderFixedStepsize integrator implemented in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegration/
 *    Version           : 3
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
 *    Date created      : 29 September, 2010
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
 *      100929    K. Kumar            File created.
 *      100929    D. Dirkx            File checked.
 *      100929    K. Kumar            Minor comment modifications.
 */

// Include statments.
#include "rungeKutta4thOrderFixedStepsize.h"

//! Default constructor.
RungeKutta4thOrderFixedStepsize::RungeKutta4thOrderFixedStepsize( )
{
}

//! Default destructor.
RungeKutta4thOrderFixedStepsize::~RungeKutta4thOrderFixedStepsize( )
{
}

//! Perform 4th-order, fixed stepsize, Runge-Kutta integration.
void RungeKutta4thOrderFixedStepsize::integrate( )
{
    // Set modified state vector to zero.
    modifiedInitialStateVector_.setZero( initialStateVector_.rows( ) );

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
void RungeKutta4thOrderFixedStepsize::computeNextStateVector( double& stepsize )
{
    // Compute next point in integration interval and set it to current point.
    integrationIntervalCurrentPoint_ += stepsize;

    // Compute k-Coefficients for 4th Order, fixed stepsize,
    // Runge-Kutta integration scheme.
    // Compute k1.
    computeStateDerivatives_( vectorOfCurrentStateVectors_.at( 0 ),
                                       stateDerivativeVector_ );
    kCoefficients_[ 0 ] = stepsize * stateDerivativeVector_;

    // Compute k2.
    modifiedInitialStateVector_ = vectorOfCurrentStateVectors_.at( 0 )
                                  + kCoefficients_[ 0 ] / 2.0;
    computeStateDerivatives_( modifiedInitialStateVector_,
                              stateDerivativeVector_ );
    kCoefficients_[ 1 ] = stepsize * stateDerivativeVector_;

    // Compute k3.
    modifiedInitialStateVector_ = vectorOfCurrentStateVectors_.at( 0 )
                                  + kCoefficients_[ 1 ] / 2.0;
    computeStateDerivatives_( modifiedInitialStateVector_,
                              stateDerivativeVector_ );
    kCoefficients_[ 2 ] = stepsize * stateDerivativeVector_;

    // Compute k4.
    modifiedInitialStateVector_ = vectorOfCurrentStateVectors_.at( 0 )
                                  + kCoefficients_[ 2 ];
    computeStateDerivatives_( modifiedInitialStateVector_,
                              stateDerivativeVector_ );
    kCoefficients_[ 3 ] = stepsize * stateDerivativeVector_;

    // Compute final state vector using 4th Order, fixed stepsize
    // Runge-Kutta algorithm
    vectorOfCurrentStateVectors_.at( 0 )
            = vectorOfCurrentStateVectors_.at( 0 )
              + ( kCoefficients_[ 0 ] + 2.0 * kCoefficients_[ 1 ]
                  + 2.0 * kCoefficients_[ 2 ]
                  + kCoefficients_[ 3 ] ) / 6.0;
}

// End of file.
