/*! \file integrator.cpp
 *    Source file that defines the base class for all integration methods
 *    included in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegrators/
 *    Version           : 10
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
 *    Date created      : 22 August, 2010
 *    Last modified     : 16 May, 2011
 *
 *    References
 *
 *    Notes
 *      Currently, the code is only setup for fixed stepsize methods.
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
 *      100908    K. Kumar          File header and footer added.
 *      100926    K. Kumar          Added computeStateDerivatives( ) function and other
 *                                  modifications for pointer-to-member-function solution.
 *      100928    K. Kumar          Added missing comments.
 *      100929    D. Dirkx          File checked.
 *      100929    K. Kumar          Minor comment modifications.
 *      110201    K. Kumar          Updated code to make use of State class; implemented adaptor
 *                                  class to replace pointer-to-member functions.
 *      110203    J. Melman         Some minor comment corrections.
 *      110204    K. Kumar          Added note about stepsize and order selection.
 *      110207    K. Kumar          Path changed.
 *      110516    K. Kumar          Changed architecture so adaptor is not used. Global functions
 *                                  can no longer be used. StateDerivativeBase used to define
 *                                  classes with state derivative function.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

// Include statements.
#include "Mathematics/basicMathematicsFunctions.h"
#include "Mathematics/NumericalIntegrators/integrator.h"

//! Tudat library namespace.
namespace tudat
{

//! Set object containing state derivative.
void Integrator::setObjectContainingStateDerivative( StateDerivativeBase*
                                                     pointerToStateDerivative )
{
    // Set pointer to object containing state derivative function.
    pointerToStateDerivative_ = pointerToStateDerivative;

    // Set boolean flag to check if state derivative is set to true.
    isStateDerivativeSet_ = true;

    // Compute internal derived integration parameters.
    computeInternalDerivedIntegrationParameters_( );
}

//! Set initial state.
void Integrator::setInitialState( State* pointerToInitialState )
{
    // Set initial state.
    pointerToInitialState_ = pointerToInitialState;

    // Set boolean flag to check if initial state is set to true.
    isInitialStateSet_ = true;

    // Compute internal derived integration parameters.
    computeInternalDerivedIntegrationParameters_( );
}

//! Set initial stepsize.
void Integrator::setInitialStepsize( double initialStepsize )
{
    // Set stepsize for integrator.
    initialStepsize_ = initialStepsize;

    // Set boolean flag to check if stepsize is set to true.
    isStepsizeSet_ = true;

    // Compute internal derived integration parameters.
    computeInternalDerivedIntegrationParameters_( );
}

//! Set start of integration interval.
void Integrator::setIntegrationIntervalStart( double integrationIntervalStart )
{
    // Set start of integration interval.
    integrationIntervalStart_ = integrationIntervalStart;

    // Set boolean flag to check if start of integration interval is set to true.
    isIntegrationIntervalStartSet_ = true;

    // Compute internal derived integration parameters.
    computeInternalDerivedIntegrationParameters_( );
}

//! Set end of integration interval.
void Integrator::setIntegrationIntervalEnd( double integrationIntervalEnd )
{
    // Set end of integration interval.
    integrationIntervalEnd_ = integrationIntervalEnd;

    // Set boolean flag to check if end of integration interval is set to true.
    isIntegrationIntervalEndSet_ = true;

    // Compute internal derived integration parameters.
    computeInternalDerivedIntegrationParameters_( );
}

//! Compute state derivative.
void Integrator::computeStateDerivative_( double& integrationIntervalCurrentPoint,
                                          State* pointerToState, State* pointerToStateDerivative )
{
    // Call state derivative function using pointer to abstract base class.
   pointerToStateDerivative_->computeStateDerivative( integrationIntervalCurrentPoint,
                                                      pointerToState, pointerToStateDerivative );
}

//! Compute internal derived integration parameters.
void Integrator::computeInternalDerivedIntegrationParameters_( )
{
    // Check if internal integration parameters can be computed.
    if ( isStepsizeSet_ == true && isInitialStateSet_ == true &&
         isIntegrationIntervalStartSet_ == true && isIntegrationIntervalEndSet_ == true &&
         isStateDerivativeSet_ == true )
    {
        // Clear vector of current states.
        vectorOfCurrentStates_.clear( );

        // Set stepsize to value of initial stepsize.
        stepsize_ = initialStepsize_;

        // Compute size of integration interval.
        integrationInterval_ = integrationIntervalEnd_ - integrationIntervalStart_;

        // Compute number of integration steps for fixed stepsize methods.
        // To prevent numerical instabilities from occuring (e.g., ceil( 4 / 2 ) !=
        // ceil( ( 4 + 1e-16 ) / 2 )), the square root of the machine precision is
        // subtracted before applying the ceil function.
        numberOfIntegrationSteps_ = static_cast < unsigned int >(
                    std::ceil( integrationInterval_ / stepsize_
                               - std::sqrt( mathematics::MACHINE_PRECISION_DOUBLES ) ) );

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
    }
}

}

// End of file.
