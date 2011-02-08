/*! \file integrator.cpp
 *    Source file that defines the base class for all integration methods
 *    included in Tudat.
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
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 22 August, 2010
 *    Last modified     : 7 February, 2011
 *
 *    References
 *
 *    Notes
 *      Currently, the code is only setup for fixed stepsize methods. In
 *      future, templates should be employed to enable selection of the
 *      stepsize and order types.
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
 *      100908    K. Kumar          File header and footer added.
 *      100926    K. Kumar          Added computeStateDerivatives( ) function
 *                                  and other modifications for pointer-to-
 *                                  member-function solution.
 *      100928    K. Kumar          Added missing comments.
 *      100929    D. Dirkx          File checked.
 *      100929    K. Kumar          Minor comment modifications.
 *      110201    K. Kumar          Updated code to make use of State class;
 *                                  implemented adaptor class to replace
 *                                  pointer-to-member functions.
 *      110203    J. Melman         Some minor comment corrections.
 *      110204    K. Kumar          Added note about stepsize and order
 *                                  selection.
 *      110207    K. Kumar          Path changed.
 */

// Include statements.
#include "integrator.h"

//! Default constructor.
Integrator::Integrator( ) : isIntegrationHistoryToBeSaved_( false ),
                            dimensionOfCurrentState_( -0 ),
                            pointerToStateTakingFunction_( NULL ),
                            pointerToIntegratorBase_( NULL ),
                            isStepsizeSet_ ( false ),
                            isInitialStateSet_( false ),
                            isIntegrationIntervalStartSet_( false ),
                            isIntegrationIntervalEndSet_( false ),
                            isStateDerivativeFunctionSet_( false )
{
}

//! Default destructor.
Integrator::~Integrator( )
{
}

//! Set function containing state derivative.
void Integrator::setStateDerivativeFunction(
        pointerToStateTakingFunction derivativeFunction )
{
    // Set pointer to state derivative function.
    pointerToStateTakingFunction_ =  derivativeFunction;

    // Set boolean flag to check if state derivative function is set, to true.
    isStateDerivativeFunctionSet_ = true;

    // Compute internal derived integration parameters.
    computeInternalDerivedIntegrationParameters_( );
}

//! Set adaptor class for Integrator.
void Integrator::setIntegratorAdaptor( IntegratorBase*
                                       pointerToIntegratorBase )
{
    // Set pointer to integrator abstract base class.
    pointerToIntegratorBase_ = pointerToIntegratorBase;

    // Set boolean flag to check if state derivative function is set, to true.
    isStateDerivativeFunctionSet_ = true;

    // Compute internal derived integration parameters.
    computeInternalDerivedIntegrationParameters_( );
}

//! Set initial state.
void Integrator::setInitialState( State* pointerToInitialState )
{
    // Set initial state.
    initialState_ = *pointerToInitialState;

    // Add pointer to initial state to container vector.
    vectorOfCurrentStatePointers_.push_back( pointerToInitialState );

    // Set boolean flag to check if initial state is set to true.
    isInitialStateSet_ = true;

    // Compute internal derived integration parameters.
    computeInternalDerivedIntegrationParameters_( );
}

//! Set initial stepsize.
void Integrator::setInitialStepsize( const double& stepsize )
{
    // Set stepsize for integrator.
    initialStepsize_ = stepsize;

    // Set boolean flag to check if stepsize is set to true.
    isStepsizeSet_ = true;

    // Compute internal derived integration parameters.
    computeInternalDerivedIntegrationParameters_( );
}

//! Set start of integration interval.
void Integrator::setIntegrationIntervalStart( const double&
                                              integrationIntervalStart )
{
    // Set start of integration interval.
    integrationIntervalStart_ = integrationIntervalStart;

    // Set boolean flag to check if start of integration interval is set to
    // true.
    isIntegrationIntervalStartSet_ = true;

    // Compute internal derived integration parameters.
    computeInternalDerivedIntegrationParameters_( );
}

//! Set end of integration interval.
void Integrator::setIntegrationIntervalEnd( const double&
                                            integrationIntervalEnd )
{
    // Set end of integration interval.
    integrationIntervalEnd_ = integrationIntervalEnd;

    // Set boolean flag to check if end of integration interval is set to true.
    isIntegrationIntervalEndSet_ = true;

    // Compute internal derived integration parameters.
    computeInternalDerivedIntegrationParameters_( );
}

//! Sets whether integration history should be saved.
void Integrator::setIsIntegrationHistoryToBeSaved(
        const bool& isIntegrationHistoryToBeSaved )
{
    // Set flag to check if integration history should be saved.
    isIntegrationHistoryToBeSaved_ = isIntegrationHistoryToBeSaved;
}


//! Get stepsize.
double Integrator::getStepsize( )
{
    // Return stepsize for integrator.
    return stepsize_;
}

//! Get number of integration steps.
unsigned int Integrator::getNumberOfIntegrationSteps( )
{
    // Return number of integration steps.
    return numberOfIntegrationSteps_;
}

//! Get start of integration interval.
double Integrator::getIntegrationIntervalStart( )
{
    // Return start of integration interval.
    return integrationIntervalStart_;
}

//! Get end of integration interval.
double Integrator::getIntegrationIntervalEnd( )
{
    // Return end of integration interval.
    return integrationIntervalEnd_;
}

//! Get pointer to initial state.
State* Integrator::getInitialState( )
{
    return &initialState_;
}

//! Get pointer to final state.
State* Integrator::getFinalState( )
{
    return &finalState_;
}

//! Get integration history.
std::map< double, State* >& Integrator::getIntegrationHistory( )
{
    // Return map integration history.
    return integrationHistory_;
}

//! Compute state derivative.
State* Integrator::computeStateDerivative_( State* pointerToState )
{
    // Check if pointer to state-taking function is set.
    if ( pointerToStateTakingFunction_ )
    {
        // Call state derivative function using pointer.
        stateDerivative_ = *pointerToStateTakingFunction_( pointerToState );
    }

    // Else check if pointer to abstract base class with derivative function is
    // set.
    else if ( pointerToIntegratorBase_ )
    {
        // Call state derivative function using pointer to abstract base class.
        stateDerivative_ = *pointerToIntegratorBase_
                           ->computeStateDerivative( pointerToState );
    }

    // Return pointer to state derivative.
    return &stateDerivative_;
}

//! Store intermediate integration result.
void Integrator::storeIntermediateIntegrationResult_( )
{
    // Pointer to intermediate state.
    State* ponterToIntermediateState_ = new State;

    // Set intermediate state to current state.
    ponterToIntermediateState_->state
            = vectorOfCurrentStatePointers_.at( 0 )->state;

    // Store intermediate integration results in integration history map.
    integrationHistory_[ integrationIntervalCurrentPoint_ ]
            = ponterToIntermediateState_;
}

//! Compute internal derived integration parameters.
void Integrator::computeInternalDerivedIntegrationParameters_( )
{
    // Check if internal integration parameters can be computed.
    if ( isStepsizeSet_ == true &&
         isInitialStateSet_ == true &&
         isIntegrationIntervalStartSet_ == true &&
         isIntegrationIntervalEndSet_ == true &&
         isStateDerivativeFunctionSet_ == true )
    {
        // Set stepsize to value of initial stepsize.
        stepsize_ = initialStepsize_;

        // Compute size of integration interval.
        integrationInterval_ = integrationIntervalEnd_
                               - integrationIntervalStart_;

        // Compute number of integration steps for fixed stepsize methods.
        numberOfIntegrationSteps_ = std::ceil( integrationInterval_
                                               / stepsize_ );

        // Set dimension of initial state.
        dimensionOfCurrentState_ = vectorOfCurrentStatePointers_.at( 0 )
                                   ->state.rows( );

        // Set dimension of final state to correspond with dimension of initial
        // state.
        finalState_.state.setZero( dimensionOfCurrentState_ );

        // Set dimension of state derivative to correpsond with dimension of
        // initial state.
        stateDerivative_.state.setZero( dimensionOfCurrentState_ );
    }
}

// End of file.
