/*! \file integrator.cpp
 *    Source file that defines the base class for all integration methods
 *    included in Tudat.
 *
 *    Path              : /Mathematics/NumericalIntegration/
 *    Version           : 6
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
 *    Date created      : 22 August, 2010
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
 *      100908    K. Kumar            File header and footer added.
 *      100926    K. Kumar            Added computeStateDerivatives( ) function
 *                                    and other modifications for pointer-to-
 *                                    member-function solution.
 *      100928    K. Kumar            Added missing comments.
 *      100929    D. Dirkx            File checked.
 *      100929    K. Kumar            Minor comment modifications.
 */

// Include statements.
#include "integrator.h"

//! Integrator class constructor.
Integrator::Integrator( ) : pointerToVectorVectorTakingFunction_( NULL ),
pointerToStateDerivativeFunction_( NULL ),
pointerToNumericalPropagator_( NULL )
{
    // Set boolean flags to check if necessary integration parameters
    // are set to false.
    isStepsizeSet_ = false;
    isInitialStateVectorSet_  = false;
    isIntegrationIntervalStartSet_ = false;
    isIntegrationIntervalEndSet_ = false;
    isStateDerivativeFunctionSet_ = false;

    // Set boolean flag to indicate if integration history should be shaved to
    // false.
    isSaveIntegrationHistory_ = false;
}

//! Integrator class destructor.
Integrator::~Integrator( )
{
}

//! Set address of function containing state derivatives.
void Integrator::setStateDerivativeFunction(
        pointerToVectorVectorTakingFunction
        pointerToDerivativeFunction )
{
    // Set pointer to state derivative function.
    pointerToVectorVectorTakingFunction_ =  pointerToDerivativeFunction;
    pointerToStateDerivativeFunction_ = NULL;
    pointerToNumericalPropagator_ = NULL;

    // Set boolean flag to check if state derivative function is set to true.
    isStateDerivativeFunctionSet_ = true;

    // Check if internal integration parameters can be computed.
    if ( isStepsizeSet_ == true &&
         isInitialStateVectorSet_ == true &&
         isIntegrationIntervalStartSet_ == true &&
         isIntegrationIntervalEndSet_ == true &&
         isStateDerivativeFunctionSet_ == true )
    {
        // Compute internally derived integration parameters
        computeInternalDerivedIntegrationParameters_( );
    }
}

//! Set address of function containing state derivatives.
void Integrator::setStateDerivativeFunction( pointerToStateDerivativeFunction
                                 derivativeFunction,
                                 NumericalPropagator*
                                 pointerToNumericalPropagator )
{
    // Set pointer to state derivative function.
    pointerToVectorVectorTakingFunction_ =  NULL;
    pointerToStateDerivativeFunction_ = derivativeFunction;
    pointerToNumericalPropagator_ = pointerToNumericalPropagator;

    // Set boolean flag to check if state derivative function is set to true.
    isStateDerivativeFunctionSet_ = true;

    // Check if internal integration parameters can be computed.
    if ( isStepsizeSet_ == true &&
         isInitialStateVectorSet_ == true &&
         isIntegrationIntervalStartSet_ == true &&
         isIntegrationIntervalEndSet_ == true &&
         isStateDerivativeFunctionSet_ == true )
    {
        // Compute internally derived integration parameters
        computeInternalDerivedIntegrationParameters_( );
    }

}

//! Set initial state vector.
void Integrator::setInitialStateVector( Vector& initialStateVector )
{
    // Set initial state vector.
    initialStateVector_ = initialStateVector;

    // Add initial state vector to container vector.
    vectorOfCurrentStateVectors_.push_back( initialStateVector );

    // Set boolean flag to check if initial state vector is set to true.
    isInitialStateVectorSet_ = true;

    // Check if internal integration parameters can be computed.
    if ( isStepsizeSet_ == true &&
         isInitialStateVectorSet_ == true &&
         isIntegrationIntervalStartSet_ == true &&
         isIntegrationIntervalEndSet_ == true &&
         isStateDerivativeFunctionSet_ == true )
    {
        // Compute internally derived integration parameters
        computeInternalDerivedIntegrationParameters_( );
    }
}

//! Set value of initial stepsize.
void Integrator::setInitialStepsize( const double& stepsize )
{
    // Set stepsize for integrator.
    initialStepsize_ = stepsize;

    // Set boolean flag to check if stepsize is set to true.
    isStepsizeSet_ = true;

    // Check if internal integration parameters can be computed.
    if ( isStepsizeSet_ == true &&
         isInitialStateVectorSet_ == true &&
         isIntegrationIntervalStartSet_ == true &&
         isIntegrationIntervalEndSet_ == true &&
         isStateDerivativeFunctionSet_ == true )
    {
        // Compute internally derived integration parameters
        computeInternalDerivedIntegrationParameters_( );
    }
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

    // Check if internal integration parameters can be computed.
    if ( isStepsizeSet_ == true &&
         isInitialStateVectorSet_ == true &&
         isIntegrationIntervalStartSet_ == true &&
         isIntegrationIntervalEndSet_ == true &&
         isStateDerivativeFunctionSet_ == true )
    {
        // Compute internally derived integration parameters
        computeInternalDerivedIntegrationParameters_( );
    }
}

//! Set end of integration interval.
void Integrator::setIntegrationIntervalEnd( const double&
                                            integrationIntervalEnd )
{
    // Set end of integration interval.
    integrationIntervalEnd_ = integrationIntervalEnd;

    // Set boolean flag to check if end of integration interval is set to true.
    isIntegrationIntervalEndSet_ = true;

    // Check if internal integration parameters can be computed.
    if ( isStepsizeSet_ == true &&
         isInitialStateVectorSet_ == true &&
         isIntegrationIntervalStartSet_ == true &&
         isIntegrationIntervalEndSet_ == true &&
         isStateDerivativeFunctionSet_ == true )
    {
        // Compute internally derived integration parameters
        computeInternalDerivedIntegrationParameters_( );
    }
}

//! Save integration history.
void Integrator::saveIntegrationHistory( const bool& isSaveIntegrationHistory )
{
    // Set flag to check if integration history should be saved.
    isSaveIntegrationHistory_ = isSaveIntegrationHistory;
}


//! Get stepsize
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

//! Get address of initial state vector.
Vector& Integrator::getInitialStateVector( )
{
    // Return initial state vector.
    return initialStateVector_;
}

//! Get address of final state vector.
Vector& Integrator::getFinalStateVector( )
{
    // Return final state vector.
    return finalStateVector_;
}

//! Get integration history.
std::map < double, Vector >& Integrator::getIntegrationHistory( )
{
    // Return map integration history.
    return integrationHistory_;
}

//! Compute state derivatives.
void Integrator::computeStateDerivatives_( Vector& stateVector,
                                           Vector& stateDerivativeVector )
{
    // Check if pointer to vector-vector-taking function is set.
    if( pointerToVectorVectorTakingFunction_ )
    {
        // Call state derivative function using pointer.
        pointerToVectorVectorTakingFunction_( stateVector,
                                              stateDerivativeVector );
    }

    // Check if pointer to state derivative function and pointer to
    // Propagator class are set.
    else if( ( pointerToStateDerivativeFunction_ )
             && ( pointerToNumericalPropagator_ ) )
    {
        // Call state derivative function in object of Propagator
        // class.
        ( pointerToNumericalPropagator_
                ->*pointerToStateDerivativeFunction_ )( stateVector,
                                                 stateDerivativeVector );
    }
}

//! Store intermediate integration result.
void Integrator::storeIntermediateIntegrationResult_( )
{
    // Store intermediate integration results in integration history map.
    integrationHistory_[ integrationIntervalCurrentPoint_ ]
            = finalStateVector_;
}

//! Compute internal derived integration parameters.
void Integrator::computeInternalDerivedIntegrationParameters_( )
{
    // Set stepsize to value of initial stepsize.
    stepsize_ = initialStepsize_;

    // Compute size of integration interval.
    integrationInterval_ = integrationIntervalEnd_
                           - integrationIntervalStart_;

    // Compute number of integration steps.
    numberOfIntegrationSteps_ =  std::ceil( integrationInterval_ / stepsize_ );

    // Set dimension of initial state vector.
    dimensionOfCurrentStateVector_ = vectorOfCurrentStateVectors_.at( 0 )
                                     .rows( );

    // Set dimension of final state vector to correspond with dimension
    // of initial state vector.
    finalStateVector_.setZero( dimensionOfCurrentStateVector_ );

    // Set dimension of state derivative vector to correpsond with
    // dimension of initial state vector.
    stateDerivativeVector_.setZero( dimensionOfCurrentStateVector_ );
};

// End of file.
