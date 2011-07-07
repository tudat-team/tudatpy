/*! \file cartesianStateNumericalPropagator.cpp
 *    Source file that implements the Cartesian state numerical propagator
 *    class included in Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 1
 *    Check status      : Unchecked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 14 May, 2011
 *    Last modified     : 14 May, 2011
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
 *      YYMMDD    Author              Comment
 *      110514    K. Kumar            File created.
 */

// Include statements.
#include "cartesianStateNumericalPropagator.h"

//! Default constructor.
CartesianStateNumericalPropagator::CartesianStateNumericalPropagator( )
    : sizeOfAssembledState_( 0 )
{
}

//! Default destructor.
CartesianStateNumericalPropagator::~CartesianStateNumericalPropagator( )
{
}

//! Propagate.
void CartesianStateNumericalPropagator::propagate( )
{
    // Declare local variables.
    // Start position for assembly of state.
    int startAssemblyPositionInState_ = 0;

    // Set integration interval start and end.
    pointerToIntegrator_
            ->setIntegrationIntervalStart( propagationIntervalStart_ );
    pointerToIntegrator_
            ->setIntegrationIntervalEnd( propagationIntervalEnd_ );

    // Loop over map of bodies to propagate.
    for ( iteratorBodiesToPropagate_ = bodiesToPropagate_.begin( );
          iteratorBodiesToPropagate_ != bodiesToPropagate_.end( );
          iteratorBodiesToPropagate_++ )
    {
        // Set propagator to this NumericalPropagator object for this body.
        iteratorBodiesToPropagate_->second.pointerToPropagator_ = this;

        // Store start of range of initial state of a given body in assembled
        // initial state.
        iteratorBodiesToPropagate_->second.stateStartIndex_
                = startAssemblyPositionInState_;

        // Assemble state using initial states of bodies to propagate.
        assembledState_.state.segment( iteratorBodiesToPropagate_
                                       ->second.stateStartIndex_,
                                       iteratorBodiesToPropagate_
                                       ->second.sizeOfState_ )
                = iteratorBodiesToPropagate_
                  ->second.pointerToInitialState_->state;

        // Start position in initial state for assembly.
        startAssemblyPositionInState_ +=
                iteratorBodiesToPropagate_->second.sizeOfState_;
    }

    // Set assembled initial state as state for integrator.
    pointerToIntegrator_->setInitialState( &assembledState_ );

    // Execute integration.
    pointerToIntegrator_->integrate( );

    // Loop over map of bodies to be propagated.
    for ( iteratorBodiesToPropagate_ = bodiesToPropagate_.begin( );
          iteratorBodiesToPropagate_ != bodiesToPropagate_.end( );
          iteratorBodiesToPropagate_++ )
    {
        // Disassemble assembled initial state by updating final state in
        // respective PropagatorDataContainer objects.
        iteratorBodiesToPropagate_->second.finalState_.state
                = pointerToIntegrator_->getFinalState( )->state
                  .segment( iteratorBodiesToPropagate_
                            ->second.stateStartIndex_,
                            iteratorBodiesToPropagate_
                            ->second.sizeOfState_ );
    }
}

//! Set initial state of body.
void CartesianStateNumericalPropagator::setInitialState(
        Body* pointerToBody, State* pointerToInitialState )
{
    // Set initial state of given body to be propagated.
    bodiesToPropagate_[ pointerToBody ]
            .pointerToInitialState_ = pointerToInitialState;

    // Set size of initial state of given body to propagate.
    bodiesToPropagate_[ pointerToBody ]
            .sizeOfState_ = pointerToInitialState->state.size( );

    // Increment size of assembled state based on size of initial state of
    // body to propagate.
    sizeOfAssembledState_ += pointerToInitialState->state.size( );

    // Set assembled state to determined size with zero values.
    assembledState_.state.setZero( sizeOfAssembledState_ );
}

//! Compute sum of state derivatives.
void CartesianStateNumericalPropagator::computeStateDerivative(
        double& independentVariable, State* pointerToAssembledState,
        State* pointerToAssembledStateDerivative )
{
    // Declare local variables.
    // State derivative.
    State stateDerivative_;

    // Vector of sum of forces.
    VectorXd sumOfForces_( 3 );

    // Loop over map of bodies to propagate.
    for ( iteratorBodiesToPropagate_ = bodiesToPropagate_.begin( );
          iteratorBodiesToPropagate_ != bodiesToPropagate_.end( );
          iteratorBodiesToPropagate_++ )
    {
        // Set state for body to propagate based on associated segment of
        // assembled state.
        iteratorBodiesToPropagate_->second.pointerToInitialState_->state
                = pointerToAssembledState->state.segment(
                        iteratorBodiesToPropagate_->second.stateStartIndex_,
                        iteratorBodiesToPropagate_->second.sizeOfState_ );

        // Reset sum of forces to zero.
        sumOfForces_.setZero( 3 );

        // Loop over container of force models for a given body to propagate.
        for ( unsigned int i = 0;
              i < iteratorBodiesToPropagate_->second
              .vectorContainerOfPointersToForceModels_.size( );
              i++ )
        {
            // Compute sum of forces for given force model.
            sumOfForces_ += iteratorBodiesToPropagate_->second
                            .vectorContainerOfPointersToForceModels_.at( i )
                            ->computeForce( iteratorBodiesToPropagate_
                                            ->second.pointerToInitialState_ );
        }

        // Set state derivative to size of state.
        stateDerivative_.state.setZero( 6 );

        // Set state derivative elements using state and computed sum of
        // forces.
        // Set derivative of position equal to the velocity.
        stateDerivative_.state.segment( 0, 3 )
                = iteratorBodiesToPropagate_->second
                  .pointerToInitialState_->state.segment( 3, 3 );

        // Set derivative of velocity equal to acceleration (specific force).
        stateDerivative_.state.segment( 3, 3 ) = sumOfForces_;

        // Add computed state derivatives to assembled state derivative.
        pointerToAssembledStateDerivative
                ->state.segment( iteratorBodiesToPropagate_
                                 ->second.stateStartIndex_, 6 )
                = stateDerivative_.state;
    }
}

// End of file.
