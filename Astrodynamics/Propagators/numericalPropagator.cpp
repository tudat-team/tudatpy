/*! \file numericalPropagator.cpp
 *    Source file that defines the numerical propagator class included in Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 8
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
 *    Date created      : 15 September, 2010
 *    Last modified     : 20 September, 2011
 *
 *    References
 *
 *    Notes
 *      The code only handles 1st-order Ordinary Differential Equations (ODEs)
 *      at the moment. In future, the code will be changed to accommodate
 *      2nd-order ODE's.
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
 *      100915    K. Kumar          File created.
 *      100928    K. Kumar          Modified propagate function; added fixed output interval
 *                                  option.
 *      100929    J. Melman         Corrected several alignments; added some comments.
 *      100929    K. Kumar          Minor comment modifications.
 *      110202    K. Kumar          Updated code to use Integrator adaptor instead of
 *                                  pointers-to-member functions; updated code to use State class.
 *      110203    J. Melman         Possibly encountered a mistake regarding the addition of the
 *                                  forces.
 *      110207    K. Kumar          Updated code with corrections.
 *      110920    K. Kumar          Corrected simple errors outlined by M. Persson.
 */

// Include statements.
#include "Astrodynamics/Propagators/numericalPropagator.h"

//! Set initial state of body.
void NumericalPropagator::setInitialState( Body* pointerToBody, State* pointerToInitialState )
{
    // Set initial state of given body to be propagated.
    bodiesToPropagate_[ pointerToBody ].pointerToInitialState = pointerToInitialState;

    // Set size of initial state of given body to propagate.
    bodiesToPropagate_[ pointerToBody ].sizeOfState = pointerToInitialState->state.size( );

    // Increment size of assembled state based on size of initial state of body to propagate.
    sizeOfAssembledState_ += pointerToInitialState->state.size( );

    // Set assembled state to determined size with zero values.
    assembledState_.state.setZero( sizeOfAssembledState_ );
}

//! Propagate.
void NumericalPropagator::propagate( )
{
    // Declare local variables.
    // Start position for assembly of state.
    int startAssemblyPositionInState_ = 0;

    // Set integration interval start and end.
    pointerToIntegrator_->setIntegrationIntervalStart( propagationIntervalStart_ );
    pointerToIntegrator_->setIntegrationIntervalEnd( propagationIntervalEnd_ );

    // Loop over map of bodies to propagate.
    for ( iteratorBodiesToPropagate_ = bodiesToPropagate_.begin( );
          iteratorBodiesToPropagate_ != bodiesToPropagate_.end( );
          iteratorBodiesToPropagate_++ )
    {
        // Store start of range of initial state of a given body in assembled
        // initial state.
        iteratorBodiesToPropagate_->second.stateStartIndex = startAssemblyPositionInState_;

        // Assemble state using initial states of bodies to propagate.
        assembledState_.state.segment( iteratorBodiesToPropagate_->second.stateStartIndex,
                                       iteratorBodiesToPropagate_->second.sizeOfState )
                = iteratorBodiesToPropagate_->second.pointerToInitialState->state;

        // Start position in initial state for assembly.
        startAssemblyPositionInState_ += iteratorBodiesToPropagate_->second.sizeOfState;
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
        iteratorBodiesToPropagate_->second.finalState.state
                = pointerToIntegrator_->getFinalState( )->state.segment(
                    iteratorBodiesToPropagate_->second.stateStartIndex,
                    iteratorBodiesToPropagate_->second.sizeOfState );
    }
}

// End of file.
