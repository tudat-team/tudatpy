/*! \file seriesPropagator.cpp
 *    Source file that defines a class that executes a series of propagation steps.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 6 May, 2011
 *    Last modified     : 20 September, 2011
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
 *      110506    K. Kumar          File created.
 *      110920    K. Kumar          Corrected simple errors outlined by M. Persson.
 */

// Include statements.
#include "Astrodynamics/Propagators/seriesPropagator.h"

//! Execute.
void SeriesPropagator::execute( )
{
    // Get map of bodies to propagate.
    bodiesToPropagate_ = pointerToPropagator_->getBodiesToPropagate( );

    // Clear and store initial states in maps for each propagated body in
    // associated propagator data containers.
    for ( BodyPropagatorDataMap::iterator iteratorPropagatedBodies_ = bodiesToPropagate_.begin( );
          iteratorPropagatedBodies_ != bodiesToPropagate_.end( ); iteratorPropagatedBodies_++ )
    {
        // Clear propagation histories.
        iteratorPropagatedBodies_->second.propagationHistory.clear( );

        // Store initial states.
        iteratorPropagatedBodies_->second.propagationHistory[ seriesPropagationStart_ ]
                = *iteratorPropagatedBodies_->second.pointerToInitialState;
    }

    for ( unsigned int i = 0; i < numberOfPropagationSteps_ - 1; i++ )
    {
        // Set start and end of propagation interval.
        pointerToPropagator_->setPropagationIntervalStart( ( i * fixedOutputInterval_ )
                                                           + seriesPropagationStart_ );
        pointerToPropagator_->setPropagationIntervalEnd( ( ( i + 1 ) * fixedOutputInterval_ )
                                                         + seriesPropagationStart_ );

        // Execute propagation.
        pointerToPropagator_->propagate( );

        // Get map of bodies to propagate.
        bodiesToPropagate_ = pointerToPropagator_->getBodiesToPropagate( );

        // Store final states in maps for each propagated body in associated
        // propagator data containers. Set initial states for each propagated
        // body to final states from current propagation segment.
        for ( BodyPropagatorDataMap::iterator iteratorPropagatedBodies_
              = bodiesToPropagate_.begin( );
              iteratorPropagatedBodies_ != bodiesToPropagate_.end( ); iteratorPropagatedBodies_++ )
        {
            // Store final state.
            iteratorPropagatedBodies_->second.propagationHistory[
                    ( i + 1 ) * fixedOutputInterval_ ]
                    = iteratorPropagatedBodies_->second.finalState;

            // Set initial state for next propagation segment to final state from
            // current segment.
            iteratorPropagatedBodies_->second.pointerToInitialState->state
                    = iteratorPropagatedBodies_->second.finalState.state;

        }

        // Store map of bodies to propagate in propagator.
        pointerToPropagator_->setBodiesToPropagate( bodiesToPropagate_ );
    }

    // Set start and end of propagation interval for last propagation step.
    pointerToPropagator_->setPropagationIntervalStart( ( ( numberOfPropagationSteps_ - 1 )
                                                         * fixedOutputInterval_ )
                                                       + seriesPropagationStart_  );
    pointerToPropagator_->setPropagationIntervalEnd( seriesPropagationEnd_ );

    // Execute propagation.
    pointerToPropagator_->propagate( );

    // Store final states in maps for each propagated body in associated
    // propagator data containers. Set initial states for each propagated
    // body to final states from current propagation segment.
    for ( BodyPropagatorDataMap::iterator iteratorPropagatedBodies_ = bodiesToPropagate_.begin( );
          iteratorPropagatedBodies_ != bodiesToPropagate_.end( ); iteratorPropagatedBodies_++ )
    {
        // Get map of bodies to propagate.
        bodiesToPropagate_ = pointerToPropagator_->getBodiesToPropagate( );

        // Store final state.
        iteratorPropagatedBodies_->second.propagationHistory[ seriesPropagationEnd_ ]
                = iteratorPropagatedBodies_->second.finalState;
    }

    // Store map of bodies to propagate in propagator.
    pointerToPropagator_->setBodiesToPropagate( bodiesToPropagate_ );
}

// End of file.
