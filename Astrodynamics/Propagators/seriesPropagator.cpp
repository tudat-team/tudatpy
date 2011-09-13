/*! \file seriesPropagator.cpp
 *    Source file that defines a class that executes a series of propagation
 *    steps.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 1
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
 *    Last modified     : 6 May, 2011
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
 */

// Include statements.
#include "seriesPropagator.h"

// Using declarations.
using std::floor;

//! Default constructor.
SeriesPropagator::SeriesPropagator( ) : fixedOutputInterval_( -0.0 )
{
}

//! Default destructor.
SeriesPropagator::~SeriesPropagator( )
{
}

//! Set start of series propagation.
void SeriesPropagator::setSeriesPropagationStart( const double&
                                                  seriesPropagationStart )
{
    seriesPropagationStart_ = seriesPropagationStart;
}

//! Set end of series propagation.
void SeriesPropagator::setSeriesPropagationEnd( const double&
                                                seriesPropagationEnd )
{
    seriesPropagationEnd_ = seriesPropagationEnd;
}

//! Set fixed output interval.
void SeriesPropagator::setFixedOutputInterval( const double&
                                               fixedOutputInterval )
{
    // Set fixed output interval for propagation output.
    fixedOutputInterval_ = fixedOutputInterval;

    // Compute number of propagation steps minus one.
    // This will only lead to a sensible result if setSeriesPropagationStart()
    // and setSeriesPropagationEnd() have been called.
    // To prevent numerical instabilities from occuring (e.g., ceil( 4 / 2 ) !=
    // ceil( ( 4 + 1e-16 ) / 2 )), the square root of the machine precision is
    // subtracted before applying the ceil function.
    numberOfPropagationSteps_
            = static_cast< unsigned int >(
                ceil( ( seriesPropagationEnd_ - seriesPropagationStart_ )
                / fixedOutputInterval_
                - sqrt( mathematics::MACHINE_PRECISION_DOUBLES ) ) );
}

//! Set initial state of body for series propagation.
void SeriesPropagator::setInitialState(
        Body* pointerToBody, State* pointerToInitialState )
{
    pointerToPropagator_->setInitialState( pointerToBody,
                                           pointerToInitialState );
}

//! Set propagator.
void SeriesPropagator::setPropagator( Propagator* pointerToPropagator )
{
    pointerToPropagator_ = pointerToPropagator;
}

//! Get fixed output interval.
double& SeriesPropagator::getFixedOutputInterval( )
{
    // Return fixed output interval.
    return fixedOutputInterval_;
}

//! Get start of series propagation.
double& SeriesPropagator::getSeriesPropagationStart( )
{
    // Return start of series propagation.
    return seriesPropagationStart_;
}

//! Get end of series propagation.
double& SeriesPropagator::getSeriesPropagationEnd( )
{
    // Return end of series propagation.
    return seriesPropagationEnd_;
}

//! Execute.
void SeriesPropagator::execute( )
{
    // Store initial states in maps for each propagated body in associated
    // propagator data containers.
    for ( iteratorPropagatedBodies_
          = pointerToPropagator_->bodiesToPropagate_.begin( );
          iteratorPropagatedBodies_
          != pointerToPropagator_->bodiesToPropagate_.end( );
          iteratorPropagatedBodies_++ )
    {
        iteratorPropagatedBodies_->second
                .propagationHistory_[ seriesPropagationStart_ ]
                = *iteratorPropagatedBodies_->second.pointerToInitialState_;
    }

    for ( unsigned int i = 0; i < numberOfPropagationSteps_ - 1; i++ )
    {
        // Set start and end of propagation interval.
        pointerToPropagator_
                ->setPropagationIntervalStart( ( i * fixedOutputInterval_ )
                                               + seriesPropagationStart_);
        pointerToPropagator_
                ->setPropagationIntervalEnd( ( ( i + 1 ) * fixedOutputInterval_ )
                                             + seriesPropagationStart_ );

        // Execute propagation.
        pointerToPropagator_->propagate( );

        // Store final states in maps for each propagated body in associated
        // propagator data containers. Set initial states for each propagated
        // body to final states from current propagation segment.
        for ( iteratorPropagatedBodies_
              = pointerToPropagator_->bodiesToPropagate_.begin( );
              iteratorPropagatedBodies_
              != pointerToPropagator_->bodiesToPropagate_.end( );
              iteratorPropagatedBodies_++ )
        {
            iteratorPropagatedBodies_->second
                    .propagationHistory_[ ( i + 1 ) * fixedOutputInterval_ ]
                    = iteratorPropagatedBodies_->second.finalState_;

            // Set initial state for next propagation segment to final state from
            // current segment.
            iteratorPropagatedBodies_->second.pointerToInitialState_->state
                    = iteratorPropagatedBodies_->second.finalState_.state;
         }
    }

    // Set start and end of propagation interval for last propagation step.
    pointerToPropagator_
            ->setPropagationIntervalStart( ( ( numberOfPropagationSteps_ - 1 )
                                           * fixedOutputInterval_ )
                                           + seriesPropagationStart_  );
    pointerToPropagator_
            ->setPropagationIntervalEnd( seriesPropagationEnd_ );

    // Execute propagation.
    pointerToPropagator_->propagate( );

    // Store final states in maps for each propagated body in associated
    // propagator data containers. Set initial states for each propagated
    // body to final states from current propagation segment.
    for ( iteratorPropagatedBodies_
          = pointerToPropagator_->bodiesToPropagate_.begin( );
          iteratorPropagatedBodies_
          != pointerToPropagator_->bodiesToPropagate_.end( );
          iteratorPropagatedBodies_++ )
    {
        iteratorPropagatedBodies_->second
                .propagationHistory_[ seriesPropagationEnd_ ]
                = iteratorPropagatedBodies_->second.finalState_;
    }
}

//! Get propagation history of body at fixed output intervals.
map< double, State >&
        SeriesPropagator::getPropagationHistoryAtFixedOutputIntervals(
                Body* pointerToBody )
{
    // Return propagation history of given body.
    return pointerToPropagator_
            ->bodiesToPropagate_[ pointerToBody ].propagationHistory_;
}

// End of file.
