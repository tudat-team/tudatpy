/*! \file propagator.cpp
 *    Source file that defines the base class for all propagators included in
 *    Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 7
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
 *    Date created      : 26 September, 2010
 *    Last modified     : 7 February, 2011
 *
 *    References
 *
 *    Notes
 *      The propagator code has now been setup to work with 1st-order ODEs. Use
 *      of 2nd-order ODEs will result in errors. This will be corrected in
 *      future versions.
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
 *      100926    K. Kumar          File created.
 *      100929    J. Melman         Deleted some superfluous comments,
 *                                  corrected several alignments.
 *      100929    K. Kumar          Minor comment modifications and obsolete
 *                                  code removed.
 *      110119    K. Kumar          Modified computeSumOfStateDerivatives_()
 *                                  to use computeForce() in ForceModel.
 *      110124    K. Kumar          Updated computeSumOfStateDerivatives_()
 *                                  and added to "Notes".
 *      110201    K. Kumar          Updated code to make use of State class;
 *                                  moved computeSumOfStateDerivatives_() to
 *                                  NumericalPropagator class.
 */

// Include statements.
#include "propagator.h"

// Using declarations.
using std::endl;

//! Default constructor.
Propagator::Propagator( ) : fixedOutputInterval_( -0.0 )
{
}

//! Default destructor.
Propagator::~Propagator( )
{
}

//! Set start of propagation interval.
void Propagator::setPropagationIntervalStart( const double&
                                              propagationIntervalStart )
{
    propagationIntervalStart_ = propagationIntervalStart;
}

//! Set end of propagation interval.
void Propagator::setPropagationIntervalEnd( const double&
                                            propagationIntervalEnd )
{
    propagationIntervalEnd_ = propagationIntervalEnd;
}

//! Add body to propagate.
void Propagator::addBody( Body* pointerToBody )
{
    // Add body as key for map, with new PropagatorDataContainer object as data.
    bodiesToPropagate_[ pointerToBody ] = new PropagatorDataContainer( );

    // Set body in new PropagatorDataContainer object as given body.
    bodiesToPropagate_[ pointerToBody ]->pointerToBody_ = pointerToBody;
}

//! Set propagator for propagation of a body.
void Propagator::setPropagator( Body* pointerToBody,
                                Propagator* pointerToPropagator )
{
    // Set pointer to object of Propagator class for given pointer to body in
    // map of bodies to be propagated.
    bodiesToPropagate_[ pointerToBody ]
            ->pointerToPropagator_ = pointerToPropagator;
}

//! Set initial state of body.
void Propagator::setInitialState( Body* pointerToBody,
                                  State* pointerToInitialState )
{
    // Set initial state of given body to be propogated.
    bodiesToPropagate_[ pointerToBody ]
            ->pointerToInitialState_ = pointerToInitialState;

    // Set size of initial state of given body to propagate.
    bodiesToPropagate_[ pointerToBody ]
            ->sizeOfState_ = pointerToInitialState->state.size( );

    // Set state in body to propagate as initial state.
    bodiesToPropagate_[ pointerToBody ]
            ->pointerToCurrentState_ = pointerToInitialState;
}

//! Set fixed output interval.
void Propagator::setFixedOutputInterval( const double& fixedOutputInterval )
{
    // Set fixed output interval for propagation output.
    fixedOutputInterval_ = fixedOutputInterval;
}

//! Get start of propagation interval.
double& Propagator::getPropagationIntervalStart( )
{
    // Return start of propagation interval.
    return propagationIntervalStart_;
}

//! Get end of propagation interval.
double& Propagator::getPropagationIntervalEnd( )
{
    // Return end of propagation interval.
    return propagationIntervalEnd_;
}

//! Get final state of body.
State* Propagator::getFinalState( Body* pointerToBody )
{
    // Return final state of given body.
    return bodiesToPropagate_[ pointerToBody ]->pointerToFinalState_;
}

//! Get propagation history of body at fixed output intervals.
std::map < double, State* >
        Propagator::getPropagationHistoryAtFixedOutputIntervals(
                Body* pointerToBody )
{
    // Return propagation history of given body.
    return bodiesToPropagate_[ pointerToBody ]->propagationHistory_;
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                          Propagator* pointerToPropagator )
{
    stream << "The start of the propagation interval is set to: " << endl;
    stream << pointerToPropagator->getPropagationIntervalStart( ) << endl;
    stream << "The end of the propagation interval is set to: " << endl;
    stream << pointerToPropagator->getPropagationIntervalEnd( ) << endl;

    // Return stream.
    return stream;
}

// End of file.
