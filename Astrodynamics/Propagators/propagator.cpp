/*! \file propagator.cpp
 *    Source file that defines the base class for all propagators included in
 *    Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 5
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
 *    Last modified     : 24 January, 2010
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
 *      YYMMDD    author              comment
 *      100926    K. Kumar            File created.
 *      100929    J. Melman           Deleted some superfluous comments,
 *                                    corrected several alignments.
 *      100929    K. Kumar            Minor comment modifications and obsolete
 *                                    code removed.
 *      110119    K. Kumar            Modified computeSumOfStateDerivatives_()
 *                                    to use computeForce() in ForceModel.
 *      110124    K. Kumar            Updated computeSumOfStateDerivatives_()
 *                                    and added to "Notes".
 */

// Include statements.
#include "propagator.h"

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

//! Add body to be propagated.
void Propagator::addBody( Body* pointerToBody )
{
    // Add body as key for map, with new BodyContainer object as data.
    bodiesToBePropagated_[ pointerToBody ] = new BodyContainer( );

    // Set body in new BodyContainer object as given body.
    bodiesToBePropagated_[ pointerToBody ]->pointerToBody_ = pointerToBody;

//    // Set pointer to Body object.
//    pointerToBody_ = pointerToBody;
}

//! Add force model for propagation of a body.
void Propagator::addForceModel( Body* pointerToBody,
                                ForceModel* pointerToForceModel )
{
    bodiesToBePropagated_[ pointerToBody ]
            ->vectorContainerOfPointersToForceModels_
                     .push_back( pointerToForceModel );
}

//! Set propagator for propagation of a body.
void Propagator::setPropagator( Body* pointerToBody,
                                Propagator* pointerToPropagator )
{
    // Set pointer to object of Propagator class for given pointer to body in
    // map of bodies to be propagated.
    bodiesToBePropagated_[ pointerToBody ]
            ->pointerToPropagator_ = pointerToPropagator;
}

//! Set initial state vector of body.
void Propagator::setInitialState( Body* pointerToBody,
                                  VectorXd& initialStateVector )
{
    // Set initial state vector of given body to be propogated.
    bodiesToBePropagated_[ pointerToBody ]
            ->initialStateVector_ = initialStateVector;

    // Set size of initial state vector of given body to be propagated.
    bodiesToBePropagated_[ pointerToBody ]
            ->sizeOfStateVector_ = initialStateVector.size( );

    // Increment size of assembled state vector based on size of initial state
    // vector of body to be propagated.
    sizeOfAssembledStateVector_ += initialStateVector.size( );

    // Set state vector in body to be propagated as initial state vector.
    bodiesToBePropagated_[ pointerToBody ]
            ->stateVector_ = initialStateVector;
}

//! Set fixed output interval.
void Propagator::setFixedOutputInterval( const double& fixedOutputInterval )
{
    // Set fixed output interval for propagation output.
    fixedOutputInterval_ = fixedOutputInterval;
}

//! Get final state of body.
VectorXd& Propagator::getFinalState( Body* pointerToBody )
{
    // Return final state of given body.
    return bodiesToBePropagated_[ pointerToBody ]->finalStateVector_;
}

//! Get propagation history of body at fixed output intervals.
std::map < double, VectorXd >
        Propagator::getPropagationHistoryAtFixedOutputIntervals(
                Body* pointerToBody )
{
    // Return propagation history of given body.
    return bodiesToBePropagated_[ pointerToBody ]->propagationHistory_;
}

void Propagator::
        computeSumOfStateDerivatives_( VectorXd& assembledStateVector,
                                       VectorXd& assembledStateDerivativeVector )
{
    // Declare local variables.
    unsigned int sizeOfStateVector_;

    // State derivative vector.
    VectorXd stateDerivativeVector_;

    // Vector of sum of forces.
    VectorXd sumOfForces_;

    // Set assembled state vector to zero.
    assembledStateDerivativeVector.setZero( assembledStateVector.rows( ) );

    // Loop over map of bodies to be propagated.
    for ( iteratorBodiesToBePropagated_ =
         bodiesToBePropagated_.begin( );
         iteratorBodiesToBePropagated_ !=
         bodiesToBePropagated_.end( );
         iteratorBodiesToBePropagated_++)
    {
        // Set state vector for body to be propagated based on associated
        // segment of assembled state vector.
        iteratorBodiesToBePropagated_->second->stateVector_ =
                assembledStateVector.segment(
                    iteratorBodiesToBePropagated_->second
                    ->stateVectorStartIndex_,
                    iteratorBodiesToBePropagated_->second
                    ->sizeOfStateVector_ );

        // Loop over container of force models for a given body to be
        // propagated.
        for ( unsigned int i = 0;
             i < iteratorBodiesToBePropagated_->second
             ->vectorContainerOfPointersToForceModels_.size( );
             i++ )
        {
            // Compute state derivatives for given force model.
            sumOfForces_ = iteratorBodiesToBePropagated_->second
                    ->vectorContainerOfPointersToForceModels_.at( i )
                    ->computeForce( iteratorBodiesToBePropagated_
                                   ->second->stateVector_ );

            // Set size of state vector.
            sizeOfStateVector_ = iteratorBodiesToBePropagated_
                    ->second->sizeOfStateVector_;

            // Set state derivative vector to size of state vector.
            stateDerivativeVector_.setZero( sizeOfStateVector_ );

            // Set state derivative vector elements using state vector and
            // computed sum of forces.
            stateDerivativeVector_.segment( 0, sizeOfStateVector_
                                           - sumOfForces_.rows( ) )
                    = iteratorBodiesToBePropagated_
                    ->second->stateVector_
                    .segment( sumOfForces_.rows( ),
                              sizeOfStateVector_
                              - sumOfForces_.rows( ) );

            stateDerivativeVector_.segment( sizeOfStateVector_
                                            - sumOfForces_.rows( ),
                                            sumOfForces_.rows( ) )
                    = sumOfForces_;

            // Add computed state derivatives to assembled state derivative
            // vector.
            assembledStateDerivativeVector
                    .segment( iteratorBodiesToBePropagated_
                              ->second->stateVectorStartIndex_,
                              sizeOfStateVector_ )
                    = stateDerivativeVector_;
        }
    }
}

// End of file.
