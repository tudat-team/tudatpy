/*! \file cartesianStateNumericalPropagator.cpp
 *    Source file that implements the Cartesian state numerical propagator
 *    class included in Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 3
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
 *    Checker           : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Date created      : 14 May, 2011
 *    Last modified     : 10 August, 2011
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
 *      YYMMDD    Author              Comment
 *      110514    K. Kumar            File created.
 *      110810    J. Leloux           Corrected doxygen documentation, unused
 *                                    variable deleted.
 *      110815    K. Kumar            Updated equations of motion, implemented mass.
 */

// Include statements.
#include "cartesianStateNumericalPropagator.h"

//! Default constructor.
CartesianStateNumericalPropagator::CartesianStateNumericalPropagator( )
{
}

//! Default destructor.
CartesianStateNumericalPropagator::~CartesianStateNumericalPropagator( )
{
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
            // Compute forces for given force model.
            iteratorBodiesToPropagate_
                    ->second.vectorContainerOfPointersToForceModels_.at( i )
                    ->computeForce( iteratorBodiesToPropagate_
                                    ->second.pointerToInitialState_ );

            // Compute sum of forces.
            sumOfForces_ += iteratorBodiesToPropagate_->second
                            .vectorContainerOfPointersToForceModels_.at( i )
                            ->getForce( );
        }

        // Set state derivative to size of state.
        stateDerivative_.state.setZero( 6 );

        // Set state derivative elements using state and computed sum of
        // forces.
        // Set derivative of position equal to the velocity.
        stateDerivative_.state.segment( 0, 3 )
                = iteratorBodiesToPropagate_->second
                  .pointerToInitialState_->state.segment( 3, 3 );

        // Set derivative of velocity equal to acceleration.
        stateDerivative_.state.segment( 3, 3 ) = sumOfForces_
                / iteratorBodiesToPropagate_->first->getMass( );

        // Add computed state derivatives to assembled state derivative.
        pointerToAssembledStateDerivative
                ->state.segment( iteratorBodiesToPropagate_
                                 ->second.stateStartIndex_, 6 )
                = stateDerivative_.state;
    }
}

// End of file.
