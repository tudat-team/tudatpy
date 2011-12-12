/*! \file cartesianStateNumericalPropagator.cpp
 *    Source file that implements the Cartesian state numerical propagator class included in Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 4
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
 *    Checker           : F.M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Date created      : 14 May, 2011
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
 *      110514    K. Kumar          File created.
 *      110810    J. Leloux         Corrected doxygen documentation, unused variable deleted.
 *      110815    K. Kumar          Updated equations of motion, implemented mass.
 *      110920    K. Kumar          Corrected simple errors outlined by M. Persson.
 */

// Macros.
#define TUDAT_UNUSED_PARAMETER( unusedParameter ) { ( void ) unusedParameter; }

// Include statements.
#include <Eigen/Core>
#include "Astrodynamics/Propagators/cartesianStateNumericalPropagator.h"

//! Tudat library namespace.
namespace tudat
{

//! Compute state derivative.
void CartesianStateNumericalPropagator::computeStateDerivative(
        double& independentVariable, State* pointerToAssembledState,
        State* pointerToAssembledStateDerivative )
{
    // Set independentVariable to unused.
    TUDAT_UNUSED_PARAMETER( independentVariable );

    // Declare local variables.
    // State derivative.
    State stateDerivative_;

    // Vector of sum of forces.
    Eigen::VectorXd sumOfForces_( 3 );

    // Loop over map of bodies to propagate.
    for ( BodyPropagatorDataMap::iterator iteratorBodiesToPropagate_ = bodiesToPropagate_.begin( );
          iteratorBodiesToPropagate_ != bodiesToPropagate_.end( ); iteratorBodiesToPropagate_++ )
    {
        // Set state for body to propagate based on associated segment of
        // assembled state.
        iteratorBodiesToPropagate_->second.pointerToInitialState->state
                = pointerToAssembledState->state.segment(
                        iteratorBodiesToPropagate_->second.stateStartIndex,
                        iteratorBodiesToPropagate_->second.sizeOfState );

        // Reset sum of forces to zero.
        sumOfForces_.setZero( 3 );

        // Loop over container of force models for a given body to propagate.
        for ( unsigned int i = 0; i < iteratorBodiesToPropagate_
              ->second.vectorContainerOfPointersToForceModels.size( ); i++ )
        {
            // Compute forces for given force model.
            iteratorBodiesToPropagate_->second.vectorContainerOfPointersToForceModels.at( i )
                    ->computeForce( iteratorBodiesToPropagate_->second.pointerToInitialState );

            // Compute sum of forces.
            sumOfForces_ += iteratorBodiesToPropagate_->second
                    .vectorContainerOfPointersToForceModels.at( i )->getForce( );
        }

        // Set state derivative to size of state.
        stateDerivative_.state.setZero( 6 );

        // Set state derivative elements using state and computed sum of forces.
        // Set derivative of position equal to the velocity.
        stateDerivative_.state.segment( 0, 3 ) = iteratorBodiesToPropagate_->second
                .pointerToInitialState->state.segment( 3, 3 );

        // Set derivative of velocity equal to acceleration.
        stateDerivative_.state.segment( 3, 3 ) = sumOfForces_
                / iteratorBodiesToPropagate_->first->getMass( );

        // Add computed state derivatives to assembled state derivative.
        pointerToAssembledStateDerivative->state.segment( iteratorBodiesToPropagate_
                                                          ->second.stateStartIndex, 6 )
                = stateDerivative_.state;
    }
}

}

// End of file.
