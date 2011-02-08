/*! \file numericalPropagator.cpp
 *    Source file that defines the numerical propagator class included in
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
 *    Date created      : 15 September, 2010
 *    Last modified     : 7 February, 2011
 *
 *    References
 *
 *    Notes
 *      The code only handles 1st-order Ordinary Differential Equations (ODEs)
 *      at the moment. In future, the code will be changed to accommodate
 *      2nd-order ODE's.
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
 *      100915    K. Kumar          File created.
 *      100928    K. Kumar          Modified propagate function; added fixed
 *                                  output interval option.
 *      100929    J. Melman         Corrected several alignments; added some
 *                                  comments.
 *      100929    K. Kumar          Minor comment modifications.
 *      110202    K. Kumar          Updated code to use Integrator adaptor
 *                                  instead of pointers-to-member functions;
 *                                  Update code to use State class.
 *      110203    J. Melman         Possibly encountered a mistake regarding
 *                                  the addition of the forces.
 *      110207    K. Kumar          Updated code with corrections.
 */

// Include statements.
#include "numericalPropagator.h"

//! Default constructor.
NumericalPropagator::NumericalPropagator( ) : sizeOfAssembledState_( 0 )
{
    // Initialize variables.
    pointerToAssembledState_ = &assembledState_;
    pointerToAssembledStateDerivative_ = &assembledStateDerivative_;
}

//! Default destructor.
NumericalPropagator::~NumericalPropagator( )
{
}

//! Set integrator for propagation.
void NumericalPropagator::setIntegrator( Integrator* pointerToIntegrator )
{
    // Set pointer to object of Integrator class.
    pointerToIntegrator_ = pointerToIntegrator;
}

//! Add force model for propagation of a body.
void NumericalPropagator::addForceModel( Body* pointerToBody,
                                         ForceModel* pointerToForceModel )
{
    // Set force model in vector container for given body.
    bodiesToPropagate_[ pointerToBody ]
            ->vectorContainerOfPointersToForceModels_
                     .push_back( pointerToForceModel );
}

//! Propagate.
void NumericalPropagator::propagate( )
{
    // Declare local variables.
    // Start position for assembly of state.
    int startAssemblyPositionInState_ = 0;

    // Declare target integration interval value, i.e., the value of
    // the independent variable at the end of the integration interval.
    double targetIntegrationIntervalValue_;

    // Set the class that contains the derivative function needed for the
    // integrator.
    integratorAdaptorForNumericalPropagator_.setClass( this );

    // Set state derivative function used by integrator.
    integratorAdaptorForNumericalPropagator_.setStateDerivativeFunction(
            &NumericalPropagator::computeSumOfStateDerivatives_ );

    // Set integration interval start and end.
    pointerToIntegrator_
            ->setIntegrationIntervalStart( propagationIntervalStart_ );
    pointerToIntegrator_
            ->setIntegrationIntervalEnd( propagationIntervalEnd_ );

    // Set the adaptor for integrator.
    pointerToIntegrator_->setIntegratorAdaptor(
            &integratorAdaptorForNumericalPropagator_ );

    // Loop over map of bodies to propagate.
    for ( iteratorBodiesToPropagate_ = bodiesToPropagate_.begin( );
          iteratorBodiesToPropagate_ != bodiesToPropagate_.end( );
          iteratorBodiesToPropagate_++ )
    {
        // Increment size of assembled state based on size of initial state of
        // body to propagate.
        sizeOfAssembledState_ += iteratorBodiesToPropagate_->second
                                 ->pointerToInitialState_->state.size( );
    }

    // Set assembled state to determined size with zero values.
    pointerToAssembledState_->state.setZero( sizeOfAssembledState_ );

    // Set assembled state derivative to determined size with zero values
    pointerToAssembledStateDerivative_->state.setZero( sizeOfAssembledState_ );

    // Loop over map of bodies to propagate.
    for ( iteratorBodiesToPropagate_ = bodiesToPropagate_.begin( );
          iteratorBodiesToPropagate_ != bodiesToPropagate_.end( );
          iteratorBodiesToPropagate_++ )
    {
        // Set propagator to this NumericalPropagator object for all bodies.
        iteratorBodiesToPropagate_->second->pointerToPropagator_ = this;

        // Store start of range of initial state of a given body in assembled
        // initial state.
        iteratorBodiesToPropagate_->second->stateStartIndex_
                = startAssemblyPositionInState_;

        // Assemble state using initial states of bodies to propagate.
        pointerToAssembledState_->state.segment( iteratorBodiesToPropagate_
                                                 ->second->stateStartIndex_,
                                                 iteratorBodiesToPropagate_
                                                 ->second->sizeOfState_ )
        = iteratorBodiesToPropagate_->second->pointerToCurrentState_->state;

        // Start position in initial state for assembly.
        startAssemblyPositionInState_ +=
                iteratorBodiesToPropagate_->second->sizeOfState_;
    }

    // Set assembled initial state as state for integrator.
    pointerToIntegrator_->setInitialState( pointerToAssembledState_ );

    // Check if fixed output interval is set.
    if ( fixedOutputInterval_ != -0.0 )
    {
        // Set boolean flag to save integration history to true.
        pointerToIntegrator_->setIsIntegrationHistoryToBeSaved( true );

        // Execute integration.
        pointerToIntegrator_->integrate( );

        // Store integration history.
        integrationHistory_ = pointerToIntegrator_->getIntegrationHistory( );

        // Compute number of output intervals.
        unsigned int numberOfOutputIntervals = std::ceil(
                ( propagationIntervalEnd_ - propagationIntervalStart_ )
                / fixedOutputInterval_ );

        // Compute values at successive output interval values.
        for ( unsigned int i = 0; i < numberOfOutputIntervals - 1; i++ )
        {
            // Compute target integration interval value for interpolation of
            // integration output.
            targetIntegrationIntervalValue_ = i * fixedOutputInterval_;

            // Compute interpolation of integration output and store in map
            // of requested integration output.
            propagationHistory_[ targetIntegrationIntervalValue_ ] =
                    mathematics::computeLinearInterpolation
                    ( integrationHistory_, targetIntegrationIntervalValue_ );
        }

        // Loop over map of propagation history.
        for ( iteratorPropagationHistory_ = propagationHistory_.begin( );
              iteratorPropagationHistory_ != propagationHistory_.end( );
              iteratorPropagationHistory_++ )
        {
            // Loop over map of bodies that were propagated.
            for ( iteratorBodiesToPropagate_ = bodiesToPropagate_.begin( );
                  iteratorBodiesToPropagate_ !=  bodiesToPropagate_.end( );
                  iteratorBodiesToPropagate_++ )
            {
                // Pointer to new State object;
                State* pointerToState_ = new State;

                // Disassemble assembled initial state by updating propagation
                // histories in respective PropagatorDataContainer objects.
                pointerToState_->state = iteratorPropagationHistory_->second
                                         ->state.segment(
                                                 iteratorBodiesToPropagate_
                                                 ->second->stateStartIndex_,
                                                 iteratorBodiesToPropagate_
                                                 ->second->sizeOfState_ );

                iteratorBodiesToPropagate_->second->propagationHistory_[
                        iteratorPropagationHistory_->first ]
                        = pointerToState_;

                // Store state as final state.
                iteratorBodiesToPropagate_->second->pointerToFinalState_->state
                        = pointerToIntegrator_->getFinalState( )->state
                          .segment( iteratorBodiesToPropagate_
                                    ->second->stateStartIndex_,
                                    iteratorBodiesToPropagate_
                                    ->second->sizeOfState_ );
            }
        }
    }

    // Perform the following if the fixed output interval is not set, i.e.,
    // the propagation history does not need to be saved, only the final state.
    else
    {
        // Execute integration.
        pointerToIntegrator_->integrate( );

        // Loop over map of bodies to be propagated.
        for ( iteratorBodiesToPropagate_ = bodiesToPropagate_.begin( );
              iteratorBodiesToPropagate_ != bodiesToPropagate_.end( );
              iteratorBodiesToPropagate_++ )
        {
            // Disassemble assembled initial state by updating state in
            // respective PropagatorDataContainer objects.
            iteratorBodiesToPropagate_->second->pointerToCurrentState_->state
                    = pointerToIntegrator_->getFinalState( )->state
                      .segment( iteratorBodiesToPropagate_
                                ->second->stateStartIndex_,
                                iteratorBodiesToPropagate_
                                ->second->sizeOfState_ );

            // Store state as final state.
            iteratorBodiesToPropagate_->second->pointerToFinalState_
                    = iteratorBodiesToPropagate_
                      ->second->pointerToCurrentState_;
        }
    }
}

//! Compute sum of state derivatives.
State* NumericalPropagator::
        computeSumOfStateDerivatives_( State* pointerToAssembledState )
{
    // Declare local variables.
    unsigned int sizeOfState_;

    // State derivative.
    State* pointerToStateDerivative_ = new State;

    // Vector of sum of forces.
    VectorXd sumOfForces_;

    // Set assembled state to zero.
    pointerToAssembledStateDerivative_
            ->state.setZero( pointerToAssembledState->state.rows( ) );

    // Loop over map of bodies to propagate.
    for ( iteratorBodiesToPropagate_ = bodiesToPropagate_.begin( );
          iteratorBodiesToPropagate_ != bodiesToPropagate_.end( );
          iteratorBodiesToPropagate_++ )
    {
        // Set state for body to propagate based on associated segment of
        // assembled state.
        iteratorBodiesToPropagate_->second->pointerToCurrentState_->state
                = pointerToAssembledState->state.segment(
                        iteratorBodiesToPropagate_->second->stateStartIndex_,
                        iteratorBodiesToPropagate_->second->sizeOfState_ );

        // Compute sum of forces for first force model to size vector.
        sumOfForces_ = iteratorBodiesToPropagate_->second
                       ->vectorContainerOfPointersToForceModels_.at( 0 )
                       ->computeForce( iteratorBodiesToPropagate_
                                       ->second->pointerToCurrentState_ );

        // Reset sized vector to zero.
        sumOfForces_.setZero( 3 );

        // Loop over container of force models for a given body to propagate.
        for ( unsigned int i = 0;
              i < iteratorBodiesToPropagate_->second
              ->vectorContainerOfPointersToForceModels_.size( );
              i++ )
        {
            // Compute sum of forces for given force model.
            sumOfForces_ += iteratorBodiesToPropagate_->second
                           ->vectorContainerOfPointersToForceModels_.at( i )
                           ->computeForce( iteratorBodiesToPropagate_
                                           ->second->pointerToCurrentState_ );
        }

            // Set size of state.
            sizeOfState_ = iteratorBodiesToPropagate_->second->sizeOfState_;

            // Set state derivative to size of state.
            pointerToStateDerivative_->state.setZero( sizeOfState_ );

            // Set state derivative elements using state and computed sum of
            // forces.
            // Set derivative of position equal to the velocity.
            pointerToStateDerivative_
                    ->state.segment( 0, sizeOfState_ - sumOfForces_.rows( ) )
                    = iteratorBodiesToPropagate_->second
                      ->pointerToCurrentState_->state.segment(
                              sumOfForces_.rows( ),
                              sizeOfState_ - sumOfForces_.rows( ) );

            // Set derivative of velocity equal to acceleration ( specific force ).
            pointerToStateDerivative_->state.segment( sizeOfState_
                                                      - sumOfForces_.rows( ),
                                                      sumOfForces_.rows( ) )
                    = sumOfForces_;

            // Add computed state derivatives to assembled state derivative.
            pointerToAssembledStateDerivative_
                    ->state.segment( iteratorBodiesToPropagate_
                                     ->second->stateStartIndex_,
                                     sizeOfState_ )
                    = pointerToStateDerivative_->state;
    }

    // De-allocate pointer to state derivative.
    delete pointerToStateDerivative_;

    // Return pointer to assembled state derivative.
    return pointerToAssembledStateDerivative_;
}

// End of file.
