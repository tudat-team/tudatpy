/*! \file numericalPropagator.cpp
 *    Source file that defines the numerical propagator class included in
 *    Tudat.
 *
 *    Path              : /Astrodynamics/Propagator/
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
 *    Date created      : 15 September, 2010
 *    Last modified     : 29 September, 2010
 *
 *    References
 *
 *    Notes
 *      Exception handling not included yet. If necessary set functions are not
 *      defined, the propagateOrbit function should throw an error; this still
 *      has to be implemented.
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
 *      100915    K. Kumar            File created.
 *      100928    K. Kumar            Modified propagate function; added fixed
 *                                    output interval option.
 *      100929    J. Melman           Corrected several alignments; added
 *                                    some comments.
 *      100929    K. Kumar            Minor comment modifications.
 */

// Include statements.
#include "numericalPropagator.h"

//! Default constructor.
NumericalPropagator::NumericalPropagator( )
{
    sizeOfAssembledStateVector_ = 0;
};

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

//! Propagate.
void NumericalPropagator::propagate( )
{
    // Declare local variables.
    // Start position for assembly of state vector.
    int startAssemblyPositionInStateVector_ = 0;

    // Declare target integration interval value, i.e., the value of
    // the independent variable at the end of the integration interval.
    double targetIntegrationIntervalValue_;

    // Set integration interval start and end.
    pointerToIntegrator_
            ->setIntegrationIntervalStart( propagationIntervalStart_ );
    pointerToIntegrator_
            ->setIntegrationIntervalEnd( propagationIntervalEnd_ );

    // Set assembled state vector to determined size with zero values.
    assembledStateVector_.setZero( sizeOfAssembledStateVector_ );

    // Set assembled state derivative vector to determined size with zero
    // values.
    assembledStateDerivativeVector_.setZero( sizeOfAssembledStateVector_ );

    // Loop over map of bodies to be propagated.
    for ( iteratorBodiesToBePropagated_ =
                  bodiesToBePropagated_.begin( );
          iteratorBodiesToBePropagated_ !=
                  bodiesToBePropagated_.end( );
          iteratorBodiesToBePropagated_++ )
    {
        // Set propagator to this NumericalPropagator object for all bodies.
        iteratorBodiesToBePropagated_
                ->second->pointerToPropagator_ = this;

        // Store start of range of initial state vector of a given
        // body in assembled initial state vector.
        iteratorBodiesToBePropagated_
                ->second->stateVectorStartIndex_ =
                        startAssemblyPositionInStateVector_;

        // Assemble state vector using initial state vectors of bodies
        // to be propagated.
        assembledStateVector_
                .segment( iteratorBodiesToBePropagated_
                                  ->second->stateVectorStartIndex_,
                          iteratorBodiesToBePropagated_
                                  ->second->sizeOfStateVector_ )
        = iteratorBodiesToBePropagated_->second->stateVector_;

        // Start Position in initial state vector for assembly.
        startAssemblyPositionInStateVector_ +=
                iteratorBodiesToBePropagated_
                        ->second->sizeOfStateVector_;
    }

    // Set assembled initial state vector as state vector for integrator.
    pointerToIntegrator_
            ->setInitialStateVector( assembledStateVector_ );

    // Set computeSumOfStateDerivatives function as state derivative function
    // for integrator.
    pointerToIntegrator_
            ->setStateDerivativeFunction(
                   &NumericalPropagator::computeSumOfStateDerivatives_, this );

    // Check if fixed output interval is set.
    if ( fixedOutputInterval_ != -0.0 )
    {
        // Set boolean flag to save integration history to true.
        pointerToIntegrator_->saveIntegrationHistory( true );

        // Execute integration.
        pointerToIntegrator_->integrate( );

        // Store integration history.
        integrationHistory_ = pointerToIntegrator_->getIntegrationHistory( );

        // Compute numerb of output intervals.
        unsigned int numberOfOutputsIntervals = std
                                       ::floor( ( propagationIntervalEnd_
                                                  - propagationIntervalStart_ )
                                                  / fixedOutputInterval_ ) + 1;

        // Compute values at successive output interval values.
        for ( unsigned int i = 0; i < numberOfOutputsIntervals; i++ )
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
            for ( iteratorBodiesToBePropagated_ =
                          bodiesToBePropagated_.begin( );
                  iteratorBodiesToBePropagated_ !=
                          bodiesToBePropagated_.end( );
                  iteratorBodiesToBePropagated_++ )
            {
                // Disassemble assembled initial state vector by updating
                // propagation histories in respective BodyContainer objects.
                // PLACEHOLDER
                iteratorBodiesToBePropagated_->second
                        ->propagationHistory_[ iteratorPropagationHistory_->first ]
                        = iteratorPropagationHistory_->second
                          .segment( iteratorBodiesToBePropagated_
                                            ->second->stateVectorStartIndex_,
                                    iteratorBodiesToBePropagated_
                                            ->second->sizeOfStateVector_ );

                // Store state vector as final state vector
                iteratorBodiesToBePropagated_->second->finalStateVector_ =
                        pointerToIntegrator_->getFinalStateVector( )
                        .segment( iteratorBodiesToBePropagated_
                                          ->second->stateVectorStartIndex_,
                                  iteratorBodiesToBePropagated_
                                          ->second->sizeOfStateVector_ );
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
        for ( iteratorBodiesToBePropagated_ =
                      bodiesToBePropagated_.begin( );
              iteratorBodiesToBePropagated_ !=
                      bodiesToBePropagated_.end( );
              iteratorBodiesToBePropagated_++ )
        {
            // Disassemble assembled initial state vector by updating state
            // vectors in respective BodyContainer objects.
            iteratorBodiesToBePropagated_->second->stateVector_ =
                    pointerToIntegrator_->getFinalStateVector( )
                    .segment( iteratorBodiesToBePropagated_
                                      ->second->stateVectorStartIndex_,
                              iteratorBodiesToBePropagated_
                                      ->second->sizeOfStateVector_ );

            // Store state vector as final state vector
            iteratorBodiesToBePropagated_->second->finalStateVector_
                    = iteratorBodiesToBePropagated_->second->stateVector_;
        }
    }
}

// End of file.
