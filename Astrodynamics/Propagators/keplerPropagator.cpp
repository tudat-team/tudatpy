/*! \file keplerPropagator.cpp
 *    Source file that defines the kepler propagator class included in
 *    Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Date created      : 3 February, 2011
 *    Last modified     : 7 February, 2011
 *
 *    References
 *
 *    Notes
 *      The code at present does not work for near-parabolic orbits
 *      ( 0.8 < eccentricity < 1.2 ). In future, this neeeds to be included
 *      and perhaps a universal method to solve Kepler's equation needs to be
 *      employed. Presently, the code will output an error if the eccentricity
 *      of the orbit to be propagated lies within this range.
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
 *      110203    K. Kumar          File created.
 *      110207    E. Iorfida        Minor changes.
 */

// Include statements.
#include "keplerPropagator.h"

// Using directives.
using mathematics::raiseToIntegerPower;

//! Default constructor.
KeplerPropagator::KeplerPropagator( ) : numberOfPropagationSteps_( -0 ),
                                        eccentricAnomaly_( -0.0 ),
                                        hyperbolicEccentricAnomaly_( -0.0 ),
                                        hyperbolicMeanAnomaly_( -0.0 ),
                                        meanAnomaly_( -0.0 ),
                                        trueAnomaly_( -0.0 )
{
}

//! Default destructor
KeplerPropagator::~KeplerPropagator( )
{
}

//! Set central body.
void KeplerPropagator::setCentralBody( Body* pointerToBody,
                                       CelestialBody* pointerToCentralBody )
{
    // Set pointer to central body for given body to propagate.
    bodiesToPropagate_[ pointerToBody ]->pointerToCentralBody_
            = pointerToCentralBody;
}

//! Set Newton-Raphson method.
void KeplerPropagator::setNewtonRaphson( NewtonRaphson*
                                         pointerToNewtonRaphson )
{
    // Set pointer to Newton-Raphson method.
    pointerToNewtonRaphson_ = pointerToNewtonRaphson;
}

//! Propagate.
void KeplerPropagator::propagate( )
{
    // Set the class that contains the functions needed for Newton-Raphson.
    newtonRaphsonAdaptorForKeplerPropagator_.setClass( this );

    // Set NewtonRaphson adaptor class.
    pointerToNewtonRaphson_->setNewtonRaphsonAdaptor(
            &newtonRaphsonAdaptorForKeplerPropagator_ );

    // Check if fixed output interval is set.
    if ( fixedOutputInterval_ != -0.0 )
    {
        // Compute number of steps required.
        numberOfPropagationSteps_ = std::floor(
                propagationIntervalEnd_ - propagationIntervalStart_ )
                / fixedOutputInterval_;
    }

    // If fixed output interval is not set, only one step is required.
    else
    {
        // Compute number of steps required.
        numberOfPropagationSteps_ = 2;

        // Set fixed output interval length of propagation interval.
        fixedOutputInterval_ = propagationIntervalEnd_
                               - propagationIntervalStart_;
    }

    // Loop over map of bodies to be propagated.
    for ( iteratorBodiesToPropagate_ =
          bodiesToPropagate_.begin( );
    iteratorBodiesToPropagate_ !=
            bodiesToPropagate_.end( );
    iteratorBodiesToPropagate_++ )
    {
        // Set propagator to this KeplerPropagator object for all bodies.
        iteratorBodiesToPropagate_
                ->second->pointerToPropagator_ = this;

        // Convert state given in Cartesian elements to Keplerian elements.
        pointerToKeplerianElements_
                = orbital_element_conversions::
                  convertCartesianToKeplerianElements(
                          static_cast< CartesianElements* > (
                                  iteratorBodiesToPropagate_
                                  ->second->pointerToInitialState_ ),
                          iteratorBodiesToPropagate_->second
                          ->pointerToCentralBody_ );

        // Check if orbit is elliptical, and not near-parabolic.
        if ( pointerToKeplerianElements_->getEccentricity( ) < 0.8 )
        {
            // Set mathematical functions.
            newtonRaphsonAdaptorForKeplerPropagator_
                    .setPointerToFunction(
                            &KeplerPropagator::
                            computeKeplerEquationForEllipticalOrbits_ );
            newtonRaphsonAdaptorForKeplerPropagator_
                    .setPointerToFirstDerivativeFunction(
                            &KeplerPropagator::
                    computeFirstDerivativeKeplerEquationForEllipticalOrbits_ );

            // Loop over all steps and store propagation history.
            for ( unsigned int i = 0; i < numberOfPropagationSteps_ - 1; i++ )
            {
                // Create pointer to CartesianElements object.
                CartesianElements* pointerToCartesianElements_
                        = new CartesianElements;

                // Compute mean anomaly using Keplerian elements for given body
                // at specific time along trajectory determined by current step
                // in for loop.
                meanAnomaly_ = sqrt( iteratorBodiesToPropagate_->second
                                     ->pointerToCentralBody_
                                     ->getGravitationalParameter( )
                                     / raiseToIntegerPower(
                                             pointerToKeplerianElements_
                                             ->getSemiMajorAxis( ), 3 ) )
                        * i * fixedOutputInterval_;

                // Set initial guess of eccentric anomaly to the computed mean
                // anomaly.
                pointerToNewtonRaphson_->setInitialGuessOfRoot( meanAnomaly_ );

                // Execute Newton-Raphon method.
                pointerToNewtonRaphson_->execute();

                // Set eccentric anomaly based on result of Newton-Raphson
                // root-finding algorithm.
                eccentricAnomaly_ = pointerToNewtonRaphson_
                                    ->getComputedRootOfFunction();

                // Compute true anomaly using eccentric anomaly.
                trueAnomaly_ = atan2( sin( eccentricAnomaly_ ) * sqrt(
                        1.0 - raiseToIntegerPower( pointerToKeplerianElements_
                                                   ->getEccentricity( ), 2 ) ),
                                      cos( eccentricAnomaly_ )
                                      - pointerToKeplerianElements_
                                      ->getEccentricity( ) );

                // Set computed true anomaly in KeplerianElements object.
                pointerToKeplerianElements_->setTrueAnomaly( trueAnomaly_ );

                // Convert state given in Keplerian elements to Cartesian
                // elements.
                pointerToCartesianElements_
                        = orbital_element_conversions::
                          convertKeplerianToCartesianElements(
                                  pointerToKeplerianElements_,
                                  iteratorBodiesToPropagate_->second
                                  ->pointerToCentralBody_ );

                // Store intermediate propagation state in propagation history
                // for given body.
                iteratorBodiesToPropagate_->second
                        ->propagationHistory_[ i * fixedOutputInterval_  ]
                        = pointerToCartesianElements_;
            }

            // Store final state in CartesianElements for given body.
            iteratorBodiesToPropagate_->second->pointerToFinalState_
                    = iteratorBodiesToPropagate_->second
                      ->propagationHistory_[ numberOfPropagationSteps_ - 1
                                             * fixedOutputInterval_ ];
        }

        // Check if orbit is near-parabolic
        else if ( pointerToKeplerianElements_->getEccentricity( ) > 0.8
                  && pointerToKeplerianElements_->getEccentricity( ) < 1.2 )
        {
            // Print cerr statement.
            std::cerr << " Orbit is near-parabolic and cannot be computed at"
                      << " present" << std::endl;
        }

        // Check if orbit is hyperbolic, and not near-parabolic
        else if ( pointerToKeplerianElements_->getEccentricity( ) > 1.2 )
        {
            // Set mathematical functions.
            newtonRaphsonAdaptorForKeplerPropagator_
                    .setPointerToFunction(
                            &KeplerPropagator::
                            computeKeplerEquationForHyperbolicOrbits_ );
            newtonRaphsonAdaptorForKeplerPropagator_
                    .setPointerToFirstDerivativeFunction(
                            &KeplerPropagator::
                    computeFirstDerivativeKeplerEquationForHyperbolicOrbits_ );

            // Loop over all steps and store propagation history.
            for ( unsigned int i = 0; i < numberOfPropagationSteps_ - 1; i++ )
            {
                // Create pointer to CartesianElements object.
                CartesianElements* pointerToCartesianElements_
                        = new CartesianElements;

                // Compute hyperbolic mean anomaly using Keplerian elements for
                // given body at specific time along trajectory determined by
                // current step in for loop.
                hyperbolicMeanAnomaly_
                        = sqrt( - iteratorBodiesToPropagate_->second
                                ->pointerToCentralBody_
                                ->getGravitationalParameter( )
                                / raiseToIntegerPower(
                                        pointerToKeplerianElements_
                                        ->getSemiMajorAxis( ), 3 ) )
                        * i * fixedOutputInterval_;

                // Check if hyperbolic mean anomaly is positive.
                if ( hyperbolicMeanAnomaly_ > 0 )
                {
                    // Set initial guess of hyperbolic eccentric anomaly.
                    pointerToNewtonRaphson_
                            ->setInitialGuessOfRoot(
                                    log( 1.8 + 2.0 * hyperbolicMeanAnomaly_
                                         / pointerToKeplerianElements_
                                         ->getEccentricity( ) ) );
                }

                else
                {
                    // Set initial guess of hyperbolic eccentric anomaly.
                    pointerToNewtonRaphson_
                            ->setInitialGuessOfRoot(
                                    -log( 1.8 - 2.0 * hyperbolicMeanAnomaly_
                                          / pointerToKeplerianElements_
                                          ->getEccentricity( ) ) );
                }

                // Execute Newton-Raphon method.
                pointerToNewtonRaphson_->execute();

                // Set hyperbolic eccentric anomaly based on result of
                // Newton-Raphson root-finding algorithm.
                hyperbolicEccentricAnomaly_ = pointerToNewtonRaphson_
                                              ->getComputedRootOfFunction();

                // Compute true anomaly using hyperbolic eccentric anomaly.
                double sineOfTrueAnomaly_
                        = sqrt( raiseToIntegerPower(
                                pointerToKeplerianElements_
                                ->getEccentricity( ), 2 ) - 1.0 )
                        * sinh( hyperbolicEccentricAnomaly_ )
                        / ( pointerToKeplerianElements_->getEccentricity( )
                            * cosh( hyperbolicEccentricAnomaly_ ) - 1.0 );

                double cosineOfTrueAnomaly_
                        = ( pointerToKeplerianElements_->getEccentricity( )
                            - cosh( hyperbolicEccentricAnomaly_ ) )
                          / ( pointerToKeplerianElements_->getEccentricity( )
                              * cosh( hyperbolicEccentricAnomaly_ ) - 1.0 );

                trueAnomaly_ = atan2( sineOfTrueAnomaly_,
                                      cosineOfTrueAnomaly_ );

                // Set computed true anomaly in KeplerianElements object.
                pointerToKeplerianElements_->setTrueAnomaly( trueAnomaly_ );

                // Convert state given in Keplerian elements to Cartesian
                // elements.
                pointerToCartesianElements_
                        = orbital_element_conversions::
                          convertKeplerianToCartesianElements(
                                  pointerToKeplerianElements_,
                                  iteratorBodiesToPropagate_->second
                                  ->pointerToCentralBody_ );
                // Store intermediate propagation state in propagation history
                // for given body.
                iteratorBodiesToPropagate_->second
                        ->propagationHistory_[ i * fixedOutputInterval_  ]
                        = pointerToCartesianElements_;
            }

            // Store final state in CartesianElements for given body.
            iteratorBodiesToPropagate_->second->pointerToFinalState_
                    = iteratorBodiesToPropagate_->second
                      ->propagationHistory_[ numberOfPropagationSteps_ - 1
                                             * fixedOutputInterval_ ];
        }
    }
}

//! Compute Kepler's equation for elliptical orbits.
double KeplerPropagator::
        computeKeplerEquationForEllipticalOrbits_( double& eccentricAnomaly )
{
    // Return value of Kepler's equation for elliptical orbits.
    return eccentricAnomaly - pointerToKeplerianElements_->getEccentricity( )
            * sin( eccentricAnomaly ) - meanAnomaly_;
}

//! Compute first-derivative of Kepler's equation for elliptical orbits.
double KeplerPropagator::
        computeFirstDerivativeKeplerEquationForEllipticalOrbits_(
                double& eccentricAnomaly )
{
    // Return value of first-derivative of Kepler's equation for elliptical
    // orbits.
    return 1.0 - pointerToKeplerianElements_->getEccentricity( )
            * cos( eccentricAnomaly );
}

//! Compute Kepler's equation for hyperbolic orbits.
double KeplerPropagator::
        computeKeplerEquationForHyperbolicOrbits_( double&
                                                 hyperbolicEccentricAnomaly )
{
    // Return value of Kepler's equation for hyperbolic orbits.
    return pointerToKeplerianElements_->getEccentricity( )
            * sinh( hyperbolicEccentricAnomaly )
            - hyperbolicEccentricAnomaly - hyperbolicMeanAnomaly_;
}

//! Compute first-derivative of Kepler's equation for hyperbolic orbits.
double KeplerPropagator::
        computeFirstDerivativeKeplerEquationForHyperbolicOrbits_(
                double& hyperbolicEccentricAnomaly )
{
    // Return value of first-derivative of Kepler's equation for hyperbolic
    // orbits.
    return pointerToKeplerianElements_->getEccentricity( )
            * cosh( hyperbolicEccentricAnomaly ) - 1.0;
}

// End of file.
