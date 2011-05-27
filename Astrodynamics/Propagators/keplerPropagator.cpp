/*! \file keplerPropagator.cpp
 *    Source file that defines the kepler propagator class included in
 *    Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 4
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
 *    Last modified     : 14 February, 2011
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
 *      110214    K. Kumar          Updated code based on new orbital
 *                                  conversion functions; optimized code.
 *      110215    E. Iorfida        Minor changes.
 */

// Include statements.
#include "keplerPropagator.h"

// Using declarations.
using std::endl;
using mathematics::raiseToIntegerPower;

//! Default constructor.
KeplerPropagator::KeplerPropagator( ) : numberOfPropagationSteps_( 0 ),
                                        trueAnomaly_( -1.0 )
{
    // Initialize variables.
    pointerToConvertMeanAnomalyToEccentricAnomaly_
            = new ConvertMeanAnomalyToEccentricAnomaly;
    pointerToConvertMeanAnomalyToHyperbolicEccentricAnomaly_
            = new ConvertMeanAnomalyToHyperbolicEccentricAnomaly;
}

//! Default destructor
KeplerPropagator::~KeplerPropagator( )
{
    // Deallocate variables.
    delete pointerToConvertMeanAnomalyToEccentricAnomaly_;
    delete pointerToConvertMeanAnomalyToHyperbolicEccentricAnomaly_;
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
    // Check if fixed output interval is set.
    if ( fixedOutputInterval_ != -0.0 )
    {
        // Compute number of steps required.
        numberOfPropagationSteps_
                = static_cast < unsigned int > ( std::floor(
                        ( propagationIntervalEnd_ - propagationIntervalStart_ )
                        / fixedOutputInterval_ ) );
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
    for ( iteratorBodiesToPropagate_ = bodiesToPropagate_.begin( );
          iteratorBodiesToPropagate_ != bodiesToPropagate_.end( );
          iteratorBodiesToPropagate_++ )
    {
        // Set propagator to this KeplerPropagator object for all bodies.
        iteratorBodiesToPropagate_
                ->second->pointerToPropagator_ = this;

        // Convert initial state given in Cartesian elements to Keplerian
        // elements.
        pointerToKeplerianElements_
                = orbital_element_conversions::
                  convertCartesianToKeplerianElements(
                          static_cast< CartesianElements* > (
                                  iteratorBodiesToPropagate_
                                  ->second->pointerToInitialState_ ),
                          iteratorBodiesToPropagate_->second
                          ->pointerToCentralBody_ );

        if ( pointerToKeplerianElements_->getEccentricity( ) < 0.8
             && pointerToKeplerianElements_->getEccentricity( ) >= 0.0 )
        {
            // Convert initial true anomaly to eccentric anomaly.
            eccentricAnomaly_ = orbital_element_conversions::
                                convertTrueAnomalyToEccentricAnomaly(
                                        pointerToKeplerianElements_
                                        ->getTrueAnomaly( ),
                                        pointerToKeplerianElements_
                                        ->getEccentricity( ) );

            // Convert initial eccentric anomaly to mean anomaly.
            meanAnomaly_ = orbital_element_conversions::
                           convertEccentricAnomalyToMeanAnomaly(
                                   eccentricAnomaly_,
                                   pointerToKeplerianElements_
                                   ->getEccentricity( ) );

            // Set Keplerian elements for mean anomaly to eccentric anomaly
            // conversion.
            pointerToConvertMeanAnomalyToEccentricAnomaly_
                    ->setEccentricity( pointerToKeplerianElements_
                                       ->getEccentricity( ) );

            // Set Newton-Raphson method.
            pointerToConvertMeanAnomalyToEccentricAnomaly_
                    ->setNewtonRaphson( pointerToNewtonRaphson_ );

            // Compute change of mean anomaly between start and end of
            // propagation step. This change is the same of each propagation
            // step since the steps are equal in time.
            meanAnomalyChange_
                    = orbital_element_conversions::
                      convertElapsedTimeToMeanAnomalyForEllipticalOrbits(
                              fixedOutputInterval_,
                              iteratorBodiesToPropagate_
                              ->second->pointerToCentralBody_,
                              pointerToKeplerianElements_
                              ->getSemiMajorAxis( ) );

            // Loop over all steps and store propagation history.
            for ( unsigned int i = 0; i < numberOfPropagationSteps_ + 1; i++ )
            {
                // Create pointer to CartesianElements object.
                CartesianElements* pointerToCartesianElements_
                        = new CartesianElements;

                // Set mean anomaly change in mean anomaly to eccentric anomaly
                // conversions
                pointerToConvertMeanAnomalyToEccentricAnomaly_
                        ->setMeanAnomaly( meanAnomaly_ );

                // Compute eccentric anomaly for mean anomaly.
                eccentricAnomaly_
                        = pointerToConvertMeanAnomalyToEccentricAnomaly_
                          ->convert( );

                // Compute true anomaly for computed eccentric anomaly.
                trueAnomaly_ = orbital_element_conversions::
                               convertEccentricAnomalyToTrueAnomaly(
                                       eccentricAnomaly_,
                                       pointerToKeplerianElements_
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
                        ->propagationHistory_[ i * fixedOutputInterval_ ]
                        = pointerToCartesianElements_;

                // Increment mean anomaly for next propagation step.
                meanAnomaly_ += meanAnomalyChange_;
            }

            // Store final state in CartesianElements for given body.
            iteratorBodiesToPropagate_->second->pointerToFinalState_
                    = iteratorBodiesToPropagate_->second
                      ->propagationHistory_[ ( numberOfPropagationSteps_ - 1 )
                                             * fixedOutputInterval_ ];
        }

        else if ( pointerToKeplerianElements_->getEccentricity( ) > 1.2 )
        {
            // Convert initial true anomaly to hyperbolic eccentric anomaly.
            hyperbolicEccentricAnomaly_
                    = orbital_element_conversions::
                      convertTrueAnomalyToHyperbolicEccentricAnomaly(
                              pointerToKeplerianElements_->getTrueAnomaly( ),
                              pointerToKeplerianElements_->getEccentricity( ) );

            // Convert initial hyperbolic eccentric anomaly to mean anomaly.
            meanAnomaly_ = orbital_element_conversions::
                           convertHyperbolicEccentricAnomalyToMeanAnomaly(
                                   hyperbolicEccentricAnomaly_,
                                   pointerToKeplerianElements_
                                   ->getEccentricity( ) );

            // Set Keplerian elements for mean anomaly to hyperbolic eccentric
            // anomaly conversion.
            pointerToConvertMeanAnomalyToHyperbolicEccentricAnomaly_
                    ->setEccentricity( pointerToKeplerianElements_
                                       ->getEccentricity( ) );

            // Set Newton-Raphson method.
            pointerToConvertMeanAnomalyToHyperbolicEccentricAnomaly_
                    ->setNewtonRaphson( pointerToNewtonRaphson_ );

            // Compute change of mean anomaly between start and end of
            // propagation step. This change is the same of each propagation
            // step since the steps are equal in time.
            meanAnomalyChange_
                    = orbital_element_conversions::
                      convertElapsedTimeToMeanAnomalyForHyperbolicOrbits(
                              fixedOutputInterval_,
                              iteratorBodiesToPropagate_
                              ->second->pointerToCentralBody_,
                              pointerToKeplerianElements_
                              ->getSemiMajorAxis( ) );

            // Loop over all steps and store propagation history.
            for ( unsigned int i = 0; i < numberOfPropagationSteps_ + 1; i++ )
            {
                // Create pointer to CartesianElements object.
                CartesianElements* pointerToCartesianElements_
                        = new CartesianElements;

                // Set mean anomaly change in mean anomaly to eccentric anomaly
                // conversions
                pointerToConvertMeanAnomalyToHyperbolicEccentricAnomaly_
                        ->setMeanAnomaly( meanAnomaly_ );

                // Compute hyperbolic eccentric anomaly for mean anomaly.
                hyperbolicEccentricAnomaly_ =
                        pointerToConvertMeanAnomalyToHyperbolicEccentricAnomaly_
                        ->convert( );

                // Compute true anomaly for computed hyperbolic eccentric anomaly.
                trueAnomaly_ = orbital_element_conversions::
                               convertHyperbolicEccentricAnomalyToTrueAnomaly(
                                       hyperbolicEccentricAnomaly_,
                                       pointerToKeplerianElements_
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
                        ->propagationHistory_[ i * fixedOutputInterval_ ]
                        = pointerToCartesianElements_;

                // Increment mean anomaly for next propagation step.
                meanAnomaly_ += meanAnomalyChange_;
            }

            // Store final state in CartesianElements for given body.
            iteratorBodiesToPropagate_->second->pointerToFinalState_
                    = iteratorBodiesToPropagate_->second
                      ->propagationHistory_[ ( numberOfPropagationSteps_ - 1 )
                                             * fixedOutputInterval_ ];
        }
    }
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream,
                          KeplerPropagator& keplerPropagator )
{
    stream << "This is a KeplerPropagator object." << endl;
    stream << "The start of the propagation interval is set to: " << endl;
    stream << keplerPropagator.getPropagationIntervalStart( ) << endl;
    stream << "The end of the propagation interval is set to: " << endl;
    stream << keplerPropagator.getPropagationIntervalEnd( ) << endl;

    // Return stream.
    return stream;
}

// End of file.
