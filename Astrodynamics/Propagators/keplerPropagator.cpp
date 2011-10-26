/*! \file keplerPropagator.cpp
 *    Source file that defines the kepler propagator class included in Tudat.
 *
 *    Path              : /Astrodynamics/Propagators/
 *    Version           : 5
 *    Check status      : Unchecked
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
 *    Last modified     : 20 September, 2011
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
 *      110203    K. Kumar          File created.
 *      110207    E. Iorfida        Minor changes.
 *      110214    K. Kumar          Updated code based on new orbital conversion functions;
 *                                  optimized code.
 *      110215    E. Iorfida        Minor changes.
 *      110920    K. Kumar          Corrected simple errors outlined by M. Persson.
 */

// Include statements.
#include "Astrodynamics/Propagators/keplerPropagator.h"
#include "Astrodynamics/States/orbitalElementConversions.h"

//! Tudat library namespace.
namespace tudat
{

//! Propagate.
void KeplerPropagator::propagate( )
{
    // Using declarations.
    using std::cerr;
    using std::endl;

    // Loop over map of bodies to be propagated.
    for ( BodyPropagatorDataMap::iterator iteratorBodiesToPropagate_ = bodiesToPropagate_.begin( );
          iteratorBodiesToPropagate_ != bodiesToPropagate_.end( ); iteratorBodiesToPropagate_++ )
    {
        // Convert initial state given in Cartesian elements to Keplerian
        // elements.
        keplerianElements_
                = orbital_element_conversions::
                convertCartesianToKeplerianElements( static_cast< CartesianElements* >
                                                     ( iteratorBodiesToPropagate_
                                                       ->second.pointerToInitialState ),
                                                     iteratorBodiesToPropagate_
                                                     ->second.pointerToCentralBody );


        if ( keplerianElements_.getEccentricity( ) < 0.8
             && keplerianElements_.getEccentricity( ) >= 0.0 )
        {
            // Convert initial true anomaly to eccentric anomaly.
            eccentricAnomaly_ = orbital_element_conversions::
                                convertTrueAnomalyToEccentricAnomaly(
                                        keplerianElements_.getTrueAnomaly( ),
                                        keplerianElements_.getEccentricity( ) );

            // Convert initial eccentric anomaly to mean anomaly.
            meanAnomaly_ = orbital_element_conversions::
                           convertEccentricAnomalyToMeanAnomaly(
                                   eccentricAnomaly_,
                                   keplerianElements_.getEccentricity( ) );

            // Set Keplerian elements for mean anomaly to eccentric anomaly
            // conversion.
            convertMeanAnomalyToEccentricAnomaly_.setEccentricity(
                    keplerianElements_.getEccentricity( ) );

            // Set Newton-Raphson method.
            convertMeanAnomalyToEccentricAnomaly_.
                    setNewtonRaphson( pointerToNewtonRaphson_ );

            // Compute change of mean anomaly between start and end of
            // propagation step.
            meanAnomalyChange_
                    = orbital_element_conversions::
                      convertElapsedTimeToMeanAnomalyForEllipticalOrbits(
                              ( propagationIntervalEnd_ - propagationIntervalStart_ ),
                              iteratorBodiesToPropagate_->second.pointerToCentralBody,
                              keplerianElements_.getSemiMajorAxis( ) );

            // Set mean anomaly change in mean anomaly to eccentric anomaly
            // conversions.
            convertMeanAnomalyToEccentricAnomaly_
                    .setMeanAnomaly( meanAnomaly_ + meanAnomalyChange_ );

            // Compute eccentric anomaly for mean anomaly.
            eccentricAnomaly_
                    = convertMeanAnomalyToEccentricAnomaly_.convert( );

            // Compute true anomaly for computed eccentric anomaly.
            trueAnomaly_ = orbital_element_conversions::
                           convertEccentricAnomalyToTrueAnomaly(
                                   eccentricAnomaly_,
                                   keplerianElements_.getEccentricity( ) );

            // Set computed true anomaly in KeplerianElements object.
            keplerianElements_.setTrueAnomaly( trueAnomaly_ );
        }

        else if ( keplerianElements_.getEccentricity( ) > 1.2 )
        {
            // Convert initial true anomaly to hyperbolic eccentric anomaly.
            hyperbolicEccentricAnomaly_
                    = orbital_element_conversions::
                      convertTrueAnomalyToHyperbolicEccentricAnomaly(
                              keplerianElements_.getTrueAnomaly( ),
                              keplerianElements_.getEccentricity( ) );

            // Convert initial hyperbolic eccentric anomaly to mean anomaly.
            meanAnomaly_ = orbital_element_conversions::
                           convertHyperbolicEccentricAnomalyToMeanAnomaly(
                                   hyperbolicEccentricAnomaly_,
                                   keplerianElements_.getEccentricity( ) );

            // Set Keplerian elements for mean anomaly to hyperbolic eccentric
            // anomaly conversion.
            convertMeanAnomalyToHyperbolicEccentricAnomaly_
                    .setEccentricity( keplerianElements_.getEccentricity( ) );

            // Set Newton-Raphson method.
            convertMeanAnomalyToHyperbolicEccentricAnomaly_
                    .setNewtonRaphson( pointerToNewtonRaphson_ );

            // Compute change of mean anomaly between start and end of
            // propagation step.
            meanAnomalyChange_
                    = orbital_element_conversions::
                      convertElapsedTimeToMeanAnomalyForHyperbolicOrbits(
                              ( propagationIntervalEnd_ - propagationIntervalStart_ ),
                              iteratorBodiesToPropagate_->second.pointerToCentralBody,
                              keplerianElements_.getSemiMajorAxis( ) );

            // Set mean anomaly change in mean anomaly to eccentric anomaly
            // conversions
            convertMeanAnomalyToHyperbolicEccentricAnomaly_
                    .setMeanAnomaly( meanAnomaly_ + meanAnomalyChange_ );

            // Compute hyperbolic eccentric anomaly for mean anomaly.
            hyperbolicEccentricAnomaly_ =
                    convertMeanAnomalyToHyperbolicEccentricAnomaly_
                    .convert( );

            // Compute true anomaly for computed hyperbolic eccentric
            // anomaly.
            trueAnomaly_ = orbital_element_conversions::
                           convertHyperbolicEccentricAnomalyToTrueAnomaly(
                                   hyperbolicEccentricAnomaly_,
                                   keplerianElements_.getEccentricity( ) );

            // Set computed true anomaly in KeplerianElements object.
            keplerianElements_.setTrueAnomaly( trueAnomaly_ );
        }

        else
        {
            cerr << "There is currently no valid implementation " << endl;
            cerr << "of the Kepler propagator for eccentricities " << endl;
            cerr << "in the range: 0.8 <= eccentricity <= 1.2" << endl;
        }

        // Convert Keplerian elements to Cartesian elements.
        cartesianElements_ = orbital_element_conversions::
                             convertKeplerianToCartesianElements( &keplerianElements_,
                                     iteratorBodiesToPropagate_->second.pointerToCentralBody );

        // Store final state in CartesianElements for given body.
        iteratorBodiesToPropagate_->second.finalState = cartesianElements_;
    }
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream, KeplerPropagator& keplerPropagator )
{
    // Using declarations.
    using std::endl;

    stream << "This is a KeplerPropagator object." << endl;
    stream << "The start of the propagation interval is set to: " << endl;
    stream << keplerPropagator.getPropagationIntervalStart( ) << endl;
    stream << "The end of the propagation interval is set to: " << endl;
    stream << keplerPropagator.getPropagationIntervalEnd( ) << endl;

    // Return stream.
    return stream;
}

}

// End of file.
