/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *      120217    K. Kumar          Updated computeModuloForSignedValues() to computeModulo() from
 *                                  Tudat Core.
 *
 *    References
 *
 */

#include <boost/exception/all.hpp>
#include <cmath>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Mathematics/BasicMathematics/basicMathematicsFunctions.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>
#include "Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToEccentricAnomaly.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphson.h"

namespace tudat
{
namespace orbital_element_conversions
{

//! Propagate Kepler orbit.
Eigen::VectorXd propagateKeplerOrbit( const Eigen::VectorXd& initialStateInKeplerianElements,
                                      const double epochOfInitialState,
                                      const double epochOfFinalState,
                                      const double centralBodyGravitationalParameter,
                                      const double newtonRaphsonConvergenceTolerance,
                                      bool useModuloOption )
{
    // Create Newton-Raphson root-finder.
    tudat::NewtonRaphson newtonRaphson_;
    newtonRaphson_.setTolerance( newtonRaphsonConvergenceTolerance );

    // Create final state in Keplerian elements.
    Eigen::VectorXd finalStateInKeplerianElements = initialStateInKeplerianElements;

    // Check if orbit is elliptical.
    if ( initialStateInKeplerianElements( eccentricityIndex ) >= 0.98
         || initialStateInKeplerianElements( eccentricityIndex ) < 0.0 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid." ) ) );
    }

    // Check if orbit is elliptical.
    if ( initialStateInKeplerianElements( eccentricityIndex ) < 0.98
         && initialStateInKeplerianElements( eccentricityIndex ) >= 0.0 )
    {
        // Set elapsed time used in the computation algorithm. This elapsed time is either
        // the complete propagation period or a fraction of an orbit.
        double elapsedTime_ = 0.0;

        // Set number of complete orbits.
        double numberOfCompleteOrbits_ = 0.0;

        // Check if modulo-option is set and subtract complete orbits from elapsed time.
        if ( useModuloOption )
        {
            // Compute orbital period of Kepler orbit.
            double orbitalPeriod_ = astrodynamics::computeKeplerOrbitalPeriod(
                        initialStateInKeplerianElements( semiMajorAxisIndex ),
                        centralBodyGravitationalParameter );

            // Determine elapsed time.
            elapsedTime_ = mathematics::computeModulo(
                        ( epochOfFinalState - epochOfInitialState ), orbitalPeriod_ );

            // Determine corresponding number of complete orbits.
            numberOfCompleteOrbits_ = std::floor(
                        ( epochOfFinalState - epochOfInitialState ) / orbitalPeriod_ );
        }

        else
        {
            elapsedTime_ = ( epochOfFinalState - epochOfInitialState );
        }

        // Convert initial true anomaly to eccentric anomaly.
        double initialEccentricAnomaly_ = convertTrueAnomalyToEccentricAnomaly(
                    initialStateInKeplerianElements( trueAnomalyIndex ),
                    initialStateInKeplerianElements( eccentricityIndex ) );

        // Convert initial eccentric anomaly to mean anomaly.
        double initialMeanAnomaly_ = convertEccentricAnomalyToMeanAnomaly(
                    initialEccentricAnomaly_,
                    initialStateInKeplerianElements( eccentricityIndex ) );

        // Set Keplerian elements for mean anomaly to eccentric anomaly
        // conversion.
        ConvertMeanAnomalyToEccentricAnomaly convertMeanAnomalyToEccentricAnomaly_;
        convertMeanAnomalyToEccentricAnomaly_.setEccentricity(
                    initialStateInKeplerianElements( eccentricityIndex ) );

        // Set Newton-Raphson method.
        convertMeanAnomalyToEccentricAnomaly_.setNewtonRaphson( &newtonRaphson_ );

        // Compute change of mean anomaly between start and end of propagation.
        double meanAnomalyChange_ = convertElapsedTimeToEllipticalMeanAnomalyChange(
                    elapsedTime_, centralBodyGravitationalParameter,
                    initialStateInKeplerianElements( semiMajorAxisIndex ) );

        // Set mean anomaly change in mean anomaly to eccentric anomaly
        // conversions.
        convertMeanAnomalyToEccentricAnomaly_.setMeanAnomaly( initialMeanAnomaly_
                                                              + meanAnomalyChange_ );

        // Compute eccentric anomaly for mean anomaly.
        double finalEccentricAnomaly_ = convertMeanAnomalyToEccentricAnomaly_.convert( );

        // Compute true anomaly for computed eccentric anomaly.
        finalStateInKeplerianElements( trueAnomalyIndex )
                = convertEccentricAnomalyToTrueAnomaly(
                    finalEccentricAnomaly_,
                    finalStateInKeplerianElements( eccentricityIndex ) );

        // If modulo option is set, add the removed completed orbits back to the final true
        // anomaly.
        if ( useModuloOption )
        {
            finalStateInKeplerianElements( trueAnomalyIndex ) += numberOfCompleteOrbits_
                    * 2.0 * mathematics::PI;
        }
    }

    return finalStateInKeplerianElements;
}

} // namespace orbital_element_conversions
} // namespace tudat

// This code has been commented out, as it was implemented in the old Kepler propagator without an
// accomanpying unit test.
// It can be re-included in the propagateKeplerOrbit() function once a unit test is written.
// The code will also have to be updated to no longer use the old Propagator architecture.

//        else if ( keplerianElements_.getEccentricity( ) > 1.2 )
//        {
//            // Convert initial true anomaly to hyperbolic eccentric anomaly.
//            hyperbolicEccentricAnomaly_
//                    = orbital_element_conversions::
//                      convertTrueAnomalyToHyperbolicEccentricAnomaly(
//                              keplerianElements_.getTrueAnomaly( ),
//                              keplerianElements_.getEccentricity( ) );

//            // Convert initial hyperbolic eccentric anomaly to mean anomaly.
//            meanAnomaly_ = orbital_element_conversions::
//                           convertHyperbolicEccentricAnomalyToMeanAnomaly(
//                                   hyperbolicEccentricAnomaly_,
//                                   keplerianElements_.getEccentricity( ) );

//            // Set Keplerian elements for mean anomaly to hyperbolic eccentric
//            // anomaly conversion.
//            convertMeanAnomalyToHyperbolicEccentricAnomaly_
//                    .setEccentricity( keplerianElements_.getEccentricity( ) );

//            // Set Newton-Raphson method.
//            convertMeanAnomalyToHyperbolicEccentricAnomaly_
//                    .setNewtonRaphson( pointerToNewtonRaphson_ );

//            // Compute change of mean anomaly between start and end of
//            // propagation step.
//            meanAnomalyChange_
//                    = orbital_element_conversions::
//                    convertElapsedTimeToHyperbolicMeanAnomalyChange(
//                        ( propagationIntervalEnd_ - propagationIntervalStart_ ),
//                        iteratorBodiesToPropagate_->second
//                        .pointerToCentralBody->getGravitationalParameter( ),
//                        keplerianElements_.getSemiMajorAxis( ) );

//            // Set mean anomaly change in mean anomaly to eccentric anomaly
//            // conversions
//            convertMeanAnomalyToHyperbolicEccentricAnomaly_
//                    .setMeanAnomaly( meanAnomaly_ + meanAnomalyChange_ );

//            // Compute hyperbolic eccentric anomaly for mean anomaly.
//            hyperbolicEccentricAnomaly_ =
//                    convertMeanAnomalyToHyperbolicEccentricAnomaly_
//                    .convert( );

//            // Compute true anomaly for computed hyperbolic eccentric
//            // anomaly.
//            trueAnomaly_ = orbital_element_conversions::
//                           convertHyperbolicEccentricAnomalyToTrueAnomaly(
//                                   hyperbolicEccentricAnomaly_,
//                                   keplerianElements_.getEccentricity( ) );

//            // Set computed true anomaly in KeplerianElements object.
//            keplerianElements_.setTrueAnomaly( trueAnomaly_ );
//        }
