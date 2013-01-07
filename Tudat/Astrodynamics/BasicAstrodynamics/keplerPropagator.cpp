/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
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
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120607    P. Musegaas       Changed interface (propagation time instead of two epochs).
 *      120713    P. Musegaas       Changed tolerance in root finder to relative tolerance.
 *      120813    P. Musegaas       Changed code to new root finding structure.
 *
 *    References
 *
 */

#include <cmath>
#include <stdexcept>

#include <boost/exception/all.hpp>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h>
#include <TudatCore/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h>
#include <TudatCore/Mathematics/BasicMathematics/basicMathematicsFunctions.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToEccentricAnomaly.h"

namespace tudat
{
namespace basic_astrodynamics
{
namespace orbital_element_conversions
{

using namespace root_finders;

//! Propagate Kepler orbit.
Eigen::VectorXd propagateKeplerOrbit( const Eigen::VectorXd& initialStateInKeplerianElements,
                                      const double propagationTime,
                                      const double centralBodyGravitationalParameter,
                                      bool useModuloOption,
                                      RootFinderPointer rootFinder )
{
    if ( !rootFinder.get( ) )
    {
        rootFinder = boost::make_shared< NewtonRaphson >(
                    boost::bind( &root_finders::termination_conditions::
                                 RootAbsoluteToleranceTerminationCondition::
                                 checkTerminationCondition,
                                 boost::make_shared< root_finders::termination_conditions::
                                 RootAbsoluteToleranceTerminationCondition >( 5.0e-15, 1000 ),
                                 _1, _2, _3, _4, _5 ) );
    }

    using tudat::orbital_element_conversions::eccentricityIndex;
    using tudat::orbital_element_conversions::semiMajorAxisIndex;
    using tudat::orbital_element_conversions::trueAnomalyIndex;

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
        double elapsedTime = 0.0;

        // Set number of complete orbits.
        double numberOfCompleteOrbits = 0.0;

        // Check if modulo-option is set and subtract complete orbits from elapsed time.
        if ( useModuloOption )
        {
            // Compute orbital period of Kepler orbit.
            const double orbitalPeriod = basic_astrodynamics::computeKeplerOrbitalPeriod(
                        initialStateInKeplerianElements( semiMajorAxisIndex ),
                        centralBodyGravitationalParameter );

            // Determine elapsed time.
            elapsedTime = basic_mathematics::computeModulo( propagationTime, orbitalPeriod );

            // Determine corresponding number of complete orbits.
            numberOfCompleteOrbits = std::floor( propagationTime / orbitalPeriod );
        }

        else
        {
            elapsedTime = ( propagationTime );
        }

        // Convert initial true anomaly to eccentric anomaly.
        const double initialEccentricAnomaly
                = tudat::orbital_element_conversions::convertTrueAnomalyToEccentricAnomaly(
                    initialStateInKeplerianElements( trueAnomalyIndex ),
                    initialStateInKeplerianElements( eccentricityIndex ) );

        // Convert initial eccentric anomaly to mean anomaly.
        const double initialMeanAnomaly
                = tudat::orbital_element_conversions::convertEccentricAnomalyToMeanAnomaly(
                    initialEccentricAnomaly,
                    initialStateInKeplerianElements( eccentricityIndex ) );

        // Compute change of mean anomaly between start and end of propagation.
        const double meanAnomalyChange
                = tudat::orbital_element_conversions::
                convertElapsedTimeToEllipticalMeanAnomalyChange(
                    elapsedTime, centralBodyGravitationalParameter,
                    initialStateInKeplerianElements( semiMajorAxisIndex ) );

        // Set Keplerian elements and Newton-Raphson root-finder for mean anomaly to eccentric
        // anomaly conversion.
        ConvertMeanAnomalyToEccentricAnomaly convertMeanAnomalyToEccentricAnomaly_(
                    initialStateInKeplerianElements( eccentricityIndex ),
                    initialMeanAnomaly + meanAnomalyChange,
                    rootFinder );

        // Compute eccentric anomaly for mean anomaly.
        const double finalEccentricAnomaly = convertMeanAnomalyToEccentricAnomaly_.convert( );

        // Compute true anomaly for computed eccentric anomaly.
        finalStateInKeplerianElements( trueAnomalyIndex )
                = tudat::orbital_element_conversions::convertEccentricAnomalyToTrueAnomaly(
                    finalEccentricAnomaly, finalStateInKeplerianElements( eccentricityIndex ) );

        // If modulo option is set, add the removed completed orbits back to the final true
        // anomaly.
        if ( useModuloOption )
        {
            finalStateInKeplerianElements( trueAnomalyIndex )
                    += numberOfCompleteOrbits * 2.0
                    * basic_mathematics::mathematical_constants::PI;
        }
    }

    return finalStateInKeplerianElements;
}

} // namespace orbital_element_conversions
} // namespace basic_astrodynamics
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
