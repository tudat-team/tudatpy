/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      120823    P. Musegaas       Added functionality for hyperbolic and near-parabolic orbits.
 *                                  Changed some parameters to const. Various small changes.
 *      120903    P. Musegaas       Removed modulo option, due to errors with it. Kepler propagator
 *                                  now simply return true anomaly in -PI to PI spectrum.
 *      121205    P. Musegaas       Updated code to final version of rootfinders.
 *
 *    References
 *
 *    Notes
 *
 */

#include <cmath>
#include <stdexcept>

#include <boost/exception/all.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/keplerPropagator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToEccentricAnomaly.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanAnomalyToHyperbolicEccentricAnomaly.h"

namespace tudat
{
\
namespace orbital_element_conversions
{

using namespace root_finders;

//! Propagate Kepler orbit.
basic_mathematics::Vector6d propagateKeplerOrbit(
        const basic_mathematics::Vector6d &initialStateInKeplerianElements,
        const double propagationTime,
        const double centralBodyGravitationalParameter,
        RootFinderPointer rootFinder )
{
    if ( !rootFinder.get( ) )
    {
        rootFinder = boost::make_shared< NewtonRaphson >(
                    boost::bind( &root_finders::termination_conditions::
                                 RootAbsoluteToleranceTerminationCondition::
                                 checkTerminationCondition,
                                 boost::make_shared< root_finders::termination_conditions::
                                 RootAbsoluteToleranceTerminationCondition >( 5.0e-14, 1000 ),
                                 _1, _2, _3, _4, _5 ) );
    }

    // Create final state in Keplerian elements.
    Eigen::VectorXd finalStateInKeplerianElements = initialStateInKeplerianElements;

    // Check if eccentricity is valid.
    if ( initialStateInKeplerianElements( eccentricityIndex ) < 0.0 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid (smaller than 0)." ) ) );
    }

    // Check if orbit is elliptical.
    else if ( initialStateInKeplerianElements( eccentricityIndex ) < 1.0 )
    {
        // Convert initial true anomaly to eccentric anomaly.
        const double initialEccentricAnomaly = convertTrueAnomalyToEccentricAnomaly(
                    initialStateInKeplerianElements( trueAnomalyIndex ),
                    initialStateInKeplerianElements( eccentricityIndex ) );

        // Convert initial eccentric anomaly to mean anomaly.
        const double initialMeanAnomaly = convertEccentricAnomalyToMeanAnomaly(
                    initialEccentricAnomaly,
                    initialStateInKeplerianElements( eccentricityIndex ) );

        // Compute change of mean anomaly between start and end of propagation.
        const double meanAnomalyChange = convertElapsedTimeToEllipticalMeanAnomalyChange(
                    propagationTime, centralBodyGravitationalParameter,
                    initialStateInKeplerianElements( semiMajorAxisIndex ) );

        // Set Keplerian elements and Newton-Raphson root-finder for mean anomaly to eccentric
        // anomaly conversion.
        ConvertMeanAnomalyToEccentricAnomaly convertMeanAnomalyToEccentricAnomaly(
                    initialStateInKeplerianElements( eccentricityIndex ),
                    initialMeanAnomaly + meanAnomalyChange,
                    rootFinder );

        // Compute eccentric anomaly for mean anomaly.
        const double finalEccentricAnomaly = convertMeanAnomalyToEccentricAnomaly.convert( );

        // Compute true anomaly for computed eccentric anomaly.
        finalStateInKeplerianElements( trueAnomalyIndex ) = convertEccentricAnomalyToTrueAnomaly(
                    finalEccentricAnomaly, finalStateInKeplerianElements( eccentricityIndex ) );
    }

    // Check if the orbit is hyperbolic.
    else if ( initialStateInKeplerianElements( eccentricityIndex ) > 1.0 )
    {
        // Convert initial true anomaly to hyperbolic eccentric anomaly.
        const double initialHyperbolicEccentricAnomaly =
                convertTrueAnomalyToHyperbolicEccentricAnomaly(
                    initialStateInKeplerianElements( trueAnomalyIndex ),
                    initialStateInKeplerianElements( eccentricityIndex ) );

        // Convert initial hyperbolic eccentric anomaly to the hyperbolic mean anomaly.
        const double initialHyperbolicMeanAnomaly =
                convertHyperbolicEccentricAnomalyToMeanAnomaly(
                    initialHyperbolicEccentricAnomaly,
                    initialStateInKeplerianElements( eccentricityIndex ) );

        // Compute change of mean anomaly because of the propagation time.
        const double hyperbolicMeanAnomalyChange =
                convertElapsedTimeToHyperbolicMeanAnomalyChange(
                    propagationTime, centralBodyGravitationalParameter,
                    initialStateInKeplerianElements( semiMajorAxisIndex ) );

        // Set Keplerian elements and Newton-Raphson root-finder for mean anomaly to eccentric
        // anomaly conversion.
        ConvertMeanAnomalyToHyperbolicEccentricAnomaly
                convertMeanAnomalyToHyperbolicEccentricAnomaly(
                    initialStateInKeplerianElements( eccentricityIndex ),
                    initialHyperbolicMeanAnomaly + hyperbolicMeanAnomalyChange,
                    rootFinder );

        // Compute hyperbolic eccentric anomaly for mean anomaly.
        const double finalHyperbolicEccentricAnomaly =
                convertMeanAnomalyToHyperbolicEccentricAnomaly.convert( );

        // Compute true anomaly for computed hyperbolic eccentric anomaly.
        finalStateInKeplerianElements( trueAnomalyIndex ) =
                convertHyperbolicEccentricAnomalyToTrueAnomaly(
                               finalHyperbolicEccentricAnomaly,
                               finalStateInKeplerianElements( eccentricityIndex ) );
    }

    // In this case the eccentricity has to be 1.0, hence the orbit is parabolic.
    else
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Parabolic orbits are not (yet) supported." ) ) );
    }

    return finalStateInKeplerianElements;
}

} // namespace orbital_element_conversions
} // namespace tudat
