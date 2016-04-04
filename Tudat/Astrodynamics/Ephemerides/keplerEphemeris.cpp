/*    Copyright (c) 2010-2014, Delft University of Technology
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
 *      150411    D. Dirkx          Migrated and updated from personal code.
 *
 *    References
 *
 *    Notes
 *
 */

#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/astrodynamicsFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/convertMeanToEccentricAnomalies.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

namespace tudat
{

namespace ephemerides
{

//! Class constructor.
KeplerEphemeris::KeplerEphemeris(
        const basic_mathematics::Vector6d& initialStateInKeplerianElements,
        const double epochOfInitialState,
        const double centralBodyGravitationalParameter,
        const std::string& referenceFrameOrigin,
        const std::string& referenceFrameOrientation,
        const double rootFinderAbsoluteTolerance,
        const double rootFinderMaximumNumberOfIterations ):
    Ephemeris( referenceFrameOrigin, referenceFrameOrientation ),
    initialStateInKeplerianElements_( initialStateInKeplerianElements ),
    epochOfInitialState_( epochOfInitialState ),
    centralBodyGravitationalParameter_( centralBodyGravitationalParameter )
{
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::root_finders;
    using namespace root_finders::termination_conditions;

    // Check whether orbit is elliptical or hyperbolic (parabola not supported).
    if( initialStateInKeplerianElements( eccentricityIndex ) < 1.0 )
    {
        isOrbitHyperbolic_ = 0;
    }
    else if( initialStateInKeplerianElements( eccentricityIndex ) > 1.0 )
    {
        isOrbitHyperbolic_ = 1;
    }
    else
    {
        boost::throw_exception( std::runtime_error( boost::str( boost::format(
            "Error, Kepler ephemeris cannot handle parabolic orbit" ) ) ) );
    }

    // Convert initial true anomaly to mean anomaly.
    if( !isOrbitHyperbolic_ )
    {
        initialMeanAnomaly_ = convertEccentricAnomalyToMeanAnomaly(
                    convertTrueAnomalyToEccentricAnomaly(
                        initialStateInKeplerianElements_( trueAnomalyIndex ),
                        initialStateInKeplerianElements_( eccentricityIndex ) ),
                    initialStateInKeplerianElements_( eccentricityIndex ) );
    }
    else
    {
        initialMeanAnomaly_ = convertHyperbolicEccentricAnomalyToMeanAnomaly(
                    convertTrueAnomalyToHyperbolicEccentricAnomaly(
                        initialStateInKeplerianElements_( trueAnomalyIndex ),
                        initialStateInKeplerianElements_( eccentricityIndex ) ),
                    initialStateInKeplerianElements_( eccentricityIndex ) );
    }


    // Calculate ancilliary variables for conversion to Cartesian elements
    eccentricity_ = initialStateInKeplerianElements_( eccentricityIndex );
    semiMajorAxis_ = initialStateInKeplerianElements_( semiMajorAxisIndex );
    semiLatusRectum_ = initialStateInKeplerianElements_( semiMajorAxisIndex ) *
            ( 1.0 - pow( initialStateInKeplerianElements_( eccentricityIndex ), 2 ) );

    rotationFromOrbitalPlane_ =
            Eigen::AngleAxisd(
                initialStateInKeplerianElements_(
                    longitudeOfAscendingNodeIndex ), Eigen::Vector3d::UnitZ( ) ) *
            Eigen::AngleAxisd(
                initialStateInKeplerianElements_(
                    inclinationIndex ), Eigen::Vector3d::UnitX( ) ) *
            Eigen::AngleAxisd(
                initialStateInKeplerianElements_(
                    argumentOfPeriapsisIndex ), Eigen::Vector3d::UnitZ( ) ) ;

    // Create root finder to be used for converting mean to eccentric anomaly
    rootFinder_ = boost::make_shared< NewtonRaphsonCore< double > >(
                boost::bind( &RootAbsoluteToleranceTerminationCondition< >::checkTerminationCondition,
                             boost::make_shared< RootAbsoluteToleranceTerminationCondition< > >(
                                 rootFinderAbsoluteTolerance, rootFinderMaximumNumberOfIterations ),
                             _1, _2, _3, _4, _5 ) );
}

//! Function to get state from ephemeris.
basic_mathematics::Vector6d KeplerEphemeris::getCartesianStateFromEphemeris(
        const double secondsSinceEpoch, const double julianDayAtEpoch )
{
    using namespace tudat::orbital_element_conversions;

    basic_mathematics::Vector6d currentCartesianState = basic_mathematics::Vector6d::Zero( );

    double propagationTime = secondsSinceEpoch - epochOfInitialState_;

    // Calculate eccentric anomaly at epoch.
    double eccentricAnomaly = TUDAT_NAN;

    // Check if the orbit is hyperbolic.
    if( !isOrbitHyperbolic_ )
    {

        // Compute change of mean anomaly between start and end of propagation.
        const double meanAnomalyChange =
                convertElapsedTimeToEllipticalMeanAnomalyChange(
                    propagationTime, centralBodyGravitationalParameter_, semiMajorAxis_ );

        // Compute eccentric anomaly for mean anomaly.
        eccentricAnomaly =
                convertMeanAnomalyToEccentricAnomaly(
                    eccentricity_,
                    initialMeanAnomaly_ + meanAnomalyChange );
    }
    else
    {
        // Compute change of mean anomaly because of the propagation time.
        const double hyperbolicMeanAnomalyChange =
                convertElapsedTimeToHyperbolicMeanAnomalyChange(
                    propagationTime, centralBodyGravitationalParameter_, semiMajorAxis_ );

        // Compute hyperbolic eccentric anomaly for mean anomaly.
        eccentricAnomaly =
                convertMeanAnomalyToHyperbolicEccentricAnomaly(
                    eccentricity_,
                    initialMeanAnomaly_ + hyperbolicMeanAnomalyChange );
    }

    // Calculate true anomaly.
    double trueAnomaly = convertEccentricAnomalyToTrueAnomaly(
                eccentricAnomaly, eccentricity_ );
    double cosineOfTrueAnomaly = cos( trueAnomaly );
    double sineOfTrueAnomaly = sin( trueAnomaly );

    // Definition of position in the perifocal coordinate system.
    currentCartesianState( 0 ) = semiLatusRectum_ * cosineOfTrueAnomaly
            / ( 1.0 + eccentricity_ * cosineOfTrueAnomaly );
    currentCartesianState( 1 ) = semiLatusRectum_ * sineOfTrueAnomaly
            / ( 1.0 + eccentricity_ * cosineOfTrueAnomaly );

    // Definition of velocity in the perifocal coordinate system.
    currentCartesianState( 3 ) =
            -sqrt( centralBodyGravitationalParameter_ / semiLatusRectum_ ) * sineOfTrueAnomaly;
    currentCartesianState( 4 ) =
            sqrt( centralBodyGravitationalParameter_ / semiLatusRectum_ )
            * ( eccentricity_ + cosineOfTrueAnomaly );

    // Rotate orbital plane to correct orientation.
    currentCartesianState.segment( 0, 3 ) = rotationFromOrbitalPlane_ *
            currentCartesianState.segment( 0, 3 );
    currentCartesianState.segment( 3, 3 ) = rotationFromOrbitalPlane_ *
            currentCartesianState.segment( 3, 3 );

    return currentCartesianState;
}

} // namespace ephemerides
} // namespace tudat
