/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/astro/basic_astro/astrodynamicsFunctions.h"
#include "tudat/astro/basic_astro/convertMeanToEccentricAnomalies.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"

namespace tudat
{

namespace ephemerides
{

//! Class constructor.
KeplerEphemeris::KeplerEphemeris(
        const Eigen::Vector6d& initialStateInKeplerianElements,
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
    if( !( centralBodyGravitationalParameter_ > 0.0 ) )
    {
        throw std::runtime_error( "Error when creating Kepler ephemeris, gravitational parameter must be larger than 0, provided value is " +
                                  std::to_string( centralBodyGravitationalParameter_ ) );
    }
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::root_finders;
    

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
        throw std::runtime_error(
            "Error, Kepler ephemeris cannot handle parabolic orbit" );
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
    rootFinder_ = root_finders::createRootFinder(
                        root_finders::newtonRaphsonRootFinderSettings(
                            TUDAT_NAN, rootFinderAbsoluteTolerance, TUDAT_NAN, rootFinderMaximumNumberOfIterations,
                    root_finders::throw_exception ) );
}

//! Function to get state from ephemeris.
Eigen::Vector6d KeplerEphemeris::getCartesianState(
        const double secondsSinceEpoch )
{
    using namespace tudat::orbital_element_conversions;

    Eigen::Vector6d currentCartesianState = Eigen::Vector6d::Zero( );

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
