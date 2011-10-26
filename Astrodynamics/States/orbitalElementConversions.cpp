/*! \file orbitalElementConversions.cpp
 *    This source file contains a namespace with orbital element conversion functions.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 18
 *    Check status      : Checked
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Date created      : 20 October, 2010
 *    Last modified     : 5 August, 2011
 *
 *    References
 *      Wertz, J. R. Mission geometry; orbit and constellation design and management.
 *      Mengali, G., Quarta, A.A. Fondamenti di meccanica del volo spaziale
 *
 *    Notes
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
 *      101020    K. Kumar          First creation of code.
 *      101025    E. Iorfida        First development of the code with conversion equations.
 *      101028    E. Iorfida        Modification of code for updating of the  Keplerian elements
 *                                  and Cartesian elements classes.
 *      101103    E. Iorfida        Additional conversion equations for extra  Keplerian elements.
 *      101104    E. Iorfida        Modification of the code (no pointers, but directly call of
 *                                  variables).
 *      101119    E. Iorfida        Removed computation for extra Keplerian elements.
 *      101130    E. Iorfida        Added different orbital cases with if-else operators.
 *      101202    J. Melman         Compile errors taken out.
 *      101203    E. Iorfida        Added gravitational parameter, and modified punctuation.
 *      101215    E. Iorfida        Added tolerance, modified punctuation, added comments, deleted
 *                                  raiseToIntegerExponent, used pow.
 *      101219    J. Melman         Suggested efficiency improvement of if-statements.
 *      110107    E. Iorfida        Written a better definition of the range in which angles are
 *                                  computed, and made some punctuation modifications.
 *      110109    J. Melman         Incorporated function computeModulo and
 *                                  determineAngleBetweenVectors. Reduced number of if-statements
 *                                  considerably and bundled all eccentricity and inclination
 *                                  checks in convertCartesianTopointerToCartesianElements_
 *      110128    K. Kumar          Changed references to pointers.
 *      110204    K. Kumar          Removed "vector" naming.
 *      110310    K. Kumar          Changed right ascension of ascending node to longitude of
 *                                  ascending node.
 *      110510    K. Kumar          Updated conversion functions to not use dynamic memory
 *                                  allocation.
 *      110805    K. Kumar          Added mean motion to semi-major axis conversion.
 */

// Include statements.
#include <cmath>
#include "Astrodynamics/States/orbitalElementConversions.h"

//! Tudat library namespace.
namespace tudat
{

// Using declarations.
using std::acos;
using std::atan2;
using std::cos;
using std::cosh;
using std::sin;
using std::sinh;
using std::sqrt;
using std::pow;
using std::fabs;
using linear_algebra::determineAngleBetweenVectors;
using tudat::mathematics::computeModulo;
using tudat::mathematics::MACHINE_PRECISION_DOUBLES;

//! Orbital element conversions namespace.
namespace orbital_element_conversions
{

//! Convert Keplerian to Cartesian orbital elements.
CartesianElements convertKeplerianToCartesianElements(
    KeplerianElements* pointerToKeplerianElements, CelestialBody* pointerToCelestialBody )
{
    // Declare local variables.
    // Declare the tolerance with which a defined double has to be equal to
    // a specific value.
    double tolerance_ = 10.0 * MACHINE_PRECISION_DOUBLES;

    // Declare CartesianElements object.
    CartesianElements cartesianElements_;

    // Pre-compute sines and cosines of involved angles for efficient
    // computation.
    double cosineOfInclination_ = cos( pointerToKeplerianElements->getInclination( ) );
    double sineOfInclination_ = sin( pointerToKeplerianElements->getInclination( ) );
    double cosineOfArgumentOfPeriapsis_ = cos(
                pointerToKeplerianElements->getArgumentOfPeriapsis( ) );
    double sineOfArgumentOfPeriapsis_ = sin(
                pointerToKeplerianElements->getArgumentOfPeriapsis( ) );
    double cosineOfLongitudeOfAscendingNode_ = cos(
                pointerToKeplerianElements->getLongitudeOfAscendingNode( ) );
    double sineOfLongitudeOfAscendingNode_ = sin(
                pointerToKeplerianElements->getLongitudeOfAscendingNode( ) );
    double cosineOfTrueAnomaly_ = cos( pointerToKeplerianElements->getTrueAnomaly( ) );
    double sineOfTrueAnomaly_ = sin( pointerToKeplerianElements->getTrueAnomaly( ) );

    // Compute semi-latus rectum in the case it is not a parabola (for which it
    // should already have been set).
    if ( pointerToKeplerianElements->getEccentricity( ) < 1.0 - tolerance_ ||
         pointerToKeplerianElements->getEccentricity( ) > 1.0 + tolerance_ )
    {
        pointerToKeplerianElements->setSemiLatusRectum(
                pointerToKeplerianElements->getSemiMajorAxis( )
                * ( 1.0 - pow( pointerToKeplerianElements->getEccentricity( ), 2.0 ) ) );
    }

    // Definition of position in the perifocal coordinate system.
    Vector2d positionPerifocal_;
    positionPerifocal_.x( ) =  pointerToKeplerianElements->getSemiLatusRectum( )
            * cosineOfTrueAnomaly_ / ( 1.0 + pointerToKeplerianElements->getEccentricity( )
                                       * cosineOfTrueAnomaly_ );
    positionPerifocal_.y( ) = pointerToKeplerianElements->getSemiLatusRectum( )
            * sineOfTrueAnomaly_ / ( 1.0 + pointerToKeplerianElements->getEccentricity( )
                                     * cosineOfTrueAnomaly_ );

    // Definition of velocity in the perifocal coordinate system.
    Vector2d velocityPerifocal_;
    velocityPerifocal_.x( ) = - sqrt( pointerToCelestialBody->getGravitationalParameter( ) /
                                      pointerToKeplerianElements->getSemiLatusRectum( ) ) *
            sineOfTrueAnomaly_;
    velocityPerifocal_.y( ) = sqrt( pointerToCelestialBody->getGravitationalParameter( ) /
                                    pointerToKeplerianElements->getSemiLatusRectum( ) ) *
            ( pointerToKeplerianElements->getEccentricity( ) + cosineOfTrueAnomaly_ );

    // Definition of the transformation matrix.
    MatrixXd transformationMatrix_;
    transformationMatrix_.setZero( 3, 2 );

    // Compute the transformation matrix.
    transformationMatrix_( 0, 0 ) = cosineOfLongitudeOfAscendingNode_ *
            cosineOfArgumentOfPeriapsis_ - sineOfLongitudeOfAscendingNode_ *
            sineOfArgumentOfPeriapsis_ * cosineOfInclination_;
    transformationMatrix_( 0, 1 ) = - cosineOfLongitudeOfAscendingNode_ *
            sineOfArgumentOfPeriapsis_ - sineOfLongitudeOfAscendingNode_ *
            cosineOfArgumentOfPeriapsis_ * cosineOfInclination_;
    transformationMatrix_( 1, 0 ) = sineOfLongitudeOfAscendingNode_ *
            cosineOfArgumentOfPeriapsis_ + cosineOfLongitudeOfAscendingNode_ *
            sineOfArgumentOfPeriapsis_ * cosineOfInclination_;
    transformationMatrix_( 1, 1 ) = - sineOfLongitudeOfAscendingNode_ *
            sineOfArgumentOfPeriapsis_ + cosineOfLongitudeOfAscendingNode_ *
            cosineOfArgumentOfPeriapsis_ * cosineOfInclination_;
    transformationMatrix_( 2, 0 ) = sineOfArgumentOfPeriapsis_ * sineOfInclination_;
    transformationMatrix_( 2, 1 ) = cosineOfArgumentOfPeriapsis_ * sineOfInclination_;

    // Compute value of position in Cartesian coordinates.
    Vector3d position_;
    position_ = ( transformationMatrix_ * positionPerifocal_ );
    cartesianElements_.setPosition( position_ );

    // Compute value of velocity in Cartesian coordinates.
    Vector3d velocity_;
    velocity_ = ( transformationMatrix_ * velocityPerifocal_ );
    cartesianElements_.setVelocity( velocity_ );

    // Return CartesianElements object.
    return cartesianElements_;
}

//! Convert Cartesian to Keplerian orbital elements.
KeplerianElements convertCartesianToKeplerianElements(
    CartesianElements* pointerToCartesianElements, CelestialBody* pointerToCelestialBody )
{
    // Declare local variables.
    // Declare the tolerance with which a computed double has to be equal to
    // a specific value.
    double tolerance_ = 10.0 * MACHINE_PRECISION_DOUBLES;

    // Declare KeplerianElements object.
    KeplerianElements keplerianElements_;

    // Norm of position in the inertial frame.
    double normOfPosition_;
    normOfPosition_ = pointerToCartesianElements->getPosition( ).norm( );

    // Norm of velocity in the inertial frame.
    double normOfVelocity_;
    normOfVelocity_ = pointerToCartesianElements->getVelocity( ).norm( );

    // Definition of orbit angular momentum.
    Vector3d orbitAngularMomentum_;
    orbitAngularMomentum_ = pointerToCartesianElements->getPosition( ).
            cross( pointerToCartesianElements->getVelocity( ) );
    double normOfOrbitAngularMomentum_ = orbitAngularMomentum_.norm( );

    // Definition of the (unit) vector to the ascending node.
    Vector3d unitVectorToAscendingNode_;
    unitVectorToAscendingNode_ = Vector3d::UnitZ( ).cross( orbitAngularMomentum_.normalized( ) );

    // Definition of eccentricity vector.
    Vector3d eccentricityVector_;
    eccentricityVector_ =
            ( ( pointerToCartesianElements->getVelocity( ).cross( orbitAngularMomentum_ ) )
              / pointerToCelestialBody->getGravitationalParameter( ) ) -
            pointerToCartesianElements->getPosition( ).normalized( );

    // Compute the total orbital energy.
    double totalOrbitalEnergy_;
    totalOrbitalEnergy_ = pow( normOfVelocity_, 2.0 ) / 2.0
            - pointerToCelestialBody->getGravitationalParameter( ) / normOfPosition_;

    // Compute the value of the eccentricity.
    keplerianElements_.setEccentricity( eccentricityVector_.norm( ) );

    // Define and compute boolean of whether orbit is circular or not.
    bool isOrbitCircular_ = keplerianElements_.getEccentricity( ) < tolerance_;

    // Compute the value of inclination.
    // Range between 0 degrees and 180 degrees.
    keplerianElements_.setInclination( acos( orbitAngularMomentum_.z( )
                                             / normOfOrbitAngularMomentum_ ) );

    // Define and compute boolean of whether orbit is equatorial or not.
    bool isOrbitEquatorial_ = unitVectorToAscendingNode_.norm( ) < tolerance_;

    // Compute the value of semi-major axis.
    // Non-parabolic orbits.
    if ( fabs( keplerianElements_.getEccentricity( ) - 1.0 ) > tolerance_ )
    {
        keplerianElements_.setSemiMajorAxis( pointerToCelestialBody->getGravitationalParameter( ) /
                                             ( -2.0 * totalOrbitalEnergy_ ) );
    }

    // Parabolic orbits.
    else
    {
        // Semi-major axis is infinite.
        keplerianElements_.setSemiMajorAxis( 1.0e100 );
    }

    // Compute value of semi-latus rectum.
    keplerianElements_.setSemiLatusRectum(
                pow( normOfOrbitAngularMomentum_, 2.0 )
                / pointerToCelestialBody->getGravitationalParameter( ) );

    // Compute the value of argument of periapsis.
    // Range between 0 degrees and 360 degrees.
    // Non-circular, inclined orbits.
    if ( !isOrbitCircular_ && !isOrbitEquatorial_ )
    {
        keplerianElements_.setArgumentOfPeriapsis(
                    determineAngleBetweenVectors( eccentricityVector_,
                                                  unitVectorToAscendingNode_ ) );

        // Quadrant check.
        if ( eccentricityVector_( 2 ) < 0.0 )
        {
            keplerianElements_.setArgumentOfPeriapsis(
                        2.0 * M_PI - keplerianElements_.getArgumentOfPeriapsis( ) );
        }
    }

    // Circular orbits.
    else if ( isOrbitCircular_ )
    {
        keplerianElements_.setArgumentOfPeriapsis( 0.0 );
    }

    // Equatorial orbits.
    // Argument of periapsis is defined as angle between eccentricity vector
    // and x-axis.
    else
    {
        keplerianElements_.setArgumentOfPeriapsis(
                    computeModulo( atan2( eccentricityVector_( 1 ),
                                          eccentricityVector_( 0 ) ), 2.0 * M_PI ) );
    }

    // Compute the value of longitude of ascending node.
    // Range between 0 degrees and 360 degrees.
    // Non-equatorial orbits.
    if ( !isOrbitEquatorial_ )
    {
        keplerianElements_.setLongitudeOfAscendingNode(
                    computeModulo( atan2( unitVectorToAscendingNode_.y( ),
                                          unitVectorToAscendingNode_.x( ) ), 2.0 * M_PI ) );
    }

    // Equatorial orbits.
    else
    {
        keplerianElements_.setLongitudeOfAscendingNode( 0.0 );
    }

    // Compute the value of true anomaly.
    // Range between 0 degrees and 360 degrees.
    // Non-circular orbits.
    if ( !isOrbitCircular_ )
    {
        keplerianElements_.setTrueAnomaly(
                    determineAngleBetweenVectors( pointerToCartesianElements->getPosition( ),
                                                  eccentricityVector_ ) );

        // Quadrant check. In the second half of the orbit, the angle
        // between position and velocity vector is larger than 90 degrees.
        if ( pointerToCartesianElements->getVelocity( ).
             dot( pointerToCartesianElements->getPosition( ) ) < 0.0 )
        {
            keplerianElements_.setTrueAnomaly(
                        2.0 * M_PI - keplerianElements_.getTrueAnomaly( ) );
        }
    }

    // Circular orbits.
    else
    {
        // Circular equatorial orbits.
        if ( isOrbitEquatorial_ )
        {
            keplerianElements_.setTrueAnomaly(
                        computeModulo(
                            atan2( pointerToCartesianElements->getCartesianElementY( ),
                                   pointerToCartesianElements->getCartesianElementX( ) ),
                            2.0 * M_PI ) );
        }

        // Circular inclined orbits.
        else
        {
            keplerianElements_.setTrueAnomaly(
                        determineAngleBetweenVectors( pointerToCartesianElements->getPosition( ),
                                                      unitVectorToAscendingNode_ ) );

            // Quadrant check. In the second half of the orbit, the body
            // will be below the xy-plane.
            if ( pointerToCartesianElements->getCartesianElementZ( ) < 0.0 )
            {
                keplerianElements_.setTrueAnomaly(
                            2.0 * M_PI - keplerianElements_.getTrueAnomaly( ) );
            }
        }
    }

    // Return KeplerianElements object.
    return keplerianElements_;
}

//! Convert true anomaly to eccentric anomaly.
double convertTrueAnomalyToEccentricAnomaly( const double& trueAnomaly,
                                             const double& eccentricity )
{
    // Declare and compute sine and cosine of eccentric anomaly.
    double sineOfEccentricAnomaly_ = sqrt( 1.0 - pow( eccentricity, 2.0 ) ) * sin( trueAnomaly )
            / ( 1.0 + eccentricity * cos( trueAnomaly ) );

    double cosineOfEccentricAnomaly_ = ( eccentricity + cos( trueAnomaly ) )
            / ( 1.0 + eccentricity * cos( trueAnomaly ) );

    // Return eccentric anomaly.
    return atan2( sineOfEccentricAnomaly_, cosineOfEccentricAnomaly_ );
}

//! Convert eccentric anomaly to true anomaly.
double convertEccentricAnomalyToTrueAnomaly( const double& eccentricAnomaly,
                                             const double& eccentricity )
{
    // Compute sine and cosine of true anomaly.
    double sineOfTrueAnomaly_ = sqrt( 1.0 - pow( eccentricity, 2.0 ) ) * sin( eccentricAnomaly )
            / ( 1.0 - eccentricity * cos( eccentricAnomaly ) );

    double cosineOfTrueAnomaly_ = ( cos( eccentricAnomaly ) - eccentricity )
            / ( 1.0 - eccentricity * cos( eccentricAnomaly ) );

    // Return true anomaly.
    return atan2( sineOfTrueAnomaly_, cosineOfTrueAnomaly_  );
}

//! Convert true anomaly to hyperbolic eccentric anomaly.
double convertTrueAnomalyToHyperbolicEccentricAnomaly( const double& trueAnomaly,
                                                       const double& eccentricity )
{
    // Compute hyperbolic sine and hyperbolic cosine of hyperbolic eccentric
    // anomaly.
    double hyperbolicSineOfHyperbolicEccentricAnomaly_ = sqrt( pow( eccentricity, 2.0 ) - 1.0 )
            * sin( trueAnomaly ) / ( 1.0 + cos( trueAnomaly ) );

    double hyperbolicCosineOfHyperbolicEccentricAnomaly_ = ( cos( trueAnomaly ) + eccentricity )
            / ( 1.0 + cos( trueAnomaly ) );

    // Return hyperbolic eccentric anomaly.
    return atanh( hyperbolicSineOfHyperbolicEccentricAnomaly_
                  / hyperbolicCosineOfHyperbolicEccentricAnomaly_ );
}

//! Convert hyperbolic eccentric anomaly to true anomaly.
double convertHyperbolicEccentricAnomalyToTrueAnomaly( const double& hyperbolicEccentricAnomaly,
                                                       const double& eccentricity )
{
    // Compute sine and cosine of true anomaly.
    double sineOfTrueAnomaly_
            = sqrt( pow( eccentricity, 2.0 ) - 1.0 ) * sinh( hyperbolicEccentricAnomaly )
            / ( eccentricity * cosh( hyperbolicEccentricAnomaly ) - 1.0 );

    double cosineOfTrueAnomaly_
            = ( eccentricity - cosh( hyperbolicEccentricAnomaly ) )
            / ( eccentricity * cosh( hyperbolicEccentricAnomaly ) - 1.0 );

    // Return true anomaly.
    return atan2( sineOfTrueAnomaly_, cosineOfTrueAnomaly_ );
}

//! Convert eccentric anomaly to mean anomaly.
double convertEccentricAnomalyToMeanAnomaly( const double& eccentricAnomaly,
                                             const double& eccentricity )
{
    return eccentricAnomaly - eccentricity * sin( eccentricAnomaly );
}

//! Convert hyperbolic eccentric anomaly to mean anomaly.
double convertHyperbolicEccentricAnomalyToMeanAnomaly(
    const double& hyperbolicEccentricAnomaly, const double& eccentricity )
{
    return eccentricity * sinh( hyperbolicEccentricAnomaly ) - hyperbolicEccentricAnomaly;
}

//! Convert elapsed time to mean anomaly for elliptical orbits.
double convertElapsedTimeToMeanAnomalyForEllipticalOrbits(
    const double& elapsedTime, CelestialBody* pointerToCentralBody, const double& semiMajorAxis )
{
    // Return mean anomaly.
    return sqrt( pointerToCentralBody->getGravitationalParameter( )
                 / pow( semiMajorAxis, 3.0 ) ) * elapsedTime;
}

//! Convert mean anomaly to elapsed time for elliptical orbits.
double convertMeanAnomalyToElapsedTimeForEllipticalOrbits(
    const double& meanAnomaly, CelestialBody* pointerToCentralBody, const double& semiMajorAxis )
{
    // Return time since last periapsis passage.
    return meanAnomaly * sqrt( pow( semiMajorAxis, 3.0 )
                               / pointerToCentralBody->getGravitationalParameter( ) );
}

//! Convert elapsed time to mean anomaly for hyperbolic orbits.
double convertElapsedTimeToMeanAnomalyForHyperbolicOrbits(
    const double& elapsedTime, CelestialBody* pointerToCentralBody, const double& semiMajorAxis )
{
    // Return mean anomaly.
    return sqrt( pointerToCentralBody->getGravitationalParameter( )
                 / pow( -semiMajorAxis, 3.0 ) ) * elapsedTime;
}

//! Convert mean anomaly to elapsed time for hyperbolic orbits.
double convertMeanAnomalyToElapsedTimeForHyperbolicOrbits(
    const double& meanAnomaly, CelestialBody* pointerToCentralBody, const double& semiMajorAxis )
{
    // Return time since last periapsis passage.
    return sqrt( pow( -semiMajorAxis, 3.0 ) / pointerToCentralBody->getGravitationalParameter( ) )
            * meanAnomaly;
}

//! Convert mean motion to semi-major axis.
double convertMeanMotionToSemiMajorAxis( const double& meanMotion,
                                         CelestialBody* pointerToCentralBody )
{
    // Return semi-major axis.
    return pow( pointerToCentralBody->getGravitationalParameter( )
                / pow( meanMotion, 2.0 ), 1.0 / 3.0 );

}

}

}

// End of file.
