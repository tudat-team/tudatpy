/*! \file orbitalElementConversions.cpp
 *    This source file contains a namespace with orbital element conversion
 *    functions.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 13
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
 *    Date created      : 20 October, 2010
 *    Last modified     : 09 January, 2011
 *
 *    References
 *      Mission geometry; orbit and constellation design and management (James R. Wertz)
 *      Fondamenti di meccanica del volo spaziale (Giovanni Mengali, Alessandro A. Quarta)
 *
 *    Notes
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
 *      YYMMDD    author        comment
 *      101020    K. Kumar      First creation of code.
 *      101025    E. Iorfida    First development of the code with conversion
 *                              equations.
 *      101028    E. Iorfida    Modification of code for updating of the
 *                              keplerian elements and cartesian elements classes.
 *      101103    E. Iorfida    Additional conversion equations for extra
 *                              keplerian elements.
 *      101104    E. Iorfida    Modification of the code (no pointers, but
 *                              directly call of variables).
 *      101119    E. Iorfida    Removed computation for extra keplerian
 *                              elements.
 *      101130    E. Iorfida    Added different orbital cases with if-else
 *                              operators.
 *      101202    J. Melman     Compile errors taken out.
 *      101203    E. Iorfida    Added gravitational parameter, and modified
 *                              punctuation.
 *      101215    E. Iorfida    Added tolerance, modified punctuation, added comments,
 *                              deleted raiseToIntegerExponent, used pow.
 *      101219    J. Melman     Suggested efficiency improvement of if-statements.
 *      110107    E. Iorfida    Written a better definition of the range in which
 *                              angles are computed, and made some punctuation
 *                              modifications.
 *      110109    J. Melman     Incorporated function determineAngleBetweenVectors
 *                              and computeModulo. Reduced the number of if-statements
 *                              considerably and bundled all eccentricity and inclination
 *                              checks in convertCartesianToKeplerianElements.
 */

// Include statements.
#include "orbitalElementConversions.h"

// Using directives.
using linear_algebra::determineAngleBetweenVectors;
using mathematics::raiseToIntegerPower;
using mathematics::computeModulo;

//! Orbital element conversions namespace.
namespace orbital_element_conversions
{

//! Function to convert Keplerian to Cartesian orbital elements.
CartesianElements convertKeplerianToCartesianElements(
        KeplerianElements& keplerianElements,
        double gravitationalParameter )
{
    // Declare local variables.
    // Declare the output of this routine.
    CartesianElements cartesianElements_;

    // Declare the tolerance with which a defined double has to be equal to
    // a specific value.
    double tolerance_ = 10.0 * mathematics::MACHINE_PRECISION_DOUBLES;

    // Precompute sines and cosines of involved angles for efficient computation.
    double cosineOfInclination_ = cos( keplerianElements.getInclination ( ) );
    double sineOfInclination_ = sin( keplerianElements.getInclination ( ) );
    double cosineOfArgumentOfPeriapsis_ = cos(
            keplerianElements.getArgumentOfPeriapsis ( ) );
    double sineOfArgumentOfPeriapsis_ = sin(
            keplerianElements.getArgumentOfPeriapsis ( ) );
    double cosineOfRightAscensionOfAscendingNode_ = cos(
            keplerianElements.getRightAscensionOfAscendingNode ( ) );
    double sineOfRightAscensionOfAscendingNode_ = sin(
            keplerianElements.getRightAscensionOfAscendingNode ( ) );
    double cosineOfTrueAnomaly_ = cos( keplerianElements.getTrueAnomaly ( ) );
    double sineOfTrueAnomaly_ = sin( keplerianElements.getTrueAnomaly ( ) );

    // Compute semi-latus rectum in the case it is not a parabola (for which it
    // should already have been set).
    if ( keplerianElements.getEccentricity( ) < 1.0 - tolerance_ ||
         keplerianElements.getEccentricity( ) > 1.0 + tolerance_ )
    {
        keplerianElements.setSemiLatusRectum(
            keplerianElements.getSemiMajorAxis( ) * ( 1.0 -
                raiseToIntegerPower( keplerianElements.getEccentricity( ),
                                     2 ) ) );
    }

    // Definition of position vector in the perifocal coordinate system.
    Vector2d positionPerifocal_;
    positionPerifocal_.x( ) =
        keplerianElements.getSemiLatusRectum( ) * cosineOfTrueAnomaly_ /
        ( 1.0 + keplerianElements.getEccentricity( ) * cosineOfTrueAnomaly_ );
    positionPerifocal_.y( ) =
        keplerianElements.getSemiLatusRectum( ) * sineOfTrueAnomaly_ /
        ( 1.0 + keplerianElements.getEccentricity( ) * cosineOfTrueAnomaly_ );

    // Definition of velocity vector in the perifocal coordinate system.
    Vector2d velocityPerifocal_;
    velocityPerifocal_.x( ) = - sqrt( gravitationalParameter /
                              keplerianElements.getSemiLatusRectum( ) ) *
                              sineOfTrueAnomaly_;
    velocityPerifocal_.y( ) = sqrt( gravitationalParameter /
                              keplerianElements.getSemiLatusRectum( ) ) *
                              ( keplerianElements.getEccentricity( ) +
                              cosineOfTrueAnomaly_ );

    // Definition of the transformation matrix.
    MatrixXd transformationMatrix_;
    transformationMatrix_.setZero( 3, 2 );

    // Compute the transformation matrix.
    transformationMatrix_( 0, 0 ) = cosineOfRightAscensionOfAscendingNode_ *
                                    cosineOfArgumentOfPeriapsis_ -
                                    sineOfRightAscensionOfAscendingNode_ *
                                    sineOfArgumentOfPeriapsis_ *
                                    cosineOfInclination_;
    transformationMatrix_( 0, 1 ) = - cosineOfRightAscensionOfAscendingNode_ *
                                    sineOfArgumentOfPeriapsis_ -
                                    sineOfRightAscensionOfAscendingNode_ *
                                    cosineOfArgumentOfPeriapsis_ *
                                    cosineOfInclination_;
    transformationMatrix_( 1, 0 ) = sineOfRightAscensionOfAscendingNode_ *
                                    cosineOfArgumentOfPeriapsis_ +
                                    cosineOfRightAscensionOfAscendingNode_ *
                                    sineOfArgumentOfPeriapsis_ *
                                    cosineOfInclination_;
    transformationMatrix_( 1, 1 ) = - sineOfRightAscensionOfAscendingNode_ *
                                    sineOfArgumentOfPeriapsis_ +
                                    cosineOfRightAscensionOfAscendingNode_ *
                                    cosineOfArgumentOfPeriapsis_ *
                                    cosineOfInclination_;
    transformationMatrix_( 2, 0 ) = sineOfArgumentOfPeriapsis_ *
                                    sineOfInclination_;
    transformationMatrix_( 2, 1 ) = cosineOfArgumentOfPeriapsis_ *
                                    sineOfInclination_;

    // Compute value of position vector in Cartesian coordinates.
    Vector3d positionVector_;
    positionVector_ = ( transformationMatrix_ *
                        positionPerifocal_ );
    cartesianElements_.setPositionVector( positionVector_ );

    // Compute value of velocity vector in Cartesian coordinates.
    Vector3d velocityVector_;
    velocityVector_ = ( transformationMatrix_ *
                        velocityPerifocal_ );
    cartesianElements_.setVelocityVector( velocityVector_ );

    return cartesianElements_;
}

//! Function to convert Cartesian to Keplerian orbital elements.
KeplerianElements convertCartesianToKeplerianElements(
        CartesianElements& cartesianElements,
        double gravitationalParameter )
{
    // Local variables.

    // Declare the output of this routine.
    KeplerianElements keplerianElements_;

    // Declare the tolerance with which a computed double has to be equal to
    // a specific value.
    double tolerance_ = 10.0 * mathematics::MACHINE_PRECISION_DOUBLES;

    // Norm of position vector in the inertial frame.
    double normOfPositionVector_;
    normOfPositionVector_ = cartesianElements.getPositionVector( ).norm( );

    // Norm of velocity vector in the inertial frame.
    double normOfVelocityVector_;
    normOfVelocityVector_ = cartesianElements.getVelocityVector( ).norm( );

    // Definition of orbit angular momentum.
    Vector3d orbitAngularMomentum_;
    orbitAngularMomentum_ = cartesianElements.getPositionVector( ).
                            cross( cartesianElements.getVelocityVector( ) );
    double normOfOrbitAngularMomentum_ = orbitAngularMomentum_.norm( );

    // Definition of the (unit) vector to the ascending node.
    Vector3d vectorToAscendingNode_;
    vectorToAscendingNode_ = Vector3d::UnitZ( ).cross( orbitAngularMomentum_ );
    Vector3d unitVectorToAscendingNode_;
    unitVectorToAscendingNode_ = vectorToAscendingNode_.normalized( );

    // Definition of eccentricity vector.
    Vector3d eccentricityVector_;
    eccentricityVector_ =
        ( ( cartesianElements.getVelocityVector( ).
            cross( orbitAngularMomentum_ ) ) / gravitationalParameter ) -
        cartesianElements.getPositionVector( ).normalized( );

    // Compute the total orbital energy.
    double totalOrbitalEnergy_;
    totalOrbitalEnergy_ =
        raiseToIntegerPower( normOfVelocityVector_, 2 ) / 2.0 -
        gravitationalParameter / normOfPositionVector_ ;

    // Compute the value of the eccentricity.
    keplerianElements_.setEccentricity( eccentricityVector_.norm( ) );

    // Define and compute boolean of whether orbit is circular or not.
    bool isOrbitCircular_ = keplerianElements_.getEccentricity( ) < tolerance_;

    // Compute the value of inclination.
    // Range between 0 degrees and 180 degrees.
    keplerianElements_.setInclination( acos( orbitAngularMomentum_.z( ) /
                                             normOfOrbitAngularMomentum_ ) );

    // Define and compute boolean of whether orbit is equatorial or not.
    bool isOrbitEquatorial_ = vectorToAscendingNode_.norm( ) < tolerance_;

    // Compute the value of semi-major axis.
    // Non-parabolic orbits.
    if ( fabs( keplerianElements_.getEccentricity( ) - 1.0 ) > tolerance_ )
    {
        keplerianElements_.setSemiMajorAxis( gravitationalParameter /
                ( -2.0 * totalOrbitalEnergy_ ) );
    }

    // Parabolic orbits.
    else
    {
        // Semi-major axis is infinite.
        keplerianElements_.setSemiMajorAxis( 1.0e100 );
    }

    // Compute value of semi-latus rectum.
    keplerianElements_.setSemiLatusRectum( raiseToIntegerPower(
            normOfOrbitAngularMomentum_, 2 ) / gravitationalParameter );

    // Compute the value of argument of periapsis.
    // Range between 0 degrees and 360 degrees.
    // Non-circular, inclined orbits.
    if ( !isOrbitCircular_ && !isOrbitEquatorial_ )
    {
        keplerianElements_.setArgumentOfPeriapsis( determineAngleBetweenVectors(
                eccentricityVector_, vectorToAscendingNode_ ) );

        // Quadrant check.
        if ( eccentricityVector_( 2 ) < 0.0 )
        {
            keplerianElements_.setArgumentOfPeriapsis( 2.0 * M_PI -
                keplerianElements_.getArgumentOfPeriapsis( ) );
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
        keplerianElements_.setArgumentOfPeriapsis( computeModulo( atan2(
                eccentricityVector_( 1 ), eccentricityVector_( 0 ) ),
                2.0 * M_PI ) );
    }

    // Compute the value of right ascension of ascending node.
    // Range between 0 degrees and 360 degrees.
    // Non-equatorial orbits.
    if ( !isOrbitEquatorial_ )
    {
        keplerianElements_.setRightAscensionOfAscendingNode( computeModulo(
                atan2( unitVectorToAscendingNode_.y( ),
                       unitVectorToAscendingNode_.x( ) ), 2.0 * M_PI ) );
    }

    // Equatorial orbits.
    else
    {
        keplerianElements_.setRightAscensionOfAscendingNode( 0.0 );
    }

    // Compute the value of true anomaly.
    // Range between 0 degrees and 360 degrees.
    // Non-circular orbits.
    if ( !isOrbitCircular_ )
    {
        keplerianElements_.setTrueAnomaly( determineAngleBetweenVectors (
                cartesianElements.getPositionVector( ),
                eccentricityVector_ ) );

        // Quadrant check. In the second half of the orbit, the angle
        // between position and velocity vector is larger than 90 degrees.
        if ( cartesianElements.getVelocityVector( ).
             dot( cartesianElements.getPositionVector( ) ) < 0.0 )
        {
            keplerianElements_.setTrueAnomaly( 2.0 * M_PI -
                keplerianElements_.getTrueAnomaly( ) );
        }
    }

    // Circular orbits.
    else
    {
        // Circular equatorial orbits.
        if ( isOrbitEquatorial_ )
        {
            keplerianElements_.setTrueAnomaly( computeModulo( atan2(
                    cartesianElements.getCartesianElementY( ),
                    cartesianElements.getCartesianElementX( ) ), 2.0 * M_PI ) );
        }

        // Circular inclined orbits.
        else
        {
            keplerianElements_.setTrueAnomaly( determineAngleBetweenVectors (
                    cartesianElements.getPositionVector( ),
                    vectorToAscendingNode_ ) );

            // Quadrant check. In the second half of the orbit, the body
            // will be below the xy-plane.
            if ( cartesianElements.getCartesianElementZ( ) < 0.0 )
            {
                keplerianElements_.setTrueAnomaly( 2.0 * M_PI -
                     keplerianElements_.getTrueAnomaly( ) );
            }
        }
    }

    return keplerianElements_;
}

}

// End of file.
