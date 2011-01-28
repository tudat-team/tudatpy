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
 *                              Keplerian elements and Cartesian elements
 *                              classes.
 *      101103    E. Iorfida    Additional conversion equations for extra
 *                              Keplerian elements.
 *      101104    E. Iorfida    Modification of the code (no pointers, but
 *                              directly call of variables).
 *      101119    E. Iorfida    Removed computation for extra Keplerian
 *                              elements.
 *      101130    E. Iorfida    Added different orbital cases with if-else
 *                              operators.
 *      101202    J. Melman     Compile errors taken out.
 *      101203    E. Iorfida    Added gravitational parameter, and modified
 *                              punctuation.
 *      101215    E. Iorfida    Added tolerance, modified punctuation, added
 *                              comments, deleted raiseToIntegerExponent, used pow.
 *      101219    J. Melman     Suggested efficiency improvement of if-statements.
 *      110107    E. Iorfida    Written a better definition of the range in which
 *                              angles are computed, and made some punctuation
 *                              modifications.
 *      110109    J. Melman     Incorporated function determineAngleBetweenVectors
 *                              and computeModulo. Reduced the number of
 *                              if-statements considerably and bundled all
 *                              eccentricity and inclination checks in
 *                              convertCartesianTopointerToCartesianElements_->
 *      110128    K. Kumar      Changed references to pointers.

 */

// Include statements.
#include "orbitalElementConversions.h"

// Using declarations.
using linear_algebra::determineAngleBetweenVectors;
using mathematics::raiseToIntegerPower;
using mathematics::computeModulo;
using mathematics::computeAbsoluteValue;

//! Orbital element conversions namespace.
namespace orbital_element_conversions
{

//! Convert Keplerian to Cartesian orbital elements.
CartesianElements* convertKeplerianToCartesianElements(
        KeplerianElements* pointerToKeplerianElements,
        CelestialBody* pointerToCelestialBody )
{
    // Declare local variables.
    // Declare the output of this routine.
    CartesianElements* pointerToCartesianElements_ = new CartesianElements;

    // Declare the tolerance with which a defined double has to be equal to
    // a specific value.
    double tolerance_ = 10.0 * mathematics::MACHINE_PRECISION_DOUBLES;

    // Pre-compute sines and cosines of involved angles for efficient
    // computation.
    double cosineOfInclination_ = cos( pointerToKeplerianElements
                                       ->getInclination ( ) );
    double sineOfInclination_ = sin( pointerToKeplerianElements
                                     ->getInclination ( ) );
    double cosineOfArgumentOfPeriapsis_ = cos(
            pointerToKeplerianElements->getArgumentOfPeriapsis ( ) );
    double sineOfArgumentOfPeriapsis_ = sin(
            pointerToKeplerianElements->getArgumentOfPeriapsis ( ) );
    double cosineOfRightAscensionOfAscendingNode_ = cos(
            pointerToKeplerianElements
            ->getRightAscensionOfAscendingNode ( ) );
    double sineOfRightAscensionOfAscendingNode_ = sin(
            pointerToKeplerianElements->getRightAscensionOfAscendingNode ( ) );
    double cosineOfTrueAnomaly_ = cos( pointerToKeplerianElements
                                       ->getTrueAnomaly ( ) );
    double sineOfTrueAnomaly_ = sin( pointerToKeplerianElements
                                     ->getTrueAnomaly ( ) );

    // Compute semi-latus rectum in the case it is not a parabola (for which it
    // should already have been set).
    if ( pointerToKeplerianElements->getEccentricity( ) < 1.0 - tolerance_ ||
         pointerToKeplerianElements->getEccentricity( ) > 1.0 + tolerance_ )
    {
        pointerToKeplerianElements->setSemiLatusRectum(
            pointerToKeplerianElements->getSemiMajorAxis( ) * ( 1.0 -
                raiseToIntegerPower( pointerToKeplerianElements
                                     ->getEccentricity( ), 2 ) ) );
    }

    // Definition of position vector in the perifocal coordinate system.
    Vector2d positionPerifocal_;
    positionPerifocal_.x( ) =
        pointerToKeplerianElements->getSemiLatusRectum( )
        * cosineOfTrueAnomaly_ /
        ( 1.0 + pointerToKeplerianElements->getEccentricity( )
          * cosineOfTrueAnomaly_ );
    positionPerifocal_.y( ) =
        pointerToKeplerianElements->getSemiLatusRectum( )
        * sineOfTrueAnomaly_ /
        ( 1.0 + pointerToKeplerianElements->getEccentricity( )
          * cosineOfTrueAnomaly_ );

    // Definition of velocity vector in the perifocal coordinate system.
    Vector2d velocityPerifocal_;
    velocityPerifocal_.x( ) = - sqrt( pointerToCelestialBody
                                      ->getGravitationalParameter( ) /
                              pointerToKeplerianElements
                              ->getSemiLatusRectum( ) ) *
                              sineOfTrueAnomaly_;
    velocityPerifocal_.y( ) = sqrt( pointerToCelestialBody
                                    ->getGravitationalParameter( ) /
                              pointerToKeplerianElements
                              ->getSemiLatusRectum( ) ) *
                              ( pointerToKeplerianElements
                                ->getEccentricity( ) +
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
    pointerToCartesianElements_->setPositionVector( positionVector_ );

    // Compute value of velocity vector in Cartesian coordinates.
    Vector3d velocityVector_;
    velocityVector_ = ( transformationMatrix_ *
                        velocityPerifocal_ );
    pointerToCartesianElements_->setVelocityVector( velocityVector_ );

    return pointerToCartesianElements_;
}

//! Convert Cartesian to Keplerian orbital elements.
KeplerianElements* convertCartesianToKeplerianElements(
        CartesianElements* pointerToCartesianElements,
        CelestialBody* pointerToCelestialBody )
{
    // Local variables.

    // Declare the output of this routine.
    KeplerianElements* pointerToKeplerianElements_ = new KeplerianElements;

    // Declare the tolerance with which a computed double has to be equal to
    // a specific value.
    double tolerance_ = 10.0 * mathematics::MACHINE_PRECISION_DOUBLES;

    // Norm of position vector in the inertial frame.
    double normOfPositionVector_;
    normOfPositionVector_ = pointerToCartesianElements
                            ->getPositionVector( ).norm( );

    // Norm of velocity vector in the inertial frame.
    double normOfVelocityVector_;
    normOfVelocityVector_ = pointerToCartesianElements
                            ->getVelocityVector( ).norm( );

    // Definition of orbit angular momentum.
    Vector3d orbitAngularMomentum_;
    orbitAngularMomentum_ = pointerToCartesianElements->getPositionVector( ).
                            cross( pointerToCartesianElements
                                   ->getVelocityVector( ) );
    double normOfOrbitAngularMomentum_ = orbitAngularMomentum_.norm( );

    // Definition of the (unit) vector to the ascending node.
    Vector3d vectorToAscendingNode_;
    vectorToAscendingNode_ = Vector3d::UnitZ( ).cross( orbitAngularMomentum_ );
    Vector3d unitVectorToAscendingNode_;
    unitVectorToAscendingNode_ = vectorToAscendingNode_.normalized( );

    // Definition of eccentricity vector.
    Vector3d eccentricityVector_;
    eccentricityVector_ =
        ( ( pointerToCartesianElements->getVelocityVector( ).
            cross( orbitAngularMomentum_ ) )
          / pointerToCelestialBody->getGravitationalParameter( ) ) -
        pointerToCartesianElements->getPositionVector( ).normalized( );

    // Compute the total orbital energy.
    double totalOrbitalEnergy_;
    totalOrbitalEnergy_ =
        raiseToIntegerPower( normOfVelocityVector_, 2 ) / 2.0 -
        pointerToCelestialBody->getGravitationalParameter( )
        / normOfPositionVector_ ;

    // Compute the value of the eccentricity.
    pointerToKeplerianElements_->setEccentricity( eccentricityVector_.norm( ) );

    // Define and compute boolean of whether orbit is circular or not.
    bool isOrbitCircular_ = pointerToKeplerianElements_->getEccentricity( )
                            < tolerance_;

    // Compute the value of inclination.
    // Range between 0 degrees and 180 degrees.
    pointerToKeplerianElements_->setInclination(
            acos( orbitAngularMomentum_.z( ) / normOfOrbitAngularMomentum_ ) );

    // Define and compute boolean of whether orbit is equatorial or not.
    bool isOrbitEquatorial_ = vectorToAscendingNode_.norm( ) < tolerance_;

    // Compute the value of semi-major axis.
    // Non-parabolic orbits.
    if ( computeAbsoluteValue( pointerToKeplerianElements_
                               ->getEccentricity( ) - 1.0 ) > tolerance_ )
    {
        pointerToKeplerianElements_->setSemiMajorAxis(
                pointerToCelestialBody->getGravitationalParameter( ) /
                ( -2.0 * totalOrbitalEnergy_ ) );
    }

    // Parabolic orbits.
    else
    {
        // Semi-major axis is infinite.
        pointerToKeplerianElements_->setSemiMajorAxis( 1.0e100 );
    }

    // Compute value of semi-latus rectum.
    pointerToKeplerianElements_->setSemiLatusRectum( raiseToIntegerPower(
            normOfOrbitAngularMomentum_, 2 ) / pointerToCelestialBody
                                           ->getGravitationalParameter( ) );

    // Compute the value of argument of periapsis.
    // Range between 0 degrees and 360 degrees.
    // Non-circular, inclined orbits.
    if ( !isOrbitCircular_ && !isOrbitEquatorial_ )
    {
        pointerToKeplerianElements_->setArgumentOfPeriapsis(
                determineAngleBetweenVectors( eccentricityVector_,
                                              vectorToAscendingNode_ ) );

        // Quadrant check.
        if ( eccentricityVector_( 2 ) < 0.0 )
        {
            pointerToKeplerianElements_->setArgumentOfPeriapsis( 2.0 * M_PI -
                pointerToKeplerianElements_->getArgumentOfPeriapsis( ) );
        }
    }

    // Circular orbits.
    else if ( isOrbitCircular_ )
    {
        pointerToKeplerianElements_->setArgumentOfPeriapsis( 0.0 );
    }

    // Equatorial orbits.
    // Argument of periapsis is defined as angle between eccentricity vector
    // and x-axis.
    else
    {
        pointerToKeplerianElements_->setArgumentOfPeriapsis( computeModulo(
                atan2( eccentricityVector_( 1 ), eccentricityVector_( 0 ) ),
                2.0 * M_PI ) );
    }

    // Compute the value of right ascension of ascending node.
    // Range between 0 degrees and 360 degrees.
    // Non-equatorial orbits.
    if ( !isOrbitEquatorial_ )
    {
        pointerToKeplerianElements_->setRightAscensionOfAscendingNode(
                computeModulo( atan2( unitVectorToAscendingNode_.y( ),
                       unitVectorToAscendingNode_.x( ) ), 2.0 * M_PI ) );
    }

    // Equatorial orbits.
    else
    {
        pointerToKeplerianElements_->setRightAscensionOfAscendingNode( 0.0 );
    }

    // Compute the value of true anomaly.
    // Range between 0 degrees and 360 degrees.
    // Non-circular orbits.
    if ( !isOrbitCircular_ )
    {
        pointerToKeplerianElements_->setTrueAnomaly(
                determineAngleBetweenVectors( pointerToCartesianElements
                                              ->getPositionVector( ),
                eccentricityVector_ ) );

        // Quadrant check. In the second half of the orbit, the angle
        // between position and velocity vector is larger than 90 degrees.
        if ( pointerToCartesianElements->getVelocityVector( ).
             dot( pointerToCartesianElements->getPositionVector( ) ) < 0.0 )
        {
            pointerToKeplerianElements_->setTrueAnomaly( 2.0 * M_PI -
                pointerToKeplerianElements_->getTrueAnomaly( ) );
        }
    }

    // Circular orbits.
    else
    {
        // Circular equatorial orbits.
        if ( isOrbitEquatorial_ )
        {
            pointerToKeplerianElements_->setTrueAnomaly( computeModulo(
                    atan2( pointerToCartesianElements->getCartesianElementY( ),
                    pointerToCartesianElements->getCartesianElementX( ) ),
                    2.0 * M_PI ) );
        }

        // Circular inclined orbits.
        else
        {
            pointerToKeplerianElements_->setTrueAnomaly(
                    determineAngleBetweenVectors( pointerToCartesianElements
                                                  ->getPositionVector( ),
                    vectorToAscendingNode_ ) );

            // Quadrant check. In the second half of the orbit, the body
            // will be below the xy-plane.
            if ( pointerToCartesianElements->getCartesianElementZ( ) < 0.0 )
            {
                pointerToKeplerianElements_->setTrueAnomaly( 2.0 * M_PI -
                     pointerToKeplerianElements_->getTrueAnomaly( ) );
            }
        }
    }

    return pointerToKeplerianElements_;
}

}

// End of file.
