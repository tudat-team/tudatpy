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
 *      120131    K. Kumar          Adapted for Tudat Core, interfaces changed to use VectorXd,
 *                                  only Keplerian <-> Cartesian conversions included.
 *      120206    K. Kumar          Added wrapper functions for orbital element conversions when
 *                                  eccentricity is not known a priori (if-statement to choose
 *                                  between elliptical and hyperbolic orbits).
 *      120422    K. Kumar          Rewrote Cartesian -> Keplerian conversion; now handles circular
 *                                  and/or equatorial solutions correctly.
 *      121205    D. Dirkx          Migrated namespace to directory-based protocol.
 *      130122    K. Kumar          Debugged true anomaly computation for limit cases in
 *                                  convertCartesianToKeplerianElements().
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *      Wertz, J. R. Mission geometry; orbit and constellation design and management.
 *      Mengali, G., Quarta, A.A. Fondamenti di meccanica del volo spaziale.
 *      Advanced Concepts Team, ESA. Keplerian Toolbox, http://sourceforge.net/projects/keptoolbox,
 *          last accessed: 21st April, 2012.
 *
 *    Notes
 *
 */

#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

#include <boost/exception/all.hpp>
#include <boost/math/special_functions/atanh.hpp>

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Convert Keplerian to Cartesian orbital elements.
basic_mathematics::Vector6d convertKeplerianToCartesianElements(
        const basic_mathematics::Vector6d& keplerianElements,
        const double centralBodyGravitationalParameter )
{
    using std::cos;
    using std::fabs;
    using std::pow;
    using std::sin;
    using std::sqrt;
    using Eigen::Vector2d;
    using Eigen::Vector3d;

    // Set tolerance.
    const double tolerance_ = std::numeric_limits< double >::epsilon( );

    // Set local keplerian elements.
    double semiMajorAxis_ = keplerianElements( semiMajorAxisIndex );
    double eccentricity_ = keplerianElements( eccentricityIndex );
    double inclination_ = keplerianElements( inclinationIndex );
    double argumentOfPeriapsis_ = keplerianElements( argumentOfPeriapsisIndex );
    double longitudeOfAscendingNode_ = keplerianElements( longitudeOfAscendingNodeIndex );
    double trueAnomaly_ = keplerianElements( trueAnomalyIndex );

    // Pre-compute sines and cosines of involved angles for efficient computation.
    double cosineOfInclination_ = cos( inclination_ );
    double sineOfInclination_ = sin( inclination_ );
    double cosineOfArgumentOfPeriapsis_ = cos( argumentOfPeriapsis_ );
    double sineOfArgumentOfPeriapsis_ = sin( argumentOfPeriapsis_ );
    double cosineOfLongitudeOfAscendingNode_ = cos( longitudeOfAscendingNode_ );
    double sineOfLongitudeOfAscendingNode_ = sin( longitudeOfAscendingNode_ );
    double cosineOfTrueAnomaly_ = cos( trueAnomaly_ );
    double sineOfTrueAnomaly_ = sin( trueAnomaly_ );

    // Declare semi-latus rectum.
    double semiLatusRectum_ = -0.0;

    // Compute semi-latus rectum in the case it is not a parabola.
    if ( fabs( eccentricity_ - 1.0 ) > tolerance_  )
    {  semiLatusRectum_ = semiMajorAxis_ * ( 1.0 - pow( eccentricity_, 2 ) ); }

    // Else set the semi-latus rectum given for a parabola as the first element in the vector
    // of Keplerian elements..
    else  { semiLatusRectum_ = semiMajorAxis_; }

    // Definition of position in the perifocal coordinate system.
    Vector2d positionPerifocal_ = Eigen::Vector2d::Zero( 2 );
    positionPerifocal_.x( ) = semiLatusRectum_ * cosineOfTrueAnomaly_
            / ( 1.0 + eccentricity_ * cosineOfTrueAnomaly_ );
    positionPerifocal_.y( ) = semiLatusRectum_ * sineOfTrueAnomaly_
            / ( 1.0 + eccentricity_ * cosineOfTrueAnomaly_ );

    // Definition of velocity in the perifocal coordinate system.
    Vector2d velocityPerifocal_(
                -sqrt( centralBodyGravitationalParameter / semiLatusRectum_ ) * sineOfTrueAnomaly_,
                sqrt( centralBodyGravitationalParameter / semiLatusRectum_ )
                * ( eccentricity_ + cosineOfTrueAnomaly_ ) );

    // Definition of the transformation matrix.
    Eigen::MatrixXd transformationMatrix_ = Eigen::MatrixXd::Zero( 3, 2 );

    // Compute the transformation matrix.
    transformationMatrix_( 0, 0 ) = cosineOfLongitudeOfAscendingNode_
            * cosineOfArgumentOfPeriapsis_ -sineOfLongitudeOfAscendingNode_
            * sineOfArgumentOfPeriapsis_ * cosineOfInclination_;
    transformationMatrix_( 0, 1 ) = -cosineOfLongitudeOfAscendingNode_
            * sineOfArgumentOfPeriapsis_ -sineOfLongitudeOfAscendingNode_
            * cosineOfArgumentOfPeriapsis_ * cosineOfInclination_;
    transformationMatrix_( 1, 0 ) = sineOfLongitudeOfAscendingNode_
            * cosineOfArgumentOfPeriapsis_ + cosineOfLongitudeOfAscendingNode_
            * sineOfArgumentOfPeriapsis_ * cosineOfInclination_;
    transformationMatrix_( 1, 1 ) = -sineOfLongitudeOfAscendingNode_
            * sineOfArgumentOfPeriapsis_ + cosineOfLongitudeOfAscendingNode_
            * cosineOfArgumentOfPeriapsis_ * cosineOfInclination_;
    transformationMatrix_( 2, 0 ) = sineOfArgumentOfPeriapsis_ * sineOfInclination_;
    transformationMatrix_( 2, 1 ) = cosineOfArgumentOfPeriapsis_ * sineOfInclination_;

    // Declare converted Cartesian elements.
    basic_mathematics::Vector6d convertedCartesianElements_ = basic_mathematics::Vector6d::Zero( );

    // Compute value of position in Cartesian coordinates.
    Vector3d position_ = transformationMatrix_ * positionPerifocal_;
    convertedCartesianElements_.segment( 0, 3 ) = position_;

    // Compute value of velocity in Cartesian coordinates.
    Vector3d velocity_ = transformationMatrix_ * velocityPerifocal_;
    convertedCartesianElements_.segment( 3, 3 ) = velocity_;

    // Return Cartesian elements.
    return convertedCartesianElements_;
}

//! Convert Cartesian to Keplerian orbital elements.
basic_mathematics::Vector6d convertCartesianToKeplerianElements(
        const basic_mathematics::Vector6d& cartesianElements,
        const double centralBodyGravitationalParameter )
{
    // Set tolerance.
    const double tolerance = 1.0e-15;

    // Declare converted Keplerian elements.
    basic_mathematics::Vector6d computedKeplerianElements_ = basic_mathematics::Vector6d::Zero( );

    // Set position and velocity vectors.
    const Eigen::Vector3d position_( cartesianElements.segment( 0, 3 ) );
    const Eigen::Vector3d velocity_( cartesianElements.segment( 3, 3 ) );

    // Compute orbital angular momentum vector.
    const Eigen::Vector3d angularMomentum_( position_.cross( velocity_ ) );

    // Compute semi-latus rectum.
    const double semiLatusRectum_ = angularMomentum_.squaredNorm( )
            / centralBodyGravitationalParameter;

    // Compute unit vector to ascending node.
    Eigen::Vector3d unitAscendingNodeVector_(
                ( Eigen::Vector3d::UnitZ( ).cross(
                      angularMomentum_.normalized( ) ) ).normalized( ) );

    // Compute eccentricity vector.
    Eigen::Vector3d eccentricityVector_(
                velocity_.cross( angularMomentum_ ) / centralBodyGravitationalParameter
                - position_.normalized( ) );

    // Store eccentricity.
    computedKeplerianElements_( eccentricityIndex ) = eccentricityVector_.norm( );

    // Compute and store semi-major axis.
    // Check if orbit is parabolic. If it is, store the semi-latus rectum instead of the
    // semi-major axis.
    if ( std::fabs( computedKeplerianElements_( eccentricityIndex ) - 1.0 ) < tolerance )
    {
        computedKeplerianElements_( semiLatusRectumIndex ) = semiLatusRectum_;
    }

    // Else the orbit is either elliptical or hyperbolic, so store the semi-major axis.
    else
    {
        computedKeplerianElements_( semiMajorAxisIndex ) = semiLatusRectum_
                / ( 1.0 - computedKeplerianElements_( eccentricityIndex )
                    * computedKeplerianElements_( eccentricityIndex ) );
    }

    // Compute and store inclination.
    computedKeplerianElements_( inclinationIndex ) = std::acos( angularMomentum_.z( )
                                                                / angularMomentum_.norm( ) );

    // Compute and store longitude of ascending node.
    // Define the quadrant condition for the argument of perigee.
    double argumentOfPeriapsisQuandrantCondition = eccentricityVector_.z( );

    // Check if the orbit is equatorial. If it is, set the vector to the line of nodes to the
    // x-axis.
    if ( std::fabs( computedKeplerianElements_( inclinationIndex ) ) < tolerance )
    {
        unitAscendingNodeVector_ = Eigen::Vector3d::UnitX( );

        // If the orbit is equatorial, eccentricityVector_.z( ) is zero, therefore the quadrant
        // condition is taken to be the y-component, eccentricityVector_.y( ).
        argumentOfPeriapsisQuandrantCondition = eccentricityVector_.y( );
    }

    // Compute and store the resulting longitude of ascending node.
    computedKeplerianElements_( longitudeOfAscendingNodeIndex )
            = acos( unitAscendingNodeVector_.x( ) );

    // Check if the quandrant is correct.
    using mathematical_constants::PI;
    if ( unitAscendingNodeVector_.y( ) < 0.0 )
    {
        computedKeplerianElements_( longitudeOfAscendingNodeIndex ) =
                2.0 * PI - computedKeplerianElements_( longitudeOfAscendingNodeIndex );
    }

    // Compute and store argument of periapsis.
    // Define the quadrant condition for the true anomaly.
    double trueAnomalyQuandrantCondition = position_.dot( velocity_ );

    // Check if the orbit is circular. If it is, set the eccentricity vector to unit vector
    // pointing to the ascending node, i.e. set the argument of periapsis to zero.
    if ( std::fabs( computedKeplerianElements_( eccentricityIndex ) ) < tolerance )
    {
        eccentricityVector_ = unitAscendingNodeVector_;

        computedKeplerianElements_( argumentOfPeriapsisIndex ) = 0.0;

        // Check if orbit is also equatorial and set true anomaly quandrant check condition
        // accordingly.
        if ( unitAscendingNodeVector_ == Eigen::Vector3d::UnitX( ) )
        {
            // If the orbit is circular, position_.dot( velocity_ ) = 0, therefore this value
            // cannot be used as a quadrant condition. Moreover, if the orbit is equatorial,
            // position_.z( ) is also zero and therefore the quadrant condition is taken to be the
            // y-component, position_.y( ).
            trueAnomalyQuandrantCondition = position_.y( );
        }

        else
        {
            // If the orbit is circular, position_.dot( velocity_ ) = 0, therefore the quadrant
            // condition is taken to be the z-component of the position, position_.z( ).
            trueAnomalyQuandrantCondition = position_.z( );
        }
    }

    // Else, compute the argument of periapsis as the angle between the eccentricity vector and
    // the unit vector to the ascending node.
    else
    {
        computedKeplerianElements_( argumentOfPeriapsisIndex )
                = std::acos( eccentricityVector_.normalized( ).dot( unitAscendingNodeVector_ ) );

        // Check if the quadrant is correct.
        using mathematical_constants::PI;
        if ( argumentOfPeriapsisQuandrantCondition < 0.0 )
        {
           computedKeplerianElements_( argumentOfPeriapsisIndex ) =
                   2.0 * PI - computedKeplerianElements_( argumentOfPeriapsisIndex );
        }
    }

    // Compute dot-product of position and eccentricity vectors.
    double dotProductPositionAndEccentricityVectors
            = position_.normalized( ).dot( eccentricityVector_.normalized( ) );

    // Check if the dot-product is one of the limiting cases: 0.0 or 1.0
    // (within prescribed tolerance).
    if ( std::fabs( 1.0 - dotProductPositionAndEccentricityVectors ) < tolerance )
    {
        dotProductPositionAndEccentricityVectors = 1.0;
    }

    if ( std::fabs( dotProductPositionAndEccentricityVectors ) < tolerance )
    {
        dotProductPositionAndEccentricityVectors  = 0.0;
    }

    // Compute and store true anomaly.
    computedKeplerianElements_( trueAnomalyIndex )
            = std::acos( dotProductPositionAndEccentricityVectors );

    // Check if the quandrant is correct.
    if ( trueAnomalyQuandrantCondition < 0.0 )
    {
        computedKeplerianElements_( trueAnomalyIndex ) =
                2.0 * PI - computedKeplerianElements_( trueAnomalyIndex );
    }

    // Return converted Keplerian elements.
    return computedKeplerianElements_;
}

//! Convert true anomaly to (elliptic) eccentric anomaly.
double convertTrueAnomalyToEllipticalEccentricAnomaly( const double trueAnomaly,
                                                       const double eccentricity )
{
    if ( eccentricity >= 1.0 || eccentricity < 0.0 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid." ) ) );
    }

    else
    {
        using std::cos;
        using std::sqrt;

        // Declare and compute sine and cosine of eccentric anomaly.
        double sineOfEccentricAnomaly_ = sqrt( 1.0 - std::pow( eccentricity, 2.0 ) )
                * std::sin( trueAnomaly ) / ( 1.0 + eccentricity * cos( trueAnomaly ) );
        double cosineOfEccentricAnomaly_ = ( eccentricity + cos( trueAnomaly ) )
                / ( 1.0 + eccentricity * cos( trueAnomaly ) );

        // Return elliptic eccentric anomaly.
        return std::atan2( sineOfEccentricAnomaly_, cosineOfEccentricAnomaly_ );
    }
}

//! Convert true anomaly to hyperbolic eccentric anomaly.
double convertTrueAnomalyToHyperbolicEccentricAnomaly( const double trueAnomaly,
                                                       const double eccentricity )
{
    if ( eccentricity <= 1.0 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid." ) ) );
    }

    else
    {
        using std::cos;

        // Compute hyperbolic sine and hyperbolic cosine of hyperbolic eccentric anomaly.
        double hyperbolicSineOfHyperbolicEccentricAnomaly_
                = std::sqrt( std::pow( eccentricity, 2.0 ) - 1.0 )
                * std::sin( trueAnomaly ) / ( 1.0 + cos( trueAnomaly ) );

        double hyperbolicCosineOfHyperbolicEccentricAnomaly_
                = ( cos( trueAnomaly ) + eccentricity ) / ( 1.0 + cos( trueAnomaly ) );

        // Return hyperbolic eccentric anomaly.
        return boost::math::atanh( hyperbolicSineOfHyperbolicEccentricAnomaly_
                      / hyperbolicCosineOfHyperbolicEccentricAnomaly_ );
    }
}

//! Convert true anomaly to eccentric anomaly.
double convertTrueAnomalyToEccentricAnomaly( const double trueAnomaly,
                                             const double eccentricity )
{
    // Declare computed eccentric anomaly.
    double eccentricAnomaly_ = -0.0;

    // Check if eccentricity is invalid and throw an error if true.
    if ( eccentricity < 0.0 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid." ) ) );
    }

    // Check if orbit is parabolic and throw an error if true.
    else if ( std::fabs( eccentricity - 1.0 ) < std::numeric_limits< double >::epsilon( ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Parabolic orbits have not yet been implemented." ) ) );
    }

    // Check if orbit is elliptical and compute eccentric anomaly.
    else if ( eccentricity >= 0.0 && eccentricity < 1.0 )
    {
        eccentricAnomaly_ = convertTrueAnomalyToEllipticalEccentricAnomaly( trueAnomaly,
                                                                            eccentricity );
    }

    else if ( eccentricity > 1.0 )
    {
        eccentricAnomaly_ = convertTrueAnomalyToHyperbolicEccentricAnomaly( trueAnomaly,
                                                                            eccentricity );
    }

    // Return computed eccentric anomaly.
    return eccentricAnomaly_;
}

//! Convert (elliptic) eccentric anomaly to true anomaly.
double convertEllipticalEccentricAnomalyToTrueAnomaly( const double ellipticEccentricAnomaly,
                                                       const double eccentricity )
{
    if ( eccentricity >= 1.0 || eccentricity < 0.0 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid." ) ) );
    }

    else
    {
        using std::cos;
        using std::sqrt;

        // Compute sine and cosine of true anomaly.
        double sineOfTrueAnomaly_ = sqrt( 1.0 - std::pow( eccentricity, 2.0 ) ) *
                std::sin( ellipticEccentricAnomaly )
                / ( 1.0 - eccentricity * cos( ellipticEccentricAnomaly ) );

        double cosineOfTrueAnomaly_ = ( cos( ellipticEccentricAnomaly ) - eccentricity )
                / ( 1.0 - eccentricity * cos( ellipticEccentricAnomaly ) );

        // Return true anomaly.
        return atan2( sineOfTrueAnomaly_, cosineOfTrueAnomaly_  );
    }
}

//! Convert hyperbolic eccentric anomaly to true anomaly.
double convertHyperbolicEccentricAnomalyToTrueAnomaly( const double hyperbolicEccentricAnomaly,
                                                       const double eccentricity )
{
    if ( eccentricity <= 1.0 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid." ) ) );
    }

    else
    {
        using std::cosh;

        // Compute sine and cosine of true anomaly.
        double sineOfTrueAnomaly_
                = std::sqrt( std::pow( eccentricity, 2.0 ) - 1.0 )
                * std::sinh( hyperbolicEccentricAnomaly )
                / ( eccentricity * cosh( hyperbolicEccentricAnomaly ) - 1.0 );

        double cosineOfTrueAnomaly_
                = ( eccentricity - cosh( hyperbolicEccentricAnomaly ) )
                / ( eccentricity * cosh( hyperbolicEccentricAnomaly ) - 1.0 );

        // Return true anomaly.
        return std::atan2( sineOfTrueAnomaly_, cosineOfTrueAnomaly_ );
    }

}

//! Convert eccentric anomaly to true anomaly.
double convertEccentricAnomalyToTrueAnomaly( const double eccentricAnomaly,
                                             const double eccentricity )
{
    // Declare computed true anomaly.
    double trueAnomaly_ = -0.0;

    // Check if eccentricity is invalid and throw an error if true.
    if ( eccentricity < 0.0 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid." ) ) );
    }

    // Check if orbit is parabolic and throw an error if true.
    else if ( std::fabs( eccentricity - 1.0 ) < std::numeric_limits< double >::epsilon( ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Parabolic orbits have not yet been implemented." ) ) );
    }

    // Check if orbit is elliptical and compute true anomaly.
    else if ( eccentricity >= 0.0 && eccentricity < 1.0 )
    {
        trueAnomaly_ = convertEllipticalEccentricAnomalyToTrueAnomaly( eccentricAnomaly,
                                                                       eccentricity );
    }

    else if ( eccentricity > 1.0 )
    {
        trueAnomaly_ = convertHyperbolicEccentricAnomalyToTrueAnomaly( eccentricAnomaly,
                                                                       eccentricity );
    }

    // Return computed true anomaly.
    return trueAnomaly_;
}

//! Convert (elliptical) eccentric anomaly to mean anomaly.
double convertEllipticalEccentricAnomalyToMeanAnomaly( const double ellipticalEccentricAnomaly,
                                                       const double eccentricity )
{
    return ellipticalEccentricAnomaly - eccentricity * std::sin( ellipticalEccentricAnomaly );
}

//! Convert hyperbolic eccentric anomaly to mean anomaly.
double convertHyperbolicEccentricAnomalyToMeanAnomaly(
    const double hyperbolicEccentricAnomaly, const double eccentricity )
{
    return eccentricity * std::sinh( hyperbolicEccentricAnomaly ) - hyperbolicEccentricAnomaly;
}

//! Convert eccentric anomaly to mean anomaly.
double convertEccentricAnomalyToMeanAnomaly( const double eccentricAnomaly,
                                             const double eccentricity )
{
    // Declare computed mean anomaly.
    double meanAnomaly_ = -0.0;

    // Check if eccentricity is invalid and throw an error if true.
    if ( eccentricity < 0.0 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid." ) ) );
    }

    // Check if orbit is parabolic and throw an error if true.
    else if ( std::fabs( eccentricity - 1.0 ) < std::numeric_limits< double >::epsilon( ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Parabolic orbits have not yet been implemented." ) ) );
    }

    // Check if orbit is elliptical and compute true anomaly.
    else if ( eccentricity >= 0.0 && eccentricity < 1.0 )
    {
        meanAnomaly_ = convertEllipticalEccentricAnomalyToMeanAnomaly( eccentricAnomaly,
                                                                       eccentricity );
    }

    else if ( eccentricity > 1.0 )
    {
        meanAnomaly_ = convertHyperbolicEccentricAnomalyToMeanAnomaly( eccentricAnomaly,
                                                                       eccentricity );
    }

    // Return computed mean anomaly.
    return meanAnomaly_;
}

//! Convert elapsed time to (elliptical) mean anomaly change.
double convertElapsedTimeToEllipticalMeanAnomalyChange(
        const double elapsedTime, const double centralBodyGravitationalParameter,
        const double semiMajorAxis )
{
    // Check if semi-major axis is invalid and throw error if true.
    if ( semiMajorAxis < 0.0 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Semi-major axis is invalid." ) ) );
    }

    // Else return elliptical mean anomaly change.
    else
    {
        return std::sqrt( centralBodyGravitationalParameter
                          / std::pow( semiMajorAxis, 3.0 ) ) * elapsedTime;
    }
}

//! Convert elapsed time to mean anomaly change for hyperbolic orbits.
double convertElapsedTimeToHyperbolicMeanAnomalyChange(
        const double elapsedTime, const double centralBodyGravitationalParameter,
        const double semiMajorAxis )
{
    // Check if semi-major axis is invalid and throw error if true.
    if ( semiMajorAxis > 0.0 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Semi-major axis is invalid." ) ) );
    }

    // Else return hyperbolic mean anomaly change.
    else
    {
        return std::sqrt( centralBodyGravitationalParameter
                          / std::pow( -semiMajorAxis, 3.0 ) ) * elapsedTime;
    }
}

//! Convert elapsed time to mean anomaly change.
double convertElapsedTimeToMeanAnomalyChange(
        const double elapsedTime, const double centralBodyGravitationalParameter,
        const double semiMajorAxis )
{
    // Declare computed mean anomaly change.
    double meanAnomalyChange_ = -0.0;

    // Check if orbit is elliptical and compute mean anomaly change.
    if ( semiMajorAxis > 0.0 )
    {
        meanAnomalyChange_ = convertElapsedTimeToEllipticalMeanAnomalyChange(
                    elapsedTime, centralBodyGravitationalParameter, semiMajorAxis );
    }

    // Else orbit is hyperbolic; compute mean anomaly change.
    else if ( semiMajorAxis < 0.0 )
    {
        meanAnomalyChange_ = convertElapsedTimeToHyperbolicMeanAnomalyChange(
                    elapsedTime, centralBodyGravitationalParameter, semiMajorAxis );
    }

    // Return computed mean anomaly change.
    return meanAnomalyChange_;
}

//! Convert (elliptical) mean anomaly change to elapsed time.
double convertEllipticalMeanAnomalyChangeToElapsedTime(
        const double ellipticalMeanAnomalyChange, const double centralBodyGravitationalParameter,
        const double semiMajorAxis )
{
    // Check if semi-major axis is invalid and throw error if true.
    if ( semiMajorAxis < 0.0 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Semi-major axis is invalid." ) ) );
    }

    // Else return elapsed time.
    else
    {
        return ellipticalMeanAnomalyChange * std::sqrt( std::pow( semiMajorAxis, 3.0 )
                                                        / centralBodyGravitationalParameter );
    }
}

//! Convert hyperbolic mean anomaly change to elapsed time.
double convertHyperbolicMeanAnomalyChangeToElapsedTime(
        const double hyperbolicMeanAnomalyChange, const double centralBodyGravitationalParameter,
        const double semiMajorAxis )
{
    // Check if semi-major axis is invalid and throw error if true.
    if ( semiMajorAxis > 0.0 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Semi-major axis is invalid." ) ) );
    }

    // Else return elapsed time.
    else
    {
        return std::sqrt( std::pow( -semiMajorAxis, 3.0 )
                          / centralBodyGravitationalParameter ) * hyperbolicMeanAnomalyChange;
    }
}

//! Convert mean anomaly change to elapsed time.
double convertMeanAnomalyChangeToElapsedTime(
        const double meanAnomalyChange, const double centralBodyGravitationalParameter,
        const double semiMajorAxis )
{
    // Declare computed elapsed time.
    double elapsedTime_ = -0.0;

    // Check if orbit is elliptical and compute elapsed time.
    if ( semiMajorAxis > 0.0 )
    {
        elapsedTime_ = convertEllipticalMeanAnomalyChangeToElapsedTime(
                    meanAnomalyChange, centralBodyGravitationalParameter, semiMajorAxis );
    }

    // Else orbit is hyperbolic; compute elapsed time.
    else if ( semiMajorAxis < 0.0 )
    {
        elapsedTime_ = convertHyperbolicMeanAnomalyChangeToElapsedTime(
                    meanAnomalyChange, centralBodyGravitationalParameter, semiMajorAxis );
    }

    // Return computed elapsed time.
    return elapsedTime_;
}

//! Convert (elliptical) mean motion to semi-major axis.
double convertEllipticalMeanMotionToSemiMajorAxis(
        const double ellipticalMeanMotion, const double centralBodyGravitationalParameter )
{
    return std::pow( centralBodyGravitationalParameter
                   / std::pow( ellipticalMeanMotion, 2.0 ), 1.0 / 3.0 );
}

//! Convert semi-major axis to elliptical mean motion.
double convertSemiMajorAxisToEllipticalMeanMotion(
        const double semiMajorAxis, const double centralBodyGravitationalParameter )
{
    // Check if semi-major axis is invalid and throw error if true.
    if ( semiMajorAxis < 0.0 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Semi-major axis is invalid." ) ) );
    }

    // Else compute and return elliptical mean motion.
    {
        return std::sqrt( centralBodyGravitationalParameter / std::pow( semiMajorAxis, 3.0 ) );
    }
}

} // namespace orbital_element_conversions

} // namespace tudat
