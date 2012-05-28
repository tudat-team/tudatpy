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
 *      101111    E. Iorfida        Creation of code.
 *      101111    E. Iorfida        Implementation of all the equations up to the Newton method.
 *      101117    E. Iorfida        Velocities computations added.
 *      101126    E. Iorfida        Get/set codes deleted.
 *      101206    E. Iorfida        LambertTargetingElements class deleted,
 *                                  added setInitialState, modified punctuation. Set single
 *                                  variables, change variables names in more understandable ones.
 *      101209    E. Iorfida        Corrected some coding errors.
 *      101213    E. Iorfida        Deleted lambertAngle, added numberOfRevolution, modified
 *                                  implementation.
 *      101214    E. Iorfida        Implementation only for the case with numberOfRevolution = 0.
 *      110113    E. Iorfida        Added necessary elements to build pointer-to-member-function
 *                                  to RootFinderAlgorithms and NewtonRaphsonMethod classes.
 *      110124    E. Iorfida        Added necessary piece of code to be able to use the last
 *                                  version of Newton-Raphson code.
 *      110126    E. Iorfida        Initialized member functions.
 *      110130    J. Melman         Simplified variable names, e.g., 'normOfdeletVelocityVector'
 *                                  became 'speed'. Requested references to specific formulas. Also
 *                                  corrected 'tangential' to 'transverse'. Simplified computation
 *                                  of radial unit vector. Corrected computation of transverse
 *                                  heliocentric velocity.
 *      110201    E. Iorfida        Added pointerToCelestialBody and modified variable names (from
 *                                  heliocentric, to inertial). Added patch for negative case of
 *                                  initialLambertGuess_. Added equations references.
 *      110206    E. Iorfida        Added unique function for Newton-Raphson method. Added
 *                                  computeAbsoluteValue to the initialLambertGuess_ for
 *                                  non-converging cases.
 *      110208    E. Iorfida        Added CartesianPositionElements objects as input and
 *                                  CartesianVelocityElements objects as output.
 *      110418    E. Iorfida        Added a new normal plane that take into account the case of two
 *                                  parallel position vector (with a relative angle of 180
 *                                  degrees). Better defined the pointers to the output
 *                                  CartesianVelocityElements.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *
 *    References        :
 *      Gooding, R.H. A procedure for the solution of Lambert's orbital boundary-value problem,
 *          Celestial Mechanics and Dynamical Astronomy, 48:145-165, 1990.
 *      Melman, J. Trajectory optimization for a mission to Neptune and Triton, MSc thesis report,
 *          Delft University of Technology, 2007.
 *      Iorfida, E. MSc thesis report, Unpublished.
 *
 *    The number of revolutions from departure to arrival body is zero by definition in this
 *    routine. This can be made user-defined later on. The resulting trajectories are in
 *    anti-clockwise direction.
 *
 *    At the moment Newton-Raphson is used for finding the root. In the future this could be
 *    replaced by a polymorphic pointer to any root finder. At the moment there is a patch applied
 *    when Newton-Raphson does not converge for negative initial guesses. This has not been stated
 *    in the paper by Gooding, but this patch has been found by trial and error.
 *
 */

#include <cmath>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>
#include <TudatCore/Mathematics/BasicMathematics/linearAlgebra.h>

#include "Tudat/Astrodynamics/MissionSegments/lambertTargeter.h"

namespace tudat
{
namespace astrodynamics
{
namespace mission_segments
{

// Using declarations.
using std::pow;
using std::sqrt;
using std::log;
using std::tan;
using std::atan;
using std::cos;
using std::sin;
using std::endl;
using mathematics::linear_algebra::computeAngleBetweenVectors;

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream, LambertTargeter& lambertTargeter )
{
    stream << "The position vector at departure is set to: "
           << lambertTargeter.cartesianPositionAtDeparture_
           << "The position vector at arrival is set to: "
           << lambertTargeter.cartesianPositionAtArrival_
           << "The velocity vector at departure is computed as: "
           << lambertTargeter.getInertialVelocityAtDeparture( )
           << "The velocity vector at departure is computed as: "
           << lambertTargeter.getInertialVelocityAtArrival( )
           << "The semi-major axis of the computed arc is: "
           << lambertTargeter.getLambertSemiMajorAxis( ) << endl;

    // Return stream.
    return stream;
}

//! Define Lambert function for positive lambertEccentricAnomaly_.
double LambertTargeter::lambertFunctionPositive( double& xParameter_ )
{
    double lambertEccentricAnomaly_ = pow( xParameter_ , 2.0 ) - 1.0;
    double yParameter_ = sqrt( fabs( lambertEccentricAnomaly_ ) );
    double zParameter_ = sqrt( 1.0 - pow( qParameter_ , 2.0 ) +
                               pow( qParameter_ * xParameter_ , 2.0 ) );
    double fParameter_ = yParameter_ * ( zParameter_ - qParameter_ * xParameter_ );
    double gParameter_ = xParameter_ * zParameter_ - qParameter_ * lambertEccentricAnomaly_;
    double dParameter_ = log( fParameter_ + gParameter_ );

    return normalizedTimeOfFlight_ -
           2.0 * ( xParameter_ - qParameter_ * zParameter_ - dParameter_ /
           yParameter_ ) / lambertEccentricAnomaly_;
}

//! Define Lambert function for negative lambertEccentricAnomaly_.
double LambertTargeter::lambertFunctionNegative( double& xParameter_ )
{
    double lambertEccentricAnomaly_ = pow( xParameter_, 2.0 ) - 1.0;
    double yParameter_ = sqrt( fabs( lambertEccentricAnomaly_ ) );
    double zParameter_ = sqrt( 1.0 - pow( qParameter_, 2.0 ) +
                               pow( qParameter_ * xParameter_, 2.0 ) );
    double fParameter_ = yParameter_ * ( zParameter_ - qParameter_ * xParameter_ );
    double gParameter_ = xParameter_ * zParameter_ - qParameter_ * lambertEccentricAnomaly_;
    double dParameter_ = atan( fParameter_ / gParameter_ );

    return normalizedTimeOfFlight_ - 2.0 * ( xParameter_ - qParameter_ * zParameter_ - dParameter_
                                             / yParameter_ ) / lambertEccentricAnomaly_;
}

//! Define first derivative of Lambert function for positive
//! lambertEccentricAnomaly_.
double LambertTargeter::lambertFirstDerivativeFunctionPositive(
         double& xParameter_ )
{
    double lambertEccentricAnomaly_ = pow( xParameter_, 2.0 ) - 1.0;
    double yParameter_ = sqrt( fabs( lambertEccentricAnomaly_ ) );
    double zParameter_ = sqrt( 1.0 - pow( qParameter_, 2.0 ) +
                               pow( qParameter_ * xParameter_, 2.0 ) );
    double fParameter_ = yParameter_ * ( zParameter_ - qParameter_ * xParameter_ );
    double gParameter_ = xParameter_ * zParameter_ - qParameter_ * lambertEccentricAnomaly_;
    double dParameter_ = log( fParameter_ + gParameter_ );
    double lambertEccentricAnomalyDerivative_ = 2.0 * xParameter_;
    double yParameterDerivative_ = xParameter_ / sqrt( pow( xParameter_ , 2.0 ) - 1.0 );
    double zParameterDerivative_ = pow( qParameter_ , 2.0 ) * xParameter_ / zParameter_;
    double fParameterDerivative_ = yParameterDerivative_  *
            ( zParameter_ - xParameter_ * qParameter_ ) + yParameter_ *
            ( zParameterDerivative_ - qParameter_ );
    double gParameterDerivative_ = zParameter_ + xParameter_ *
            zParameterDerivative_ - qParameter_ * lambertEccentricAnomalyDerivative_;
    double dParameterDerivative_ = ( fParameterDerivative_ + gParameterDerivative_ ) /
            ( fParameter_ + gParameter_ );

    return - 2.0 * ( lambertEccentricAnomaly_ * ( 1.0 - qParameter_ *
           zParameterDerivative_ - ( ( dParameterDerivative_ * yParameter_ -
           yParameterDerivative_ * dParameter_ ) / pow( yParameter_, 2.0 ) ) ) - ( ( xParameter_ -
           qParameter_ * zParameter_ - ( dParameter_ / yParameter_ ) ) *
           lambertEccentricAnomalyDerivative_ ) ) / ( pow( lambertEccentricAnomaly_, 2.0 ) );
}

//! Define first derivative of Lambert function for negative lambertEccentricAnomaly_.
double LambertTargeter::lambertFirstDerivativeFunctionNegative( double& xParameter_ )
{
    double lambertEccentricAnomaly_ = pow( xParameter_, 2.0 ) - 1.0;
    double yParameter_ = sqrt( fabs( lambertEccentricAnomaly_ ) );
    double zParameter_ = sqrt( 1.0 - pow( qParameter_, 2.0 ) +
                               pow( qParameter_ * xParameter_, 2.0 ) );
    double fParameter_ = yParameter_ * ( zParameter_ - qParameter_ * xParameter_ );
    double gParameter_ = xParameter_ * zParameter_ - qParameter_ * lambertEccentricAnomaly_;
    double dParameter_ = atan( fParameter_ / gParameter_ );
    double lambertEccentricAnomalyDerivative_ = 2.0 * xParameter_;
    double yParameterDerivative_ = -1.0 * xParameter_ / sqrt( 1.0 - pow( xParameter_, 2.0 ) );
    double zParameterDerivative_ = pow( qParameter_, 2.0 ) * xParameter_ / zParameter_;
    double fParameterDerivative_ = yParameterDerivative_ *
            ( zParameter_ - xParameter_ * qParameter_ ) + yParameter_ *
            ( zParameterDerivative_ - qParameter_ );
    double gParameterDerivative_ = zParameter_ + xParameter_ *
            zParameterDerivative_ - qParameter_ *  lambertEccentricAnomalyDerivative_;
    double dParameterDerivative_ = ( fParameterDerivative_ * gParameter_ -
            fParameter_ * gParameterDerivative_ ) /
            ( pow( fParameter_ , 2.0 ) + pow( gParameter_ , 2.0 ) );

    return - 2.0 * ( lambertEccentricAnomaly_ * ( 1.0 - qParameter_ *
           zParameterDerivative_ - ( ( dParameterDerivative_ * yParameter_ -
           yParameterDerivative_ * dParameter_ ) /
           pow( yParameter_, 2.0 ) ) ) - ( ( xParameter_ -
           qParameter_ * zParameter_ - ( dParameter_ / yParameter_ ) ) *
           lambertEccentricAnomalyDerivative_ ) ) / ( pow( lambertEccentricAnomaly_, 2.0 ) );
}

//! Define general Lambert function.
double LambertTargeter::lambertFunction( double &xParameter_ )
{
    double lambertEccentricAnomaly_ = pow( xParameter_, 2.0 ) - 1.0;

    if ( lambertEccentricAnomaly_ > 0.0 )
    {
        return lambertFunctionPositive( xParameter_);
    }

    else
    {
        return lambertFunctionNegative( xParameter_);
    }
}

//! Define first derivative of general Lambert function.
double LambertTargeter::lambertFirstDerivativeFunction( double &xParameter_ )
{
    double lambertEccentricAnomaly_ = pow( xParameter_, 2 ) - 1.0;

    if ( lambertEccentricAnomaly_ > 0.0 )
    {
        return lambertFirstDerivativeFunctionPositive( xParameter_ );
    }

    else
    {
        return lambertFirstDerivativeFunctionNegative( xParameter_ );
    }
}

//! Execute Lambert targeting solver.
void LambertTargeter::execute( )
{
    // Normalize positions.
    double radiusAtDeparture_;
    double radiusAtArrival_;
    radiusAtDeparture_ = cartesianPositionAtDeparture_.norm( );
    radiusAtArrival_  = cartesianPositionAtArrival_.norm( );

    // Compute angle between positions.
    double reducedLambertAngle_;

    Eigen::Vector3d planeNormalPosition_;
    planeNormalPosition_ = cartesianPositionAtDeparture_.
                  cross( cartesianPositionAtArrival_ ).normalized( );

    reducedLambertAngle_ = computeAngleBetweenVectors( cartesianPositionAtDeparture_,
                                                       cartesianPositionAtArrival_ );

    if ( planeNormalPosition_.z( ) < 0.0 )
    {
        reducedLambertAngle_ = 2.0 * mathematics::PI - reducedLambertAngle_;
    }

    // Compute chord.
    double chord_;
    chord_ = sqrt( pow( radiusAtDeparture_, 2.0 ) + pow( radiusAtArrival_, 2.0 ) -
                   2.0 * radiusAtDeparture_ * radiusAtArrival_ * cos( reducedLambertAngle_ ) );

    // Compute semi-perimeter.
    double semiPerimeter_;
    semiPerimeter_ = ( radiusAtDeparture_ + radiusAtArrival_ + chord_ ) / 2.0;

    // Compute normalized time of flight.
    // Formula (7) [1].
    normalizedTimeOfFlight_ = sqrt( 8.0 * celestialBody_->getGravityFieldModel( )->
            getGravitationalParameter( ) / pow( semiPerimeter_, 3.0 ) ) * timeOfFlight_;

    // Compute q-parameter.
    // Formula (5) [1].
    qParameter_ = ( sqrt( radiusAtDeparture_ * radiusAtArrival_ ) /
                    semiPerimeter_ ) * cos( reducedLambertAngle_ / 2.0 ) ;

    // Declare initial guess of x-parameter.
    double initialLambertGuess_;

    // Value of T(x) for x=0.
    double tFunctionInitial_;

    // Compute values of parameters needed to compute the initial guess of
    // x-parameter.
    double zParameterInitial_;
    zParameterInitial_ = sqrt( 1.0 - pow( qParameter_, 2.0 ) );

    // Formula (6.11) [2].
    double lambertEccentricAnomalyInitial_ ;
    lambertEccentricAnomalyInitial_ = -1.0;

    // Formula (6.12) [2].
    double yParameterInitial_;
    yParameterInitial_ = sqrt( fabs( lambertEccentricAnomalyInitial_ ) );

    // Formula (6.14) [2].
    double fParameterInitial_;
    fParameterInitial_ = yParameterInitial_ * zParameterInitial_;

    // Formula (6.15) [2].
    double gParameterInitial_;
    gParameterInitial_= -1.0 * qParameter_ * lambertEccentricAnomalyInitial_;

    // Page 47 [2].
    double dParameterInitial_;
    dParameterInitial_ = atan( fParameterInitial_ / gParameterInitial_ );

    // Formula (6.9) [2].
    tFunctionInitial_ = -2.0 * ( qParameter_ * zParameterInitial_ +
            dParameterInitial_ / yParameterInitial_ ) / lambertEccentricAnomalyInitial_;

    // Determine initial Lambert guess.
    if ( tFunctionInitial_ > normalizedTimeOfFlight_ )
    {
        // Formula (11) [1].
        initialLambertGuess_ = tFunctionInitial_ *
                ( tFunctionInitial_ - normalizedTimeOfFlight_ ) /
                ( 4.0 * normalizedTimeOfFlight_ );
    }
    else
    {
        // Formula (13) [1].
        double x01 = -1.0 * ( normalizedTimeOfFlight_ - tFunctionInitial_ ) /
                     ( normalizedTimeOfFlight_ - tFunctionInitial_ + 4.0 );

        // Formula (15) [1].
        double x02 = -1.0 * sqrt( ( normalizedTimeOfFlight_ -
                     tFunctionInitial_ ) / ( normalizedTimeOfFlight_ +
                     0.5 * tFunctionInitial_ ) );

        // Formula (20) [1].
        double phi = 2.0 * atan2( ( 1.0 - pow( qParameter_, 2.0 ) ), 2.0 * qParameter_ );

        // Formula (16) [1].
        double W = x01 + 1.7 * sqrt( 2.0 - phi / mathematics::PI );

        // Formula (17) [1].
        double x03;
        if ( W >= 0.0 )
        {
           x03 = x01;
        }
        else
        {
           x03 = x01 + pow( -W , 1/16 ) * ( x02 - x01 );
        }

        // Formula (19) [1].
        double lambdax;
        lambdax = 1 + 0.5 * x03 * ( 1 + x01 ) - 0.03 * pow( x03 , 2.0 ) * sqrt( 1.0 + x01 );

        // Formula (18) [1].
        initialLambertGuess_ = lambdax * x03;
    }

    // Newton-Raphson method implementation.
    // Set the class that contains the functions needed for Newton-Raphson.
    newtonRaphsonAdaptorForLambertTargeter_.setClass( this );

    // Set the functions needed for Newton-Raphson method.
    newtonRaphsonAdaptorForLambertTargeter_.setPointerToFunction(
            &LambertTargeter::lambertFunction );
    newtonRaphsonAdaptorForLambertTargeter_.setPointerToFirstDerivativeFunction(
            &LambertTargeter::lambertFirstDerivativeFunction );

    // Set initial guess of the variable computed in Newton-Rapshon method.
    // A patch for negative initialLambertGuess_ is applied.
    // This has not been stated in the paper by Gooding, but this patch
    // has been found by trial and error.
    if ( pow( initialLambertGuess_, 2.0 ) - 1.0 < 0.0 )
    {
        newtonRaphson_->setInitialGuessOfRoot( fabs( initialLambertGuess_ ) );
    }

    else
    {
        newtonRaphson_->setInitialGuessOfRoot( initialLambertGuess_ );
    }

    // Set tolerance for Newton-Raphson method.
    // Page 155 [1].
    newtonRaphson_->setTolerance( 1.0e-13 );

    // Set maximum number of iterations that can be computed in Newton-Raphson.
    newtonRaphson_->setMaximumNumberOfIterations( 100 );

    // Set the adaptor for Newton-Raphson method.
    newtonRaphson_->setNewtonRaphsonAdaptor(
            &newtonRaphsonAdaptorForLambertTargeter_ );

    // Execute Newton-Raphson method.
    newtonRaphson_->execute( );

    // Define xParameter_ as the output value of Newton-Raphson method.
    xParameter_ = newtonRaphson_->getComputedRootOfFunction( );

    // Compute semi-major axis of the transfer orbit.
    lambertSemiMajorAxis_ = semiPerimeter_ / ( 2.0 * ( 1.0 - pow( xParameter_ , 2.0 ) ) );

    // Compute velocities at departure and at arrival.

    // Compute gamma, rho and sigma parameters, needed to compute the velocities.
    double lambertGamma_;
    double lambertRho_;
    double lambertSigma_;
    double zParameter_;

    // Formula (8) [1].
    zParameter_ = sqrt( 1.0 - pow( qParameter_, 2.0 ) +
            pow( qParameter_, 2.0 ) * pow( xParameter_, 2.0 ) );

    // Formula (6.23) [2].
    lambertGamma_ = sqrt(
            celestialBody_->getGravityFieldModel( )->getGravitationalParameter( ) *
                semiPerimeter_ / 2.0 );

    // Formula (6.24) [2].
    lambertRho_ = ( radiusAtDeparture_ -
            radiusAtArrival_ ) / chord_;

    // Formula (6.25) [2].
    lambertSigma_ = 2.0 *( sqrt(
            radiusAtDeparture_ * radiusAtArrival_ / pow( chord_, 2.0 ) ) ) *
            sin( reducedLambertAngle_ / 2.0 );

    // Compute radial speeds at departure and at arrival.
    // Formula (23) [1].
    radialSpeedAtDeparture_ =
            lambertGamma_ * ( qParameter_ * zParameter_ - xParameter_ -
            lambertRho_ * ( qParameter_ * zParameter_ + xParameter_ ) ) /
            radiusAtDeparture_;

    // Formula (24) [1].
    radialSpeedAtArrival_ =
            -lambertGamma_ * ( qParameter_ * zParameter_ - xParameter_ +
            lambertRho_ * ( qParameter_ * zParameter_ + xParameter_ ) ) /
            radiusAtArrival_;

    // Compute transverse speeds at departure and at arrival.
    // Compute large part of formulas (25) and (26).
    double transverseSpeedHelper_ = lambertGamma_ * lambertSigma_ *
            ( zParameter_ + qParameter_ * xParameter_ );

    // Formula (25) [1].
    transverseSpeedAtDeparture_ = transverseSpeedHelper_ / radiusAtDeparture_;

    // Formula (26) [1].
    transverseSpeedAtArrival_ = transverseSpeedHelper_ / radiusAtArrival_;

    // Compute inertial velocities (those velocity are computed respect
    // the central body, so they can be either heliocentric or
    // planetocentric).

    // Compute radial unit vectors at departure and at arrival.
    Eigen::Vector3d radialUnitVectorAtDeparture_ = cartesianPositionAtDeparture_.normalized( );
    Eigen::Vector3d radialUnitVectorAtArrival_ = cartesianPositionAtArrival_.normalized( );

    Eigen::Vector3d unitZVector_;
    unitZVector_.x( ) = 0.0;
    unitZVector_.y( ) = 0.0;
    unitZVector_.z( ) = 1.0;

    Eigen::Vector3d planeNormal_;
    planeNormal_ =
            ( ( cartesianPositionAtDeparture_.cross( unitZVector_ ) ).cross(
                    cartesianPositionAtDeparture_ ) ).normalized( );

    // Compute unit vector that is normal to the plane in which the
    // trajectory takes place, and points in the positive z-direction.
    if ( planeNormal_.z( ) < 0.0 )
    {
        planeNormal_ = -planeNormal_;
    }

    // Compute transverse unit vectors at departure and at arrival.
    Eigen::Vector3d transverseUnitVectorAtDeparture_ = planeNormal_.cross(
            cartesianPositionAtDeparture_.normalized( ) );
    Eigen::Vector3d transverseUnitVectorAtArrival_   = planeNormal_.cross(
            cartesianPositionAtArrival_.normalized( ) );

    // Compute radial inertial velocities at departure and
    // at arrival.
    Eigen::Vector3d radialInertialVelocityAtDeparture_ =
            radialSpeedAtDeparture_ * radialUnitVectorAtDeparture_;
    Eigen::Vector3d radialInertialVelocityAtArrival_ =
            radialSpeedAtArrival_ * radialUnitVectorAtArrival_;

    // Compute transverse heliocentric velocities at departure and
    // at arrival.
    Eigen::Vector3d transverseInertialVelocityAtDeparture_ =
            transverseSpeedAtDeparture_ * transverseUnitVectorAtDeparture_;
    Eigen::Vector3d transverseInertialVelocityAtArrival_ =
            transverseSpeedAtArrival_ * transverseUnitVectorAtArrival_;

    // Compute heliocentric velocities at departure and at
    // arrival.

    // Define output velocities.
    cartesianVelocityAtDeparture_ = radialInertialVelocityAtDeparture_ +
                           transverseInertialVelocityAtDeparture_;
    cartesianVelocityAtArrival_ = radialInertialVelocityAtArrival_ +
                         transverseInertialVelocityAtArrival_;
}

} // namespace mission_segments
} // namespace astrodynamics
} // namespace tudat
