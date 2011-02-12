/*! \file lambertTargeter.cpp
 *    Source file of the Lambert targeting solver implemented in Tudat.
 *
 *    Path              : /Astrodynamics/Trajectory/LambertTargeter/
 *    Version           : 19
 *    Check status      : Checked
 *
 *    Author            : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 11 November, 2010
 *    Last modified     : 08 February, 2011
 *
 *    References
 *      [1] R.H. Gooding, "A procedure for the solution of Lambert's orbital
 *      boundary-value problem", Celestial Mechanics and Dynamical Astronomy
 *      48:145-165, 1990
 *      [2] J. Melman, "Trajectory optimization for a mission to Neptune and
 *      Triton", MSc thesis report, TU Delft, 2007
 *      [3] E. Iorfida, "Unknown name", MSc thesis report, TU Delft, 2011,
 *      (Unpublished)
 *
 *    Notes
 *      The number of revolutions from departure to arrival body is zero
 *      by definition in this routine. This can be made user-defined later on.
 *      The resulting trajectories are in anti-clockwise direction.
 *      At the moment Newton-Raphson is used for finding the root. In the future
 *      this could be replaced by a polymorphic pointer to any root finder.
 *      At the moment there is a patch applied when Newton-Raphson does not
 *      converge for negative initial guesses. This has not been stated in the
 *      paper by Gooding, but this patch has been found by trial and error.
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
 *      101111    E. Iorfida    First creation of code.
 *      101111    E. Iorfida    Implementation of all the equations up to the
 *                              Newton method.
 *      101117    E. Iorfida    Velocities computations added.
 *      101126    E. Iorfida    Get/set codes deleted.
 *      101206    E. Iorfida    LambertTargetingElements class deleted,
 *                              added setInitialState, modified punctuation.
 *      101207    E. Iorfida    Set single variables, change variables names
 *                              in more understandable ones.
 *      101209    E. Iorfida    Corrected some coding errors.
 *      101213    E. Iorfida    Deleted lambertAngle, added numberOfRevolution,
 *                              modified implementation.
 *      101214    E. Iorfida    Implementation only for the case with
 *                              numberOfRevolution = 0.
 *      110113    E. Iorfida    Added necessary elements to build
 *                              pointer-to-member-function to RootFinderAlgorithms
 *                              and NewtonRaphsonMethod classes.
 *      110124    E. Iorfida    Added necessary piece of code to be able to use
 *                              the last version of Newton-Raphson code.
 *      110126    E. Iorfida    Initialized member functions.
 *      110130    J. Melman     Simplified variable names, e.g.,
 *                              'normOfVelocityVector' became 'speed'.
 *                              Requested references to specific formulas.
 *                              Also corrected 'tangential' to 'transverse'.
 *                              Simplified computation of radial unit vector.
 *                              Corrected computation of transverse heliocentric
 *                              velocity.
 *      110201    E. Iorfida    Added pointerToCelestialBody and modified
 *                              variable names (from heliocentric, to inertial).
 *                              Added patch for negative case of
 *                              initialLambertGuess_.
 *                              Added equations references.
 *      110206    E. Iorfida    Added unique function for Newton-Raphson method.
 *                              Added computeAbsoluteValue to the
 *                              initialLambertGuess_ for non-converging cases.
 *      110208    E. Iorfida    Added CartesianPositionElements objects as input
 *                              and CartesianVelocityElements objects as output.
 */

// Include statementes.
# include "lambertTargeter.h"

// Using directives.
using linear_algebra::determineAngleBetweenVectors;
using mathematics::raiseToIntegerPower;
using mathematics::computeAbsoluteValue;

//! Default constructor.
LambertTargeter::LambertTargeter( ):
        xParameter_( -0.0 ),
        normalizedTimeOfFlight_( -0.0 ),
        qParameter_( -0.0 )
{
    pointerToCartesianVelocityAtDeparture_ = &cartesianVelocityAtDeparture_;
    pointerToCartesianVelocityAtArrival_ = &cartesianVelocityAtArrival_;
}

//! Default destructor.
LambertTargeter::~LambertTargeter( )
{
}

//! Set position at departure.
void LambertTargeter::setPositionAtDeparture(
        CartesianPositionElements *pointerToCartesianPositionAtDeparture )
{
    pointerToCartesianPositionAtDeparture_ =
            pointerToCartesianPositionAtDeparture;
}

//! Set position at arrival.
void LambertTargeter::setPositionAtArrival(
        CartesianPositionElements *pointerToCartesianPositionAtArrival )
{
    pointerToCartesianPositionAtArrival_ = pointerToCartesianPositionAtArrival;
}

//! Set number of revolutions.
void LambertTargeter::setNumberOfRevolutions(
        const int &numberOfRevolutions )
{
    numberOfRevolutions_ = numberOfRevolutions;
}

//! Set time of flight.
void LambertTargeter::setTimeOfFlight( const double &timeOfFlight )
{
    timeOfFlight_ = timeOfFlight;
}

//! Set central body.
void LambertTargeter::setCentralBody(
        CelestialBody *pointerToCelestialBody )
{
    pointerToCelestialBody_ = pointerToCelestialBody;
}

//! Set pointer to Newton-Raphson method for Lambert targeting algorithm.
void LambertTargeter::setNewtonRaphsonMethod(
        NewtonRaphson *pointerToNewtonRaphson )
{
    // Set pointer to object of NewtonRaphson class.
    pointerToNewtonRaphson_ = pointerToNewtonRaphson;
}

//! Get semi-major axis.
double& LambertTargeter::getLambertSemiMajorAxis( )
{
    return lambertSemiMajorAxis_;
}

//! Get radial speed at departure.
double& LambertTargeter::getRadialSpeedAtDeparture( )
{
    return radialSpeedAtDeparture_;
}

//! Get radial speed at arrival.
double& LambertTargeter::getRadialSpeedAtArrival( )
{
    return radialSpeedAtArrival_;
}

//! Get transverse speed at departure.
double& LambertTargeter::
        getTransverseSpeedAtDeparture( )
{
    return transverseSpeedAtDeparture_;
}

//! Get transverse speed at arrival.
double& LambertTargeter::getTransverseSpeedAtArrival( )
{
    return transverseSpeedAtArrival_;
}

//! Get inertial velocity at departure.
CartesianVelocityElements* LambertTargeter::getInertialVelocityAtDeparture( )
{
    return pointerToCartesianVelocityAtDeparture_;
}

//! Get inertial velocity at arrival.
CartesianVelocityElements* LambertTargeter::getInertialVelocityAtArrival( )
{
    return pointerToCartesianVelocityAtArrival_;
}

//! Define Lambert function for positive lambertEccentricAnomaly_.
// Page 47 [2].
double LambertTargeter::lambertFunctionPositive( double& xParameter_ )
{
    double lambertEccentricAnomaly_ =
            raiseToIntegerPower( xParameter_ , 2 ) - 1.0;
    double yParameter_ = sqrt( computeAbsoluteValue(
            lambertEccentricAnomaly_ ) );
    double zParameter_ = sqrt( 1.0 - raiseToIntegerPower( qParameter_ , 2 ) +
            raiseToIntegerPower( qParameter_ * xParameter_ , 2 ) );
    double fParameter_ = yParameter_ * ( zParameter_ - qParameter_ *
            xParameter_ );
    double gParameter_ = xParameter_ * zParameter_ - qParameter_ *
            lambertEccentricAnomaly_;
    double dParameter_ = log( fParameter_ + gParameter_ );

    return normalizedTimeOfFlight_ -
           2.0 * ( xParameter_ - qParameter_ * zParameter_ - dParameter_ /
           yParameter_ ) / lambertEccentricAnomaly_;
}

//! Define Lambert function for negative lambertEccentricAnomaly_.
// Page 47 [2].
double LambertTargeter::lambertFunctionNegative( double& xParameter_ )
{
    double lambertEccentricAnomaly_ =
            raiseToIntegerPower( xParameter_ , 2 ) - 1.0;
    double yParameter_ = sqrt( computeAbsoluteValue(
            lambertEccentricAnomaly_ ) );
    double zParameter_ = sqrt( 1.0 - raiseToIntegerPower( qParameter_ , 2 ) +
            raiseToIntegerPower( qParameter_ * xParameter_ , 2 ) );
    double fParameter_ = yParameter_ * ( zParameter_ - qParameter_ *
            xParameter_ );
    double gParameter_ = xParameter_ * zParameter_ - qParameter_ *
            lambertEccentricAnomaly_;
    double dParameter_ = atan( fParameter_ / gParameter_ );

    return normalizedTimeOfFlight_ -
           2.0 * ( xParameter_ - qParameter_ * zParameter_ - dParameter_ /
           yParameter_ ) / lambertEccentricAnomaly_;
}

//! Define first derivative of Lambert function for positive
//! lambertEccentricAnomaly_.
// Page xx [3].
double LambertTargeter::lambertFirstDerivativeFunctionPositive(
         double& xParameter_ )
{
    double lambertEccentricAnomaly_ =
            raiseToIntegerPower( xParameter_ , 2 ) - 1.0;
    double yParameter_ = sqrt( computeAbsoluteValue(
            lambertEccentricAnomaly_ ) );
    double zParameter_ = sqrt( 1.0 - raiseToIntegerPower( qParameter_ , 2 ) +
            raiseToIntegerPower( qParameter_ * xParameter_ , 2 ) );
    double fParameter_ = yParameter_ * ( zParameter_ - qParameter_ *
            xParameter_ );
    double gParameter_ = xParameter_ * zParameter_ - qParameter_ *
            lambertEccentricAnomaly_;
    double dParameter_ = log( fParameter_ + gParameter_ );
    double lambertEccentricAnomalyDerivative_ = 2.0 * xParameter_;
    double yParameterDerivative_ = xParameter_ /
            sqrt( raiseToIntegerPower( xParameter_ , 2 ) - 1.0 );
    double zParameterDerivative_ = raiseToIntegerPower( qParameter_ , 2 ) *
            xParameter_ / zParameter_;
    double fParameterDerivative_ = yParameterDerivative_ * ( zParameter_ -
            xParameter_ * qParameter_ ) + yParameter_ *
            ( zParameterDerivative_ - qParameter_ );
    double gParameterDerivative_ = zParameter_ + xParameter_ *
            zParameterDerivative_ - qParameter_ *
            lambertEccentricAnomalyDerivative_;
    double dParameterDerivative_ = ( fParameterDerivative_ +
            gParameterDerivative_ ) / ( fParameter_ + gParameter_ );

    return - 2.0 * ( lambertEccentricAnomaly_ * ( 1.0 - qParameter_ *
           zParameterDerivative_ - ( ( dParameterDerivative_ * yParameter_ -
           yParameterDerivative_ * dParameter_ ) /
           raiseToIntegerPower( yParameter_ , 2 ) ) ) - ( ( xParameter_ -
           qParameter_ * zParameter_ - ( dParameter_ / yParameter_ ) ) *
           lambertEccentricAnomalyDerivative_ ) ) /
           ( raiseToIntegerPower( lambertEccentricAnomaly_ , 2 ) );
}

//! Define first derivative of Lambert function for negative
//! lambertEccentricAnomaly_.
// Page xx [3].
double LambertTargeter::lambertFirstDerivativeFunctionNegative(
         double& xParameter_ )
{
    double lambertEccentricAnomaly_ =
            raiseToIntegerPower( xParameter_ , 2 ) - 1.0;
    double yParameter_ = sqrt( computeAbsoluteValue(
            lambertEccentricAnomaly_ ) );
    double zParameter_ = sqrt( 1.0 - raiseToIntegerPower( qParameter_ , 2 ) +
            raiseToIntegerPower( qParameter_ * xParameter_ , 2 ) );
    double fParameter_ = yParameter_ * ( zParameter_ - qParameter_ *
            xParameter_ );
    double gParameter_ = xParameter_ * zParameter_ - qParameter_ *
            lambertEccentricAnomaly_;
    double dParameter_ = atan( fParameter_ / gParameter_ );
    double lambertEccentricAnomalyDerivative_ = 2.0 * xParameter_;
    double yParameterDerivative_ = -1.0 * xParameter_ /
            sqrt( 1.0 - raiseToIntegerPower( xParameter_ , 2 ) );
    double zParameterDerivative_ = raiseToIntegerPower( qParameter_ , 2 ) *
            xParameter_ / zParameter_;
    double fParameterDerivative_ = yParameterDerivative_ * ( zParameter_ -
            xParameter_ * qParameter_ ) + yParameter_ *
            ( zParameterDerivative_ - qParameter_ );
    double gParameterDerivative_ = zParameter_ + xParameter_ *
            zParameterDerivative_ - qParameter_ *
            lambertEccentricAnomalyDerivative_;
    double dParameterDerivative_ = ( fParameterDerivative_ * gParameter_ -
            fParameter_ * gParameterDerivative_ ) /
            ( raiseToIntegerPower( fParameter_ , 2 ) +
              raiseToIntegerPower( gParameter_ , 2 ) );

    return - 2.0 * ( lambertEccentricAnomaly_ * ( 1.0 - qParameter_ *
           zParameterDerivative_ - ( ( dParameterDerivative_ * yParameter_ -
           yParameterDerivative_ * dParameter_ ) /
           raiseToIntegerPower( yParameter_ , 2 ) ) ) - ( ( xParameter_ -
           qParameter_ * zParameter_ - ( dParameter_ / yParameter_ ) ) *
           lambertEccentricAnomalyDerivative_ ) ) /
           ( raiseToIntegerPower( lambertEccentricAnomaly_ , 2 ) );
}

//! Define general Lambert function.
double LambertTargeter::lambertFunction( double &xParameter_ )
{
    double lambertEccentricAnomaly_ =
            raiseToIntegerPower( xParameter_ , 2 ) - 1.0;

    if ( lambertEccentricAnomaly_ > 0.0 )
    {
        return LambertTargeter::lambertFunctionPositive( xParameter_);
    }
    else
    {
        return LambertTargeter::lambertFunctionNegative( xParameter_);
    }

}

//! Define first derivative of general Lambert function.
double LambertTargeter::lambertFirstDerivativeFunction( double &xParameter_ )
{
    double lambertEccentricAnomaly_ =
            raiseToIntegerPower( xParameter_ , 2 ) - 1.0;

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
    // Implement Lambert targeting method.

    // Define local position vectors.
    Vector3d positionAtDeparture_;
    positionAtDeparture_.x( ) = pointerToCartesianPositionAtDeparture_->
                                getCartesianElementX( );
    positionAtDeparture_.y( ) = pointerToCartesianPositionAtDeparture_->
                                getCartesianElementY( );
    positionAtDeparture_.z( ) = pointerToCartesianPositionAtDeparture_->
                                getCartesianElementZ( );

    Vector3d positionAtArrival_;
    positionAtArrival_.x( ) = pointerToCartesianPositionAtArrival_->
                              getCartesianElementX( );
    positionAtArrival_.y( ) = pointerToCartesianPositionAtArrival_->
                              getCartesianElementY( );
    positionAtArrival_.z( ) = pointerToCartesianPositionAtArrival_->
                              getCartesianElementZ( );

    // Normalize positions.
    double radiusAtDeparture_;
    double radiusAtArrival_;
    radiusAtDeparture_ = positionAtDeparture_.norm( );
    radiusAtArrival_  = positionAtArrival_.norm( );

    // Compute angle between positions.
    double reducedLambertAngle_;
    Vector3d planeNormal_ = positionAtDeparture_.
                            cross( positionAtArrival_ ).normalized( );
    reducedLambertAngle_ = determineAngleBetweenVectors( positionAtDeparture_,
                                                         positionAtArrival_ );

    if ( planeNormal_.z( ) < 0.0 )
    {
        reducedLambertAngle_ = 2.0 * M_PI - reducedLambertAngle_;
    }

    // Compute chord.
    double chord_;
    chord_ = sqrt( raiseToIntegerPower( radiusAtDeparture_, 2 ) +
                   raiseToIntegerPower( radiusAtArrival_, 2 ) -
                   2.0 * radiusAtDeparture_ * radiusAtArrival_ *
                   cos( reducedLambertAngle_ ) );

    // Compute semi-perimeter.
    double semiPerimeter_;
    semiPerimeter_ = ( radiusAtDeparture_ + radiusAtArrival_ + chord_ ) / 2.0;

    // Compute normalized time of flight.
    // Formula (7) [1].
    normalizedTimeOfFlight_ = sqrt( 8.0 * pointerToCelestialBody_->
            getGravitationalParameter( ) /
            raiseToIntegerPower( semiPerimeter_ , 3 ) ) * timeOfFlight_;

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
    zParameterInitial_ = sqrt( 1.0 - raiseToIntegerPower( qParameter_ , 2 ) );

    // Formula (6.11) [2].
    double lambertEccentricAnomalyInitial_ ;
    lambertEccentricAnomalyInitial_ = -1.0;

    // Formula (6.12) [2].
    double yParameterInitial_;
    yParameterInitial_ = sqrt(
            computeAbsoluteValue( lambertEccentricAnomalyInitial_ ) );

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
            dParameterInitial_ / yParameterInitial_ ) /
            lambertEccentricAnomalyInitial_;

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
        double phi = 2.0 * atan2( ( 1.0 -
                     raiseToIntegerPower( qParameter_ , 2 ) ) ,
                     2.0 * qParameter_ );

        // Formula (16) [1].
        double W = x01 + 1.7 * sqrt( 2.0 - phi / M_PI );

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
        lambdax = 1 + 0.5 * x03 * ( 1 + x01 ) - 0.03 *
                  raiseToIntegerPower( x03 , 2 ) * sqrt( 1.0 + x01 );

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
    pointerToNewtonRaphson_->setInitialGuessOfRoot( initialLambertGuess_ );

    // Set tolerance for Newton-Raphson method.
    // Page 155 [1].
    pointerToNewtonRaphson_->setTolerance( 1.0e-13 );

    // Set maximum number of iterations that can be computed in Newton-Raphson.
    pointerToNewtonRaphson_->setMaximumNumberOfIterations( 100 );

    // Set the adaptor for Newton-Raphson method.
    pointerToNewtonRaphson_->setNewtonRaphsonAdaptor(
            &newtonRaphsonAdaptorForLambertTargeter_ );

    // Execute Newton-Raphson method.
    pointerToNewtonRaphson_->execute( );

    // Check convergence Newton-Raphson method.
    if ( pointerToNewtonRaphson_->getComputedRootOfFunction( ) == -0.0 )
    {
        // Apparently Newton-Raphson has not converged.
        // Print message about second attempt.
        std::cerr << "A second attempt to converge is made." << std::endl;

        // Set initial guess as absolute value of initialLambertGuess_.
        // This has not been stated in the paper by Gooding, but this patch
        // has been found by trial and error.
        pointerToNewtonRaphson_->setInitialGuessOfRoot( computeAbsoluteValue(
                initialLambertGuess_ ) );
        pointerToNewtonRaphson_->execute( );

        // Check convergence Newton-Raphson method again.
        if ( pointerToNewtonRaphson_->getComputedRootOfFunction( ) != -0.0 )
        {
            // Apparently Newton-Raphson has converged this time.
            std::cerr << "The second attempt was successful." << std::endl;
        }
    }

    // Define xParameter_ as the output value of Newton-Raphson method.
    xParameter_ = pointerToNewtonRaphson_->getComputedRootOfFunction( );

    // Compute semi-major axis of the transfer orbit.
    lambertSemiMajorAxis_ = semiPerimeter_ / ( 2.0 *
            ( 1.0 - raiseToIntegerPower( xParameter_ , 2 ) ) );

    // Compute velocities at departure and at arrival.

    // Compute gamma, rho and sigma parameters, needed to compute the velocities.
    double lambertGamma_;
    double lambertRho_;
    double lambertSigma_;
    double zParameter_;

    // Formula (8) [1].
    zParameter_ = sqrt( 1.0 - raiseToIntegerPower( qParameter_, 2 ) +
            raiseToIntegerPower( qParameter_ , 2 ) *
            raiseToIntegerPower( xParameter_ , 2 ) );

    // Formula (6.23) [2].
    lambertGamma_ = sqrt(
            pointerToCelestialBody_->getGravitationalParameter( ) *
            semiPerimeter_ / 2.0 );

    // Formula (6.24) [2].
    lambertRho_ = ( radiusAtDeparture_ -
            radiusAtArrival_ ) / chord_;

    // Formula (6.25) [2].
    lambertSigma_ = 2.0 *( sqrt(
            radiusAtDeparture_ * radiusAtArrival_ /
            raiseToIntegerPower( chord_ , 2 ) ) ) *
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
    Vector3d radialUnitVectorAtDeparture_ = positionAtDeparture_.normalized( );
    Vector3d radialUnitVectorAtArrival_ = positionAtArrival_.normalized( );

    // Compute unit vector that is normal to the plane in which the
    // trajectory takes place, and points in the positive z-direction.
    if ( planeNormal_.z( ) < 0 )
    {
        planeNormal_ = -planeNormal_;
    }

    // Compute transverse unit vectors at departure and at arrival.
    Vector3d transverseUnitVectorAtDeparture_ = planeNormal_.cross(
            positionAtDeparture_.normalized( ) );
    Vector3d transverseUnitVectorAtArrival_   = planeNormal_.cross(
            positionAtArrival_.normalized( ) );

    // Compute radial inertial velocities at departure and
    // at arrival.
    Vector3d radialInertialVelocityAtDeparture_ =
            radialSpeedAtDeparture_ * radialUnitVectorAtDeparture_;
    Vector3d radialInertialVelocityAtArrival_ =
            radialSpeedAtArrival_ * radialUnitVectorAtArrival_;

    // Compute transverse heliocentric velocities at departure and
    // at arrival.
    Vector3d transverseInertialVelocityAtDeparture_ =
            transverseSpeedAtDeparture_ *
            transverseUnitVectorAtDeparture_;
    Vector3d transverseInertialVelocityAtArrival_ =
            transverseSpeedAtArrival_ *
            transverseUnitVectorAtArrival_;

    // Compute heliocentric velocities at departure and at
    // arrival.

    // Define local velocities.
    Vector3d velocityAtDeparture_;
    Vector3d velocityAtArrival_;
    velocityAtDeparture_ = radialInertialVelocityAtDeparture_ +
                           transverseInertialVelocityAtDeparture_;
    velocityAtArrival_ = radialInertialVelocityAtArrival_ +
                         transverseInertialVelocityAtArrival_;

    pointerToCartesianVelocityAtDeparture_->
            setCartesianElementXDot( velocityAtDeparture_.x( ) );
    pointerToCartesianVelocityAtDeparture_->
            setCartesianElementXDot( velocityAtDeparture_.y( ) );
    pointerToCartesianVelocityAtDeparture_->
            setCartesianElementXDot( velocityAtDeparture_.z( ) );

    pointerToCartesianVelocityAtArrival_->
            setCartesianElementXDot( velocityAtArrival_.x( ) );
    pointerToCartesianVelocityAtArrival_->
            setCartesianElementXDot( velocityAtArrival_.y( ) );
    pointerToCartesianVelocityAtArrival_->
            setCartesianElementXDot( velocityAtArrival_.z( ) );


}

// End of file.
