/*! \file gravityAssist.cpp
 *    This source file contains the code of the gravity assist method.
 *
 *    Path              : /Astrodynamics/MissionSegments/GravityAssist/
 *    Version           : 7
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
 *    Date created      : 17 January, 2011
 *    Last modified     : 14 February, 2011
 *
 *    References
 *      Melman J. Trajectory optimization for a mission to Neptune and
 *          Triton, MSc thesis report, Delft University of Technology, 2007.
 *
 *    Notes
 *      Gravity assist and swing-by are different words for the same thing.
 *      The delta-V that is computed for a powered swing-by has not been
 *      proven to be the optimum (lowest) to achieve the desired geometry
 *      of incoming and outgoing hyperbolic legs. Some literature research
 *      will have to be done to look at the alternatives.
 *      For the moment in this code the smallestPeriapsisDistanceFactor
 *      is given as external input by the user, but in the future it should
 *      be part of the CelestialBody object.
 *      Also, the velocity of the central body will need to be computed by
 *      the ephemeris code.
 *      At the moment the shape of the central body is a sphere segment,
 *      and the radius of the planet is set externally by the user.
 *      In the future it should be possible to get the radius of each planet
 *      directly from the CelestialBody class, by a link to GeometricShape
 *      class.
 *      At the moment, this code uses a Newton-Raphson root finder by default.
 *      In the future it should be possible to apply for example the Halley
 *      method by using polymorphism.
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
 *      YYMMDD    Author            Comment
 *      110117    E. Iorfida        First creation of code.
 *      110128    E. Iorfida        Added boolean variable that sets the
 *                                  necessity of Newton-Raphson method,
 *                                  particularly in unit test.
 *      110202    J. Melman         Renamed certain parameters and added
 *                                  comments to clarify the code more.
 *      110203    E. Iorfida        Changed some variables names and modified
 *                                  punctuation.
 *      110205    J. Melman         Removed the trailing underscores in some
 *                                  public variables. Some comment rephrasing.
 *                                  Changed and added some notes.
 *      110212    J. Melman         Added a reference to my own thesis.
 *      110214    E. Iorfida        Deleted temporary centralBodyRadius,
 *                                  replaced by an element of GeometricShapes.
 */

// Include statements.
#include "gravityAssist.h"

// Using directives.
using linear_algebra::determineAngleBetweenVectors;
using mathematics::raiseToIntegerPower;
using mathematics::computeAbsoluteValue;
using unit_conversions::convertRadiansToDegrees;

//! Default constructor.
GravityAssist::GravityAssist( ) :
        incomingEccentricity( -0.0 ),
        outgoingEccentricity( -0.0 ),
        isRootFinderRequiredForChecking( false ),
        deltaV_( -0.0 ),
        bendingAngle_( -0.0 ),
        incomingSemiMajorAxis_( -0.0 ),
        outgoingSemiMajorAxis_( -0.0 ),
        bendingEffectDeltaV_( -0.0 ),
        velocityEffectDeltaV_( -0.0 )
{
}

//! Default destructor.
GravityAssist::~GravityAssist( )
{
}

//! Set central body of the swing-by.
void GravityAssist::setCentralBody( CelestialBody *pointerToCentralBody )
{
    pointerToCentralBody_ = pointerToCentralBody;
}

// Its input will come from the ephemeris class.
//! Set velocity of the swing-by central body.
void GravityAssist::setCentralBodyVelocity( Vector3d centralBodyVelocity )
{
    centralBodyVelocity_ = centralBodyVelocity;
}

// TEMPORARY!! Needs to be part of CelestialBody object.
//! Set smallest periapsis distance factor.
void GravityAssist::setSmallestPeriapsisDistanceFactor(
        const double &smallestPeriapsisDistanceFactor )
{
    smallestPeriapsisDistanceFactor_ = smallestPeriapsisDistanceFactor;
}

 //! Set pointer to incoming velocity of the satellite.
void GravityAssist::setPointerToIncomingVelocity(
        CartesianVelocityElements *pointerToIncomingVelocity )
{
    pointerToIncomingVelocity_ = pointerToIncomingVelocity;
}

//! Set pointer to outgoing velocity of the satellite.
void GravityAssist::setPointerToOutgoingVelocity(
        CartesianVelocityElements *pointerToOutgoingVelocity )
{
    pointerToOutgoingVelocity_ = pointerToOutgoingVelocity;
}

//! Set pointer to Newton-Raphson method for gravity assist algorithm.
void GravityAssist::setNewtonRaphsonMethod(
        NewtonRaphson *pointerToNewtonRaphson )
{
    // Set pointer to object of NewtonRaphson class.
    pointerToNewtonRaphson_ = pointerToNewtonRaphson;
}

//! Define root-finder function for the velocity-effect delta-V.
double GravityAssist::velocityEffectFunction( double &incomingEccentricity_ )
{
    return asin( 1.0 /incomingEccentricity_ ) + asin( 1.0 / ( 1.0 -
           incomingSemiMajorAxis_ / outgoingSemiMajorAxis_ *
           ( 1.0 -incomingEccentricity_ ) ) ) - bendingAngle_ ;
}

//! Define root-finder first-derivative function for the velocity-effect
//! delta-V.
double GravityAssist::firstDerivativeVelocityEffectFunction(
        double &incomingEccentricity_ )
{
    double eccentricitySquareMinusOne_ =
           raiseToIntegerPower( incomingEccentricity_, 2 ) - 1.0;
    double semiMajorAxisRatio_ = incomingSemiMajorAxis_ /
                                 outgoingSemiMajorAxis_ ;
    double bParameter_ = 1.0 - semiMajorAxisRatio_ * ( 1.0 -
            incomingEccentricity_ );

    return -1.0 / ( incomingEccentricity_ *
           sqrt( eccentricitySquareMinusOne_ ) ) -
           semiMajorAxisRatio_ / ( bParameter_ * sqrt(
           raiseToIntegerPower( bParameter_ , 2 ) - 1.0 ) );
}

//! Compute the delta-V of the powered swing-by.
const double& GravityAssist::computeDeltaV( )
{
    // Get shape model of central body.
    pointerToCentralBodySphere_ = static_cast< SphereSegment* >
                                  ( pointerToCentralBody_->getShapeModel( ) );

    // Define local velocity vectors.
    Vector3d incomingVelocity_;
    incomingVelocity_.x( ) =
            pointerToIncomingVelocity_->getCartesianElementXDot( );
    incomingVelocity_.y( ) =
            pointerToIncomingVelocity_->getCartesianElementYDot( );
    incomingVelocity_.z( ) =
            pointerToIncomingVelocity_->getCartesianElementZDot( );

    Vector3d outgoingVelocity_;
    outgoingVelocity_.x( ) =
            pointerToOutgoingVelocity_->getCartesianElementXDot( );
    outgoingVelocity_.y( ) =
            pointerToOutgoingVelocity_->getCartesianElementYDot( );
    outgoingVelocity_.z( ) =
            pointerToOutgoingVelocity_->getCartesianElementZDot( );

    // Compute incoming hyperbolic excess velocity.
    incomingHyperbolicExcessVelocity_ = incomingVelocity_ -
                                        centralBodyVelocity_;

    // Compute outgoing hyperbolic excess velocity.
    outgoingHyperbolicExcessVelocity_ = outgoingVelocity_ -
                                        centralBodyVelocity_;

    // Compute bending angle.
    bendingAngle_ = determineAngleBetweenVectors(
            incomingHyperbolicExcessVelocity_ ,
            outgoingHyperbolicExcessVelocity_ );

    // Compute smallest allowable periapsis distance.
    double smallestPeriapsisDistance;
    smallestPeriapsisDistance = pointerToCentralBodySphere_->getRadius( ) *
                                smallestPeriapsisDistanceFactor_;

    // Compute maximum achievable bending angle.
    double maximumBendingAngle_;
    maximumBendingAngle_ =
            asin( 1.0 / ( 1.0 + ( smallestPeriapsisDistance *
            incomingHyperbolicExcessVelocity_.squaredNorm( ) /
            pointerToCentralBody_->getGravitationalParameter( ) ) ) ) +
            asin( 1.0 / ( 1.0 + ( smallestPeriapsisDistance *
            outgoingHyperbolicExcessVelocity_.squaredNorm( ) /
            pointerToCentralBody_->getGravitationalParameter( ) ) ) );

    // Verify necessity to apply an extra swing-by delta-V due to the
    // incapability of the body's gravity to bend the trajectory a
    // sufficient amount given the incoming and outgoing velocities.
    // Compute required extra bending angle.
    double extraBendingAngle_;
    if ( bendingAngle_ > maximumBendingAngle_ )
    {
        extraBendingAngle_ = bendingAngle_ - maximumBendingAngle_;
    }
    else
    {
        extraBendingAngle_ = 0.0;
    }

    // Compute hyperbolic excess speeds.
    double incomingHyperbolicExcessSpeed_ =
            incomingHyperbolicExcessVelocity_.norm( );
    double outgoingHyperbolicExcessSpeed_ =
            outgoingHyperbolicExcessVelocity_.norm( );

    // Compute necessary delta-V due to bending-effect.
    bendingEffectDeltaV_ = 2.0 * std::min(
            incomingHyperbolicExcessSpeed_ ,
            outgoingHyperbolicExcessSpeed_ ) *
            sin( extraBendingAngle_ / 2.0 );

    // An excess speed difference of less than 10 cm/s is deemed to be
    // negligible.
    double speedTolerance_ = 1.0e-01;

    // Verify necessity to apply a swing-by delta-V due to the effect of
    // a difference in excess speeds. For checking purposes one might want
    // to let the root finder do his work anyway; setting the boolean below
    // equal to true serves that purpose.
    if ( computeAbsoluteValue( incomingHyperbolicExcessSpeed_ -
         outgoingHyperbolicExcessSpeed_ ) <= speedTolerance_ &&
         !isRootFinderRequiredForChecking )
    {
        // Set delta-V due to velocity effect equal to zero.
        velocityEffectDeltaV_ = 0.0;
    }
    else
    {
        // Compute semi-major axis of hyperbolic legs.
        incomingSemiMajorAxis_ =
                -1.0 * pointerToCentralBody_->getGravitationalParameter( ) /
                incomingHyperbolicExcessVelocity_.squaredNorm( );
        outgoingSemiMajorAxis_ =
                -1.0 * pointerToCentralBody_->getGravitationalParameter( ) /
                outgoingHyperbolicExcessVelocity_.squaredNorm( );

        // Newton-Raphson method implementation.
        // Set the class that contains the functions needed for Newton-Raphson.
        newtonRaphsonAdaptorForGravityAssist_.setClass( this );

        // Set the functions needed for Newton-Raphson method.
        newtonRaphsonAdaptorForGravityAssist_.setPointerToFunction(
                &GravityAssist::velocityEffectFunction );
        newtonRaphsonAdaptorForGravityAssist_.
                setPointerToFirstDerivativeFunction(
                &GravityAssist::firstDerivativeVelocityEffectFunction );

        // Set initial guess of the variable computed in Newton-Rapshon method.
        pointerToNewtonRaphson_->setInitialGuessOfRoot( 1.01 );

        // Set maximum number of iterations that can be computed in
        // Newton-Raphson.
        pointerToNewtonRaphson_->setMaximumNumberOfIterations( 100 );

        // Set tolerance for Newton-Raphson method.
        pointerToNewtonRaphson_->setTolerance( 10.0 * mathematics::
                                               MACHINE_PRECISION_DOUBLES );

        // Set the adaptor for Newton-Raphson method.
        pointerToNewtonRaphson_->setNewtonRaphsonAdaptor(
                &newtonRaphsonAdaptorForGravityAssist_ );

        // Execute Newton-Raphson method.
        pointerToNewtonRaphson_->execute( );

        // Define incoming hyperbolic leg eccentricity as the output value of
        // Newton-Raphson method.
        incomingEccentricity =
                pointerToNewtonRaphson_->getComputedRootOfFunction( );

        // Compute outgoing hyperbolic leg eccentricity.
        outgoingEccentricity = 1.0 - ( incomingSemiMajorAxis_ /
                outgoingSemiMajorAxis_ ) * ( 1.0 - incomingEccentricity );

        // Compute incoming and outgoing velocities at periapsis.
        double incomingVelocityAtPeriapsis_;
        double outgoingVelocityAtPeriapsis_;

        incomingVelocityAtPeriapsis_ = incomingHyperbolicExcessSpeed_ *
                sqrt( ( incomingEccentricity + 1.0 ) /
                      ( incomingEccentricity - 1.0 ) );
        outgoingVelocityAtPeriapsis_ = outgoingHyperbolicExcessSpeed_ *
                sqrt( ( outgoingEccentricity + 1.0 ) /
                      ( outgoingEccentricity - 1.0 ) );

        // Compute necessary delta-V due to velocity-effect.
        velocityEffectDeltaV_ =
                computeAbsoluteValue( outgoingVelocityAtPeriapsis_ -
                                      incomingVelocityAtPeriapsis_ );
    }

    // Compute total delta-V.
    deltaV_ = bendingEffectDeltaV_ + velocityEffectDeltaV_;

    // Return DeltaV.
    return deltaV_;
}

// End of file.
