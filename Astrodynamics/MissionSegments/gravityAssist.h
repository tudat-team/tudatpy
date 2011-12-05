/*! \file gravityAssist.h
 *    This header file contains a base class for the gravity assist method.
 *
 *    Path              : /Astrodynamics/MissionSegments/
 *    Version           : 8
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
 *
 *    Notes
 *      Gravity assist and swing-by are different words for the same thing.
 *      The delta-V that is computed for a powered swing-by has not been
 *      proven to be the optimum (lowest) to achieve the desired geometry
 *      of incoming and outgoing hyperbolic legs.
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
 *      In the future it should be possible to apply, for example, the Halley
 *      method by using polymorphism.
 *
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
 *      110117    E. Iorfida        First creation of code.
 *      110128    E. Iorfida        Added member variables necessary for the unit tests.
 *      110202    J. Melman         Changed several names and suggested minor restructuring of
 *                                  class.
 *      110203    E. Iorfida        Changed some variable names and modified punctuation.
 *      110205    J. Melman         Removed the trailing underscores in some public variables.
 *                                  Changed and added some notes.
 *      110208    E. Iorfida        Added CartesianVelocityElements objects as input. Deleted
 *                                  inheritance from TrajectoryDesignMethod.
 *      110212    J. Melman         Made delta-V private. getDeltaV changed into computeDeltaV.
 *      110214    E. Iorfida        Deleted temporary centralBodyRadius, replaced by an element
 *                                  of GeometricShapes.
 */

#ifndef GRAVITYASSIST_H
#define GRAVITYASSIST_H

// Include statements.
#include <iostream>
#include "Astrodynamics/Bodies/celestialBody.h"
#include "Astrodynamics/States/cartesianVelocityElements.h"
#include "Mathematics/GeometricShapes/sphereSegment.h"
#include "Mathematics/RootFindingMethods/newtonRaphson.h"
#include "Mathematics/RootFindingMethods/newtonRaphsonAdaptor.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! Gravity assist method class.
/*!
 * Gravity assist method class.
 */
class GravityAssist
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    GravityAssist( ) : pointerToCentralBody_( NULL ), pointerToCentralBodySphere_( NULL ),
        centralBodyVelocity_( Vector3d::Zero( ) ), smallestPeriapsisDistanceFactor_( -0.0 ),
        pointerToIncomingVelocity_( NULL ), pointerToOutgoingVelocity_( NULL ),
        incomingHyperbolicExcessVelocity_ ( Vector3d::Zero( ) ),
        outgoingHyperbolicExcessVelocity_ ( Vector3d::Zero( ) ), deltaV_( -0.0 ),
        bendingAngle_( -0.0 ), incomingEccentricity_( -0.0 ), outgoingEccentricity_( -0.0 ),
        incomingSemiMajorAxis_( -0.0 ), outgoingSemiMajorAxis_( -0.0 ),
        bendingEffectDeltaV_( -0.0 ), velocityEffectDeltaV_( -0.0 ),
        pointerToNewtonRaphson_( NULL ) { }

    //! Set central body of the swing-by.
    /*!
     * Sets pointer to central body of the swing-by.
     * \param pointerToCentralBody Central body of the swing-by.
     */
    void setCentralBody( CelestialBody *pointerToCentralBody )
    { pointerToCentralBody_ = pointerToCentralBody; }

    // Its input will come from the ephemeris class.
    //! Set velocity of the swing-by central body.
    /*!
     * Sets the velocity of the central body involved in the swing-by,
     * which will come from the Ephemeris class.
     * \param centralBodyVelocity Velocity of swing-by central body.
     */
    void setCentralBodyVelocity( Vector3d centralBodyVelocity )
    { centralBodyVelocity_ = centralBodyVelocity; }

    // TEMPORARY!! Needs to be part of CelestialBody object.
    //! Set smallest periapsis distance factor.
    /*!
     * Sets the smallest allowable periapsis distance factor that has to be
     * multiplied with the central body radius to find the smallest allowable
     * periapsis distance.
     * \param smallestPeriapsisDistanceFactor Smallest periapsis distance factor.
     */
    void setSmallestPeriapsisDistanceFactor( double smallestPeriapsisDistanceFactor )
    { smallestPeriapsisDistanceFactor_ = smallestPeriapsisDistanceFactor; }

    //! Set pointer to incoming velocity of the satellite.
    /*!
     * Sets pointer to the incoming velocity vector of the satellite.
     * \param pointerToIncomingVelocity Pointer to incoming velocity of the
     *          satellite.
     */
    void setPointerToIncomingVelocity( CartesianVelocityElements* pointerToIncomingVelocity )
    { pointerToIncomingVelocity_ = pointerToIncomingVelocity; }

    //! Set pointer to outgoing velocity of the satellite.
    /*!
     * Sets pointer to the outgoing velocity of the satellite.
     * \param pointerToOutgoingVelocity Pointer to outgoing velocity of the satellite.
     */
    void setPointerToOutgoingVelocity( CartesianVelocityElements* pointerToOutgoingVelocity )
    { pointerToOutgoingVelocity_ = pointerToOutgoingVelocity; }

    //! Set pointer to Newton-Raphson method for gravity assist algorithm.
    /*!
     * Sets a pointer to the Newton-Raphson method for gravity assist
     * algorithm.
     * \param pointerToNewtonRaphson Pointer to NewtonRaphson object.
     */
    void setNewtonRaphsonMethod( NewtonRaphson *pointerToNewtonRaphson )
    { pointerToNewtonRaphson_ = pointerToNewtonRaphson; }

    //! Compute the delta-V of a powered swing-by.
    /*!
     * Computes the necessary delta-V to perform a powered swing-by.
     * \return deltaV Necessary delta-V of swing-by.
     */
    double computeDeltaV( );

    //! Overload ostream to print class information.
    /*!
     * Overloads ostream to print class information.
     * \param stream Stream object.
     * \param gravityAssist Gravity Assist.
     * \return Stream object.
     */
    friend std::ostream& operator<<( std::ostream& stream, GravityAssist& gravityAssist );

protected:

private:

    //! Pointer to CelestialBody class for swing-by.
    /*!
     * Pointer to CelestialBody class for swing-by.
     */
    CelestialBody* pointerToCentralBody_;

    //! Pointer to SphereSegment class for central body.
    /*!
     * Pointer to SphereSegment class for central body.
     */
    SphereSegment* pointerToCentralBodySphere_;

    //! Velocity of the swing-by central body.
    /*!
     * Velocity vector of the central body involved in the swing-by.
     */
    Vector3d centralBodyVelocity_;

    //! Smallest periapsisDistanceFactor.
    /*!
     * The smallest allowable periapsis distance factor
     * that is to be multiplied with the central body radius to find the
     * smallest allowable periapsis distance.
     */
    double smallestPeriapsisDistanceFactor_;

    //! Pointer to CartesianVelocityElements object.
    /*!
     * Pointer to CartesianVelocityElements object.
     */
    CartesianVelocityElements* pointerToIncomingVelocity_;

    //! Pointer to CartesianVelocityElements object.
    /*!
     * Pointer to CartesianVelocityElements object.
     */
    CartesianVelocityElements* pointerToOutgoingVelocity_;

    //! Hyperbolic excess velocity of the incoming leg.
    /*!
     * Hyperbolic excess velocity of the incoming leg.
     */
    Vector3d incomingHyperbolicExcessVelocity_;

    //! Hyperbolic excess velocity of the outgoing leg.
    /*!
     * Hyperbolic excess velocity of the outgoing leg.
     */
    Vector3d outgoingHyperbolicExcessVelocity_;

    //! Delta-V of powered gravity assist.
    /*!
     * Necessary delta-V to perform a powered swing-by.
     */
    double deltaV_;

    //! Bending angle between the excess velocities.
    /*!
     * Bending angle between the excess velocities.
     */
    double bendingAngle_;

    //! Eccentricity of the incoming hyperbolic leg.
    /*!
     * Eccentricity of the incoming hyperbolic leg.
     */
    double incomingEccentricity_;

    //! Eccentricity of the outgoing hyperbolic leg.
    /*!
     * Eccentricity of the outgoing hyperbolic leg.
     */
    double outgoingEccentricity_;

    //! Semi-major axis of the incoming hyperbolic leg.
    /*!
     *  Semi-major axis of the incoming hyperbolic leg.
     */
    double incomingSemiMajorAxis_;

    //! Semi-major axis of the outgoing hyperbolic leg.
    /*!
     * Semi-major axis of the outgoing hyperbolic leg.
     */
    double outgoingSemiMajorAxis_;

    //! Delta-V due to the effect of large bending angle.
    /*!
     * Delta-V due to the effect of large bending angle between excess
     * velocities.
     */
    double bendingEffectDeltaV_;

    //! Delta-V due to the effect of different excess speeds.
    /*!
     * Delta-V due to the effect of different excess speeds.
     */
    double velocityEffectDeltaV_;

    //! Pointer to object of NewtonRaphson class.
    /*!
     * Pointer to object of NewtonRaphson class.
     */
    NewtonRaphson* pointerToNewtonRaphson_;

    //! Pointer to adaptor object of NewtonRaphsonAdaptor class.
    /*!
     * Pointer to adaptor object of NewtonRaphsonAdaptor class. The template
     * parameter passed is this class.
     */
    NewtonRaphsonAdaptor< GravityAssist > newtonRaphsonAdaptorForGravityAssist_;

    //! Define root-finder function.
    /*!
     * Defines root-finder function for the velocity-effect delta-V.
     * \param incomingEccentricity Incoming eccentricity.
     * \return Value of root-finder function at defined eccentricity.
     */
    double velocityEffectFunction( double& incomingEccentricity );

    //! Define first derivative of root-finder function.
    /*!
     * Defines first derivative of root-finder function for the velocity-effect
     * delta-V.
     * \param incomingEccentricity_ Incoming eccentricity.
     * \return Value of first derivative of root-finder function at defined
     * eccentricity.
     */
    double firstDerivativeVelocityEffectFunction( double& incomingEccentricity );
};

}

#endif // GRAVITYASSIST_H

// End of file.
