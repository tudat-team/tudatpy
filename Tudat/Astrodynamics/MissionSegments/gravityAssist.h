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
 *      110117    E. Iorfida        Creation of code.
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
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120417    T. Secretin       Moved set functions to constructor.
 *      120508    P. Musegaas       The gravitational parameter is now passed as double, made some
 *                                  variables constant.
 *
 *    References
 *
 *    Gravity assist and swing-by are different words for the same thing. The delta-V that is
 *    computed for a powered swing-by has not been proven to be the optimum (lowest) to achieve the
 *    desired geometry of incoming and outgoing hyperbolic legs. Some literature research will have
 *    to be done to look at the alternatives.
 *
 *    For the moment in this code the smallestPeriapsisDistanceFactor is given as external input by
 *    the user, but in the future it should be part of the CelestialBody object.
 *
 *    Also, the velocity of the central body will need to be computed by the ephemeris code.
 *    At the moment the shape of the central body is a sphere segment, and the radius of the planet
 *    is set externally by the user. In the future it should be possible to get the radius of each
 *    planet directly from the CelestialBody class, by a link to GeometricShape class.
 *
 *    At the moment, this code uses a Newton-Raphson root finder by default. In the future it
 *    should be possible to apply for example the Halley method by using polymorphism.
 *
 */

#ifndef TUDAT_GRAVITY_ASSIST_H
#define TUDAT_GRAVITY_ASSIST_H

#include <iostream>

#include <Eigen/Core>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Bodies/celestialBody.h"
#include "Tudat/Astrodynamics/States/cartesianVelocityElements.h"
#include "Tudat/Mathematics/GeometricShapes/sphereSegment.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphson.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphsonAdaptor.h"

namespace tudat
{
namespace astrodynamics
{
namespace mission_segments
{

//! Gravity assist method class.
/*!
 * Gravity assist method class.
 */
class GravityAssist
{
public:

    //! Typedef for shared pointer to Newton-Raphson method.
    /*!
     * Typedef for shared pointer to Newton-Raphson method.
     */
    typedef boost::shared_ptr< NewtonRaphson > NewtonRaphsonPointer;

    //! Default constructor.
    /*!
     * Default constructor.
     */
    GravityAssist( const double centralBodyGravitationalParameter,
                   const double smallestPeriapsisDistance,
                   const Eigen::Vector3d& centralBodyVelocity,
                   const Eigen::Vector3d& incomingVelocity,
                   const Eigen::Vector3d& outgoingVelocity,
                   boost::shared_ptr< NewtonRaphson > newtonRaphson )
        : centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
          smallestPeriapsisDistance_( smallestPeriapsisDistance ),
          centralBodyVelocity_( centralBodyVelocity ),
          incomingVelocity_( incomingVelocity ),
          outgoingVelocity_( outgoingVelocity ),
          incomingHyperbolicExcessVelocity_( Eigen::Vector3d::Zero( ) ),
          outgoingHyperbolicExcessVelocity_( Eigen::Vector3d::Zero( ) ),
          newtonRaphson_( newtonRaphson ),
          deltaV_( TUDAT_NAN ),
          bendingAngle_( TUDAT_NAN ),
          incomingEccentricity_( TUDAT_NAN ),
          outgoingEccentricity_( TUDAT_NAN ),
          incomingSemiMajorAxis_( TUDAT_NAN ),
          outgoingSemiMajorAxis_( TUDAT_NAN ),
          bendingEffectDeltaV_( TUDAT_NAN ),
          velocityEffectDeltaV_( TUDAT_NAN )
    { }

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

    //! The gravitational parameter of the central body.
    /*!
     * The gravitational parameter of the central body involved in the swingby [m^3 s^-2].
     */
    const double centralBodyGravitationalParameter_;

    //! Smallest periapsisDistance.
    /*!
     * The smallest allowable periapsis distance. This is the radius of closest possible approach to
     * the planet. For maximal swing-by energy, this is the central body radius.
     */
    const double smallestPeriapsisDistance_;

    //! Velocity of the swing-by central body.
    /*!
     * Velocity vector of the central body involved in the swing-by.
     */
    const Eigen::Vector3d centralBodyVelocity_;

    //! Incoming velocity of object.
    /*!
     * Incoming velocity of object.
     */
    const Eigen::Vector3d incomingVelocity_;

    //! Outgoing velocity of object.
    /*!
     * Outgoing velocity of object.
     */
    const Eigen::Vector3d outgoingVelocity_;

    //! Hyperbolic excess velocity of the incoming leg.
    /*!
     * Hyperbolic excess velocity of the incoming leg.
     */
    Eigen::Vector3d incomingHyperbolicExcessVelocity_;

    //! Hyperbolic excess velocity of the outgoing leg.
    /*!
     * Hyperbolic excess velocity of the outgoing leg.
     */
    Eigen::Vector3d outgoingHyperbolicExcessVelocity_;

    //! Shared pointer to object of NewtonRaphson class.
    /*!
     * Shared pointer to object of NewtonRaphson class.
     */
    NewtonRaphsonPointer newtonRaphson_;

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
     * Defines first derivative of root-finder function for the velocity-effect delta-V.
     * \param incomingEccentricity Incoming eccentricity.
     * \return Value of first derivative of root-finder function at defined eccentricity.
     */
    double firstDerivativeVelocityEffectFunction( double& incomingEccentricity );
};

} // namespace mission_segments
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_GRAVITY_ASSIST_H
