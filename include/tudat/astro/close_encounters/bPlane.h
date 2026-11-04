/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Kizner, W. A method of describing miss distances for lunar and interplanetary trajectories.
 *          Planetary and Space Science, 7, 125-131, 1961.
 *      Valsecchi, G.B., Milani, A., Gronchi, G.F., Chesley, S.R. Resonant returns to close approaches:
 *          Analytical theory. Astronomy & Astrophysics, 408, 1179-1196, 2003.
 *
 *    Notes
 *      All functions describe the encounter with the osculating two-body hyperbola that matches the
 *      given relative state. The input state and all output vectors are expressed in the same
 *      (inertial) frame orientation.
 */

#ifndef TUDAT_B_PLANE_H
#define TUDAT_B_PLANE_H

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"

namespace tudat
{
namespace close_encounters
{

//! Data of the B-plane (target plane) of a hyperbolic encounter: coordinates and reference frame.
/*!
 *  The B-plane passes through the centre of the encountered body and is perpendicular to the incoming
 *  asymptote S. Its axes are T = (V x S) / |V x S| and R = S x T, with V a reference direction, so that
 *  (T, R, S) is a right-handed orthonormal triad. The B vector points from the centre of the encountered
 *  body to the point where the incoming asymptote crosses the B-plane.
 */
struct BPlaneData {
    //! Asymptotic velocity v_inf (magnitude of the velocity at infinity, or hyperbolic excess velocity) [m/s].
    double asymptoticVelocity;

    //! Impact parameter |B| [m].
    double impactParameter;

    //! B vector [m].
    Eigen::Vector3d bVector;

    //! Direction of the incoming asymptote (direction of the incoming velocity at infinity).
    Eigen::Vector3d sUnitVector;

    //! First axis of the B-plane, perpendicular to the reference direction and to S.
    Eigen::Vector3d tUnitVector;

    //! Second axis of the B-plane, R = S x T.
    Eigen::Vector3d rUnitVector;

    //! Component of the B vector along T [m].
    double bT;

    //! Component of the B vector along R [m].
    double bR;

    //! Rotation matrix from the frame of the input state to the B-plane frame.
    /*!
     *  The rows of the matrix are (T, R, S), so that a vector x in the frame of the input state has
     *  components (x.T, x.R, x.S) in the B-plane frame. This can be used to project the B vector of a
     *  different trajectory on the axes of this B-plane.
     *  \return Rotation matrix from the frame of the input state to the B-plane frame.
     */
    Eigen::Matrix3d getRotationMatrixToBPlaneFrame( ) const;
};

//! Compute the direction of the incoming asymptote of a hyperbolic relative state.
/*!
 *  \param relativeState Cartesian state of the small body relative to the encountered body [m, m/s].
 *  \param centralBodyGravitationalParameter Gravitational parameter of the encountered body [m^3/s^2].
 *  \return Unit vector along the incoming velocity at infinity.
 */
Eigen::Vector3d computeIncomingAsymptoteDirection( const Eigen::Vector6d& relativeState, const double centralBodyGravitationalParameter );

//! Compute the direction of the outgoing asymptote of a hyperbolic relative state.
/*!
 *  \param relativeState Cartesian state of the small body relative to the encountered body [m, m/s].
 *  \param centralBodyGravitationalParameter Gravitational parameter of the encountered body [m^3/s^2].
 *  \return Unit vector along the outgoing velocity at infinity.
 */
Eigen::Vector3d computeOutgoingAsymptoteDirection( const Eigen::Vector6d& relativeState, const double centralBodyGravitationalParameter );

//! Compute the B-plane data of a hyperbolic relative state.
/*!
 *  The default reference direction is the -z axis of the frame of the input state. When the input state is
 *  given with ICRF (J2000 equatorial) axes, T then lies in the equatorial plane, as in the B-plane definition
 *  of Kizner (also used by GMAT in the Earth MJ2000 equatorial frame).
 *  \param relativeState Cartesian state of the small body relative to the encountered body [m, m/s].
 *  \param centralBodyGravitationalParameter Gravitational parameter of the encountered body [m^3/s^2].
 *  \param referenceDirection Reference direction V that defines the B-plane axes (any norm). It must not be
 *  parallel to the incoming asymptote.
 *  \return B-plane data (coordinates and reference frame) of the encounter.
 */
BPlaneData computeBPlaneCoordinates( const Eigen::Vector6d& relativeState,
                                     const double centralBodyGravitationalParameter,
                                     const Eigen::Vector3d& referenceDirection = -Eigen::Vector3d::UnitZ( ) );

//! Compute the coordinates (xi, zeta) of the Opik theory on the B-plane of a hyperbolic relative state.
/*!
 *  The (xi, eta, zeta) frame of Valsecchi et al. (2003) is right-handed, with eta along the incoming asymptote
 *  and zeta opposite to the projection on the B-plane of the direction of motion of the encountered body. That
 *  direction is the Y axis of the planetocentric Opik frame (X from the Sun to the encountered body, Y in the
 *  direction of motion and perpendicular to X, Z = X x Y), i.e. the transverse component of the heliocentric
 *  velocity of the encountered body, v - (v.r_hat) r_hat. For a circular orbit of the encountered body it equals
 *  the full heliocentric velocity; for an eccentric orbit the two differ by the flight path angle (up to ~e rad).
 *  Then xi = B.T and zeta = -B.R, with T and R computed with this transverse direction as reference direction.
 *  \param relativeState Cartesian state of the small body relative to the encountered body [m, m/s].
 *  \param centralBodyGravitationalParameter Gravitational parameter of the encountered body [m^3/s^2].
 *  \param encounteredBodyHeliocentricState Heliocentric Cartesian state of the encountered body, in the frame
 *  orientation of the relative state [m, m/s].
 *  \return Opik coordinates (xi, zeta) [m].
 */
Eigen::Vector2d computeOpikCoordinates( const Eigen::Vector6d& relativeState,
                                        const double centralBodyGravitationalParameter,
                                        const Eigen::Vector6d& encounteredBodyHeliocentricState );

//! Compute the direction of motion of the encountered body that defines the Opik zeta axis.
/*!
 *  \param encounteredBodyHeliocentricState Heliocentric Cartesian state of the encountered body [m, m/s].
 *  \return Unit vector along the transverse component of the heliocentric velocity, v - (v.r_hat) r_hat.
 */
Eigen::Vector3d computeOpikReferenceDirection( const Eigen::Vector6d& encounteredBodyHeliocentricState );

//! Compute the radius of the capture cross section of a body on the B-plane.
/*!
 *  Due to gravitational focusing, any trajectory with impact parameter smaller than this radius impacts the body.
 *  \param bodyRadius Radius of the encountered body [m].
 *  \param centralBodyGravitationalParameter Gravitational parameter of the encountered body [m^3/s^2].
 *  \param asymptoticVelocity Asymptotic velocity v_inf of the encounter [m/s].
 *  \return Radius of the capture cross section [m].
 */
double computeCaptureRadius( const double bodyRadius, const double centralBodyGravitationalParameter, const double asymptoticVelocity );

}  // namespace close_encounters

}  // namespace tudat

#endif  // TUDAT_B_PLANE_H
