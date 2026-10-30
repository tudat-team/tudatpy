/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_close_encounters.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <tudat/astro/close_encounters/bPlane.h>

namespace py = pybind11;
namespace tce = tudat::close_encounters;

namespace tudatpy
{
namespace astro
{
namespace close_encounters
{

void expose_close_encounters( py::module& m )
{
    py::class_< tce::BPlaneData >( m, "BPlaneData", R"doc(

 Data of the B-plane (target plane) of a hyperbolic encounter: coordinates and reference frame.

 The B-plane passes through the centre of the encountered body and is perpendicular to the incoming
 asymptote :math:`\hat{S}`. Its axes are :math:`\hat{T} = (\hat{V} \times \hat{S}) / |\hat{V} \times \hat{S}|`
 and :math:`\hat{R} = \hat{S} \times \hat{T}`, with :math:`\hat{V}` a reference direction, so that
 :math:`(\hat{T}, \hat{R}, \hat{S})` is a right-handed orthonormal triad. The B vector points from the centre of
 the encountered body to the point where the incoming asymptote crosses the B-plane.

 All vectors are expressed in the frame of the relative state from which the B-plane was computed. Objects of
 this class are created with :func:`compute_b_plane_coordinates`.

     )doc" )
            .def_readonly( "asymptotic_velocity", &tce::BPlaneData::asymptoticVelocity, R"doc(

 **read-only**

 Asymptotic velocity :math:`v_\infty` (magnitude of the velocity at infinity, or hyperbolic excess velocity) [m/s].

 :type: float
     )doc" )
            .def_readonly( "impact_parameter", &tce::BPlaneData::impactParameter, R"doc(

 **read-only**

 Impact parameter :math:`B = |\vec{B}|` [m].

 :type: float
     )doc" )
            .def_readonly( "b_vector", &tce::BPlaneData::bVector, R"doc(

 **read-only**

 B vector [m].

 :type: numpy.ndarray
     )doc" )
            .def_readonly( "s_unit_vector", &tce::BPlaneData::sUnitVector, R"doc(

 **read-only**

 Direction :math:`\hat{S}` of the incoming asymptote (direction of the incoming velocity at infinity).

 :type: numpy.ndarray
     )doc" )
            .def_readonly( "t_unit_vector", &tce::BPlaneData::tUnitVector, R"doc(

 **read-only**

 First axis :math:`\hat{T}` of the B-plane, perpendicular to the reference direction and to :math:`\hat{S}`.

 :type: numpy.ndarray
     )doc" )
            .def_readonly( "r_unit_vector", &tce::BPlaneData::rUnitVector, R"doc(

 **read-only**

 Second axis :math:`\hat{R} = \hat{S} \times \hat{T}` of the B-plane.

 :type: numpy.ndarray
     )doc" )
            .def_readonly( "b_t", &tce::BPlaneData::bT, R"doc(

 **read-only**

 Component :math:`B_T = \vec{B} \cdot \hat{T}` of the B vector [m].

 :type: float
     )doc" )
            .def_readonly( "b_r", &tce::BPlaneData::bR, R"doc(

 **read-only**

 Component :math:`B_R = \vec{B} \cdot \hat{R}` of the B vector [m].

 :type: float
     )doc" )
            .def_property_readonly( "rotation_matrix_to_b_plane_frame", &tce::BPlaneData::getRotationMatrixToBPlaneFrame, R"doc(

 **read-only**

 Rotation matrix from the frame of the relative state to the B-plane frame.

 The rows of the matrix are :math:`(\hat{T}, \hat{R}, \hat{S})`, so that a vector :math:`\vec{x}` has components
 :math:`(\vec{x} \cdot \hat{T}, \vec{x} \cdot \hat{R}, \vec{x} \cdot \hat{S})` in the B-plane frame. This can be
 used to project the B vector of a different trajectory on the axes of this B-plane.

 :type: numpy.ndarray
     )doc" );

    m.def( "compute_b_plane_coordinates",
           &tce::computeBPlaneCoordinates,
           py::arg( "relative_state" ),
           py::arg( "gravitational_parameter" ),
           py::arg( "reference_direction" ) = Eigen::Vector3d( -Eigen::Vector3d::UnitZ( ) ),
           R"doc(

 Compute the B-plane data of a hyperbolic relative state.

 The B-plane describes the encounter with the osculating two-body hyperbola that matches the relative state. The
 default reference direction is the :math:`-z` axis of the frame of the relative state. When the relative state is
 given with ICRF (J2000 equatorial) axes, :math:`\hat{T}` then lies in the equatorial plane, as in the B-plane
 definition of Kizner.

 Parameters
 ----------
 relative_state : numpy.ndarray
     Cartesian state of the small body relative to the encountered body [m, m/s].
 gravitational_parameter : float
     Gravitational parameter of the encountered body [m^3/s^2].
 reference_direction : numpy.ndarray, default = [0, 0, -1]
     Reference direction :math:`\vec{V}` that defines the B-plane axes, in the frame of the relative state (any
     norm). It must not be parallel to the incoming asymptote.

 Returns
 -------
 BPlaneData
     B-plane data (coordinates and reference frame) of the encounter.

 Raises
 ------
 RuntimeError
     If the relative orbit is not hyperbolic, or if the reference direction is zero or (almost) parallel to the
     incoming asymptote.

     )doc" );

    m.def( "compute_opik_coordinates",
           &tce::computeOpikCoordinates,
           py::arg( "relative_state" ),
           py::arg( "gravitational_parameter" ),
           py::arg( "encountered_body_heliocentric_state" ),
           R"doc(

 Compute the coordinates :math:`(\xi, \zeta)` of the Öpik theory on the B-plane of a hyperbolic relative state.

 The :math:`(\xi, \eta, \zeta)` frame of Valsecchi et al. (2003) is right-handed, with :math:`\eta` along the
 incoming asymptote and :math:`\zeta` opposite to the projection on the B-plane of the direction of motion of the
 encountered body. That direction is the :math:`Y` axis of the planetocentric Öpik frame (:math:`X` from the Sun to
 the encountered body, :math:`Y` in the direction of motion and perpendicular to :math:`X`,
 :math:`Z = X \times Y`), i.e. the transverse component
 :math:`\vec{v} - (\vec{v} \cdot \hat{r}) \hat{r}` of the heliocentric velocity of the encountered body
 (see :func:`compute_opik_reference_direction`). For an eccentric orbit of the encountered body this differs from
 the full heliocentric velocity by the flight path angle. Then :math:`\xi = \vec{B} \cdot \hat{T}` and
 :math:`\zeta = -\vec{B} \cdot \hat{R}`, with :math:`\hat{T}` and :math:`\hat{R}` computed with this
 direction as reference direction (see :func:`compute_b_plane_coordinates`).

 Parameters
 ----------
 relative_state : numpy.ndarray
     Cartesian state of the small body relative to the encountered body [m, m/s].
 gravitational_parameter : float
     Gravitational parameter of the encountered body [m^3/s^2].
 encountered_body_heliocentric_state : numpy.ndarray
     Heliocentric Cartesian state of the encountered body, in the frame orientation of the relative state [m, m/s].

 Returns
 -------
 numpy.ndarray
     Öpik coordinates :math:`(\xi, \zeta)` [m].

 Raises
 ------
 RuntimeError
     If the relative orbit is not hyperbolic, or if the heliocentric velocity of the encountered body has no
     component perpendicular to its position.

     )doc" );

    m.def( "compute_opik_reference_direction",
           &tce::computeOpikReferenceDirection,
           py::arg( "encountered_body_heliocentric_state" ),
           R"doc(

 Compute the direction of motion of the encountered body that defines the Öpik :math:`\zeta` axis.

 Parameters
 ----------
 encountered_body_heliocentric_state : numpy.ndarray
     Heliocentric Cartesian state of the encountered body [m, m/s].

 Returns
 -------
 numpy.ndarray
     Unit vector along the transverse component :math:`\vec{v} - (\vec{v} \cdot \hat{r}) \hat{r}` of the
     heliocentric velocity of the encountered body.

     )doc" );

    m.def( "compute_incoming_asymptote_direction",
           &tce::computeIncomingAsymptoteDirection,
           py::arg( "relative_state" ),
           py::arg( "gravitational_parameter" ),
           R"doc(

 Compute the direction of the incoming asymptote of a hyperbolic relative state.

 Parameters
 ----------
 relative_state : numpy.ndarray
     Cartesian state of the small body relative to the encountered body [m, m/s].
 gravitational_parameter : float
     Gravitational parameter of the encountered body [m^3/s^2].

 Returns
 -------
 numpy.ndarray
     Unit vector along the incoming velocity at infinity.

     )doc" );

    m.def( "compute_outgoing_asymptote_direction",
           &tce::computeOutgoingAsymptoteDirection,
           py::arg( "relative_state" ),
           py::arg( "gravitational_parameter" ),
           R"doc(

 Compute the direction of the outgoing asymptote of a hyperbolic relative state.

 Parameters
 ----------
 relative_state : numpy.ndarray
     Cartesian state of the small body relative to the encountered body [m, m/s].
 gravitational_parameter : float
     Gravitational parameter of the encountered body [m^3/s^2].

 Returns
 -------
 numpy.ndarray
     Unit vector along the outgoing velocity at infinity.

     )doc" );

    m.def( "compute_capture_radius",
           &tce::computeCaptureRadius,
           py::arg( "body_radius" ),
           py::arg( "gravitational_parameter" ),
           py::arg( "asymptotic_velocity" ),
           R"doc(

 Compute the radius of the capture cross section of a body on the B-plane.

 Due to gravitational focusing, any trajectory with impact parameter smaller than
 :math:`b = R \sqrt{1 + 2\mu / (R v_\infty^2)}` impacts the body.

 Parameters
 ----------
 body_radius : float
     Radius :math:`R` of the encountered body [m].
 gravitational_parameter : float
     Gravitational parameter :math:`\mu` of the encountered body [m^3/s^2].
 asymptotic_velocity : float
     Asymptotic velocity :math:`v_\infty` of the encounter [m/s].

 Returns
 -------
 float
     Radius of the capture cross section [m].

     )doc" );
}

}  // namespace close_encounters
}  // namespace astro
}  // namespace tudatpy
