/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define PYBIND11_DETAILED_ERROR_MESSAGES

#include "expose_two_body_dynamics.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <tudat/astro/basic_astro.h>
#include <tudat/astro/mission_segments.h>

namespace py = pybind11;
namespace tms = tudat::mission_segments;
namespace trf = tudat::root_finders;
namespace toec = tudat::orbital_element_conversions;

namespace tudatpy
{

namespace trampoline
{

class PyLambertTargeter : public tms::LambertTargeter
{
public:
    /* Inherit the constructors */
    using tms::LambertTargeter::LambertTargeter;

    void execute( ) override
    {
        PYBIND11_OVERLOAD_PURE( void,                 /* Return type */
                                tms::LambertTargeter, /* Parent class */
                                execute,              /* Name of function in C++ (must
                                                         match Python name) */
                                                      /* Argument(s) */
        );
    }

    /* Trampoline (need one for each virtual function) */
    Eigen::Vector3d get_departure_velocity( )
    {
        PYBIND11_OVERLOAD( Eigen::Vector3d,                /* Return type */
                           tms::LambertTargeter,           /* Parent class */
                           getInertialVelocityAtDeparture, /* Name of function in C++
                                                              (must match Python name)
                                                            */
                                                           /* Argument(s) */
        );
    }
};

}  // namespace trampoline

namespace astro
{
namespace two_body_dynamics
{
void expose_two_body_dynamics( py::module& m )
{
    m.def( "compute_escape_or_capture_delta_v",
           &tms::computeEscapeOrCaptureDeltaV,
           py::arg( "gravitational_param" ),
           py::arg( "semi_major_axis" ),
           py::arg( "eccentricity" ),
           py::arg( "excess_velocity" ),
           R"doc(

 Compute the escape or capture delta-v budget for a spacecraft.

 This function calculates the required change in velocity (delta-v) for a spacecraft to escape from or
 be captured by the gravitational influence of a central body. The calculation is based on the pericenter
 of the orbit, the orbital parameters, and the excess velocity of the spacecraft. It is commonly used in
 mission design for estimating propulsion requirements in orbital transfers or interplanetary trajectories.

 Parameters
 ----------
 gravitational_parameter : float
     Gravitational parameter of the central body, defined as the product of the gravitational constant (G)
     and the mass of the body (M).
 semi_major_axis : float
     Semi-major axis of the spacecraft's orbit, representing the average distance from the central body.
 eccentricity : float
     Eccentricity of the spacecraft's orbit, which defines its shape. Must be valid for elliptical
     or hyperbolic orbits (e.g., 0 <= eccentricity < 1 for elliptical orbits).
 excess_velocity : float
     Excess velocity of the spacecraft, representing its velocity relative to the central body
     at infinity.

 Returns
 -------
 deltaV : float
     The delta-v required for the escape or capture maneuver. This is the difference between the velocity
     needed to achieve the specified excess velocity at infinity and the current orbital velocity at the
     pericenter.


     )doc" );

    py::class_< tms::PericenterFindingFunctions, std::shared_ptr< tms::PericenterFindingFunctions > >( m, "PericenterFindingFunctions" )
            .def( py::init< const double, const double, const double >( ),
                  py::arg( "absolute_incoming_semi_major_axis" ),
                  py::arg( "absolute_outgoing_semi_major_axis" ),
                  py::arg( "bending_angle" ) )
            .def( "compute_pericenter_radius_fn", &tms::PericenterFindingFunctions::computePericenterRadiusFunction )
            .def( "compute_derivative_pericenter_radius_fn",
                  &tms::PericenterFindingFunctions::computeFirstDerivativePericenterRadiusFunction );

    py::class_< tms::EccentricityFindingFunctions, std::shared_ptr< tms::EccentricityFindingFunctions > >( m,
                                                                                                           "EccentricityFindingFunctions",
                                                                                                           R"doc(
Class containing functions for finding the eccentricity during a gravity assist.

This class provides the objective function and its derivative for solving the incoming 
eccentricity in a gravity assist maneuver, given the incoming and outgoing hyperbolic 
semi-major axes and the bending angle. Used as input for root-finding algorithms.

)doc" )
            .def( py::init< const double, const double, const double >( ),
                  py::arg( "absolute_incoming_semi_major_axis" ),
                  py::arg( "absolute_outgoing_semi_major_axis" ),
                  py::arg( "bending_angle" ),
                  R"doc(
Constructor for EccentricityFindingFunctions.

Parameters
----------
absolute_incoming_semi_major_axis : float
    Absolute value of the semi-major axis of the incoming hyperbolic trajectory [m].
absolute_outgoing_semi_major_axis : float
    Absolute value of the semi-major axis of the outgoing hyperbolic trajectory [m].
absolute_bending_angle : float
    Bending angle of the gravity assist [rad].

)doc" )
            .def( "compute_incoming_eccentricity_fn",
                  &tms::EccentricityFindingFunctions::computeIncomingEccentricityFunction,
                  py::arg( "incoming_eccentricity" ),
                  R"doc(
Compute the incoming eccentricity function value.

Parameters
----------
incoming_eccentricity : float
    Incoming eccentricity for which to evaluate the function [-].

Returns
-------
float
    Function value.

)doc" )
            .def( "compute_derivative_incoming_eccentricity_fn",
                  &tms::EccentricityFindingFunctions::computeFirstDerivativeIncomingEccentricityFunction,
                  py::arg( "incoming_eccentricity" ),
                  R"doc(
Compute the first derivative of the incoming eccentricity function.

Parameters
----------
incoming_eccentricity : float
    Incoming eccentricity for which to evaluate the derivative [-].

Returns
-------
float
    First derivative of the function.

)doc" );

    py::class_< tms::LambertTargeter, std::shared_ptr< tms::LambertTargeter >, trampoline::PyLambertTargeter >( m,
                                                                                                                "LambertTargeter",
                                                                                                                R"doc(
Base class for Lambert targeting algorithms.

This abstract base class defines the interface for Lambert problem solvers. The Lambert 
problem consists of finding the orbit that connects two position vectors in a specified 
time of flight. Derived classes implement specific solution algorithms (e.g., Gooding, Izzo).

)doc" )
            // This class required a trampoline class to inherit from
            // due to the virtual ~execute member.
            .def( py::init< const Eigen::Vector3d&, const Eigen::Vector3d&, const double, const double >( ),
                  py::arg( "departure_position" ),
                  py::arg( "arrival_position" ),
                  py::arg( "time_of_flight" ),
                  py::arg( "gravitational_parameter" ),
                  R"doc(
Constructor for LambertTargeter.

Parameters
----------
departure_position : numpy.ndarray
    Cartesian position vector at departure [m].
arrival_position : numpy.ndarray
    Cartesian position vector at arrival [m].
time_of_flight : float
    Time of flight between departure and arrival [s].
gravitational_parameter : float
    Gravitational parameter of the central body [m^3/s^2].

)doc" )
            .def( "get_departure_velocity",
                  &tms::LambertTargeter::getInertialVelocityAtDeparture,
                  R"doc(
Get the inertial velocity at departure.

Returns
-------
numpy.ndarray
    Cartesian velocity vector at departure [m/s].

)doc" )
            .def( "get_arrival_velocity",
                  &tms::LambertTargeter::getInertialVelocityAtArrival,
                  R"doc(
Get the inertial velocity at arrival.

Returns
-------
numpy.ndarray
    Cartesian velocity vector at arrival [m/s].

)doc" )
            .def( "get_velocity_vectors",
                  &tms::LambertTargeter::getInertialVelocityVectors,
                  R"doc(
Get both velocity vectors as a pair.

Returns
-------
tuple[numpy.ndarray, numpy.ndarray]
    Tuple containing the departure and arrival velocity vectors [m/s].

)doc" );

    //////////////////////////////////////////////////////////////////////
    //  lambertTargeterGooding.cpp (complete)
    //////////////////////////////////////////////////////////////////////
    py::class_< tms::LambertTargeterGooding, std::shared_ptr< tms::LambertTargeterGooding >, tms::LambertTargeter >(
            m, "LambertTargeterGooding" )
            .def( py::init< const Eigen::Vector3d&, const Eigen::Vector3d&, const double, const double, trf::RootFinderPointer >( ),
                  py::arg( "departure_position" ),
                  py::arg( "arrival_position" ),
                  py::arg( "time_of_flight" ),
                  py::arg( "gravitational_parameter" ),
                  py::arg( "root_finder" ) = trf::RootFinderPointer( ),
                  R"doc(
Constructor for LambertTargeterGooding.

Parameters
----------
departure_position : numpy.ndarray
    Cartesian position vector at departure [m].
arrival_position : numpy.ndarray
    Cartesian position vector at arrival [m].
time_of_flight : float
    Time of flight between departure and arrival [s].
gravitational_parameter : float
    Gravitational parameter of the central body [m^3/s^2].
root_finder : RootFinder, default=None
    Root finder to use for solving the Lambert equation. If None, a default Newton-Raphson 
    solver with 1000 iterations and 1e-12 relative tolerance is used.

)doc" )
            .def( "get_radial_departure_velocity",
                  &tms::LambertTargeterGooding::getRadialVelocityAtDeparture,
                  R"doc(
Get the radial velocity component at departure.

Returns
-------
float
    Radial velocity at departure [m/s].

)doc" )
            .def( "get_radial_arrival_velocity",
                  &tms::LambertTargeterGooding::getRadialVelocityAtArrival,
                  R"doc(
Get the radial velocity component at arrival.

Returns
-------
float
    Radial velocity at arrival [m/s].

)doc" )
            .def( "get_transverse_departure_velocity",
                  &tms::LambertTargeterGooding::getTransverseVelocityAtDeparture,
                  R"doc(
Get the transverse velocity component at departure.

Returns
-------
float
    Transverse velocity at departure [m/s].

)doc" )
            .def( "get_transverse_arrival_velocity",
                  &tms::LambertTargeterGooding::getTransverseVelocityAtArrival,
                  R"doc(
Get the transverse velocity component at arrival.

Returns
-------
float
    Transverse velocity at arrival [m/s].

)doc" )
            .def( "get_semi_major_axis",
                  &tms::LambertTargeterGooding::getSemiMajorAxis,
                  R"doc(
Get the semi-major axis of the transfer orbit.

Returns
-------
float
    Semi-major axis of the conic section connecting departure and arrival [m].

)doc" );

    //////////////////////////////////////////////////////////////////////
    //  lambertTargeterIzzo.cpp (complete)
    //////////////////////////////////////////////////////////////////////
    py::class_< tms::LambertTargeterIzzo, std::shared_ptr< tms::LambertTargeterIzzo >, tms::LambertTargeter >( m,
                                                                                                               "LambertTargeterIzzo",
                                                                                                               R"doc(
Lambert targeter using Izzo's algorithm.

Implementation of Izzo's Lambert targeting algorithm. This method is particularly robust 
for near-pi transfers and does not suffer from singularities that affect other methods. 
It supports both prograde and retrograde orbits.

References
----------
Izzo, D., "Revisiting Lambert's problem", Celestial Mechanics and Dynamical Astronomy, 
Vol. 121, 2015.

)doc" )
            .def( py::init< const Eigen::Vector3d&,
                            const Eigen::Vector3d&,
                            const double,
                            const double,
                            const bool,
                            const double,
                            const int >( ),
                  py::arg( "departure_position" ),
                  py::arg( "arrival_position" ),
                  py::arg( "time_of_flight" ),
                  py::arg( "gravitational_parameter" ),
                  py::arg( "is_retrograde" ) = false,
                  py::arg( "tolerance" ) = 1e-9,
                  py::arg( "max_iter" ) = 50,
                  R"doc(
Constructor for LambertTargeterIzzo.

Parameters
----------
departure_position : numpy.ndarray
    Cartesian position vector at departure [m].
arrival_position : numpy.ndarray
    Cartesian position vector at arrival [m].
time_of_flight : float
    Time of flight between departure and arrival [s].
gravitational_parameter : float
    Gravitational parameter of the central body [m^3/s^2].
is_retrograde : bool, default=False
    If True, computes retrograde orbit; if False, computes prograde orbit.
tolerance : float, default=1e-9
    Convergence tolerance for the iterative solution.
max_iter : int, default=50
    Maximum number of iterations for the solution procedure.

)doc" )
            .def( "get_radial_departure_velocity",
                  &tms::LambertTargeterIzzo::getRadialVelocityAtDeparture,
                  R"doc(
Get the radial velocity component at departure.

Returns
-------
float
    Radial velocity at departure [m/s].

)doc" )
            .def( "get_radial_arrival_velocity",
                  &tms::LambertTargeterIzzo::getRadialVelocityAtArrival,
                  R"doc(
Get the radial velocity component at arrival.

Returns
-------
float
    Radial velocity at arrival [m/s].

)doc" )
            .def( "get_transverse_departure_velocity",
                  &tms::LambertTargeterIzzo::getTransverseVelocityAtDeparture,
                  R"doc(
Get the transverse velocity component at departure.

Returns
-------
float
    Transverse velocity at departure [m/s].

)doc" )
            .def( "get_transverse_arrival_velocity",
                  &tms::LambertTargeterIzzo::getTransverseVelocityAtArrival,
                  R"doc(
Get the transverse velocity component at arrival.

Returns
-------
float
    Transverse velocity at arrival [m/s].

)doc" )
            .def( "get_semi_major_axis",
                  &tms::LambertTargeterIzzo::getSemiMajorAxis,
                  R"doc(
Get the semi-major axis of the transfer orbit.

Returns
-------
float
    Semi-major axis of the conic section connecting departure and arrival [m].

)doc" );

    //////////////////////////////////////////////////////////////////////
    //  zeroRevolutionLambertTargeterIzzo.cpp (complete)
    //////////////////////////////////////////////////////////////////////
    py::class_< tms::ZeroRevolutionLambertTargeterIzzo, std::shared_ptr< tms::ZeroRevolutionLambertTargeterIzzo >, tms::LambertTargeter >(
            m,
            "ZeroRevolutionLambertTargeterIzzo",
            R"doc(
Zero-revolution Lambert targeter using Izzo's algorithm.

Specialized implementation of Izzo's algorithm for zero-revolution transfers (direct transfers 
without completing full orbits). This is a more focused version that handles the most common 
case efficiently.

)doc" )
            .def( py::init< const Eigen::Vector3d&,
                            const Eigen::Vector3d&,
                            const double,
                            const double,
                            const bool,
                            const double,
                            const int >( ),
                  py::arg( "departure_position" ),
                  py::arg( "arrival_position" ),
                  py::arg( "time_of_flight" ),
                  py::arg( "gravitational_parameter" ),
                  py::arg( "is_retrograde" ) = false,
                  py::arg( "tolerance" ) = 1e-9,
                  py::arg( "max_iter" ) = 50,
                  R"doc(
Constructor for ZeroRevolutionLambertTargeterIzzo.

Parameters
----------
departure_position : numpy.ndarray
    Cartesian position vector at departure [m].
arrival_position : numpy.ndarray
    Cartesian position vector at arrival [m].
time_of_flight : float
    Time of flight between departure and arrival [s].
gravitational_parameter : float
    Gravitational parameter of the central body [m^3/s^2].
is_retrograde : bool, default=False
    If True, computes retrograde orbit; if False, computes prograde orbit.
tolerance : float, default=1e-9
    Convergence tolerance for the iterative solution.
max_iter : int, default=50
    Maximum number of iterations for the solution procedure.

)doc" )
            .def( "get_radial_departure_velocity",
                  &tms::ZeroRevolutionLambertTargeterIzzo::getRadialVelocityAtDeparture,
                  R"doc(
Get the radial velocity component at departure.

Returns
-------
float
    Radial velocity at departure [m/s].

)doc" )
            .def( "get_radial_arrival_velocity",
                  &tms::ZeroRevolutionLambertTargeterIzzo::getRadialVelocityAtArrival,
                  R"doc(
Get the radial velocity component at arrival.

Returns
-------
float
    Radial velocity at arrival [m/s].

)doc" )
            .def( "get_transverse_departure_velocity",
                  &tms::ZeroRevolutionLambertTargeterIzzo::getTransverseVelocityAtDeparture,
                  R"doc(
Get the transverse velocity component at departure.

Returns
-------
float
    Transverse velocity at departure [m/s].

)doc" )
            .def( "get_transverse_arrival_velocity",
                  &tms::ZeroRevolutionLambertTargeterIzzo::getTransverseVelocityAtArrival,
                  R"doc(
Get the transverse velocity component at arrival.

Returns
-------
float
    Transverse velocity at arrival [m/s].

)doc" )
            .def( "get_semi_major_axis",
                  &tms::ZeroRevolutionLambertTargeterIzzo::getSemiMajorAxis,
                  R"doc(
Get the semi-major axis of the transfer orbit.

Returns
-------
float
    Semi-major axis of the conic section connecting departure and arrival [m].

)doc" );

    //////////////////////////////////////////////////////////////////////
    //  multiRevolutionLambertTargeterIzzo.cpp (complete)
    //////////////////////////////////////////////////////////////////////
    py::class_< tms::MultiRevolutionLambertTargeterIzzo,
                std::shared_ptr< tms::MultiRevolutionLambertTargeterIzzo >,
                tms::ZeroRevolutionLambertTargeterIzzo >( m,
                                                          "MultiRevolutionLambertTargeterIzzo",
                                                          R"doc(
Multi-revolution Lambert targeter using Izzo's algorithm.

Extension of Izzo's algorithm to handle multi-revolution transfers. Supports computing 
solutions for transfers that complete one or more full orbits before arrival, with both 
left-branch and right-branch solutions available.

)doc" )
            .def( py::init< const Eigen::Vector3d&,
                            const Eigen::Vector3d&,
                            const double,
                            const double,
                            const int,
                            const bool,
                            const bool,
                            const double,
                            const int >( ),
                  py::arg( "departure_position" ),
                  py::arg( "arrival_position" ),
                  py::arg( "time_of_flight" ),
                  py::arg( "gravitational_parameter" ),
                  py::arg( "n_revolutions" ) = 0,
                  py::arg( "is_right_branch" ) = false,
                  py::arg( "is_retrograde" ) = false,
                  py::arg( "tolerance" ) = 1e-9,
                  py::arg( "max_iter" ) = 50,
                  R"doc(
Constructor for MultiRevolutionLambertTargeterIzzo.

Parameters
----------
departure_position : numpy.ndarray
    Cartesian position vector at departure [m].
arrival_position : numpy.ndarray
    Cartesian position vector at arrival [m].
time_of_flight : float
    Time of flight between departure and arrival [s].
gravitational_parameter : float
    Gravitational parameter of the central body [m^3/s^2].
n_revolutions : int, default=0
    Number of complete revolutions before arrival.
is_right_branch : bool, default=False
    If True, uses right branch solution; if False, uses left branch solution.
is_retrograde : bool, default=False
    If True, computes retrograde orbit; if False, computes prograde orbit.
tolerance : float, default=1e-9
    Convergence tolerance for the iterative solution.
max_iter : int, default=50
    Maximum number of iterations for the solution procedure.

)doc" )
            //        .def("NO_MAXIMUM_REVOLUTIONS",
            //        &tms::MultiRevolutionLambertTargeterIzzo::NO_MAXIMUM_REVOLUTIONS)
            .def( "compute_for_revolutions_and_branch",
                  &tms::MultiRevolutionLambertTargeterIzzo::computeForRevolutionsAndBranch,
                  py::arg( "n_revolutions" ),
                  py::arg( "is_right_branch" ),
                  R"doc(
Compute the Lambert solution for specified revolutions and branch.

Parameters
----------
n_revolutions : int
    Number of complete revolutions.
is_right_branch : bool
    If True, uses right branch solution; if False, uses left branch solution.

)doc" )
            .def( "get_maximum_number_of_revolutions",
                  &tms::MultiRevolutionLambertTargeterIzzo::getMaximumNumberOfRevolutions,
                  R"doc(
Get the maximum number of revolutions for which a solution exists.

Returns
-------
int
    Maximum number of revolutions possible for the given geometry and time of flight.

)doc" );

    //////////////////////////////////////////////////////////////////////
    //  keplerPropagator.h (complete)
    /////////////////////////////////////////////////////////////////

    m.def( "propagate_kepler_orbit",
           &toec::propagateKeplerOrbit< double >,
           py::arg( "initial_kepler_elements" ),
           py::arg( "propagation_time" ),
           py::arg( "gravitational_parameter" ),
           py::arg( "root_finder" ) = trf::RootFinderPointer( ),
           R"doc(

 Function to propagate Keplerian elements to a later epoch, assuming an unperturbed system.

 Function to propagate Keplerian elements to a later epoch, assuming an unperturbed system. This function will
 take the initial Keplerian elements, and propagate the true anomaly in time as per the requested input. This
 is done by converting true anomaly to mean anomaly, apply the constant rate in mean motion for the requested
 time, and converting the result back to true anomaly. Currently both elliptic and hyperbolic orbits are supported.
 Parabolic orbits are not supported and will result in an error message.


 Parameters
 ----------
 initial_kepler_elements : numpy.ndarray
     Keplerian elements that are to be propagated (see :ref:`element_conversion` for order)
 propagation_time : astro.time_representation.Time
     Time object for which the elements are to be propagated w.r.t. the initial elements
 gravitational_parameter : float
     Gravitational parameter of central body used for propagation
 root_finder : RootFinder, default = None
     Root finder used to solve Kepler's equation when converting mean to eccentric anomaly. When no root finder is specified, the default option of the mean to eccentric anomaly function is used (see :func:`~tudatpy.astro.element_conversion.mean_to_eccentric_anomaly`).
 Returns
 -------
 numpy.ndarray
     Keplerian elements, propagated in time from initial elements assuming unperturbed dynamics. Note that the true anomaly is returned within the -PI to PI spectrum. If the user desires a different spectrum (possibly including the number of revolutions), these should be added by the user a posteriori.






     )doc" );
}
}  // namespace two_body_dynamics
}  // namespace astro
}  // namespace tudatpy
