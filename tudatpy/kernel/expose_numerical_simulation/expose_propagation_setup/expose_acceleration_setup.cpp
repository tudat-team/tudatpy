/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_acceleration_setup.h"
// #include "kernel/expose_numerical_simulation/deprecation_support.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/basics/deprecationWarnings.h>
#include <tudat/simulation/propagation_setup.h>

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;
namespace tmrf = tudat::root_finders;

namespace tudat
{
namespace simulation_setup
{

inline std::shared_ptr< AccelerationSettings > customAccelerationSettingsDeprecated(
        const std::function< Eigen::Vector3d( const double ) > accelerationFunction )
{
    static bool isWarningPrinted = false;
    if( isWarningPrinted == false )
    {
        tudat::utilities::printDeprecationWarning(
                "tudatpy.numerical_simulation.propagation_setup."
                "acceleration.custom",
                "tudatpy.numerical_simulation.propagation_setup."
                "acceleration.custom_acceleration" );
        isWarningPrinted = true;
    }

    return customAccelerationSettings( accelerationFunction );
}

inline std::shared_ptr< AccelerationSettings > thrustAccelerationRemoved1(
        const std::shared_ptr< tss::ThrustDirectionSettings > thrustDirectionSettings,
        const std::shared_ptr< tss::ThrustMagnitudeSettings > thrustMagnitudeSettings )
{
    tudat::utilities::printDeprecationError(
            "tudatpy.numerical_simulation.propagation_setup.acceleration."
            "thrust_from_direction_and_magnitude",
            "https://docs.tudat.space/en/stable/_src_user_guide/"
            "state_propagation/environment_setup/thrust_refactor/"
            "thrust_refactor.html#thrust-acceleration" );
    return nullptr;
}

inline std::shared_ptr< AccelerationSettings > thrustAccelerationRemoved2(
        const std::function< Eigen::Vector3d( const double ) > thrustForceFunction,
        const std::function< double( const double ) > specificImpulseFunction,
        const ThrustFrames thrustFrame = unspecified_thrust_frame,
        const std::string centralBody = "" )
{
    tudat::utilities::printDeprecationError(
            "tudatpy.numerical_simulation.propagation_setup.acceleration."
            "thrust_and_isp_from_custom_function",
            "https://docs.tudat.space/en/stable/_src_user_guide/"
            "state_propagation/environment_setup/thrust_refactor/"
            "thrust_refactor.html#thrust-acceleration" );
    return nullptr;
}

inline std::shared_ptr< AccelerationSettings > thrustAccelerationRemoved3(
        const std::function< Eigen::Vector3d( const double ) > thrustForceFunction,
        const double constantSpecificImpulse,
        const ThrustFrames thrustFrame = unspecified_thrust_frame,
        const std::string centralBody = "" )
{
    tudat::utilities::printDeprecationError(
            "tudatpy.numerical_simulation.propagation_setup.acceleration."
            "thrust_from_custom_function",
            "https://docs.tudat.space/en/stable/_src_user_guide/"
            "state_propagation/environment_setup/thrust_refactor/"
            "thrust_refactor.html#thrust-acceleration" );
    return nullptr;
}

//! @get_docstring(customAccelerationSettings)
inline std::shared_ptr< AccelerationSettings > customAccelerationSettings(
        const std::function< Eigen::Vector3d( const double ) > accelerationFunction )
{
    return std::make_shared< CustomAccelerationSettings >( accelerationFunction );
}

}  // namespace simulation_setup

}  // namespace tudat

namespace tudatpy
{
namespace numerical_simulation
{
namespace propagation_setup
{
namespace acceleration
{

void expose_acceleration_setup( py::module &m )
{
    /*
     * This contains the addition of IntegratorSettings and
     * AvailableIntegrators and AvailableAccelerations which
     * should be relocated in the tudat source.
     */

    py::enum_< tba::AvailableAcceleration >( m, "AvailableAcceleration", R"doc(

        Enumeration of available acceleration types.

        Enumeration of acceleration types supported by tudat.





     )doc" )
            .value( "undefined_acceleration_type",
                    tba::AvailableAcceleration::undefined_acceleration,
                    R"doc(
     )doc" )
            .value( "point_mass_gravity_type",
                    tba::AvailableAcceleration::point_mass_gravity,
                    R"doc(
     )doc" )
            .value( "aerodynamic_type", tba::AvailableAcceleration::aerodynamic, R"doc(
     )doc" )
            .value( "cannonball_radiation_pressure_type",
                    tba::AvailableAcceleration::cannon_ball_radiation_pressure,
                    R"doc(
     )doc" )
            .value( "spherical_harmonic_gravity_type",
                    tba::AvailableAcceleration::spherical_harmonic_gravity,
                    R"doc(
     )doc" )
            .value( "mutual_spherical_harmonic_gravity_type",
                    tba::AvailableAcceleration::mutual_spherical_harmonic_gravity,
                    R"doc(
     )doc" )
            .value( "polyhedron_gravity_type",
                    tba::AvailableAcceleration::polyhedron_gravity,
                    R"doc(
     )doc" )
            .value( "ring_gravity_type", tba::AvailableAcceleration::ring_gravity, R"doc(
     )doc" )
            .value( "thrust_acceleration_type",
                    tba::AvailableAcceleration::thrust_acceleration,
                    R"doc(
     )doc" )
            .value( "relativistic_correction_acceleration_type",
                    tba::AvailableAcceleration::relativistic_correction_acceleration,
                    R"doc(
     )doc" )
            .value( "empirical_acceleration_type",
                    tba::AvailableAcceleration::empirical_acceleration,
                    R"doc(
     )doc" )
            .value( "direct_tidal_dissipation_in_central_body_"
                    "acceleration_type",
                    tba::AvailableAcceleration::direct_tidal_dissipation_in_central_body_acceleration,
                    R"doc(
     )doc" )
            .value( "direct_tidal_dissipation_in_orbiting_body_"
                    "acceleration_type",
                    tba::AvailableAcceleration::direct_tidal_dissipation_in_orbiting_body_acceleration,
                    R"doc(
     )doc" )
            .value( "quasi_impulsive_shots_acceleration_type",
                    tba::AvailableAcceleration::momentum_wheel_desaturation_acceleration,
                    R"doc(
     )doc" )
            .value( "custom_acceleration_type",
                    tba::AvailableAcceleration::custom_acceleration,
                    R"doc(
     )doc" )
            .value( "radiation_pressure_type",
                    tba::AvailableAcceleration::radiation_pressure,
                    R"doc(
     )doc" )
            .export_values( );

    //////////////////////////////////////////////////////////////////////////////
    // accelerationSettings.h
    //////////////////////////////////////////////////////////////////////////////

    py::class_< tss::AccelerationSettings, std::shared_ptr< tss::AccelerationSettings > >( m,
                                                                                           "AccelerationSettings",
                                                                                           R"doc(

        Functional base class to define settings for accelerations.

        Class for providing settings for acceleration model. This class is a functional (base) class for
        settings of acceleration models that  require no information in addition to their type.
        Classes defining settings for acceleration models requiring additional information must be derived from this class.
        Bodies exerting and undergoing acceleration are set externally from this class.
        This class can be used for the easy setup of acceleration models
        (see createAccelerationModels.h), but users may also chose to do so manually.
        (Derived) Class members are all public, for ease of access and modification.





     )doc" );
    //            .def(py::init<const
    //            tudat::basic_astrodynamics::AvailableAcceleration>(),
    //                 py::arg("acceleration_type"));

    py::class_< tss::SphericalHarmonicAccelerationSettings,
                std::shared_ptr< tss::SphericalHarmonicAccelerationSettings >,
                tss::AccelerationSettings >( m,
                                             "SphericalHarmonicAccelerationSettings",
                                             R"doc(

        `AccelerationSettings`-derived class to define settings for the spherical harmonic acceleration.

        Class for providing settings for spherical harmonics acceleration model,
        including the maximum degree and order up to which the field is to be expanded. Note that
        the minimum degree and order are currently always set to zero.





     )doc" );
    //            .def(py::init<const int, const int>(),
    //            py::arg("maximum_degree"),
    //                 py::arg("maximum_order"));

    py::class_< tss::MutualSphericalHarmonicAccelerationSettings,
                std::shared_ptr< tss::MutualSphericalHarmonicAccelerationSettings >,
                tss::AccelerationSettings >( m,
                                             "MutualSphericalHarmonicAccelerationSettings",
                                             R"doc(

        `AccelerationSettings`-derived class to define settings for the mutual spherical harmonic acceleration.

        Class for providing settings for the mutual spherical harmonics acceleration model,
        including the maximum degree and order up to which the fields of the bodies are to be expanded. Note that
        the minimum degree and order are currently always set to zero.





     )doc" );

    py::class_< tss::EmpiricalAccelerationSettings, std::shared_ptr< tss::EmpiricalAccelerationSettings >, tss::AccelerationSettings >(
            m,
            "EmpiricalAccelerationSettings",
            R"doc(

        `AccelerationSettings`-derived class to define settings for the empirical acceleration.

        Class to provide settings for empirical accelerations. These are expressed in the
        RSW frame, for which the magnitude is determined empirically (typically during an orbit determination process).
        The acceleration components are defined according to Montenbruck and Gill (2000), with a total of 9 components:
        a constant, sine and cosine term (with true anomaly as argument) for each of the three independent directions of
        the RSW frame.





     )doc" );

    py::class_< tss::RelativisticAccelerationCorrectionSettings,
                std::shared_ptr< tss::RelativisticAccelerationCorrectionSettings >,
                tss::AccelerationSettings >( m,
                                             "RelativisticAccelerationCorrectionSettings",
                                             R"doc(

        `AccelerationSettings`-derived class to define settings for the relativistic acceleration correction.

        Class to provide settings for typical relativistic corrections to the dynamics of an orbiter: the
        Schwarzschild, Lense-Thirring and de Sitter terms (see 'General relativity and Space Geodesy' by L. Combrinck,
        2012).





     )doc" );

    py::class_< tss::CustomAccelerationSettings, std::shared_ptr< tss::CustomAccelerationSettings >, tss::AccelerationSettings >(
            m,
            "CustomAccelerationSettings",
            R"doc(

        `AccelerationSettings`-derived class to define settings for custom acceleration.

        Class to provide settings for custom accelerations. This is done by means of a function and, if necessary,
        an associated scaling function.





     )doc" );

    py::class_< tss::DirectTidalDissipationAccelerationSettings,
                std::shared_ptr< tss::DirectTidalDissipationAccelerationSettings >,
                tss::AccelerationSettings >( m,
                                             "DirectTidalDissipationAccelerationSettings",
                                             R"doc(

        `AccelerationSettings`-derived class to define settings for direct tidal dissipation acceleration.

        Class to provide settings for direct tidal dissipation accelerations. Creates settings for tidal accelerations.
        The direct of tidal effects in a satellite system is applied directly as an acceleration
        (as opposed to a modification of spherical harmonic coefficients).





     )doc" );

    py::class_< tss::MomentumWheelDesaturationAccelerationSettings,
                std::shared_ptr< tss::MomentumWheelDesaturationAccelerationSettings >,
                tss::AccelerationSettings >( m,
                                             "MomentumWheelDesaturationAccelerationSettings",
                                             R"doc(

        `AccelerationSettings`-derived class to define settings for momentum wheel desaturation acceleration.

        Class to provide settings for momentum wheel desaturation acceleration. Settings for the direction and magnitude
        of the thrust are included.





     )doc" );

    py::class_< tss::ThrustAccelerationSettings, std::shared_ptr< tss::ThrustAccelerationSettings >, tss::AccelerationSettings >(
            m,
            "ThrustAccelerationSettings",
            R"doc(

        `AccelerationSettings`-derived class to define settings for thrust acceleration, listing the engine models that are to be used

        Class to provide settings for thrust acceleration, listing the engine models that are to be used





     )doc" )
            .def_property_readonly(
                    "direction_settings",
                    &tss::ThrustAccelerationSettings::printDeprecationError< std::shared_ptr< tss::ThrustDirectionSettings > > )
            .def_property_readonly(
                    "magnitude_settings",
                    &tss::ThrustAccelerationSettings::printDeprecationError< std::shared_ptr< tss::ThrustMagnitudeSettings > > );

    // Unified interface functions for acceleration settings
    //  m.def("acceleration", &tss::acceleration,
    //  py::arg("acceleration_type"));
    m.def( "point_mass_gravity",
           &tss::pointMassGravityAcceleration,
           R"doc(

Creates settings for the point-mass gravity acceleration.

Creates settings for the point-mass gravity acceleration. The direct acceleration (acceleration w.r.t. an inertial frame) is computed from:

.. math::
   \mathbf{a}=\frac{\mu}{{r}^{2}}\hat{\mathbf{r}}

with :math:`\mathbf{r}` the position vector measured from the center of mass of the body exerting the acceleration.

The body exerting the acceleration needs to have a gravity field model (:ref:`\`\`gravity_field\`\`` module) defined to use this acceleration.

Depending on the body undergoing the acceleration :math:`A`, the body exerting the acceleration :math:`B`, and the central body of propagation :math:`C`, choosing this option may create a direct point-mass attraction (:math:`\mu=\mu_{B}`), a central point-mass attraction (:math:`\mu=\mu_{B}+\mu_{A}`) or a third-body point-mass attraction (see `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/third_body_acceleration.html>`_ for more details).

Returns
-------
AccelerationSettings
    Acceleration settings object.





Examples
--------
In this example, we define the point mass gravity acceleration exerted by the Earth on the vehicle:

.. code-block:: python

   # Create acceleration dict
   accelerations_acting_on_vehicle = dict()
   # Add aerodynamic acceleration exerted by Earth
   accelerations_acting_on_vehicle["Earth"] = [propagation_setup.acceleration.point_mass_gravity()]


    )doc" );

    m.def( "aerodynamic",
           &tss::aerodynamicAcceleration,
           R"doc(

Creates settings for the aerodynamic acceleration.

Creates settings for the aerodynamic acceleration. The acceleration is computed from:

.. math::
   \mathbf{a}=-\frac{1}{m}\mathbf{R}^{(I/\text{Aero})}\left(\frac{1}{2}\rho v_{\text{air}}^{2}S_{ref}\begin{pmatrix} C_{D} \\ C_{S} \\ C_{L}\end{pmatrix}\right)

with :math:`\mathbf{R}^{(I/\text{Aero})}` the rotation matrix from the aerodynamic frame of the body undergoing acceleration to the inertial frame (computed from the body's current state, and the rotation of the body exerting the acceleration), :math:`\rho` the local freestream atmospheric density, :math:`v_{\text{air}}` the airspeed,  :math:`C_{D,S,L}` the drag, side and lift coefficients (which may depend on any number of properties of the body/environment) with reference area :math:`S_{ref}`, and :math:`m` the mass of the body undergoing acceleration
The body exerting the acceleration needs to have an
atmosphere (:ref:`\`\`atmosphere\`\`` module), shape (:ref:`\`\`shape\`\`` module) and rotation model (:ref:`\`\`rotation_model\`\`` module) defined. The body undergoing the acceleration needs to have aerodynamic coefficients (:ref:`\`\`aerodynamic_coefficients\`\`` module) defined.

Returns
-------
AccelerationSettings
    Acceleration settings object.





Examples
--------
In this example, we define the aerodynamic acceleration exerted by the Earth on the vehicle:

.. code-block:: python

   # Create acceleration dict
   accelerations_acting_on_vehicle = dict()
   # Add aerodynamic acceleration exerted by Earth
   accelerations_acting_on_vehicle["Earth"] = [propagation_setup.acceleration.aerodynamic()]


    )doc" );

    m.def( "cannonball_radiation_pressure",
           &tss::cannonBallRadiationPressureAcceleration,
           R"doc(

Creates settings for the cannonball radiation pressure acceleration.

Creates settings for the radiation pressure acceleration, for which a cannonball model is used. The acceleration is computed from:

.. math::

   \mathbf{a}=\left(\frac{P}{4\pi c}\right)\left(\frac{C_{r}S_{ref}}{m}\right)\frac{\hat{\mathbf{r}}}{r^{2}}

with :math:`P` the total emitted radiation power for the body exerting the acceleration, :math:`C_{r}` the radiation pressure coefficient with reference area :math:`S_{ref}`, :math:`\mathbf{r}` the vector from the body exerting the acceleration to the body undergoing the acceleration, and :math:`m` the mass of the body undergoing acceleration

In this model,
the effective acceleration is colinear with the vector connecting the source of radiation and the target.
The body undergoing the acceleration needs to have a radiation pressure model defined, while the body emitting
radiation needs to have radiative properties defined.

Returns
-------
AccelerationSettings
    Acceleration settings object.





Examples
--------
In this example, we define the aerodynamic acceleration exerted by the Sun on the vehicle:

.. code-block:: python

   # Create acceleration dict
   accelerations_acting_on_vehicle = dict()
   # Add cannonball radiation pressure acceleration exerted by Sun
   accelerations_acting_on_vehicle["Sun"] = [propagation_setup.acceleration.cannonball_radiation_pressure()]


    )doc" );

    m.def( "radiation_pressure",
           &tss::radiationPressureAcceleration,
           py::arg( "target_type" ) = tss::undefined_target,
           R"doc(

Creates settings for the radiation pressure acceleration.

Creates settings for the radiation pressure acceleration. This function takes the source model of the body exerting the acceleration, and the target model of the
body undergoing the acceleration, and links these models to set up the specific acceleration model, regardless of how the source and target models are defined.
For more extensive details on how this is done, check out the
`radiation pressure acceleration page <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/radiation_pressure_acceleration.html>`_.

Returns
-------
AccelerationSettings
    Acceleration settings object.






    )doc" );

    m.def( "spherical_harmonic_gravity",
           &tss::sphericalHarmonicAcceleration,
           py::arg( "maximum_degree" ),
           py::arg( "maximum_order" ),
           R"doc(

Creates settings for the spherical harmonic gravity acceleration.

Creates settings for the spherical harmonic gravity acceleration, accounting for a finite (given as input) number
of degrees and orders. The direct acceleration (acceleration w.r.t. an inertial origin) is computed from:

.. math::
   \mathbf{a}=\mathbf{R}^{(I/B)}\nabla^{(B)}U(\mathbf{r})

with :math:`\mathbf{r}` the position vector measured from the center of mass of the body exerting the acceleration, :math:`\mathbf{R}^{(I/B)}` the rotation matrix from body-fixed to inertial frame, and :math:`\nabla^{(B)}` the gradient operator in a body-fixed frame, and :math:`U` the spherical harmonic gravitational potential, expanded up to the provided ``maximum_degree`` and ``maximum_order``.

The body exerting the acceleration needs to have a spherical harmonic gravity field model (see :class:`~tudatpy.numerical_simulation.environment_setup.gravity_field.spherical_harmonic`) and a rotation model (:ref:`\`\`rotation_model\`\`` module) defined.

Depending on the body undergoing the acceleration :math:`A`, the body exerting the acceleration :math:`B`, and the central body of propagation :math:`C`, choosing this option may create a direct spherical harmonic attraction (:math:`\mu=\mu_{B}`), a central spherical harmonic attraction (:math:`\mu=\mu_{B}+\mu_{A}`) or a third-body spherical harmonic attraction (see `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/third_body_acceleration.html>`_ for more details).


Parameters
----------
maximum_degree : int
    Maximum degree of the spherical harmonic expansion.
maximum_order : int
    Maximum order of the spherical harmonic expansion.
Returns
-------
SphericalHarmonicAccelerationSettings
    Spherical harmonic acceleration settings object.





Examples
--------
In this example, we define the spherical harmonic gravity acceleration (where the gravity field is expanded
up to degree 12 and order 6) exerted by the Earth on the vehicle:

.. code-block:: python

   # Define the maximum degree and order
   maximum_degree = 12
   maximum_order = 6
   # Create acceleration dict
   accelerations_acting_on_vehicle = dict()
   # Add aerodynamic acceleration exerted by Earth
   accelerations_acting_on_vehicle["Earth"] = [propagation_setup.acceleration.spherical_harmonic_gravity(
        maximum_degree,
        maximum_order)]


    )doc" );

    m.def( "mutual_spherical_harmonic_gravity",
           &tss::mutualSphericalHarmonicAcceleration,
           py::arg( "maximum_degree_body_exerting" ),
           py::arg( "maximum_order_body_exerting" ),
           py::arg( "maximum_degree_body_undergoing" ),
           py::arg( "maximum_order_body_undergoing" ),
           py::arg( "maximum_degree_central_body" ) = 0,
           py::arg( "maximum_order_central_body" ) = 0,
           R"doc(

Creates settings for the mutual spherical harmonic gravity acceleration.

Creates settings for the mutual spherical harmonic gravity acceleration. This model computes the total spherical harmonic acceleration exerted by a body :math:`B` on a body :math:`A`, where the influence of the gravity field coefficients of body :math:`A` itself has been included. The model includes couplings between the mass of each body, and the gravity field coefficients of the other body. It does not include the 'figure-figure' interactions (coupling between the two-bodies' gravity field coefficients). It corresponds to the model presented by Lainey et al. (2004); Dirkx et al. (2016).
The model combines the spherical harmonic accelerations of the two bodies (see :func:`~tudatpy.numerical_simulation.propagation_setup.acceleration.spherical_harmonic_gravity`) on each other. The direct acceleration (acceleration w.r.t. an inertial origin) is computed from:

.. math::

   \mathbf{a}={-\frac{\mu_{_{B}}}{{r}^{2}}\hat{\mathbf{r}}}+{\mathbf{R}^{(I/B)}\nabla^{(B)}U_{\hat{B}}(\mathbf{r})}-{\frac{\mu_{_{B}}}{\mu_{_{A}}}\mathbf{R}^{(I/A)}\nabla^{(A)}U_{\hat{A}}(-\mathbf{r})}

where :math:`U_{\hat{B}}` and :math:`U_{\hat{A}}` denote the spherical harmonic gravity fields a degree :math:`>=1` of bodies :math:`B` and :math:`A`, respectively.
Both the body exerting the acceleration and the body undergoing it need to
have spherical harmonic gravity field and rotation models defined.

Depending on the body undergoing the acceleration :math:`A`, the body exerting the acceleration :math:`B`, and the central body of propagation :math:`C`, choosing this option may create a direct spherical harmonic attraction (as above), a central spherical harmonic attraction (:math:`\mu_{B}\rightarrow\mu_{B}+\mu_{A}`, in the above equation and in :math:`U_{\hat{B}}`) or a third-body spherical harmonic attraction (see `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/third_body_acceleration.html>`_ for more details).

For the case where a third-body mutual spherical harmonic acceleration,
additional parameters have to be provided that denote the expansion degree/order of the central body (``maximum_degree_central_body`` and ``maximum_order_central_body``)


Parameters
----------
maximum_degree_body_exerting : int
    Maximum degree of the spherical harmonic expansion for the body exerting the acceleration.
maximum_order_body_exerting : int
    Maximum order of the spherical harmonic expansion for the body exerting the acceleration.
maximum_degree_body_undergoing : int
    Maximum degree of the spherical harmonic expansion for the body undergoing the acceleration.
maximum_order_body_undergoing : int
    Maximum order of the spherical harmonic expansion for the body undergoing the acceleration.
maximum_degree_central_body : int, default=0
    Maximum degree of the spherical harmonic expansion for the central body, if needed.
maximum_order_central_body : int, default=0
    Maximum order of the spherical harmonic expansion for the central body, if needed.
Returns
-------
MutualSphericalHarmonicAccelerationSettings
    Spherical harmonic acceleration settings object.





Examples
--------
In this example, we define the spherical harmonic gravity accelerations exerted on Io by Jupiter and vice-versa.

.. code-block:: python

   # Define the maximum degree and order for both bodies
   maximum_degree_of_io = 12
   maximum_order_of_io = 12
   maximum_degree_of_jupiter = 4
   maximum_order_of_jupiter = 4
   # Create acceleration dict
   acceleration_settings_on_io = dict()
   # Add aerodynamic acceleration exerted by Earth
   acceleration_settings_on_io["Jupiter"] = [propagation_setup.acceleration.mutual_spherical_harmonic_gravity(
        maximum_degree_of_jupiter,
        maximum_order_of_jupiter,
        maximum_degree_of_io,
        maximum_order_of_io)]

For the case where the mutual spherical harmonic acceleration is a third-body acceleration,
additional parameters have to be provided to denote the expansion of the spherical harmonics of the central body.
In the following example, we consider the spherical harmonic gravity acceleration mutually exerted between
Ganymede and Io when propagating w.r.t. Jupiter:

.. code-block:: python

   # Define the maximum degree and order for both bodies
   maximum_degree_of_jupiter = 4
   maximum_order_of_jupiter = 4
   maximum_degree_of_ganymede = 4
   maximum_order_of_ganymede = 4
   maximum_degree_of_io = 12
   maximum_order_of_io = 12
   # Create acceleration dict
   acceleration_settings_on_io = dict()
   # Add the acceleration to the dict
   acceleration_settings_on_io["Jupiter"] = [propagation_setup.acceleration.mutual_spherical_harmonic_gravity(
        maximum_degree_of_jupiter,
        maximum_order_of_jupiter,
        maximum_degree_of_ganymede,
        maximum_order_of_ganymede,
        maximum_degree_of_io,
        maximum_order_of_io)]


    )doc" );

    m.def( "polyhedron_gravity",
           &tss::polyhedronAcceleration,
           R"doc(

Creates settings for the polyhedron gravity acceleration.

Creates settings for the polyhedron gravity acceleration, which follows from defining a body to have polyhedral gravity.

Returns
-------
AccelerationSettings
    Acceleration settings object.






    )doc" );

    m.def( "ring_gravity",
           &tss::ringAcceleration,
           R"doc(

Creates settings for the ring gravity acceleration.

Creates settings for the ring gravity acceleration, which follows from defining a body to have ring gravity.

Returns
-------
AccelerationSettings
    Acceleration settings object.






    )doc" );

    m.def( "relativistic_correction",
           &tss::relativisticAccelerationCorrection,
           py::arg( "use_schwarzschild" ) = false,
           py::arg( "use_lense_thirring" ) = false,
           py::arg( "use_de_sitter" ) = false,
           py::arg( "de_sitter_central_body" ) = "",
           py::arg( "lense_thirring_angular_momentum" ) = Eigen::Vector3d::Zero( ),
           R"doc(

Creates settings for the relativistic acceleration correction.

Creates settings for typical relativistic acceleration corrections: the Schwarzschild, Lense-Thirring and de
Sitter terms, where each of the three terms can be toggled on or of (see 'General relativity and Space Geodesy' by L. Combrinck, 2012). It implements the model of
2010 Conventions (chapter 10, section 3). Here, the ‘primary body’ for a planetary orbiter should always be set
as the Sun (only relevant for de Sitter correction). The angular momentum vector of the orbited body is only
relevant for Lense-Thirring correction.


Parameters
----------
use_schwarzschild : bool, default=False
    Boolean defining whether or not to use the Schwarzschild contribution to the acceleration correction
use_lense_thirring : bool
    Boolean defining whether or not to use the Lense-Thirring contribution to the acceleration correction
use_de_sitter : bool, default=False
    Boolean defining whether or not to use the de Sitter contribution to the acceleration correction
de_sitter_central_body : str, default=""
    Body used as 'perturbed' in the calculation of the de Sitter acceleration. For the case of an Earth-orbiting satellite, this would be the Sun
lense_thirring_angular_momentum : numpy.ndarray, default=numpy.array([0, 0, 0])
    Angular momentum vector per unit mass (in global frame) that is to be used for the calculation of the Lense-Thirring acceleration
Returns
-------
RelativisticAccelerationCorrectionSettings
    Relativistic acceleration correction settings object.





Examples
--------
In this example, we define the relativistic correction acceleration for a Mars orbiter:

.. code-block:: python

   # Select terms to be used
   use_schwarzschild = True
   use_lense_thirring = True
   use_de_sitter = True
   # Define central body for De-Sitter term
   de_sitter_central_body = "Sun"
   # Define angular momentum for the Lense-Thirring term
   lense_thirring_angular_momentum = ... # numpy.ndarray 3D vector
   # Create acceleration dict
   acceleration_settings_on_vehicle = dict()
   # Add the acceleration to the dict
   acceleration_settings_on_vehicle["Mars"] = [propagation_setup.acceleration.relativistic_correction(
      use_schwarzschild,
      use_lense_thirring,
      use_de_sitter,
      de_sitter_central_body,
      lense_thirring_angular_momentum)]


    )doc" );

    m.def( "empirical",
           &tss::empiricalAcceleration,
           py::arg( "constant_acceleration" ) = Eigen::Vector3d::Zero( ),
           py::arg( "sine_acceleration" ) = Eigen::Vector3d::Zero( ),
           py::arg( "cosine_acceleration" ) = Eigen::Vector3d::Zero( ),
           R"doc(

Creates settings for empirical acceleration.

Creates settings for empirical accelerations. These are expressed in the
RSW frame, for which the magnitude is determined empirically (typically during an orbit determination process).
The acceleration components are defined according to Montenbruck and Gill (2000), with a total of 9 components:
a constant, sine and cosine term (with true anomaly as argument) for each of the three independent directions of
the RSW frame. The empirical acceleration is calculated as:

 .. math::

    \mathbf{a}=R^{I/RSW}\left(\mathbf{a}_{\text{const.}}+\mathbf{a}_{\sin}\sin\theta+\mathbf{a}_{\cos}\cos\theta \right)

Here, :math:`R^{I/RSW}` is the rotation matrix from the RSW frame (of the body undergoing the acceleration w.r.t. the
body exerting the acceleration), :math:`theta` is the true anomaly, and the three constituent acceleration vectors are
the inputs provided in the above code block. The body 'exerting' the acceleration is considered to be the
central body, w.r.t. which the true anomaly is calculated.


Parameters
----------
constant_acceleration : numpy.ndarray, default=numpy.array([0, 0, 0])
    Constant term, defined in the RSW frame.
sine_acceleration : numpy.ndarray, default=numpy.array([0, 0, 0])
    Sine term (function of the true anomaly), defined in the RSW frame..
cosine_acceleration : numpy.ndarray, default=numpy.array([0, 0, 0])
    Cosine term (function of the true anomaly), defined in the RSW frame..
Returns
-------
EmpiricalAccelerationSettings
    Empirical acceleration settings object.





Examples
--------
In this example, we define the empirical acceleration acting on the vehicle. The body that 'exerts' the acceleration
is here used to determine the body w.r.t. which the true anomaly has o be calculated when determining the sine/cosine
contributions. This central body must be endowed with a gravity field (so that it possesses a gravitational parameter
for the Cartesian to Keplerian conversion)

.. code-block:: python

   # Define the nine terms (constant, sine and cosine)
   constant_acceleration = ... # 3D numpy.ndarray vector
   sine_acceleration = ... # 3D numpy.ndarray vector
   cosine_acceleration = ... # 3D numpy.ndarray vector
   # Create acceleration dict
   acceleration_settings_on_vehicle = dict()
   # Add the acceleration to the dict
   acceleration_settings_on_vehicle["Sun"] = [propagation_setup.acceleration.empirical(
       constant_acceleration,
       sine_acceleration,
       cosine_acceleration)]


    )doc" );

    m.def( "yarkovsky",
           &tss::yarkovskyAcceleration,
           py::arg( "yarkovsky_parameter" ),
           R"doc(

Creates settings for the Yarkovsky acceleration.

Creates settings for the Yarkovsky acceleration, which is calculated based on Pérez-Hernández & Benet (2022). The acceleration is only
considered in the tangential direction and is proportional to:

.. math::

   \mathbf{a}=A_{2} \cdot (\frac{r_{0}}{r_{S}})^{2}

where :math:`A_{2}` is the Yarkovsky parameter, :math:`r_{0} = 1` AU and :math:`r_{S}` is the heliocentric distance in AU.


Parameters
----------
yarkovsky_parameter : float
    Value of the Yarkovsky parameter.
Returns
-------
AccelerationSettings
    Acceleration settings object.






    )doc" );

    m.def( "custom", &tss::customAccelerationSettingsDeprecated, py::arg( "acceleration_function" ) );

    m.def( "custom_acceleration",
           py::overload_cast< std::function< Eigen::Vector3d( const double ) > >( &tss::customAccelerationSettings ),
           py::arg( "acceleration_function" ),
           R"doc(

Creates settings for custom acceleration.

Creates settings for a custom accelerations, this acceleration must be parameterized as a function of time,
and expressed with an inertial orientation.


Parameters
----------
acceleration_function : callable[[float], list]
    Custom acceleration function with time as an independent variable, returning the acceleration in an inertial frame (*e.g.* with global frame orientation) as a function of time.
Returns
-------
CustomAccelerationSettings
    Custom acceleration settings object.





Examples
--------
In this example, we define a simple, standalone, custom acceleration (depending only on time),
with the following (arbitrary, for the sake of example) mathematical definition:

 .. math::

    \mathbf{a}=\begin{pmatrix}C\\0\\0 \end{pmatrix}\sin\left(\frac{t-t_{0}}{T}\right)

with :math:`C=10^{-8}`, :math:`t_{0}=0` and :math:`T=86400`.
More complex custom acceleration functions can be created by, for instance, extracting the
custom function from a user-defined class, which may in turn have other simulation objects
(*e.g.* :class:`SystemOfBodies`) as members, allowing the custom acceleration to depend on
the current simulation state/environment

.. code-block:: python

   # Define custom function
   def custom_function( current_time ):
       period = 86400.0
       reference_time = 0.0
       phase = np.pi * ( current_time - reference_time ) / period
       return np.array([1.0E-8, 0.0, 0.0]) * np.sin(phase)

   acceleration_settings_on_vehicle["Vehicle"] = [propagation_setup.acceleration.custom_acceleration(
       custom_function)]


    )doc" );

    m.def( "direct_tidal_dissipation_acceleration",
           &tss::directTidalDissipationAcceleration,
           py::arg( "k2_love_number" ),
           py::arg( "time_lag" ),
           py::arg( "include_direct_radial_component" ) = true,
           py::arg( "use_tide_raised_on_planet" ) = true,
           py::arg( "explicit_libraional_tide_on_satellite" ) = false,
           R"doc(

Creates settings for tidal acceleration.

Creates settings for tidal accelerations. The direct of tidal effects in a satellite system is applied directly as
an acceleration (as opposed to a modification of spherical harmonic coefficients).
The model is based on Lainey et al. (2007, 2012). It can compute the acceleration due to tides, and in
particular tidal dissipation, on a planetary satellite. The acceleration computed can account for either the
effect of tide raised on the satellite by the planet or on the planet by the satellite. The satellite is assumed
to be tidally locked to the planet.


Parameters
----------
k2_love_number : float
    Value of the k2 Love number.
time_lag : float
    Value of the tidal time lag.
include_direct_radial_component : bool, default=True
    It denotes whether the term independent of the time lag is to be computed.
use_tide_raised_on_planet : bool, default=True
    It denotes whether the tide raised on the planet is to be modelled (if true) or the tide raised on the satellite (if false).
Returns
-------
DirectTidalDissipationAccelerationSettings
    Direct tidal dissipation acceleration settings object.





Examples
--------
In this example, we define the tidal dissipation exerted by Jupiter on Io directly, instead of computing it
through the spherical harmonic gravity:

.. code-block:: python

   # Define parameters
   love_number = 0.1
   time_lag = 100.0
   # Add entry to acceleration settings dict
   acceleration_settings_on_io["Jupiter"] = [propagation_setup.acceleration.direct_tidal_dissipation(
      love_number,
      time_lag,
      False,
      False)]


    )doc" );

    m.def( "direct_tidal_dissipation_acceleration",
           &tss::directTidalDissipationAccelerationFromInvQ,
           py::arg( "k2_love_number" ),
           py::arg( "inverse_tidal_quality_factor" ),
           py::arg( "tidal_period" ),
           py::arg( "include_direct_radial_component" ) = true,
           py::arg( "use_tide_raised_on_planet" ) = true,
           py::arg( "explicit_libraional_tide_on_satellite" ) = false,
           R"doc(

Creates settings for tidal acceleration.

Creates settings for tidal accelerations. The direct of tidal effects in a satellite system is applied directly as
an acceleration (as opposed to a modification of spherical harmonic coefficients).
The model is based on Lainey et al. (2007, 2012). It can compute the acceleration due to tides, and in
particular tidal dissipation, on a planetary satellite. The acceleration computed can account for either the
effect of tide raised on the satellite by the planet or on the planet by the satellite. The satellite is assumed
to be tidally locked to the planet.


Parameters
----------
k2_love_number : float
    Value of the k2 Love number.
time_lag : float
    Value of the tidal time lag.
include_direct_radial_component : bool, default=True
    It denotes whether the term independent of the time lag is to be computed.
use_tide_raised_on_planet : bool, default=True
    It denotes whether the tide raised on the planet is to be modelled (if true) or the tide raised on the satellite (if false).
Returns
-------
DirectTidalDissipationAccelerationSettings
    Direct tidal dissipation acceleration settings object.





Examples
--------
In this example, we define the tidal dissipation exerted by Jupiter on Io directly, instead of computing it
through the spherical harmonic gravity:

.. code-block:: python

   # Define parameters
   love_number = 0.1
   time_lag = 100.0
   # Add entry to acceleration settings dict
   acceleration_settings_on_io["Jupiter"] = [propagation_setup.acceleration.direct_tidal_dissipation(
      love_number,
      time_lag,
      False,
      False)]


    )doc" );

    m.def( "quasi_impulsive_shots_acceleration",
           &tss::momentumWheelDesaturationAcceleration,
           py::arg( "thrust_mid_times" ),
           py::arg( "delta_v_values" ),
           py::arg( "total_maneuver_time" ),
           py::arg( "maneuver_rise_time" ),
           R"doc(

Creates settings for incorporating quasi-impulsive shots into the acceleration.

The acceleration model is purpose-built to represent short bursts of thrust, such as a momentum wheel desaturation.
A typical use case is precise orbit determination, but the functionality can be used just as well in propagation
(for instance to model an impulsive manuever in a continuous manner when going from preliminary modelling to
full modelling). The thrust is modelled similarly to Fig. 3 of Alessi et al. (2012), with the main difference
being that a third-order polynomial to go from zero acceleration to the maximum acceleration level is employed.
By using a 3rd-order polynomial and imposing continuity in the value and first derivative of the acceleration,
defining the rise time (time it takes acceleration to go from 0 to its maximum level), the total time where
there is non-zero thrust (total maneuver time), and the total Delta V exerted by a single maneuver,
the acceleration profile is fully defined.


Parameters
----------
thrust_mid_times : list[float]
    Set of middle point in times in the maneuver denoting the epoch of each maneuver.
delta_v_values : list[numpy.ndarray]
    Set of delta V, one for each maneuver.
total_maneuver_time : float
    Total duration of every maneuver.
maneuver_rise_time : float
    Time taken by the acceleration to go from zero to its maximum level.
Returns
-------
MomentumWheelDesaturationAccelerationSettings
    Momentum wheel desaturation acceleration settings object.





Examples
--------
In this example, we define an acceleration model to represent two quasi-impulsive shots, with a total duration of
30 seconds, and a rise time of 5 seconds. The maneuvers are to be done at :math:`t=86400` and :math:`t=2*86400`.
The first maneuver is exclusively in :math:`x`-direction, and the second one exclusively in :math:`y`-direction,
with both maneuvers having a magnitude of 1 mm/s

.. code-block:: python

   # Define the quasi-impulsive shot settings
   mid_times = [ 86400.0, 2.0 * 86400.0]
   delta_v_values = [ np.array([1.0E-3, 0.0, 0.0]), np.array([0.0, 1.0E-3, 0.0]) ]
   maneuver_duration = 30
   maneuver_duration = 5

   # Create acceleration dict
   acceleration_settings_on_spacecraft = dict()

   # Add quasi-impulsive acceleration exerted by Spacecraft itself (!)
   acceleration_settings_on_spacecraft["Spacecraft"] = [propagation_setup.acceleration.quasi_impulsive_shots_acceleration(
        mid_times,
        delta_v_values,
        maneuver_duration,
        maneuver_duration)]


    )doc" );

    m.def( "thrust_from_engines",
           &tss::thrustAcceleration,
           py::arg( "engine_names" ),
           R"doc(

Creates settings for thrust acceleration using a list of engine models.

Creates settings for thrust acceleration using a list of engine models. See the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/thrust_models.html>`_
for more details on the definition of a thrust model in Tudat.


Parameters
----------
engine_names : List[str]
    List of engine names to use when computing thrust.
Returns
-------
ThrustAccelerationSettings
    Thrust acceleration settings object.






    )doc" );

    m.def( "thrust_from_engine",
           &tss::thrustAccelerationFromSingleEngine,
           py::arg( "engine_name" ),
           R"doc(

Creates settings for thrust acceleration using a single engine models.

Creates settings for thrust acceleration using a single engine models. See the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/thrust_models.html>`_
for more details on the definition of a thrust model in Tudat.


Parameters
----------
engine_name : str
    Name of engine to use when computing thrust.
Returns
-------
ThrustAccelerationSettings
    Thrust acceleration settings object.






    )doc" );

    m.def( "thrust_from_all_engines",
           &tss::thrustAccelerationFromAllEngines,
           R"doc(

Creates settings for thrust acceleration using a single engine models.

Creates settings for thrust acceleration by combining thrust from all engines defined in the body. See the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/translational/thrust_models.html>`_
for more details on the definition of a thrust model in Tudat.

Returns
-------
ThrustAccelerationSettings
    Thrust acceleration settings object.






    )doc" );

    m.def( "thrust_from_direction_and_magnitude",
           &tss::thrustAccelerationRemoved1,
           py::arg( "thrust_direction_settings" ),
           py::arg( "thrust_magnitude_settings" ) );

    m.def( "thrust_from_custom_function",
           &tss::thrustAccelerationRemoved2,
           py::arg( "thrust_force_function" ),
           py::arg( "specific_impulse_function" ),
           py::arg( "thrust_frame" ) = tss::ThrustFrames::inertial_thrust_frame,
           py::arg( "central_body" ) = "" );

    m.def( "thrust_and_isp_from_custom_function",
           &tss::thrustAccelerationRemoved3,
           py::arg( "thrust_force_function" ),
           py::arg( "constant_specific_impulse" ),
           py::arg( "thrust_frame" ) = tss::ThrustFrames::inertial_thrust_frame,
           py::arg( "central_body" ) = "",
           R"doc(No documentation found.)doc" );
}

}  // namespace acceleration
}  // namespace propagation_setup
}  // namespace numerical_simulation
}  // namespace tudatpy
