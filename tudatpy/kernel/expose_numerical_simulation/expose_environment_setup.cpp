/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_environment_setup.h"

#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/astro/reference_frames/referenceFrameTransformations.h>
#include <tudat/simulation/environment_setup.h>

#include "expose_environment_setup/expose_aerodynamic_coefficient_setup.h"
#include "expose_environment_setup/expose_atmosphere_setup.h"
#include "expose_environment_setup/expose_ephemeris_setup.h"
#include "expose_environment_setup/expose_gravity_field_setup.h"
#include "expose_environment_setup/expose_gravity_field_variation_setup.h"
#include "expose_environment_setup/expose_ground_station_setup.h"
#include "expose_environment_setup/expose_radiation_pressure_setup.h"
#include "expose_environment_setup/expose_rigid_body_setup.h"
#include "expose_environment_setup/expose_rotation_model_setup.h"
#include "expose_environment_setup/expose_shape_deformation_setup.h"
#include "expose_environment_setup/expose_shape_setup.h"
#include "expose_environment_setup/expose_vehicle_systems_setup.h"
#include "scalarTypes.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace te = tudat::ephemerides;
namespace ti = tudat::interpolators;
namespace tba = tudat::basic_astrodynamics;
namespace ta = tudat::aerodynamics;
namespace trf = tudat::reference_frames;
namespace tg = tudat::gravitation;
namespace tcc = tudat::coordinate_conversions;
namespace tp = tudat::propagators;

namespace tudatpy
{
namespace numerical_simulation
{
namespace environment_setup
{

void expose_environment_setup( py::module &m )
{
    //        m.def("get_body_gravitational_parameter",
    //              &tss::getBodyGravitationalParameter,
    //              py::arg("body_collection"),
    //              py::arg("body_name"));

    py::class_< tss::BodySettings, std::shared_ptr< tss::BodySettings > >( m, "BodySettings", R"doc(

        Class for defining settings for the creation of a single body.

        Class for defining settings for the creation of a single body, this object is typically stored inside a
        :class:`BodyListSettings` object.

     )doc" )
            .def_readwrite( "constant_mass", &tss::BodySettings::constantMass, R"doc(

        Mass that gets assigned to the vehicle. This mass does *not* automatically define a gravity field
        model, but is instead used for the calculation of non-conservative forces only. When creating a body with a gravity field,
        leave this entry empty. NOTE: this option is a shorthand for assigning a mass-only
        :func:`~tudatpy.numerical_simulation.environment_setup.rigid_body.constant_rigid_body_properties` to ``mass_property_settings``, and will be deprecated.

void expose_environment_setup( py::module &m )
{
    //        m.def("get_body_gravitational_parameter",
    //              &tss::getBodyGravitationalParameter,
    //              py::arg("body_collection"), py::arg("body_name"));

        :type: float
     )doc" )
            .def_readwrite( "atmosphere_settings",
                            &tss::BodySettings::atmosphereSettings,
                            R"doc(

        Object that defines the settings of the atmosphere model that is to be created. Note that wind model settings
        may be defined inside this object. A variable of this type is typically assigned by using a function from the
        :ref:`\`\`atmosphere\`\`` module.

    m.def( "get_default_body_settings_time_limited",
           py::overload_cast< const std::vector< std::string > &,
                              const double,
                              const double,
                              const std::string,
                              const std::string,
                              const double >( &tss::getDefaultBodySettings ),
           py::arg( "bodies" ),
           py::arg( "initial_time" ),
           py::arg( "final_time" ),
           py::arg( "base_frame_origin" ) = "SSB",
           py::arg( "base_frame_orientation" ) = "ECLIPJ2000",
           py::arg( "time_step" ) = 300.0,
           get_docstring( "get_default_body_settings_time_limited" ).c_str( ) );

        :type: AtmosphereSettings
     )doc" )
            .def_readwrite( "ephemeris_settings", &tss::BodySettings::ephemerisSettings, R"doc(

        Object that defines the settings of the ephemeris model that is to be created. A variable of this type is typically
        assigned by using a function from the :ref:`\`\`ephemeris\`\`` module.

    m.def( "get_default_single_body_settings_time_limited",
           py::overload_cast< const std::string &, const double, const double, const std::string &, const double >(
                   &tss::getDefaultSingleBodySettings ),
           py::arg( "body_name" ),
           py::arg( "initial_time" ),
           py::arg( "final_time" ),
           py::arg( "base_frame_orientation" ) = "ECLIPJ2000",
           py::arg( "time_step" ) = 300.0,
           get_docstring( "get_default_single_body_settings_time_limited" ).c_str( ) );

        :type: EphemerisSettings
     )doc" )
            .def_readwrite( "gravity_field_settings",
                            &tss::BodySettings::gravityFieldSettings,
                            R"doc(

        Object that defines the settings of the gravity field model that is to be created. A variable of this type is typically
        assigned by using a function from the :ref:`\`\`gravity_field\`\`` module.


        :type: GravityFieldSettings
     )doc" )
            .def_readwrite( "rotation_model_settings",
                            &tss::BodySettings::rotationModelSettings,
                            R"doc(

        Object that defines the settings of the rotation model that is to be created. A variable of this type is typically
        assigned by using a function from the :ref:`\`\`rotation_model\`\`` module.


        :type: RotationModelSettings
     )doc" )
            .def_readwrite( "shape_settings",
                            &tss::BodySettings::shapeModelSettings,
                            R"doc(

        Object that defines the settings of the shape model that is to be created. A variable of this type is typically
        assigned by using a function from the :ref:`\`\`shape\`\`` module.


        :type: BodyShapeSettings
     )doc" )
            .def_readwrite( "aerodynamic_coefficient_settings",
                            &tss::BodySettings::aerodynamicCoefficientSettings,
                            R"doc(

        Object that defines the settings of the aerodynamic coefficient model that is to be created. A variable of this type is typically
        assigned by using a function from the :ref:`\`\`aerodynamic_coefficients\`\`` module.


        :type: AerodynamicCoefficientSettings
     )doc" )
            .def_readwrite( "gravity_field_variation_settings",
                            &tss::BodySettings::gravityFieldVariationSettings,
                            R"doc(

        List of objects that define the settings of time variations of the gravity field variation models that are to be created. Variables in this list are typically
        assigned by using a function from the :ref:`\`\`gravity_field_variations\`\`` module.


        :type: list[GravityFieldVariationSettings]
     )doc" )
            .def_readwrite( "shape_deformation_settings",
                            &tss::BodySettings::bodyDeformationSettings,
                            R"doc(

        List of objects that define the settings of time variations of the exterior shape of natural bodies are to be created. Variables in this list are typically
        assigned by using a function from the :ref:`\`\`shape_deformation\`\`` module.


        :type: list[BodyDeformationSettings]
     )doc" )
            .def_readwrite( "ground_station_settings", &tss::BodySettings::groundStationSettings, R"doc(No documentation found.)doc" )
            .def_readwrite( "rigid_body_settings",
                            &tss::BodySettings::rigidBodyPropertiesSettings,
                            R"doc(

        Object that defines the settings of the body rigid body (mass, center of mass, inertia) properties that are to be created. A variable of this type is typically
        assigned by using a function from the :ref:`\`\`rigid_body\`\`` module. Note that this setting does *not* define
        the gravity field, but rather only the mass, center of mass and inertia tensor.


        :type: RigidBodyPropertiesSettings
     )doc" )
            .def_readwrite( "radiation_pressure_target_settings",
                            &tss::BodySettings::radiationPressureTargetModelSettings,
                            R"doc(

        Object that defines the settings of the radiation pressure target model that is to be created. A variable of this type is typically
        assigned by using a function from the :ref:`\`\`radiation_pressure\`\`` module. 


        :type: RadiationPressureTargetModelSettings
     )doc" )
            .def_readwrite( "radiation_source_settings",
                            &tss::BodySettings::radiationSourceModelSettings,
                            R"doc(

        Object that defines the settings of the radiation source model that is to be created. A variable of this type is typically
        assigned by using a function from the :ref:`\`\`radiation_pressure\`\`` module. 


        :type: RadiationSourceModelSettings
     )doc" )
            .def_readwrite( "vehicle_shape_settings", &tss::BodySettings::bodyExteriorPanelSettings_, R"doc(

        Object that defines the settings of an exterior panelled vehicle shape that is to be created. A variable of this type is typically
        assigned by using a function from the :ref:`\`\`vehicle_systems\`\`` module.


        :type: FullPanelledBodySettings
     )doc" )
            .def_readwrite( "radiation_pressure_settings",
                            &tss::BodySettings::radiationPressureSettings,
                            R"doc(

        .. warning::

            This interface is deprecated and will be removed in a future release. Use :attr:`~tudatpy.numerical_simulation.environment_setup.BodySettings.radiation_source_settings` and :attr:`~tudatpy.numerical_simulation.environment_setup.BodySettings.radiation_pressure_target_settings` instead.


     )doc" );

    py::class_< tss::BodyListSettings, std::shared_ptr< tss::BodyListSettings > >( m, "BodyListSettings", R"doc(

        Class for defining settings for the creation of a system of bodies.

        Class for defining settings for the creation of a system of bodies. This object is typically created from default settings, and
        then adapted to the user's specific needs.





     )doc" )
            .def( py::init< const std::string, const std::string >( ), py::arg( "frame_origin" ), py::arg( "frame_orientation" ), R"doc(

        Class initialization method.

        Class method to initialize an empty BodyListSettings object.

        .. note::

            When creating BodyListSettings from this method, the settings for each body will have to be added manually.
            It is typically more convenient to use the :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings` function to create a BodyListSettings object with default settings for all bodies, and then modify the settings for specific bodies as needed.


        Parameters
        ----------
        frame_origin : str
            Definition of the global frame origin for the bodies.
        frame_orientation : str
            Definition of the global frame orientation for the bodies.


     )doc" )
            .def( "get", &tss::BodyListSettings::get, py::arg( "body_name" ), R"doc(

        This function extracts a single BodySettings object.


        Parameters
        ----------
        body_name : str
            Name of the body for which settings are to be retrieved


        Returns
        -------
        BodySettings
            Settings for the requested body


    )doc" )
            .def( "add_settings",
                  py::overload_cast< std::shared_ptr< tss::BodySettings >, const std::string >( &tss::BodyListSettings::addSettings ),
                  py::arg( "settings_to_add" ),
                  py::arg( "body_name" ),
                  R"doc(

        Add a single :class:`BodySettings` object to the :class:`BodyListSettings` instance.

        .. warning::

            This method is rarely called by the user, as :class:`BodySettings` objects cannot be created directly but only be extracted from a BodyListSettings instance.
            Instead, users are recommended to use the :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings` to create settings for major celestial bodies, and the :func:`~tudatpy.numerical_simulation.environment_setup.BodyListSettings.add_empty_settings` function to create settings for custom bodies.
            See the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/creation_celestial_body_settings.html>`_ for more information.


        Parameters
        ----------
        settings_to_add : BodySettings
            Settings to be added
        body_name : str
            Name of the body for which settings are added




    )doc" )
            .def( "add_empty_settings",
                  py::overload_cast< const std::string >( &tss::BodyListSettings::addSettings ),
                  py::arg( "body_name" ),
                  R"doc(

        This method adds empty settings to the :class:`BodyListSettings` instance.

        Adds empty settings to the :class:`BodyListSettings` instance. This is typically used to add settings for custom bodies, for which no default settings are available.
        See the `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/creation_celestial_body_settings.html>`_ for more information. 

        Parameters
        ----------
        body_name : str
            Name of the body for which settings are added




    )doc" )
            .def_property_readonly( "frame_origin",
                                    &tss::BodyListSettings::getFrameOrigin,
                                    R"doc(

        **read-only**

        Definition of the global frame origin for the bodies

        :type: str
     )doc" )
            .def_property_readonly( "frame_orientation", &tss::BodyListSettings::getFrameOrientation, R"doc(

        **read-only**

        Definition of the global frame orientation for the bodies

        :type: str
     )doc" );

    m.def( "get_default_body_settings",
           py::overload_cast< const std::vector< std::string > &, const std::string, const std::string >( &tss::getDefaultBodySettings ),
           py::arg( "bodies" ),
           py::arg( "base_frame_origin" ) = "SSB",
           py::arg( "base_frame_orientation" ) = "ECLIPJ2000",
           R"doc(

Function that retrieves the default settings for the given set of input bodies.

Function that retrieves the default settings for the given set of input bodies. Default settings are described in
detail `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/default_env_models.html>`_ .
Note that if a body is provided as input for which default settings do not exist, an exception is thrown. In addition
to settings for each separate body, this function returns an object that defines the global frame origin and orientation,


Parameters
----------
bodies : list[str]
    List of name of bodies for which default settings are to be retrieved and created.
base_frame_origin : str, default = 'SSB'
    Base frame origin of the set of bodies that is to be created. It defaults to the solar system barycenter (SSB), but it can by any of the bodies in `bodies_to_create` (provided it has an ephemeris defined).
base_frame_orientation : str, default = 'ECLIPJ2000'
    Base frame orientation of the set of bodies that is to be created. It can be either ECLIPJ2000 (default) or J2000.
Returns
-------
BodyListSettings
    Object containing the settings for the SystemOfBodies that are to be created






    )doc" );

    m.def( "get_default_body_settings_time_limited",
           py::overload_cast< const std::vector< std::string > &,
                              const double,
                              const double,
                              const std::string,
                              const std::string,
                              const double >( &tss::getDefaultBodySettings ),
           py::arg( "bodies" ),
           py::arg( "initial_time" ),
           py::arg( "final_time" ),
           py::arg( "base_frame_origin" ) = "SSB",
           py::arg( "base_frame_orientation" ) = "ECLIPJ2000",
           py::arg( "time_step" ) = 300.0,
           R"doc(

Function that retrieves the default settings for the given set of input bodies, with a limited valid time interval.

Same as :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings`, but with body settings valid over a limited time interval. This makes the
the extraction of states from ephemerides more computationally efficient, at the expense of more RAM usage, and a
constrained time interval over which the ephemerides are valid. See `this page <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/environment_setup/default_env_models/default_bodies_limited_time_range.html>`_ for more details.


Parameters
----------
bodies : list[str]
    List of name of bodies for which default settings are to be retrieved and created.
initial_time : float
    Start time from which the environment settings should be created.
final_time : float
    End time up to which the environment settings should be created.
base_frame_origin : str
    Base frame origin of the set of bodies that is to be created.
base_frame_orientation : str
    Base frame orientation of the set of bodies that is to be created.
time_step : float, default = 300.0
    Time step to be used for the tabulated ephemeris.
Returns
-------
BodyListSettings
    Object containing the settings for the SystemOfBodies that are to be created






    )doc" );

    m.def( "get_default_single_body_settings",
           py::overload_cast< const std::string &, const std::string & >( &tss::getDefaultSingleBodySettings ),
           py::arg( "body_name" ),
           py::arg( "base_frame_orientation" ) = "ECLIPJ2000",
           R"doc(

Function that retrieves the default settings for a single body.

As :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings`, but for retrieving default settings of only a single body


Parameters
----------
body_name : str
    Name of body for which default settings are to be retrieved and created.
base_frame_orientation : str, default = 'ECLIPJ2000'
    Base frame orientation of the body settings. It can be either ECLIPJ2000 (default) or J2000.
Returns
-------
BodySettings
    Object containing the settings for the body that is to be created






    )doc" );

    m.def( "get_default_single_body_settings_time_limited",
           py::overload_cast< const std::string &, const double, const double, const std::string &, const double >(
                   &tss::getDefaultSingleBodySettings ),
           py::arg( "body_name" ),
           py::arg( "initial_time" ),
           py::arg( "final_time" ),
           py::arg( "base_frame_orientation" ) = "ECLIPJ2000",
           py::arg( "time_step" ) = 300.0,
           R"doc(

Function that retrieves the default settings for a single body, with a limited valid time interval.

As :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings_time_limited`, but for retrieving default settings of only a single body


Parameters
----------
body_name : str
    Name of body for which default settings are to be retrieved and created.
initial_time : float
    Start time from which the environment settings should be created.
final_time : float
    End time up to which the environment settings should be created.
base_frame_orientation : str, default = 'ECLIPJ2000'
    Base frame orientation of the body settings. It can be either ECLIPJ2000 (default) or J2000.
time_step : float, default = 300.0
    Time step to be used for the tabulated ephemeris.
Returns
-------
BodySettings
    Object containing the settings for the body that is to be created






    )doc" );

    m.def( "get_default_single_alternate_body_settings",
           py::overload_cast< const std::string &, const std::string &, const std::string & >(
                   &tss::getDefaultSingleAlternateNameBodySettings ),
           py::arg( "body_name" ),
           py::arg( "source_body_name" ),
           py::arg( "base_frame_orientation" ) = "ECLIPJ2000",
           R"doc(

Function that retrieves the default settings for a single body, and assigns them to another body.

As :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings`, but for retrieving default settings of only a single body,
where the default settings of body with name ``source_body_name`` are retrieved and assigned to a body with name ``body_name``.
For instance, if ``source_body_name`` is set to "Mars", and ````body_name`` is set to "Earth" body name Earth will be created, with all the properties
of Mars


Parameters
----------
body_name : str
    Name of body for which default settings are to be created.
source_body_name : str
    Name of body for which default settings are to be retrieved, and assigned to a body with name ``body_name``.
base_frame_orientation : str, default = 'ECLIPJ2000'
    Base frame orientation of the body settings. It can be either ECLIPJ2000 (default) or J2000.
Returns
-------
BodySettings
    Object containing the settings for the body that is to be created






    )doc" );

    m.def( "get_default_single_alternate_body_settings_time_limited",
           py::overload_cast< const std::string &, const std::string &, const double, const double, const std::string &, const double >(
                   &tss::getDefaultSingleAlternateNameBodySettings ),
           py::arg( "body_name" ),
           py::arg( "source_body_name" ),
           py::arg( "initial_time" ),
           py::arg( "final_time" ),
           py::arg( "base_frame_orientation" ) = "ECLIPJ2000",
           py::arg( "time_step" ) = 300.0,
           R"doc(

Function that retrieves the default settings for a single body, with a limited valid time interval.

As :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings_time_limited`, but for retrieving default settings of only a single body,
where the default settings of body with name ``source_body_name`` are retrieved and assigned to a body with name ``body_name``.
For instance, if ``source_body_name`` is set to "Mars", and ````body_name`` is set to "Earth" body name Earth will be created, with all the properties
of Mars


Parameters
----------
body_name : str
    Name of body for which default settings are to be retrieved.
source_body_name : str
    Name of body for which default settings are to be retrieved, and assigned to a body with name ``body_name``.
initial_time : float
    Start time from which the environment settings should be created.
final_time : float
    End time up to which the environment settings should be created.
base_frame_orientation : str, default = 'ECLIPJ2000'
    Base frame orientation of the body settings. It can be either ECLIPJ2000 (default) or J2000.
time_step : float, default = 300.0
    Time step to be used for the tabulated ephemeris.
Returns
-------
BodySettings
    Object containing the settings for the body that is to be created






    )doc" );

    m.def( "create_simplified_system_of_bodies",
           &tss::createSimplifiedSystemOfBodies,
           py::arg( "initial_time" ) = 0,
           R"doc(

Function that creates a simplified System of bodies.

Function that creates a simplified system of bodies. The following bodies are created in this system: the Sun, all planets of the Solar system, and Pluto.
All bodies in this system use Gtop ephemerides and point mass gravity. The Earth is setup with a spherical shape model and a simple rotation model.
The reference frame used to setup this simplified system of bodies has its origin at the SSB, and has an ECLIPJ2000 orientation.


Parameters
----------
initial_time : float, optional, default=0
    Initial system time in seconds since J2000.
Returns
-------
:class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
    Object containing the objects for bodies and environment models constituting the physical environment






    )doc" );

    m.def( "create_system_of_bodies",
           &tss::createSystemOfBodies< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "body_settings" ),
           R"doc(

Function that creates a System of bodies from associated settings.

Function that creates a System of bodies from associated settings. This function creates the separate :class:`~tudatpy.numerical_simulation.Body`
objects and stores them in a :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` object. This object represents the full
physical environment in the simulation.


Parameters
----------
body_settings : BodyListSettings
    Object defining the physical environment, with all properties of artificial and natural bodies.
Returns
-------
:class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
    Object containing the objects for bodies and environment models constituting the physical environment






    )doc" );

    m.def( "add_empty_tabulated_ephemeris",
           &tp::addEmptyTabulatedEphemeris< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "bodies" ),
           py::arg( "body_name" ),
           py::arg( "ephemeris_origin" ) = "",
           py::arg( "is_part_of_multi_arc" ) = false,
           R"doc(No documentation found.)doc" );

    m.def( "create_tabulated_ephemeris_from_spice",
           &tss::createTabulatedEphemerisFromSpice< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "body" ),
           py::arg( "initial_time" ),
           py::arg( "end_time" ),
           py::arg( "time_step" ),
           py::arg( "observer_name" ),
           py::arg( "reference_frame_name" ),
           py::arg( "interpolator_settings" ) = std::make_shared< tudat::interpolators::LagrangeInterpolatorSettings >( 8 ) );

    m.def( "create_body_ephemeris",
           &tss::createBodyEphemeris< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "ephemeris_settings" ),
           py::arg( "body_name" ),
           R"doc(

Function that creates an Ephemeris object.

Function that creates an :class:`~tudatpy.numerical_simulation.environment.Ephemeris` object, but does *not*
associate it with any specific body (e.g., it does not go into the environment, but can be used independently of it)


Parameters
----------
ephemeris_settings : EphemerisSettings
    Object defining the ephemeris settings.
body_name : str
    Name of body for which the ephemeris is created. Note that this input is only relevant for some ephemeris settings (for instance, a spice ephemeris setting), and it does *not* imply that the ephemeris object is associated with a Body object of this name.
Returns
-------
:class:`~tudatpy.numerical_simulation.environment.Ephemeris`
    Ephemeris object, created according to the provided settings






    )doc" );

    m.def( "create_ground_station_ephemeris",
           py::overload_cast< const std::shared_ptr< tss::Body >, const std::string &, const tss::SystemOfBodies & >(
                   &tss::createReferencePointEphemerisFromId< TIME_TYPE, STATE_SCALAR_TYPE > ),
           "body_with_ground_station",
           "station_name" );

    m.def( "get_safe_interpolation_interval", &tss::getSafeInterpolationInterval, py::arg( "ephemeris_model" ) );

    m.def( "add_aerodynamic_coefficient_interface",
           &tss::addAerodynamicCoefficientInterface,
           py::arg( "bodies" ),
           py::arg( "body_name" ),
           py::arg( "coefficient_settings" ),
           R"doc(

Function that creates an aerodynamic coefficient interface from settings, and adds it to an existing body.

This function can be used to add an aerodynamic coefficient interface to an existing body. It requires
settings for the aerodynamic coefficients, created using one of the functions from the `~tudatpy.numerical_simulation_environment_setup.aerodynamic_coefficient` module.
This function creates the actual coefficient interface from these settings, and assigns it to the
selected body. In addition to the identifier for the body to which it is assigned, this function
requires the full :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` as input, to facilitate
inter-body dependencies in the coefficient interface

    m.def( "add_aerodynamic_coefficient_interface",
           &tss::addAerodynamicCoefficientInterface,
           py::arg( "bodies" ),
           py::arg( "body_name" ),
           py::arg( "coefficient_settings" ),
           get_docstring( "add_aerodynamic_coefficient_interface" ).c_str( ) );

Parameters
----------
bodies : SystemOfBodies
    Object defining the physical environment, with all properties of artificial and natural bodies.
body_name : str
    Name of the body to which the aerodynamic coefficients are to be assigned
coefficient_settings : AerodynamicCoefficientSettings
    Settings defining the coefficient interface that is to be created.





    )doc" );

    m.def( "create_aerodynamic_coefficient_interface",
           &tss::createAerodynamicCoefficientInterfaceDeprecated,
           py::arg( "coefficient_settings" ),
           py::arg( "body" ) );

    m.def( "create_aerodynamic_coefficient_interface",
           &tss::createAerodynamicCoefficientInterface,
           py::arg( "coefficient_settings" ),
           py::arg( "body" ),
           py::arg( "bodies" ),
           R"doc(No documentation found.)doc" );

    m.def( "add_radiation_pressure_interface",
           &tss::addRadiationPressureInterface,
           py::arg( "bodies" ),
           py::arg( "body_name" ),
           py::arg( "radiation_pressure_settings" ) );

    m.def( "add_radiation_pressure_target_model",
           &tss::addRadiationPressureTargetModel,
           py::arg( "bodies" ),
           py::arg( "body_name" ),
           py::arg( "radiation_pressure_target_settings" ),
           R"doc(

Function that creates an radiation pressure interface from settings, and adds it to an existing body.

This function can be used to add an radiation pressure interface to an existing body. It requires
settings for the radiation pressure interface, created using one of the functions from the :ref:`\`\`radiation_pressure\`\`` module.
This function creates the actual coefficient interface from these settings, and assigns it to the
selected body. In addition to the identifier for the body to which it is assigned, this function
requires the full :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` as input, to facilitate
inter-body dependencies in the radiation pressure interface


Parameters
----------
bodies : SystemOfBodies
    Object defining the physical environment, with all properties of artificial and natural bodies.
body_name : str
    Name of the body to which the radiation pressure interface is to be assigned
radiation_pressure_settings : RadiationPressureInterfaceSettings
    Settings defining the radiation pressure interface that is to be created.





    )doc" );

    m.def( "add_rotation_model",
           &tss::addRotationModel,
           py::arg( "bodies" ),
           py::arg( "body_name" ),
           py::arg( "rotation_model_settings" ),
           R"doc(

Function that creates a rotation model, and adds it to an existing body.

This function can be used to add  a :class:`~tudatpy.numerical_simulation.environment.RotationalEphemeris` object to an existing body.
Typically, the ``RotationalEphemeris`` is created along with the :class:`~tudatpy.numerical_simulation.environment.Body` itself. However, in some cases it may be useful
to create a rotation model after the Body objects have been created. This function requires
settings for the rotation model, created using one of the functions from the :ref:`rotation_model` module.
This function creates the actual coefficient interface from these settings, and assigns it to the
selected body. In addition to the identifier for the body to which it is assigned, this function
requires the full :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` as input, to facilitate
inter-body dependencies in the radiation model


Parameters
----------
bodies : SystemOfBodies
    Object defining the physical environment, with all properties of artificial and natural bodies.
body_name : str
    Name of the body to which the rotation model is to be assigned
rotation_model_settings
    Settings defining the rotation model that is to be created.





    )doc" );

    m.def( "add_gravity_field_model",
           &tss::addGravityFieldModel,
           py::arg( "bodies" ),
           py::arg( "body_name" ),
           py::arg( "gravity_field_settings" ),
           py::arg( "gravity_field_variation_settings" ) = std::vector< std::shared_ptr< tss::GravityFieldVariationSettings > >( ),
           R"doc(No documentation found.)doc" );

    m.def( "add_mass_properties_model",
           &tss::addRigidBodyProperties,
           py::arg( "bodies" ),
           py::arg( "body_name" ),
           py::arg( "mass_property_settings" ) );

    m.def( "add_rigid_body_properties",
           &tss::addRigidBodyProperties,
           py::arg( "bodies" ),
           py::arg( "body_name" ),
           py::arg( "rigid_body_property_settings" ),
           R"doc(
        
Function that creates a rigid body property model, and adds it to an existing body.

This function can be used to add a :class:`~tudatpy.numerical_simulation.environment.RigidBodyProperties` object to an existing body.
Typically, the ``RigidBodyProperties`` are created along with the :class:`~tudatpy.numerical_simulation.environment.Body` itself. However, in some cases it may be useful
to create body mass properties after the Body objects have been created. This function requires
settings for the rigid body properties, created using one of the functions from the :ref:`rigid_body` module.
This function creates the actual rigid body properties from these settings, and assigns it to the
selected body.


Parameters
----------
bodies : SystemOfBodies
    Object defining the physical environment, with all properties of artificial and natural bodies.
body_name : str
    Name of the body to which the model is to be assigned
rigid_body_property_settings : RigidBodyPropertiesSettings
    Settings defining the rigid body properties model that is to be created.





    )doc" );

    m.def( "add_engine_model",
           &tss::addEngineModel,
           py::arg( "body_name" ),
           py::arg( "engine_name" ),
           py::arg( "thrust_magnitude_settings" ),
           py::arg( "bodies" ),
           py::arg( "body_fixed_thrust_direction" ) = Eigen::Vector3d::UnitX( ),
           R"doc(

Function that creates an engine model (to be used for thrust calculations), and adds it to an existing body.

Function that creates an engine model (to be used for thrust calculations), and adds it to an existing body. It creates and
object of class :class:`~tudatpy.numerical_simulation.environment.EngineModel`, and adds it to an existing body. Properties
assigned to this engine model are:

* The (constant) direction in body-fixed frame in which the engine is pointing (e.g. the body-fixed thrust direction when the engine is on)
* Settings for computing the thrust magnitude (as a function of time and/or other parameters), using a suitable function from the :ref:`\`\`thrust\`\`` submodule


Parameters
----------
body_name : str
    Name of the body to which the engine is to be added.
engine_name : str
    Name (e.g. unique identifier) of the engine that is to be added to the body
thrust_magnitude_settings : ThrustMagnitudeSettings
    Settings for computing the thrust magnitude (and specific impulse) as a function of time
bodies : SystemOfBodies
    Object defining the physical environment, with all properties of artificial and natural bodies.
body_fixed_thrust_direction : numpy.ndarray[numpy.float64[3, 1]], default = [1,0,0]
    Unit vector along which the thrust from the engine will point in a body-fixed frame





    )doc" );

    m.def( "add_variable_direction_engine_model",
           &tss::addVariableDirectionEngineModel,
           py::arg( "body_name" ),
           py::arg( "engine_name" ),
           py::arg( "thrust_magnitude_settings" ),
           py::arg( "bodies" ),
           py::arg( "body_fixed_thrust_direction_function" ),
           R"doc(

Function that creates an engine model (to be used for thrust calculations), and adds it to an existing body.

Same as :func:`add_engine_model`, but with a time-variable body-fixed thrust direction


Parameters
----------
body_name : str
    Name of the body to which the engine is to be added.
engine_name : str
    Name (e.g. unique identifier) of the engine that is to be added to the body
thrust_magnitude_settings : ThrustMagnitudeSettings
    Settings for computing the thrust magnitude (and specific impulse) as a function of time
bodies : SystemOfBodies
    Object defining the physical environment, with all properties of artificial and natural bodies.
body_fixed_thrust_direction_function : Callable[[float], numpy.ndarray[numpy.float64[3, 1]]]
    Function returning a unit vector, as a function of time, along which the thrust from the engine will point in a body-fixed frame





    )doc" );

    m.def( "add_flight_conditions",
           &tss::addFlightConditions,
           py::arg( "bodies" ),
           py::arg( "body_name" ),
           py::arg( "central_body_name" ),
           R"doc(

Function that creates a flight conditions, and adds it to an existing body.

This function can be used to add  a :class:`~tudatpy.numerical_simulation.environment.FlightConditions` object to an existing body.
Typically, the ``FlightConditions`` are created automatically when they are required (for the calculation of an
aerodynamic acceleration, or the saving of certain dependent variables). However, in some cases it may be useful
to manually trigger their creation, which is done through this function. If the ``central_body_name`` input
denotes a body that is endowed with an :class:`~tudatpy.numerical_simulation.environment.AtmosphereModel`, this function
automatically creates an :class:`~tudatpy.numerical_simulation.environment.AtmosphericFlightConditions` object (capable of
calculating density, speed of sound, etc.), instead of the more basic :class:`~tudatpy.numerical_simulation.environment.FlightConditions`
(which is limited to properties such as altitude, latitude, etc.)


Parameters
----------
bodies : SystemOfBodies
    Object defining the physical environment, with all properties of artificial and natural bodies.
body_name : str
    Name of the body for which the flight conditions are to be created
central_body_name : str
    Name of the central body w.r.t. which the flight conditions are to be created (typically, but not necessarily, the central body of propagation)/





    )doc" );

    m.def( "convert_ground_station_state_between_itrf_frames",
           &trf::convertGroundStationStateBetweenItrfFrames,
           py::arg( "ground_station_state" ),
           py::arg( "epoch" ),
           py::arg( "base_frame" ),
           py::arg( "target_frame" ),
           R"doc(No documentation found.)doc" );

    m.def( "add_ground_station",
           py::overload_cast< const std::shared_ptr< tss::Body >,
                              const std::string,
                              const Eigen::Vector3d,
                              const tcc::PositionElementTypes,
                              const std::vector< std::shared_ptr< tss::GroundStationMotionSettings > > >( &tss::createGroundStation ),
           py::arg( "body" ),
           py::arg( "ground_station_name" ),
           py::arg( "ground_station_position" ),
           py::arg( "position_type" ) = tcc::cartesian_position,
           py::arg( "station_motion_settings" ) = std::vector< std::shared_ptr< tss::GroundStationMotionSettings > >( ) );

    m.def( "add_ground_station",
           py::overload_cast< const std::shared_ptr< tss::Body >, const std::shared_ptr< tss::GroundStationSettings > >(
                   &tss::createGroundStation ),
           py::arg( "body" ),
           py::arg( "ground_station_settings" ),
           R"doc(No documentation found.)doc" );

    m.def( "create_radiation_pressure_interface",
           &tss::createRadiationPressureInterface,
           py::arg( "radiationPressureInterfaceSettings" ),
           py::arg( "body_name" ),
           py::arg( "body_dict" ) );

    m.def( "get_ground_station_list", &tss::getGroundStationsLinkEndList, py::arg( "body" ) );

    //        m.def("get_target_elevation_angles",
    //              &tss::getTargetElevationAngles,
    //              py::arg( "observing_body" ),
    //              py::arg( "target_body" ),
    //              py::arg( "station_name" ),
    //              py::arg( "times" ) );

    auto aerodynamic_coefficient_setup = m.def_submodule( "aerodynamic_coefficients" );
    aerodynamic_coefficients::expose_aerodynamic_coefficient_setup( aerodynamic_coefficient_setup );

    auto radiation_pressure_setup = m.def_submodule( "radiation_pressure" );
    radiation_pressure::expose_radiation_pressure_setup( radiation_pressure_setup );

    auto rotation_model_setup = m.def_submodule( "rotation_model" );
    rotation_model::expose_rotation_model_setup( rotation_model_setup );

    auto gravity_field_setup = m.def_submodule( "gravity_field" );
    gravity_field::expose_gravity_field_setup( gravity_field_setup );

    auto ephemeris_setup = m.def_submodule( "ephemeris" );
    ephemeris::expose_ephemeris_setup( ephemeris_setup );

    auto atmosphere_setup = m.def_submodule( "atmosphere" );
    atmosphere::expose_atmosphere_setup( atmosphere_setup );

    auto shape_setup = m.def_submodule( "shape" );
    shape::expose_shape_setup( shape_setup );

    auto gravity_variation_setup = m.def_submodule( "gravity_field_variation" );
    gravity_field_variation::expose_gravity_field_variation_setup( gravity_variation_setup );

    auto shape_deformation_setup = m.def_submodule( "shape_deformation" );
    shape_deformation::expose_shape_deformation_setup( shape_deformation_setup );

    auto ground_station_setup = m.def_submodule( "ground_station" );
    ground_station::expose_ground_station_setup( ground_station_setup );

    auto rigid_body_setup = m.def_submodule( "rigid_body" );
    rigid_body::expose_rigid_body_setup( rigid_body_setup );

    auto vehicle_systems_setup = m.def_submodule( "vehicle_systems" );
    vehicle_systems::expose_vehicle_systems_setup( vehicle_systems_setup );

    //        auto system_model_setup =
    //        m.def_submodule("system_models");
    //        gravity_field_variation::expose_system_model_setup(system_model_setup);

    // Function removed; error is shown
    m.def( "set_aerodynamic_guidance",
           py::overload_cast< const std::shared_ptr< ta::AerodynamicGuidance >, const std::shared_ptr< tss::Body >, const bool >(
                   &tss::setGuidanceAnglesFunctions ),
           py::arg( "aerodynamic_guidance" ),
           py::arg( "body" ),
           py::arg( "silence_warnings" ) = false );

    // Function removed; error is shown
    m.def( "set_aerodynamic_orientation_functions",
           &tss::setAerodynamicOrientationFunctions,
           py::arg( "body" ),
           py::arg( "angle_of_attack_function" ) = std::function< double( ) >( ),
           py::arg( "sideslip_angle_function" ) = std::function< double( ) >( ),
           py::arg( "bank_angle_function" ) = std::function< double( ) >( ),
           py::arg( "update_function" ) = std::function< void( const double ) >( ) );

    // Function removed; error is shown
    m.def( "set_constant_aerodynamic_orientation",
           &tss::setConstantAerodynamicOrientation,
           py::arg( "body" ),
           py::arg( "angle_of_attack" ),
           py::arg( "sideslip_angle" ),
           py::arg( "bank_angle" ),
           py::arg( "silence_warnings" ) = false );
}

}  // namespace environment_setup
}  // namespace numerical_simulation
}  // namespace tudatpy
