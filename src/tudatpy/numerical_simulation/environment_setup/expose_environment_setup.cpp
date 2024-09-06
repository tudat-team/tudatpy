/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/astro/reference_frames/referenceFrameTransformations.h>
#include <tudat/simulation/environment_setup.h>

#include "tudatpy/scalarTypes.h"

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


namespace tudatpy {
    namespace numerical_simulation {
        namespace environment_setup {

            PYBIND11_MODULE(expose_environment_setup, m) {
                py::module_::import("tudatpy.math.interpolators");
                py::module_::import(
                    "tudatpy.numerical_simulation.environment_setup.gravity_"
                    "field_variation");
                py::module_::import(
                    "tudatpy.numerical_simulation.environment_setup."
                    "ground_station");
                py::module_::import("tudatpy.astro.element_conversion");
                //        m.def("get_body_gravitational_parameter",
                //              &tss::getBodyGravitationalParameter,
                //              py::arg("body_collection"),
                //              py::arg("body_name"));


                py::class_<tss::BodySettings,
                           std::shared_ptr<tss::BodySettings>>(
                    m, "BodySettings",
                    R"doc(Class for defining settings for the creation of a single body.

	Class for defining settings for the creation of a single body, this object is typically stored inside a
	:class:`BodyListSettings`, object.

)doc")
                    .def_readwrite(
                        "constant_mass", &tss::BodySettings::constantMass,
                        R"doc(Mass that gets assigned to the vehicle. Note that this mass does *not* automatically define a gravity field
model, but is instead used for the calculation of non-conservative forces only. When creating a body with a gravity field,
leave this entry empty.

	)doc")
                    .def_readwrite(
                        "atmosphere_settings",
                        &tss::BodySettings::atmosphereSettings,
                        R"doc(Object that defines the settings of the atmosphere model that is to be created. Note that wind model settings
may be defined inside this object. A variable of this type is typically assigned by using a factory function from the
:ref:`\`\`atmosphere\`\`` module.

	)doc")
                    .def_readwrite(
                        "ephemeris_settings",
                        &tss::BodySettings::ephemerisSettings,
                        R"doc(Object that defines the settings of the ephemeris model that is to be created. A variable of this type is typically
assigned by using a factory function from the :ref:`\`\`ephemeris\`\`` module.

	)doc")
                    .def_readwrite(
                        "gravity_field_settings",
                        &tss::BodySettings::gravityFieldSettings,
                        R"doc(Object that defines the settings of the gravity field model that is to be created. A variable of this type is typically
assigned by using a factory function from the :ref:`\`\`gravity_field\`\`` module.

	)doc")
                    .def_readwrite(
                        "rotation_model_settings",
                        &tss::BodySettings::rotationModelSettings,
                        R"doc(Object that defines the settings of the rotation model that is to be created. A variable of this type is typically
assigned by using a factory function from the :ref:`\`\`rotation_model\`\`` module.

	)doc")
                    .def_readwrite(
                        "shape_settings",
                        &tss::BodySettings::shapeModelSettings,
                        R"doc(Object that defines the settings of the shape model that is to be created. A variable of this type is typically
assigned by using a factory function from the :ref:`\`\`shape\`\`` module.

	)doc")
                    .def_readwrite(
                        "aerodynamic_coefficient_settings",
                        &tss::BodySettings::aerodynamicCoefficientSettings,
"")
                    .def_readwrite(
                        "gravity_field_variation_settings",
                        &tss::BodySettings::gravityFieldVariationSettings,
"")
                    .def_readwrite(
                        "shape_deformation_settings",
                        &tss::BodySettings::bodyDeformationSettings,
"")
                    .def_readwrite(
                        "ground_station_settings",
                        &tss::BodySettings::groundStationSettings,
"")
                    .def_readwrite(
                        "rigid_body_settings",
                        &tss::BodySettings::rigidBodyPropertiesSettings,
"")
                    .def_readwrite(
                        "radiation_pressure_target_settings",
                        &tss::BodySettings::
                            radiationPressureTargetModelSettings,
"")
                    .def_readwrite(
                        "radiation_source_settings",
                        &tss::BodySettings::radiationSourceModelSettings,
"")
                    .def_readwrite(
                        "vehicle_shape_settings",
                        &tss::BodySettings::bodyExteriorPanelSettings_,
"")
                    .def_readwrite(
                        "radiation_pressure_settings",
                        &tss::BodySettings::radiationPressureSettings,
"");


                py::class_<tss::BodyListSettings,
                           std::shared_ptr<tss::BodyListSettings>>(
                    m, "BodyListSettings",
                    R"doc(Class for defining settings for the creation of a system of bodies.

	Class for defining settings for the creation of a system of bodies. This object is typically created from default settings, and
	then adapted to the user's specific needs.

)doc")
                    .def(py::init<const std::string, const std::string>(),
                         py::arg("frame_origin"), py::arg("frame_orientation"))
                    .def(
                        "get", &tss::BodyListSettings::get,
                        R"doc(This function extracts a single BodySettings object .

	:param body_name:
		Name of the body for which settings are to be retrieved

)doc")
                    .def("add_settings",
                         py::overload_cast<std::shared_ptr<tss::BodySettings>,
                                           const std::string>(
                             &tss::BodyListSettings::addSettings),
                         py::arg("settings_to_add"), py::arg("body_name"))
                    .def("add_empty_settings",
                         py::overload_cast<const std::string>(
                             &tss::BodyListSettings::addSettings),
                         py::arg("body_name"))
                    .def_property_readonly(
                        "frame_origin", &tss::BodyListSettings::getFrameOrigin,
                        R"doc(Definition of the global frame origin for the bodies
	)doc")
                    .def_property_readonly(
                        "frame_orientation",
                        &tss::BodyListSettings::getFrameOrientation,
                        R"doc(Definition of the global frame orientation for the bodies
	)doc");

                m.def(
                    "get_default_body_settings",
                    py::overload_cast<const std::vector<std::string> &,
                                      const std::string, const std::string>(
                        &tss::getDefaultBodySettings),
                    py::arg("bodies"), py::arg("base_frame_origin") = "SSB",
                    py::arg("base_frame_orientation") = "ECLIPJ2000",
                    R"doc(Function that retrieves the default settings for the given set of input bodies.

	Function that retrieves the default settings for the given set of input bodies. Default settings are described in
	detail `here <https://tudat-space.readthedocs.io/en/latest/_src_user_guide/state_propagation/environment_setup/create_bodies/default_settings.html>`_ .
	Note that if a body is provided as input for which default settings do not exist, an exception is thrown. In addition
	to settings for each separate body, this function returns an object that defines the global frame origin and orientation,


	:param bodies:
		List of name of bodies for which default settings are to be retrieved.
	:param base_frame_origin:
		Base frame origin of the set of bodies that is to be created. It defaults to the solar system barycenter (SSB), but it can by any of the bodies in `bodies_to_create` (provided it has an ephemeris defined).
	:param base_frame_orientation:
		Base frame orientation of the set of bodies that is to be created. It can be either ECLIPJ2000 (default) or J2000.
	:return:
		Object containing the settings for the SystemOfBodies that are to be created
)doc");

                m.def(
                    "get_default_body_settings_time_limited",
                    py::overload_cast<const std::vector<std::string> &,
                                      const double, const double,
                                      const std::string, const std::string,
                                      const double>(
                        &tss::getDefaultBodySettings),
                    py::arg("bodies"), py::arg("initial_time"),
                    py::arg("final_time"), py::arg("base_frame_origin") = "SSB",
                    py::arg("base_frame_orientation") = "ECLIPJ2000",
                    py::arg("time_step") = 300.0,
                    R"doc(Function that retrieves the default settings for the given set of input bodies, with a limited valid time interval.

	Same as :func:`~tudatpy.numerical_simulation.environment_setup.get_default_body_settings`, but with body settings valid over a limited time interval. This makes the
	the extraction of states from ephemerides more computationally efficient, at the expense of more RAM usage, and a
	constrained time interval over which the ephemerides are valid. See `this page <https://tudat-space.readthedocs.io/en/latest/_src_user_guide/state_propagation/environment_setup/valid_time_range.html>`_ for more details.


	:param bodies:
		List of name of bodies for which default settings are to be retrieved.
	:param initial_time:
		Start time from which the environment settings should be created.
	:param final_time:
		End time up to which the environment settings should be created.
	:param base_frame_origin:
		Base frame origin of the set of bodies that is to be created.
	:param base_frame_orientation:
		Base frame orientation of the set of bodies that is to be created.
	:param time_step:
		Time step to be used for the tabulated ephemeris.
	:return:
		Object containing the settings for the SystemOfBodies that are to be created
)doc");

                m.def(
                    "get_default_single_body_settings",
                    py::overload_cast<const std::string &, const std::string &>(
                        &tss::getDefaultSingleBodySettings),
                    py::arg("body_name"),
                    py::arg("base_frame_orientation") = "ECLIPJ2000",
"");

                m.def("get_default_single_body_settings_time_limited",
                      py::overload_cast<const std::string &, const double,
                                        const double, const std::string &,
                                        const double>(
                          &tss::getDefaultSingleBodySettings),
                      py::arg("body_name"), py::arg("initial_time"),
                      py::arg("final_time"),
                      py::arg("base_frame_orientation") = "ECLIPJ2000",
                      py::arg("time_step") = 300.0,
"");

                m.def(
                    "get_default_single_alternate_body_settings",
                    py::overload_cast<const std::string &, const std::string &,
                                      const std::string &>(
                        &tss::getDefaultSingleAlternateNameBodySettings),
                    py::arg("body_name"), py::arg("source_body_name"),
                    py::arg("base_frame_orientation") = "ECLIPJ2000",
"");

                m.def(
                    "get_default_single_alternate_body_settings_time_limited",
                    py::overload_cast<const std::string &, const std::string &,
                                      const double, const double,
                                      const std::string &, const double>(
                        &tss::getDefaultSingleAlternateNameBodySettings),
                    py::arg("body_name"), py::arg("source_body_name"),
                    py::arg("initial_time"), py::arg("final_time"),
                    py::arg("base_frame_orientation") = "ECLIPJ2000",
                    py::arg("time_step") = 300.0,
"");

                m.def("create_simplified_system_of_bodies",
                      &tss::createSimplifiedSystemOfBodies,
                      py::arg("initial_time") = 0,
                      R"doc(Function that creates a simplified System of bodies.

	Function that creates a simplified system of bodies. The following bodies are created in this system: the Sun, all planets of the Solar system, and Pluto.
	All bodies in this system use Gtop ephemerides and point mass gravity. The Earth is setup with a spherical shape model and a simple rotation model.
	The reference frame used to setup this simplified system of bodies has its origin at the SSB, and has an ECLIPJ2000 orientation.


	:param initial_time:
		Initial system time in seconds since J2000.
	:return:
		Object containing the objects for bodies and environment models constituting the physical environment
)doc");

                m.def(
                    "create_system_of_bodies",
                    &tss::createSystemOfBodies<double, TIME_TYPE>,
                    py::arg("body_settings"),
                    R"doc(Function that creates a System of bodies from associated settings.

	Function that creates a System of bodies from associated settings. This function creates the separate :class:`~tudatpy.numerical_simulation.Body`
	objects and stores them in a :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` object. This object represents the full
	physical environment in the simulation.


	:param body_settings:
		Object defining the physical environment, with all properties of artificial and natural bodies.
	:return:
		Object containing the objects for bodies and environment models constituting the physical environment
)doc");

                m.def("add_empty_tabulated_ephemeris",
                      &tp::addEmptyTabulatedEphemeris<double, TIME_TYPE>,
                      py::arg("bodies"), py::arg("body_name"),
                      py::arg("ephemeris_origin") = "",
                      py::arg("is_part_of_multi_arc") = false,
"");

                m.def(
                    "create_tabulated_ephemeris_from_spice",
                    &tss::createTabulatedEphemerisFromSpice<double, TIME_TYPE>,
                    py::arg("body"), py::arg("initial_time"),
                    py::arg("end_time"), py::arg("time_step"),
                    py::arg("observer_name"), py::arg("reference_frame_name"),
                    py::arg("interpolator_settings") = std::make_shared<
                        tudat::interpolators::LagrangeInterpolatorSettings>(8));

                m.def("create_body_ephemeris",
                      &tss::createBodyEphemeris<double, TIME_TYPE>,
                      py::arg("ephemeris_settings"), py::arg("body_name"),
"");

                m.def(
                    "create_ground_station_ephemeris",
                    py::overload_cast<const std::shared_ptr<tss::Body>,
                                      const std::string &>(
                        &tss::createReferencePointEphemeris<TIME_TYPE, double>),
                    py::arg("body"), py::arg("station_name"));

                m.def("get_safe_interpolation_interval",
                      &tss::getSafeInterpolationInterval,
                      py::arg("ephemeris_model"));


                m.def(
                    "add_aerodynamic_coefficient_interface",
                    &tss::addAerodynamicCoefficientInterface, py::arg("bodies"),
                    py::arg("body_name"), py::arg("coefficient_settings"),
                    R"doc(Function that creates an aerodynamic coefficient interface from settings, and adds it to an existing body.

	This function can be used to add an aerodynamic coefficient interface to an existing body. It requires
	settings for the aerodynamic coefficients, created using one of the factory functions from the `~tudatpy.numerical_simulation_environment_setup.aerodynamic_coefficient` module.
	This function creates the actual coefficient interface from these settings, and assigns it to the
	selected body. In addition to the identifier for the body to which it is assigned, this function
	requires the full :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` as input, to facilitate
	inter-body dependencies in the coefficient interface


	:param bodies:
		Object defining the physical environment, with all properties of artificial and natural bodies.
	:param body_name:
		Name of the body to which the aerodynamic coefficients are to be assigned
	:param coefficient_settings:
		Settings defining the coefficient interface that is to be created.
)doc");

                m.def("create_aerodynamic_coefficient_interface",
                      &tss::createAerodynamicCoefficientInterfaceDeprecated,
                      py::arg("coefficient_settings"), py::arg("body"));

                m.def("create_aerodynamic_coefficient_interface",
                      &tss::createAerodynamicCoefficientInterface,
                      py::arg("coefficient_settings"), py::arg("body"),
                      py::arg("bodies"),
"");

                m.def("add_radiation_pressure_interface",
                      &tss::addRadiationPressureInterface, py::arg("bodies"),
                      py::arg("body_name"),
                      py::arg("radiation_pressure_settings"));


                m.def(
                    "add_radiation_pressure_target_model",
                    &tss::addRadiationPressureTargetModel, py::arg("bodies"),
                    py::arg("body_name"),
                    py::arg("radiation_pressure_target_settings"),
                    R"doc(Function that creates an radiation pressure interface from settings, and adds it to an existing body.

	This function can be used to add an radiation pressure interface to an existing body. It requires
	settings for the radiation pressure interface, created using one of the factory functions from the :ref:`\`\`radiation_pressure\`\`` module.
	This function creates the actual coefficient interface from these settings, and assigns it to the
	selected body. In addition to the identifier for the body to which it is assigned, this function
	requires the full :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` as input, to facilitate
	inter-body dependencies in the radiation pressure interface


	:param bodies:
		Object defining the physical environment, with all properties of artificial and natural bodies.
	:param body_name:
		Name of the body to which the radiation pressure interface is to be assigned
	:param radiation_pressure_settings:
		Settings defining the radiation pressure interface that is to be created.
)doc");

                m.def(
                    "add_rotation_model", &tss::addRotationModel,
                    py::arg("bodies"), py::arg("body_name"),
                    py::arg("rotation_model_settings"),
                    R"doc(Function that creates a rotation model, and adds it to an existing body.

	This function can be used to add  a :class:`~tudatpy.numerical_simulation.environment.RotationalEphemeris` object to an existing body.
	Typically, the ``RotationalEphemeris`` is created along with the `~tudatpy.numerical_simulation.environment.Body` itself However, in some cases it may be useful
	to create a rotation model after the Body objects have been created. This function requires
	settings for the rotation model, created using one of the factory functions from the :ref:`~tudatpy.numerical_simulation_environment_setup.rotation_model` module.
	This function creates the actual coefficient interface from these settings, and assigns it to the
	selected body. In addition to the identifier for the body to which it is assigned, this function
	requires the full :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies` as input, to facilitate
	inter-body dependencies in the radiation model


	:param bodies:
		Object defining the physical environment, with all properties of artificial and natural bodies.
	:param body_name:
		Name of the body to which the rotation model is to be assigned
	:param rotation_model_settings:
		Settings defining the rotation model that is to be created.
)doc");

                m.def(
                    "add_gravity_field_model", &tss::addGravityFieldModel,
                    py::arg("bodies"), py::arg("body_name"),
                    py::arg("gravity_field_settings"),
                    py::arg("gravity_field_variation_settings") = std::vector<
                        std::shared_ptr<tss::GravityFieldVariationSettings>>(),
"");

                m.def("add_mass_properties_model", &tss::addRigidBodyProperties,
                      py::arg("bodies"), py::arg("body_name"),
                      py::arg("mass_property_settings"));

                m.def("add_rigid_body_properties", &tss::addRigidBodyProperties,
                      py::arg("bodies"), py::arg("body_name"),
                      py::arg("rigid_body_property_settings"),
"");


                m.def(
                    "add_engine_model", &tss::addEngineModel,
                    py::arg("body_name"), py::arg("engine_name"),
                    py::arg("thrust_magnitude_settings"), py::arg("bodies"),
                    py::arg("body_fixed_thrust_direction") =
                        Eigen::Vector3d::UnitX(),
                    R"doc(Function that creates an engine model (to be used for thrust calculations), and adds it to an existing body.

	Function that creates an engine model (to be used for thrust calculations), and adds it to an existing body. It creates and
	object of class :class:`~tudatpy.numerical_simulation.environment.EngineModel`, and adds it to an existing body. Properties
	assigned to this engine model are:
	* The (constant) direction in body-fixed frame in which the engine is pointing (e.g. the body-fixed thrust direction when the engine is on)
	* Settings for computing the thrust magnitude (as a function of time and/or other parameters), using a suitable function from the :ref:`\`\`thrust\`\`` submodule


	:param body_name:
		Name of the body to which the engine is to be added.
	:param engine_name:
		Name (e.g. unique identifier) of the engine that is to be added to the body
	:param thrust_magnitude_settings:
		Settings for computing the thrust magnitude (and specific impulse) as a function of time
	:param bodies:
		Object defining the physical environment, with all properties of artificial and natural bodies.
	:param body_name:
		Name of the body to which the rotation model is to be assigned
	:param body_fixed_thrust_direction:
		Unit vector along which the thrust from the engine will point in a body-fixed frame
)doc");

                m.def(
                    "add_variable_direction_engine_model",
                    &tss::addVariableDirectionEngineModel, py::arg("body_name"),
                    py::arg("engine_name"),
                    py::arg("thrust_magnitude_settings"), py::arg("bodies"),
                    py::arg("body_fixed_thrust_direction_function"),
                    R"doc(Function that creates an engine model (to be used for thrust calculations), and adds it to an existing body.

	Same as :func:`add_engine_model`, but with a time-variable body-fixed thrust direction


	:param body_name:
		Name of the body to which the engine is to be added.
	:param engine_name:
		Name (e.g. unique identifier) of the engine that is to be added to the body
	:param thrust_magnitude_settings:
		Settings for computing the thrust magnitude (and specific impulse) as a function of time
	:param bodies:
		Object defining the physical environment, with all properties of artificial and natural bodies.
	:param body_name:
		Name of the body to which the rotation model is to be assigned
	:param body_fixed_thrust_direction_function:
		Function returning a unit vector, as a function of time, along which the thrust from the engine will point in a body-fixed frame
)doc");

                m.def(
                    "add_flight_conditions", &tss::addFlightConditions,
                    py::arg("bodies"), py::arg("body_name"),
                    py::arg("central_body_name"),
                    R"doc(Function that creates a flight conditions, and adds it to an existing body.

	This function can be used to add  a :class:`~tudatpy.numerical_simulation.environment.FlightConditions` object to an existing body.
	Typically, the ``FlightConditions`` are created automatically when they are required (for the calulcation of an
	aerodynamic acceleration, or the saving of certain dependent variables). However, in some cases it may be useful
	to manually trigger their creation, which is done through this function. If the ``central_body_name`` input
	denotes a body that is endowed with an :class:`~tudatpy.numerical_simulation.environment.AtmosphereModel`, this function
	automically creates an :class:`~tudatpy.numerical_simulation.environment.AtmosphericFlightConditions` object (capable of
	calculating density, speed of sound, etc.), instead of the more basic :class:`~tudatpy.numerical_simulation.environment.FlightConditions`
	(which is limited to properties such as altitude, latitude, etc.)


	:param bodies:
		Object defining the physical environment, with all properties of artificial and natural bodies.
	:param body_name:
		Name of the body for which the flight conditions are to be created
	:param central_body_name:
		Name of the cenral body w.r.t. which the flight conditions are to be created (typically, but not necesarilly, the central body of propagation)/
)doc");

                m.def(
                    "convert_ground_station_state_between_itrf_frames",
                    &trf::convertGroundStationStateBetweenItrfFrames,
                    py::arg("ground_station_state"), py::arg("epoch"),
                    py::arg("base_frame"), py::arg("target_frame"),
"");

                m.def(
                    "add_ground_station",
                    py::overload_cast<
                        const std::shared_ptr<tss::Body>, const std::string,
                        const Eigen::Vector3d, const tcc::PositionElementTypes,
                        const std::vector<
                            std::shared_ptr<tss::GroundStationMotionSettings>>>(
                        &tss::createGroundStation),
                    py::arg("body"), py::arg("ground_station_name"),
                    py::arg("ground_station_position"),
                    py::arg("position_type") = tcc::cartesian_position,
                    py::arg("station_motion_settings") = std::vector<
                        std::shared_ptr<tss::GroundStationMotionSettings>>());

                m.def("add_ground_station",
                      py::overload_cast<
                          const std::shared_ptr<tss::Body>,
                          const std::shared_ptr<tss::GroundStationSettings>>(
                          &tss::createGroundStation),
                      py::arg("body"), py::arg("ground_station_settings"),
"");

                m.def("create_radiation_pressure_interface",
                      &tss::createRadiationPressureInterface,
                      py::arg("radiationPressureInterfaceSettings"),
                      py::arg("body_name"), py::arg("body_dict"));

                m.def("get_ground_station_list",
                      &tss::getGroundStationsLinkEndList, py::arg("body"));

                m.def("get_target_elevation_angles",
                      &tss::getTargetElevationAngles, py::arg("observing_body"),
                      py::arg("target_body"), py::arg("station_name"),
                      py::arg("times"));

                // Function removed; error is shown
                m.def("set_aerodynamic_guidance",
                      py::overload_cast<
                          const std::shared_ptr<ta::AerodynamicGuidance>,
                          const std::shared_ptr<tss::Body>, const bool>(
                          &tss::setGuidanceAnglesFunctions),
                      py::arg("aerodynamic_guidance"), py::arg("body"),
                      py::arg("silence_warnings") = false);

                // Function removed; error is shown
                m.def(
                    "set_aerodynamic_orientation_functions",
                    &tss::setAerodynamicOrientationFunctions, py::arg("body"),
                    py::arg("angle_of_attack_function") =
                        std::function<double()>(),
                    py::arg("sideslip_angle_function") =
                        std::function<double()>(),
                    py::arg("bank_angle_function") = std::function<double()>(),
                    py::arg("update_function") =
                        std::function<void(const double)>());

                // Function removed; error is shown
                m.def("set_constant_aerodynamic_orientation",
                      &tss::setConstantAerodynamicOrientation, py::arg("body"),
                      py::arg("angle_of_attack"), py::arg("sideslip_angle"),
                      py::arg("bank_angle"),
                      py::arg("silence_warnings") = false);
            }

        }  // namespace environment_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
