/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/basics/deprecationWarnings.h>
#include <tudat/simulation/propagation_setup.h>

#include "tudatpy/docstrings.h"

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;
namespace tmrf = tudat::root_finders;

namespace tudat {
    namespace simulation_setup {
        inline std::shared_ptr<TorqueSettings> customTorqueSettingsDeprecated(
            const std::function<Eigen::Vector3d(const double)> torqueFunction,
            const std::function<double(const double)> scalingFunction =
                nullptr) {
            static bool isWarningPrinted = false;
            if(isWarningPrinted == false) {
                tudat::utilities::printDeprecationWarning(
                    "tudatpy.numerical_simulation.propagation_setup."
                    "acceleration.custom",
                    "tudatpy.numerical_simulation.propagation_setup."
                    "acceleration.custom_torque");
                isWarningPrinted = true;
            }

            return customTorqueSettings(torqueFunction, scalingFunction);
        }
    }  // namespace simulation_setup
}  // namespace tudat
namespace tudatpy {
    namespace numerical_simulation {
        namespace propagation_setup {
            namespace torque {

                PYBIND11_MODULE(expose_torque, m) {
                    py::enum_<tba::AvailableTorque>(
                        m, "AvailableTorque",
R"doc(Enumeration of available torque types.

	Enumeration of torque types supported by tudat.
	

	:member torque_free_type:
	:member undefined_torque_type:
	:member second_order_gravitational_torque_type:
	:member aerodynamic_torque_type:
	:member spherical_harmonic_gravitational_torque_type:
	:member inertial_torque_type:
	:member dissipative_torque_type:
	:member custom_torque_type:
)doc")
                        .value("torque_free_type",
                               tba::AvailableTorque::torque_free,
get_docstring("AvailableTorque.torque_free_type").c_str())
                        .value("underfined_type",
                               tba::AvailableTorque::underfined_torque,
get_docstring("AvailableTorque.underfined_type").c_str())
                        .value("second_order_gravitational_type",
                               tba::AvailableTorque::
                                   second_order_gravitational_torque,
get_docstring("AvailableTorque.second_order_gravitational_type").c_str())
                        .value("aerodynamic_type",
                               tba::AvailableTorque::aerodynamic_torque,
get_docstring("AvailableTorque.aerodynamic_type").c_str())
                        .value("radiation_pressure_torque_type",
                               tba::AvailableTorque::radiation_pressure_torque,
get_docstring("AvailableTorque.radiation_pressure_torque_type").c_str())
                        .value("spherical_harmonic_gravitational_type",
                               tba::AvailableTorque::
                                   spherical_harmonic_gravitational_torque,
get_docstring("AvailableTorque.spherical_harmonic_gravitational_type").c_str())
                        .value("inertial_type",
                               tba::AvailableTorque::inertial_torque,
get_docstring("AvailableTorque.inertial_type").c_str())
                        .value("dissipative_type",
                               tba::AvailableTorque::dissipative_torque,
get_docstring("AvailableTorque.dissipative_type").c_str())
                        .export_values();

                    py::class_<tss::TorqueSettings,
                               std::shared_ptr<tss::TorqueSettings>>(
                        m, "TorqueSettings",
R"doc(Functional base class to define settings for torques.

	This is a functional base class to define settings for torques that require no information in addition to their type.
	Classes defining settings for torque models requiring additional information must be
	derived from this class.
	Bodies exerting and undergoing torque are set outside of this class.
	This class can be used for the easy setup of torque models
	(see createTorqueModels.h), but users may also chose to do so manually.
	(Derived) Class members are all public, for ease of access and modification.
	
)doc");

                    py::class_<
                        tss::SphericalHarmonicTorqueSettings,
                        std::shared_ptr<tss::SphericalHarmonicTorqueSettings>,
                        tss::TorqueSettings>(
                        m, "SphericalHarmonicTorqueSettings",
R"doc(`TorqueSettings`-derived class to define settings for torques caused by spherical harmonic gravity.

	`TorqueSettings`-derived class to define settings for torques caused by spherical harmonic gravity.
	
)doc");


                    m.def("aerodynamic", &tss::aerodynamicTorque,
R"doc(Creates the settings for the aerodynamic torque.

	Creates the settings for the aerodynamic torque exerted by a body with an atmosphere model and shape model on
	another body. The body exerting the torque needs to have both an atmosphere model and a shape model defined.
	Furthermore, the body undergoing the torque needs to have the aerodynamic coefficient interface and its moment
	coefficients defined. In the case that the aerodynamic coefficients are defined as a function of the vehicle
	orientation (e.g. angle of attack and sideslip angle), these angles can be manually or automatically defined.
	
	:return:
		Torque settings object.
)doc");

                    m.def("radiation_pressure_torque",
                          &tss::radiationPressureTorque,
get_docstring("radiation_pressure_torque").c_str());

                    m.def("second_degree_gravitational",
                          &tss::secondDegreeGravitationalTorque,
R"doc(Creates the settings for the second-degree gravitational torque.

	Torque exerted by a point mass on a body with a degree two spherical harmonics mass distribution.
	A degree two spherical harmonics mass distribution can be represented by an inertia tensor; thus,
	for this torque model, the body undergoing the torque needs to have an inertia tensor defined.
	The body exerting the torque only needs to have a gravitational model defined (either point-mass or spherical
	harmonics).
	
	:return:
		Torque settings object.
)doc");

                    m.def("spherical_harmonic_gravitational",
                          &tss::sphericalHarmonicGravitationalTorque,
                          py::arg("maximum_degree"), py::arg("maximum_order"),
R"doc(Creates the settings for the spherical harmonic torque.

	Torque exerted by a point mass on a body with an arbitrary degree/order spherical harmonics mass distribution.
	The body exerting the torque only needs to have a gravitational model defined (point-mass or spherical harmonic),
	while the body undergoing the torque needs to have a spherical harmonic gravity field defined.
	

	:param maximum_degree:
		Maximum degree of the spherical harmonic expansion.
	:param maximum_order:
		Maximum order of the spherical harmonic expansion.
	:return:
		Torque settings object.
)doc");

                    m.def("custom_torque", &tss::customTorqueSettings,
                          py::arg("torque_function"),
                          py::arg("scaling_function") = nullptr,
get_docstring("custom").c_str());

                    m.def("custom", &tss::customTorqueSettingsDeprecated,
                          py::arg("torque_function"),
                          py::arg("scaling_function") = nullptr);

                    // NOTE: the only unexposed torque model is
                    // dissipativeTorque, but it is probably obsolete
                }

            }  // namespace torque
        }  // namespace propagation_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
