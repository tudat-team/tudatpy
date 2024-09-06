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

namespace tudatpy {
    namespace numerical_simulation {
        namespace propagation_setup {

            PYBIND11_MODULE(expose_propagation_setup, m) {
                m.def(
                    "create_acceleration_models",
                    py::overload_cast<const tss::SystemOfBodies &,
                                      const tss::SelectedAccelerationMap &,
                                      const std::vector<std::string> &,
                                      const std::vector<std::string> &>(
                        &tss::createAccelerationModelsMap),
                    py::arg("body_system"),
                    py::arg("selected_acceleration_per_body"),
                    py::arg("bodies_to_propagate"), py::arg("central_bodies"),
                    R"doc(Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.

	Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
	bodies and central bodies are provided as two separate lists with the same order.


	:param body_system:
		System of bodies to be used in the propagation.
	:param selected_acceleration_per_body:
		Key-value container, with key denoting the body undergoing the acceleration, and the value containing an additional key-value container, with the body exerting acceleration, and list of acceleration settings exerted by this body.
	:param bodies_to_propagate:
		List of bodies to propagate.
	:param central_bodies:
		List of central bodies, each referred to each propagated body in the same order.
	:return:
		Set of accelerations acting on the bodies to propagate, provided as dual key-value container, similar to the acceleration settings input, but now with ``AccelerationModel`` lists as inner value
)doc");

                m.def(
                    "create_torque_models", &tss::createTorqueModelsMap,
                    py::arg("body_system"), py::arg("selected_torque_per_body"),
                    py::arg("bodies_to_propagate"),
                    R"doc(Function to create a set of acceleration models from a dictionary of bodies linked to acceleration model types.

	Function to create a set of acceleration models from a map of bodies and acceleration model types. The propagated
	bodies is provided as a list.


	:param body_system:
		System of bodies to be used in the propagation.
	:param selected_torque_per_body:
		Key-value container, with key denoting the body undergoing the torque, and the value containing an additional key-value container, with the body exerting torque, and list of torque settings exerted by this body.
	:param bodies_to_propagate:
		List of bodies to propagate.
	:return:
		Set of torques acting on the bodies to propagate, provided as dual key-value container, similar to the torque settings input, but now with ``TorqueModel`` lists as inner value
)doc");

                m.def(
                    "create_mass_rate_models", &tss::createMassRateModelsMap,
                    py::arg("body_system"),
                    py::arg("selected_mass_rates_per_body"),
                    py::arg("acceleration_models") = nullptr,
                    R"doc(Function to create a set of mass-rate models from associated settings.

	Function to create a set of mass-rate models from a map of bodies and mass-rate model types.
	If the mass-rate depends on any acceleration models (e.g. thrust), the acceleration
	models must be provided as an input.


	:param body_system:
		System of bodies to be used in the propagation.
	:param selected_mass_rates_per_body:
		Key-value container, with key denoting the body with changing mass, and the value containing a list of mass rate settings (in most cases, this list will have only a single entry)
	:param acceleration_models:
		Sorted list of acceleration models, as created by :func:`create_acceleration_models`
	:return:
		Set of mass-rate models, as key-value container, same as the settings input, with the difference that the rate settings objects have been processed into the associated objects calculating the actual mass-rate changes.
)doc");
            }
        }  // namespace propagation_setup
    }  // namespace numerical_simulation
}  // namespace tudatpy
